#include "ElectronRepulsionPrimRecSISF.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_sisf(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_sisf,
                                  size_t idx_eri_0_sgsf,
                                  size_t idx_eri_1_sgsf,
                                  size_t idx_eri_1_shsd,
                                  size_t idx_eri_0_shsf,
                                  size_t idx_eri_1_shsf,
                                  CSimdArray<double>& factors,
                                  const size_t idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double a_exp,
                                  const double b_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(WP) distances

    auto wp_x = factors.data(idx_wp);

    auto wp_y = factors.data(idx_wp + 1);

    auto wp_z = factors.data(idx_wp + 2);

    // set up R(PB) distances

    const auto xyz = r_pb.coordinates();

    const auto pb_x = xyz[0];

    const auto pb_y = xyz[1];

    const auto pb_z = xyz[2];

    /// Set up components of auxilary buffer : SGSF

    auto g_0_xxxx_0_xxx_0 = pbuffer.data(idx_eri_0_sgsf);

    auto g_0_xxxx_0_xxy_0 = pbuffer.data(idx_eri_0_sgsf + 1);

    auto g_0_xxxx_0_xxz_0 = pbuffer.data(idx_eri_0_sgsf + 2);

    auto g_0_xxxx_0_xyy_0 = pbuffer.data(idx_eri_0_sgsf + 3);

    auto g_0_xxxx_0_xyz_0 = pbuffer.data(idx_eri_0_sgsf + 4);

    auto g_0_xxxx_0_xzz_0 = pbuffer.data(idx_eri_0_sgsf + 5);

    auto g_0_xxxx_0_yyy_0 = pbuffer.data(idx_eri_0_sgsf + 6);

    auto g_0_xxxx_0_yyz_0 = pbuffer.data(idx_eri_0_sgsf + 7);

    auto g_0_xxxx_0_yzz_0 = pbuffer.data(idx_eri_0_sgsf + 8);

    auto g_0_xxxx_0_zzz_0 = pbuffer.data(idx_eri_0_sgsf + 9);

    auto g_0_xxxy_0_xxx_0 = pbuffer.data(idx_eri_0_sgsf + 10);

    auto g_0_xxxy_0_xxz_0 = pbuffer.data(idx_eri_0_sgsf + 12);

    auto g_0_xxxy_0_xzz_0 = pbuffer.data(idx_eri_0_sgsf + 15);

    auto g_0_xxxz_0_xxx_0 = pbuffer.data(idx_eri_0_sgsf + 20);

    auto g_0_xxxz_0_xxy_0 = pbuffer.data(idx_eri_0_sgsf + 21);

    auto g_0_xxxz_0_xyy_0 = pbuffer.data(idx_eri_0_sgsf + 23);

    auto g_0_xxyy_0_xxx_0 = pbuffer.data(idx_eri_0_sgsf + 30);

    auto g_0_xxyy_0_xxy_0 = pbuffer.data(idx_eri_0_sgsf + 31);

    auto g_0_xxyy_0_xxz_0 = pbuffer.data(idx_eri_0_sgsf + 32);

    auto g_0_xxyy_0_xyy_0 = pbuffer.data(idx_eri_0_sgsf + 33);

    auto g_0_xxyy_0_xyz_0 = pbuffer.data(idx_eri_0_sgsf + 34);

    auto g_0_xxyy_0_xzz_0 = pbuffer.data(idx_eri_0_sgsf + 35);

    auto g_0_xxyy_0_yyy_0 = pbuffer.data(idx_eri_0_sgsf + 36);

    auto g_0_xxyy_0_yyz_0 = pbuffer.data(idx_eri_0_sgsf + 37);

    auto g_0_xxyy_0_yzz_0 = pbuffer.data(idx_eri_0_sgsf + 38);

    auto g_0_xxyy_0_zzz_0 = pbuffer.data(idx_eri_0_sgsf + 39);

    auto g_0_xxzz_0_xxx_0 = pbuffer.data(idx_eri_0_sgsf + 50);

    auto g_0_xxzz_0_xxy_0 = pbuffer.data(idx_eri_0_sgsf + 51);

    auto g_0_xxzz_0_xxz_0 = pbuffer.data(idx_eri_0_sgsf + 52);

    auto g_0_xxzz_0_xyy_0 = pbuffer.data(idx_eri_0_sgsf + 53);

    auto g_0_xxzz_0_xyz_0 = pbuffer.data(idx_eri_0_sgsf + 54);

    auto g_0_xxzz_0_xzz_0 = pbuffer.data(idx_eri_0_sgsf + 55);

    auto g_0_xxzz_0_yyy_0 = pbuffer.data(idx_eri_0_sgsf + 56);

    auto g_0_xxzz_0_yyz_0 = pbuffer.data(idx_eri_0_sgsf + 57);

    auto g_0_xxzz_0_yzz_0 = pbuffer.data(idx_eri_0_sgsf + 58);

    auto g_0_xxzz_0_zzz_0 = pbuffer.data(idx_eri_0_sgsf + 59);

    auto g_0_xyyy_0_xxy_0 = pbuffer.data(idx_eri_0_sgsf + 61);

    auto g_0_xyyy_0_xyy_0 = pbuffer.data(idx_eri_0_sgsf + 63);

    auto g_0_xyyy_0_xyz_0 = pbuffer.data(idx_eri_0_sgsf + 64);

    auto g_0_xyyy_0_yyy_0 = pbuffer.data(idx_eri_0_sgsf + 66);

    auto g_0_xyyy_0_yyz_0 = pbuffer.data(idx_eri_0_sgsf + 67);

    auto g_0_xyyy_0_yzz_0 = pbuffer.data(idx_eri_0_sgsf + 68);

    auto g_0_xyyy_0_zzz_0 = pbuffer.data(idx_eri_0_sgsf + 69);

    auto g_0_xzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sgsf + 92);

    auto g_0_xzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sgsf + 94);

    auto g_0_xzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sgsf + 95);

    auto g_0_xzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sgsf + 96);

    auto g_0_xzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sgsf + 97);

    auto g_0_xzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sgsf + 98);

    auto g_0_xzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sgsf + 99);

    auto g_0_yyyy_0_xxx_0 = pbuffer.data(idx_eri_0_sgsf + 100);

    auto g_0_yyyy_0_xxy_0 = pbuffer.data(idx_eri_0_sgsf + 101);

    auto g_0_yyyy_0_xxz_0 = pbuffer.data(idx_eri_0_sgsf + 102);

    auto g_0_yyyy_0_xyy_0 = pbuffer.data(idx_eri_0_sgsf + 103);

    auto g_0_yyyy_0_xyz_0 = pbuffer.data(idx_eri_0_sgsf + 104);

    auto g_0_yyyy_0_xzz_0 = pbuffer.data(idx_eri_0_sgsf + 105);

    auto g_0_yyyy_0_yyy_0 = pbuffer.data(idx_eri_0_sgsf + 106);

    auto g_0_yyyy_0_yyz_0 = pbuffer.data(idx_eri_0_sgsf + 107);

    auto g_0_yyyy_0_yzz_0 = pbuffer.data(idx_eri_0_sgsf + 108);

    auto g_0_yyyy_0_zzz_0 = pbuffer.data(idx_eri_0_sgsf + 109);

    auto g_0_yyyz_0_xxy_0 = pbuffer.data(idx_eri_0_sgsf + 111);

    auto g_0_yyyz_0_xyy_0 = pbuffer.data(idx_eri_0_sgsf + 113);

    auto g_0_yyyz_0_yyy_0 = pbuffer.data(idx_eri_0_sgsf + 116);

    auto g_0_yyzz_0_xxx_0 = pbuffer.data(idx_eri_0_sgsf + 120);

    auto g_0_yyzz_0_xxy_0 = pbuffer.data(idx_eri_0_sgsf + 121);

    auto g_0_yyzz_0_xxz_0 = pbuffer.data(idx_eri_0_sgsf + 122);

    auto g_0_yyzz_0_xyy_0 = pbuffer.data(idx_eri_0_sgsf + 123);

    auto g_0_yyzz_0_xyz_0 = pbuffer.data(idx_eri_0_sgsf + 124);

    auto g_0_yyzz_0_xzz_0 = pbuffer.data(idx_eri_0_sgsf + 125);

    auto g_0_yyzz_0_yyy_0 = pbuffer.data(idx_eri_0_sgsf + 126);

    auto g_0_yyzz_0_yyz_0 = pbuffer.data(idx_eri_0_sgsf + 127);

    auto g_0_yyzz_0_yzz_0 = pbuffer.data(idx_eri_0_sgsf + 128);

    auto g_0_yyzz_0_zzz_0 = pbuffer.data(idx_eri_0_sgsf + 129);

    auto g_0_yzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sgsf + 130);

    auto g_0_yzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sgsf + 132);

    auto g_0_yzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sgsf + 134);

    auto g_0_yzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sgsf + 135);

    auto g_0_yzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sgsf + 137);

    auto g_0_yzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sgsf + 138);

    auto g_0_yzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sgsf + 139);

    auto g_0_zzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sgsf + 140);

    auto g_0_zzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sgsf + 141);

    auto g_0_zzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sgsf + 142);

    auto g_0_zzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sgsf + 143);

    auto g_0_zzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sgsf + 144);

    auto g_0_zzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sgsf + 145);

    auto g_0_zzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sgsf + 146);

    auto g_0_zzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sgsf + 147);

    auto g_0_zzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sgsf + 148);

    auto g_0_zzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sgsf + 149);

    /// Set up components of auxilary buffer : SGSF

    auto g_0_xxxx_0_xxx_1 = pbuffer.data(idx_eri_1_sgsf);

    auto g_0_xxxx_0_xxy_1 = pbuffer.data(idx_eri_1_sgsf + 1);

    auto g_0_xxxx_0_xxz_1 = pbuffer.data(idx_eri_1_sgsf + 2);

    auto g_0_xxxx_0_xyy_1 = pbuffer.data(idx_eri_1_sgsf + 3);

    auto g_0_xxxx_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 4);

    auto g_0_xxxx_0_xzz_1 = pbuffer.data(idx_eri_1_sgsf + 5);

    auto g_0_xxxx_0_yyy_1 = pbuffer.data(idx_eri_1_sgsf + 6);

    auto g_0_xxxx_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 7);

    auto g_0_xxxx_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 8);

    auto g_0_xxxx_0_zzz_1 = pbuffer.data(idx_eri_1_sgsf + 9);

    auto g_0_xxxy_0_xxx_1 = pbuffer.data(idx_eri_1_sgsf + 10);

    auto g_0_xxxy_0_xxz_1 = pbuffer.data(idx_eri_1_sgsf + 12);

    auto g_0_xxxy_0_xzz_1 = pbuffer.data(idx_eri_1_sgsf + 15);

    auto g_0_xxxz_0_xxx_1 = pbuffer.data(idx_eri_1_sgsf + 20);

    auto g_0_xxxz_0_xxy_1 = pbuffer.data(idx_eri_1_sgsf + 21);

    auto g_0_xxxz_0_xyy_1 = pbuffer.data(idx_eri_1_sgsf + 23);

    auto g_0_xxyy_0_xxx_1 = pbuffer.data(idx_eri_1_sgsf + 30);

    auto g_0_xxyy_0_xxy_1 = pbuffer.data(idx_eri_1_sgsf + 31);

    auto g_0_xxyy_0_xxz_1 = pbuffer.data(idx_eri_1_sgsf + 32);

    auto g_0_xxyy_0_xyy_1 = pbuffer.data(idx_eri_1_sgsf + 33);

    auto g_0_xxyy_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 34);

    auto g_0_xxyy_0_xzz_1 = pbuffer.data(idx_eri_1_sgsf + 35);

    auto g_0_xxyy_0_yyy_1 = pbuffer.data(idx_eri_1_sgsf + 36);

    auto g_0_xxyy_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 37);

    auto g_0_xxyy_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 38);

    auto g_0_xxyy_0_zzz_1 = pbuffer.data(idx_eri_1_sgsf + 39);

    auto g_0_xxzz_0_xxx_1 = pbuffer.data(idx_eri_1_sgsf + 50);

    auto g_0_xxzz_0_xxy_1 = pbuffer.data(idx_eri_1_sgsf + 51);

    auto g_0_xxzz_0_xxz_1 = pbuffer.data(idx_eri_1_sgsf + 52);

    auto g_0_xxzz_0_xyy_1 = pbuffer.data(idx_eri_1_sgsf + 53);

    auto g_0_xxzz_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 54);

    auto g_0_xxzz_0_xzz_1 = pbuffer.data(idx_eri_1_sgsf + 55);

    auto g_0_xxzz_0_yyy_1 = pbuffer.data(idx_eri_1_sgsf + 56);

    auto g_0_xxzz_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 57);

    auto g_0_xxzz_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 58);

    auto g_0_xxzz_0_zzz_1 = pbuffer.data(idx_eri_1_sgsf + 59);

    auto g_0_xyyy_0_xxy_1 = pbuffer.data(idx_eri_1_sgsf + 61);

    auto g_0_xyyy_0_xyy_1 = pbuffer.data(idx_eri_1_sgsf + 63);

    auto g_0_xyyy_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 64);

    auto g_0_xyyy_0_yyy_1 = pbuffer.data(idx_eri_1_sgsf + 66);

    auto g_0_xyyy_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 67);

    auto g_0_xyyy_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 68);

    auto g_0_xyyy_0_zzz_1 = pbuffer.data(idx_eri_1_sgsf + 69);

    auto g_0_xzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sgsf + 92);

    auto g_0_xzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 94);

    auto g_0_xzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sgsf + 95);

    auto g_0_xzzz_0_yyy_1 = pbuffer.data(idx_eri_1_sgsf + 96);

    auto g_0_xzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 97);

    auto g_0_xzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 98);

    auto g_0_xzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sgsf + 99);

    auto g_0_yyyy_0_xxx_1 = pbuffer.data(idx_eri_1_sgsf + 100);

    auto g_0_yyyy_0_xxy_1 = pbuffer.data(idx_eri_1_sgsf + 101);

    auto g_0_yyyy_0_xxz_1 = pbuffer.data(idx_eri_1_sgsf + 102);

    auto g_0_yyyy_0_xyy_1 = pbuffer.data(idx_eri_1_sgsf + 103);

    auto g_0_yyyy_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 104);

    auto g_0_yyyy_0_xzz_1 = pbuffer.data(idx_eri_1_sgsf + 105);

    auto g_0_yyyy_0_yyy_1 = pbuffer.data(idx_eri_1_sgsf + 106);

    auto g_0_yyyy_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 107);

    auto g_0_yyyy_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 108);

    auto g_0_yyyy_0_zzz_1 = pbuffer.data(idx_eri_1_sgsf + 109);

    auto g_0_yyyz_0_xxy_1 = pbuffer.data(idx_eri_1_sgsf + 111);

    auto g_0_yyyz_0_xyy_1 = pbuffer.data(idx_eri_1_sgsf + 113);

    auto g_0_yyyz_0_yyy_1 = pbuffer.data(idx_eri_1_sgsf + 116);

    auto g_0_yyzz_0_xxx_1 = pbuffer.data(idx_eri_1_sgsf + 120);

    auto g_0_yyzz_0_xxy_1 = pbuffer.data(idx_eri_1_sgsf + 121);

    auto g_0_yyzz_0_xxz_1 = pbuffer.data(idx_eri_1_sgsf + 122);

    auto g_0_yyzz_0_xyy_1 = pbuffer.data(idx_eri_1_sgsf + 123);

    auto g_0_yyzz_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 124);

    auto g_0_yyzz_0_xzz_1 = pbuffer.data(idx_eri_1_sgsf + 125);

    auto g_0_yyzz_0_yyy_1 = pbuffer.data(idx_eri_1_sgsf + 126);

    auto g_0_yyzz_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 127);

    auto g_0_yyzz_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 128);

    auto g_0_yyzz_0_zzz_1 = pbuffer.data(idx_eri_1_sgsf + 129);

    auto g_0_yzzz_0_xxx_1 = pbuffer.data(idx_eri_1_sgsf + 130);

    auto g_0_yzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sgsf + 132);

    auto g_0_yzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 134);

    auto g_0_yzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sgsf + 135);

    auto g_0_yzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 137);

    auto g_0_yzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 138);

    auto g_0_yzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sgsf + 139);

    auto g_0_zzzz_0_xxx_1 = pbuffer.data(idx_eri_1_sgsf + 140);

    auto g_0_zzzz_0_xxy_1 = pbuffer.data(idx_eri_1_sgsf + 141);

    auto g_0_zzzz_0_xxz_1 = pbuffer.data(idx_eri_1_sgsf + 142);

    auto g_0_zzzz_0_xyy_1 = pbuffer.data(idx_eri_1_sgsf + 143);

    auto g_0_zzzz_0_xyz_1 = pbuffer.data(idx_eri_1_sgsf + 144);

    auto g_0_zzzz_0_xzz_1 = pbuffer.data(idx_eri_1_sgsf + 145);

    auto g_0_zzzz_0_yyy_1 = pbuffer.data(idx_eri_1_sgsf + 146);

    auto g_0_zzzz_0_yyz_1 = pbuffer.data(idx_eri_1_sgsf + 147);

    auto g_0_zzzz_0_yzz_1 = pbuffer.data(idx_eri_1_sgsf + 148);

    auto g_0_zzzz_0_zzz_1 = pbuffer.data(idx_eri_1_sgsf + 149);

    /// Set up components of auxilary buffer : SHSD

    auto g_0_xxxxx_0_xx_1 = pbuffer.data(idx_eri_1_shsd);

    auto g_0_xxxxx_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 1);

    auto g_0_xxxxx_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 2);

    auto g_0_xxxxx_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 3);

    auto g_0_xxxxx_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 4);

    auto g_0_xxxxx_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 5);

    auto g_0_xxxxz_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 14);

    auto g_0_xxxxz_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 16);

    auto g_0_xxxxz_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 17);

    auto g_0_xxxyy_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 18);

    auto g_0_xxxyy_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 19);

    auto g_0_xxxyy_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 20);

    auto g_0_xxxyy_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 21);

    auto g_0_xxxyy_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 22);

    auto g_0_xxxyy_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 23);

    auto g_0_xxxzz_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 30);

    auto g_0_xxxzz_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 31);

    auto g_0_xxxzz_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 32);

    auto g_0_xxxzz_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 33);

    auto g_0_xxxzz_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 34);

    auto g_0_xxxzz_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 35);

    auto g_0_xxyyy_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 36);

    auto g_0_xxyyy_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 37);

    auto g_0_xxyyy_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 38);

    auto g_0_xxyyy_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 39);

    auto g_0_xxyyy_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 40);

    auto g_0_xxyyy_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 41);

    auto g_0_xxzzz_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 54);

    auto g_0_xxzzz_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 55);

    auto g_0_xxzzz_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 56);

    auto g_0_xxzzz_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 57);

    auto g_0_xxzzz_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 58);

    auto g_0_xxzzz_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 59);

    auto g_0_xyyyy_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 61);

    auto g_0_xyyyy_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 63);

    auto g_0_xyyyy_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 64);

    auto g_0_xyyzz_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 76);

    auto g_0_xzzzz_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 86);

    auto g_0_xzzzz_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 88);

    auto g_0_xzzzz_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 89);

    auto g_0_yyyyy_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 90);

    auto g_0_yyyyy_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 91);

    auto g_0_yyyyy_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 92);

    auto g_0_yyyyy_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 93);

    auto g_0_yyyyy_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 94);

    auto g_0_yyyyy_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 95);

    auto g_0_yyyyz_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 98);

    auto g_0_yyyyz_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 100);

    auto g_0_yyyyz_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 101);

    auto g_0_yyyzz_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 102);

    auto g_0_yyyzz_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 103);

    auto g_0_yyyzz_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 104);

    auto g_0_yyyzz_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 105);

    auto g_0_yyyzz_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 106);

    auto g_0_yyyzz_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 107);

    auto g_0_yyzzz_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 108);

    auto g_0_yyzzz_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 109);

    auto g_0_yyzzz_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 110);

    auto g_0_yyzzz_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 111);

    auto g_0_yyzzz_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 112);

    auto g_0_yyzzz_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 113);

    auto g_0_yzzzz_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 115);

    auto g_0_yzzzz_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 116);

    auto g_0_yzzzz_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 117);

    auto g_0_yzzzz_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 118);

    auto g_0_yzzzz_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 119);

    auto g_0_zzzzz_0_xx_1 = pbuffer.data(idx_eri_1_shsd + 120);

    auto g_0_zzzzz_0_xy_1 = pbuffer.data(idx_eri_1_shsd + 121);

    auto g_0_zzzzz_0_xz_1 = pbuffer.data(idx_eri_1_shsd + 122);

    auto g_0_zzzzz_0_yy_1 = pbuffer.data(idx_eri_1_shsd + 123);

    auto g_0_zzzzz_0_yz_1 = pbuffer.data(idx_eri_1_shsd + 124);

    auto g_0_zzzzz_0_zz_1 = pbuffer.data(idx_eri_1_shsd + 125);

    /// Set up components of auxilary buffer : SHSF

    auto g_0_xxxxx_0_xxx_0 = pbuffer.data(idx_eri_0_shsf);

    auto g_0_xxxxx_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 1);

    auto g_0_xxxxx_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 2);

    auto g_0_xxxxx_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 3);

    auto g_0_xxxxx_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 4);

    auto g_0_xxxxx_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 5);

    auto g_0_xxxxx_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 6);

    auto g_0_xxxxx_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 7);

    auto g_0_xxxxx_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 8);

    auto g_0_xxxxx_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 9);

    auto g_0_xxxxy_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 10);

    auto g_0_xxxxy_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 11);

    auto g_0_xxxxy_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 12);

    auto g_0_xxxxy_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 13);

    auto g_0_xxxxy_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 15);

    auto g_0_xxxxy_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 16);

    auto g_0_xxxxz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 20);

    auto g_0_xxxxz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 21);

    auto g_0_xxxxz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 22);

    auto g_0_xxxxz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 23);

    auto g_0_xxxxz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 24);

    auto g_0_xxxxz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 25);

    auto g_0_xxxxz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 27);

    auto g_0_xxxxz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 28);

    auto g_0_xxxxz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 29);

    auto g_0_xxxyy_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 30);

    auto g_0_xxxyy_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 31);

    auto g_0_xxxyy_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 32);

    auto g_0_xxxyy_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 33);

    auto g_0_xxxyy_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 34);

    auto g_0_xxxyy_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 35);

    auto g_0_xxxyy_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 36);

    auto g_0_xxxyy_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 37);

    auto g_0_xxxyy_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 38);

    auto g_0_xxxyy_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 39);

    auto g_0_xxxzz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 50);

    auto g_0_xxxzz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 51);

    auto g_0_xxxzz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 52);

    auto g_0_xxxzz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 53);

    auto g_0_xxxzz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 54);

    auto g_0_xxxzz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 55);

    auto g_0_xxxzz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 56);

    auto g_0_xxxzz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 57);

    auto g_0_xxxzz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 58);

    auto g_0_xxxzz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 59);

    auto g_0_xxyyy_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 60);

    auto g_0_xxyyy_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 61);

    auto g_0_xxyyy_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 62);

    auto g_0_xxyyy_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 63);

    auto g_0_xxyyy_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 64);

    auto g_0_xxyyy_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 65);

    auto g_0_xxyyy_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 66);

    auto g_0_xxyyy_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 67);

    auto g_0_xxyyy_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 68);

    auto g_0_xxyyy_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 69);

    auto g_0_xxyyz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 71);

    auto g_0_xxyyz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 73);

    auto g_0_xxyzz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 80);

    auto g_0_xxyzz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 82);

    auto g_0_xxyzz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 85);

    auto g_0_xxzzz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 90);

    auto g_0_xxzzz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 91);

    auto g_0_xxzzz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 92);

    auto g_0_xxzzz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 93);

    auto g_0_xxzzz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 94);

    auto g_0_xxzzz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 95);

    auto g_0_xxzzz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 96);

    auto g_0_xxzzz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 97);

    auto g_0_xxzzz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 98);

    auto g_0_xxzzz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 99);

    auto g_0_xyyyy_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 100);

    auto g_0_xyyyy_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 101);

    auto g_0_xyyyy_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 103);

    auto g_0_xyyyy_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 104);

    auto g_0_xyyyy_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 106);

    auto g_0_xyyyy_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 107);

    auto g_0_xyyyy_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 108);

    auto g_0_xyyyy_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 109);

    auto g_0_xyyzz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 124);

    auto g_0_xyyzz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 126);

    auto g_0_xyyzz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 127);

    auto g_0_xyyzz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 128);

    auto g_0_xyyzz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 129);

    auto g_0_xzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 140);

    auto g_0_xzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 142);

    auto g_0_xzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 144);

    auto g_0_xzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 145);

    auto g_0_xzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 146);

    auto g_0_xzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 147);

    auto g_0_xzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 148);

    auto g_0_xzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 149);

    auto g_0_yyyyy_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 150);

    auto g_0_yyyyy_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 151);

    auto g_0_yyyyy_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 152);

    auto g_0_yyyyy_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 153);

    auto g_0_yyyyy_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 154);

    auto g_0_yyyyy_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 155);

    auto g_0_yyyyy_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 156);

    auto g_0_yyyyy_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 157);

    auto g_0_yyyyy_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 158);

    auto g_0_yyyyy_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 159);

    auto g_0_yyyyz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 161);

    auto g_0_yyyyz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 162);

    auto g_0_yyyyz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 163);

    auto g_0_yyyyz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 164);

    auto g_0_yyyyz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 165);

    auto g_0_yyyyz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 166);

    auto g_0_yyyyz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 167);

    auto g_0_yyyyz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 168);

    auto g_0_yyyyz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 169);

    auto g_0_yyyzz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 170);

    auto g_0_yyyzz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 171);

    auto g_0_yyyzz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 172);

    auto g_0_yyyzz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 173);

    auto g_0_yyyzz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 174);

    auto g_0_yyyzz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 175);

    auto g_0_yyyzz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 176);

    auto g_0_yyyzz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 177);

    auto g_0_yyyzz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 178);

    auto g_0_yyyzz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 179);

    auto g_0_yyzzz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 180);

    auto g_0_yyzzz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 181);

    auto g_0_yyzzz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 182);

    auto g_0_yyzzz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 183);

    auto g_0_yyzzz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 184);

    auto g_0_yyzzz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 185);

    auto g_0_yyzzz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 186);

    auto g_0_yyzzz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 187);

    auto g_0_yyzzz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 188);

    auto g_0_yyzzz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 189);

    auto g_0_yzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 190);

    auto g_0_yzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 191);

    auto g_0_yzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 192);

    auto g_0_yzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 193);

    auto g_0_yzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 194);

    auto g_0_yzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 195);

    auto g_0_yzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 196);

    auto g_0_yzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 197);

    auto g_0_yzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 198);

    auto g_0_yzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 199);

    auto g_0_zzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_shsf + 200);

    auto g_0_zzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_shsf + 201);

    auto g_0_zzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_shsf + 202);

    auto g_0_zzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_shsf + 203);

    auto g_0_zzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_shsf + 204);

    auto g_0_zzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_shsf + 205);

    auto g_0_zzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_shsf + 206);

    auto g_0_zzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_shsf + 207);

    auto g_0_zzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_shsf + 208);

    auto g_0_zzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_shsf + 209);

    /// Set up components of auxilary buffer : SHSF

    auto g_0_xxxxx_0_xxx_1 = pbuffer.data(idx_eri_1_shsf);

    auto g_0_xxxxx_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 1);

    auto g_0_xxxxx_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 2);

    auto g_0_xxxxx_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 3);

    auto g_0_xxxxx_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 4);

    auto g_0_xxxxx_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 5);

    auto g_0_xxxxx_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 6);

    auto g_0_xxxxx_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 7);

    auto g_0_xxxxx_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 8);

    auto g_0_xxxxx_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 9);

    auto g_0_xxxxy_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 10);

    auto g_0_xxxxy_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 11);

    auto g_0_xxxxy_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 12);

    auto g_0_xxxxy_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 13);

    auto g_0_xxxxy_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 15);

    auto g_0_xxxxy_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 16);

    auto g_0_xxxxz_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 20);

    auto g_0_xxxxz_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 21);

    auto g_0_xxxxz_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 22);

    auto g_0_xxxxz_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 23);

    auto g_0_xxxxz_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 24);

    auto g_0_xxxxz_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 25);

    auto g_0_xxxxz_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 27);

    auto g_0_xxxxz_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 28);

    auto g_0_xxxxz_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 29);

    auto g_0_xxxyy_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 30);

    auto g_0_xxxyy_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 31);

    auto g_0_xxxyy_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 32);

    auto g_0_xxxyy_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 33);

    auto g_0_xxxyy_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 34);

    auto g_0_xxxyy_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 35);

    auto g_0_xxxyy_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 36);

    auto g_0_xxxyy_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 37);

    auto g_0_xxxyy_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 38);

    auto g_0_xxxyy_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 39);

    auto g_0_xxxzz_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 50);

    auto g_0_xxxzz_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 51);

    auto g_0_xxxzz_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 52);

    auto g_0_xxxzz_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 53);

    auto g_0_xxxzz_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 54);

    auto g_0_xxxzz_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 55);

    auto g_0_xxxzz_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 56);

    auto g_0_xxxzz_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 57);

    auto g_0_xxxzz_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 58);

    auto g_0_xxxzz_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 59);

    auto g_0_xxyyy_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 60);

    auto g_0_xxyyy_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 61);

    auto g_0_xxyyy_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 62);

    auto g_0_xxyyy_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 63);

    auto g_0_xxyyy_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 64);

    auto g_0_xxyyy_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 65);

    auto g_0_xxyyy_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 66);

    auto g_0_xxyyy_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 67);

    auto g_0_xxyyy_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 68);

    auto g_0_xxyyy_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 69);

    auto g_0_xxyyz_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 71);

    auto g_0_xxyyz_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 73);

    auto g_0_xxyzz_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 80);

    auto g_0_xxyzz_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 82);

    auto g_0_xxyzz_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 85);

    auto g_0_xxzzz_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 90);

    auto g_0_xxzzz_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 91);

    auto g_0_xxzzz_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 92);

    auto g_0_xxzzz_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 93);

    auto g_0_xxzzz_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 94);

    auto g_0_xxzzz_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 95);

    auto g_0_xxzzz_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 96);

    auto g_0_xxzzz_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 97);

    auto g_0_xxzzz_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 98);

    auto g_0_xxzzz_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 99);

    auto g_0_xyyyy_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 100);

    auto g_0_xyyyy_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 101);

    auto g_0_xyyyy_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 103);

    auto g_0_xyyyy_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 104);

    auto g_0_xyyyy_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 106);

    auto g_0_xyyyy_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 107);

    auto g_0_xyyyy_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 108);

    auto g_0_xyyyy_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 109);

    auto g_0_xyyzz_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 124);

    auto g_0_xyyzz_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 126);

    auto g_0_xyyzz_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 127);

    auto g_0_xyyzz_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 128);

    auto g_0_xyyzz_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 129);

    auto g_0_xzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 140);

    auto g_0_xzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 142);

    auto g_0_xzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 144);

    auto g_0_xzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 145);

    auto g_0_xzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 146);

    auto g_0_xzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 147);

    auto g_0_xzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 148);

    auto g_0_xzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 149);

    auto g_0_yyyyy_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 150);

    auto g_0_yyyyy_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 151);

    auto g_0_yyyyy_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 152);

    auto g_0_yyyyy_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 153);

    auto g_0_yyyyy_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 154);

    auto g_0_yyyyy_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 155);

    auto g_0_yyyyy_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 156);

    auto g_0_yyyyy_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 157);

    auto g_0_yyyyy_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 158);

    auto g_0_yyyyy_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 159);

    auto g_0_yyyyz_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 161);

    auto g_0_yyyyz_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 162);

    auto g_0_yyyyz_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 163);

    auto g_0_yyyyz_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 164);

    auto g_0_yyyyz_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 165);

    auto g_0_yyyyz_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 166);

    auto g_0_yyyyz_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 167);

    auto g_0_yyyyz_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 168);

    auto g_0_yyyyz_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 169);

    auto g_0_yyyzz_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 170);

    auto g_0_yyyzz_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 171);

    auto g_0_yyyzz_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 172);

    auto g_0_yyyzz_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 173);

    auto g_0_yyyzz_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 174);

    auto g_0_yyyzz_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 175);

    auto g_0_yyyzz_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 176);

    auto g_0_yyyzz_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 177);

    auto g_0_yyyzz_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 178);

    auto g_0_yyyzz_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 179);

    auto g_0_yyzzz_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 180);

    auto g_0_yyzzz_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 181);

    auto g_0_yyzzz_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 182);

    auto g_0_yyzzz_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 183);

    auto g_0_yyzzz_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 184);

    auto g_0_yyzzz_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 185);

    auto g_0_yyzzz_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 186);

    auto g_0_yyzzz_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 187);

    auto g_0_yyzzz_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 188);

    auto g_0_yyzzz_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 189);

    auto g_0_yzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 190);

    auto g_0_yzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 191);

    auto g_0_yzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 192);

    auto g_0_yzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 193);

    auto g_0_yzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 194);

    auto g_0_yzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 195);

    auto g_0_yzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 196);

    auto g_0_yzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 197);

    auto g_0_yzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 198);

    auto g_0_yzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 199);

    auto g_0_zzzzz_0_xxx_1 = pbuffer.data(idx_eri_1_shsf + 200);

    auto g_0_zzzzz_0_xxy_1 = pbuffer.data(idx_eri_1_shsf + 201);

    auto g_0_zzzzz_0_xxz_1 = pbuffer.data(idx_eri_1_shsf + 202);

    auto g_0_zzzzz_0_xyy_1 = pbuffer.data(idx_eri_1_shsf + 203);

    auto g_0_zzzzz_0_xyz_1 = pbuffer.data(idx_eri_1_shsf + 204);

    auto g_0_zzzzz_0_xzz_1 = pbuffer.data(idx_eri_1_shsf + 205);

    auto g_0_zzzzz_0_yyy_1 = pbuffer.data(idx_eri_1_shsf + 206);

    auto g_0_zzzzz_0_yyz_1 = pbuffer.data(idx_eri_1_shsf + 207);

    auto g_0_zzzzz_0_yzz_1 = pbuffer.data(idx_eri_1_shsf + 208);

    auto g_0_zzzzz_0_zzz_1 = pbuffer.data(idx_eri_1_shsf + 209);

    /// Set up 0-10 components of targeted buffer : SISF

    auto g_0_xxxxxx_0_xxx_0 = pbuffer.data(idx_eri_0_sisf);

    auto g_0_xxxxxx_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 1);

    auto g_0_xxxxxx_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 2);

    auto g_0_xxxxxx_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 3);

    auto g_0_xxxxxx_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 4);

    auto g_0_xxxxxx_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 5);

    auto g_0_xxxxxx_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 6);

    auto g_0_xxxxxx_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 7);

    auto g_0_xxxxxx_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 8);

    auto g_0_xxxxxx_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 9);

    #pragma omp simd aligned(g_0_xxxx_0_xxx_0, g_0_xxxx_0_xxx_1, g_0_xxxx_0_xxy_0, g_0_xxxx_0_xxy_1, g_0_xxxx_0_xxz_0, g_0_xxxx_0_xxz_1, g_0_xxxx_0_xyy_0, g_0_xxxx_0_xyy_1, g_0_xxxx_0_xyz_0, g_0_xxxx_0_xyz_1, g_0_xxxx_0_xzz_0, g_0_xxxx_0_xzz_1, g_0_xxxx_0_yyy_0, g_0_xxxx_0_yyy_1, g_0_xxxx_0_yyz_0, g_0_xxxx_0_yyz_1, g_0_xxxx_0_yzz_0, g_0_xxxx_0_yzz_1, g_0_xxxx_0_zzz_0, g_0_xxxx_0_zzz_1, g_0_xxxxx_0_xx_1, g_0_xxxxx_0_xxx_0, g_0_xxxxx_0_xxx_1, g_0_xxxxx_0_xxy_0, g_0_xxxxx_0_xxy_1, g_0_xxxxx_0_xxz_0, g_0_xxxxx_0_xxz_1, g_0_xxxxx_0_xy_1, g_0_xxxxx_0_xyy_0, g_0_xxxxx_0_xyy_1, g_0_xxxxx_0_xyz_0, g_0_xxxxx_0_xyz_1, g_0_xxxxx_0_xz_1, g_0_xxxxx_0_xzz_0, g_0_xxxxx_0_xzz_1, g_0_xxxxx_0_yy_1, g_0_xxxxx_0_yyy_0, g_0_xxxxx_0_yyy_1, g_0_xxxxx_0_yyz_0, g_0_xxxxx_0_yyz_1, g_0_xxxxx_0_yz_1, g_0_xxxxx_0_yzz_0, g_0_xxxxx_0_yzz_1, g_0_xxxxx_0_zz_1, g_0_xxxxx_0_zzz_0, g_0_xxxxx_0_zzz_1, g_0_xxxxxx_0_xxx_0, g_0_xxxxxx_0_xxy_0, g_0_xxxxxx_0_xxz_0, g_0_xxxxxx_0_xyy_0, g_0_xxxxxx_0_xyz_0, g_0_xxxxxx_0_xzz_0, g_0_xxxxxx_0_yyy_0, g_0_xxxxxx_0_yyz_0, g_0_xxxxxx_0_yzz_0, g_0_xxxxxx_0_zzz_0, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxx_0_xxx_0[i] = 5.0 * g_0_xxxx_0_xxx_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxx_1[i] * fti_ab_0 + 3.0 * g_0_xxxxx_0_xx_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxx_0[i] * pb_x + g_0_xxxxx_0_xxx_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxy_0[i] = 5.0 * g_0_xxxx_0_xxy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxy_1[i] * fti_ab_0 + 2.0 * g_0_xxxxx_0_xy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxy_0[i] * pb_x + g_0_xxxxx_0_xxy_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxz_0[i] = 5.0 * g_0_xxxx_0_xxz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxx_0_xz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxz_0[i] * pb_x + g_0_xxxxx_0_xxz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xyy_0[i] = 5.0 * g_0_xxxx_0_xyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xyy_1[i] * fti_ab_0 + g_0_xxxxx_0_yy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyy_0[i] * pb_x + g_0_xxxxx_0_xyy_1[i] * wp_x[i];

        g_0_xxxxxx_0_xyz_0[i] = 5.0 * g_0_xxxx_0_xyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xyz_1[i] * fti_ab_0 + g_0_xxxxx_0_yz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyz_0[i] * pb_x + g_0_xxxxx_0_xyz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xzz_0[i] = 5.0 * g_0_xxxx_0_xzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xzz_1[i] * fti_ab_0 + g_0_xxxxx_0_zz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xzz_0[i] * pb_x + g_0_xxxxx_0_xzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_yyy_0[i] = 5.0 * g_0_xxxx_0_yyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yyy_1[i] * fti_ab_0 + g_0_xxxxx_0_yyy_0[i] * pb_x + g_0_xxxxx_0_yyy_1[i] * wp_x[i];

        g_0_xxxxxx_0_yyz_0[i] = 5.0 * g_0_xxxx_0_yyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yyz_1[i] * fti_ab_0 + g_0_xxxxx_0_yyz_0[i] * pb_x + g_0_xxxxx_0_yyz_1[i] * wp_x[i];

        g_0_xxxxxx_0_yzz_0[i] = 5.0 * g_0_xxxx_0_yzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yzz_1[i] * fti_ab_0 + g_0_xxxxx_0_yzz_0[i] * pb_x + g_0_xxxxx_0_yzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_zzz_0[i] = 5.0 * g_0_xxxx_0_zzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_zzz_1[i] * fti_ab_0 + g_0_xxxxx_0_zzz_0[i] * pb_x + g_0_xxxxx_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 10-20 components of targeted buffer : SISF

    auto g_0_xxxxxy_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 10);

    auto g_0_xxxxxy_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 11);

    auto g_0_xxxxxy_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 12);

    auto g_0_xxxxxy_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 13);

    auto g_0_xxxxxy_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 14);

    auto g_0_xxxxxy_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 15);

    auto g_0_xxxxxy_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 16);

    auto g_0_xxxxxy_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 17);

    auto g_0_xxxxxy_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 18);

    auto g_0_xxxxxy_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 19);

    #pragma omp simd aligned(g_0_xxxxx_0_xx_1, g_0_xxxxx_0_xxx_0, g_0_xxxxx_0_xxx_1, g_0_xxxxx_0_xxy_0, g_0_xxxxx_0_xxy_1, g_0_xxxxx_0_xxz_0, g_0_xxxxx_0_xxz_1, g_0_xxxxx_0_xy_1, g_0_xxxxx_0_xyy_0, g_0_xxxxx_0_xyy_1, g_0_xxxxx_0_xyz_0, g_0_xxxxx_0_xyz_1, g_0_xxxxx_0_xz_1, g_0_xxxxx_0_xzz_0, g_0_xxxxx_0_xzz_1, g_0_xxxxx_0_yy_1, g_0_xxxxx_0_yyy_0, g_0_xxxxx_0_yyy_1, g_0_xxxxx_0_yyz_0, g_0_xxxxx_0_yyz_1, g_0_xxxxx_0_yz_1, g_0_xxxxx_0_yzz_0, g_0_xxxxx_0_yzz_1, g_0_xxxxx_0_zz_1, g_0_xxxxx_0_zzz_0, g_0_xxxxx_0_zzz_1, g_0_xxxxxy_0_xxx_0, g_0_xxxxxy_0_xxy_0, g_0_xxxxxy_0_xxz_0, g_0_xxxxxy_0_xyy_0, g_0_xxxxxy_0_xyz_0, g_0_xxxxxy_0_xzz_0, g_0_xxxxxy_0_yyy_0, g_0_xxxxxy_0_yyz_0, g_0_xxxxxy_0_yzz_0, g_0_xxxxxy_0_zzz_0, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxy_0_xxx_0[i] = g_0_xxxxx_0_xxx_0[i] * pb_y + g_0_xxxxx_0_xxx_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxy_0[i] = g_0_xxxxx_0_xx_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxy_0[i] * pb_y + g_0_xxxxx_0_xxy_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxz_0[i] = g_0_xxxxx_0_xxz_0[i] * pb_y + g_0_xxxxx_0_xxz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xyy_0[i] = 2.0 * g_0_xxxxx_0_xy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyy_0[i] * pb_y + g_0_xxxxx_0_xyy_1[i] * wp_y[i];

        g_0_xxxxxy_0_xyz_0[i] = g_0_xxxxx_0_xz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyz_0[i] * pb_y + g_0_xxxxx_0_xyz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xzz_0[i] = g_0_xxxxx_0_xzz_0[i] * pb_y + g_0_xxxxx_0_xzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_yyy_0[i] = 3.0 * g_0_xxxxx_0_yy_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyy_0[i] * pb_y + g_0_xxxxx_0_yyy_1[i] * wp_y[i];

        g_0_xxxxxy_0_yyz_0[i] = 2.0 * g_0_xxxxx_0_yz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyz_0[i] * pb_y + g_0_xxxxx_0_yyz_1[i] * wp_y[i];

        g_0_xxxxxy_0_yzz_0[i] = g_0_xxxxx_0_zz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yzz_0[i] * pb_y + g_0_xxxxx_0_yzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_zzz_0[i] = g_0_xxxxx_0_zzz_0[i] * pb_y + g_0_xxxxx_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 20-30 components of targeted buffer : SISF

    auto g_0_xxxxxz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 20);

    auto g_0_xxxxxz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 21);

    auto g_0_xxxxxz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 22);

    auto g_0_xxxxxz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 23);

    auto g_0_xxxxxz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 24);

    auto g_0_xxxxxz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 25);

    auto g_0_xxxxxz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 26);

    auto g_0_xxxxxz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 27);

    auto g_0_xxxxxz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 28);

    auto g_0_xxxxxz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 29);

    #pragma omp simd aligned(g_0_xxxxx_0_xx_1, g_0_xxxxx_0_xxx_0, g_0_xxxxx_0_xxx_1, g_0_xxxxx_0_xxy_0, g_0_xxxxx_0_xxy_1, g_0_xxxxx_0_xxz_0, g_0_xxxxx_0_xxz_1, g_0_xxxxx_0_xy_1, g_0_xxxxx_0_xyy_0, g_0_xxxxx_0_xyy_1, g_0_xxxxx_0_xyz_0, g_0_xxxxx_0_xyz_1, g_0_xxxxx_0_xz_1, g_0_xxxxx_0_xzz_0, g_0_xxxxx_0_xzz_1, g_0_xxxxx_0_yy_1, g_0_xxxxx_0_yyy_0, g_0_xxxxx_0_yyy_1, g_0_xxxxx_0_yyz_0, g_0_xxxxx_0_yyz_1, g_0_xxxxx_0_yz_1, g_0_xxxxx_0_yzz_0, g_0_xxxxx_0_yzz_1, g_0_xxxxx_0_zz_1, g_0_xxxxx_0_zzz_0, g_0_xxxxx_0_zzz_1, g_0_xxxxxz_0_xxx_0, g_0_xxxxxz_0_xxy_0, g_0_xxxxxz_0_xxz_0, g_0_xxxxxz_0_xyy_0, g_0_xxxxxz_0_xyz_0, g_0_xxxxxz_0_xzz_0, g_0_xxxxxz_0_yyy_0, g_0_xxxxxz_0_yyz_0, g_0_xxxxxz_0_yzz_0, g_0_xxxxxz_0_zzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxz_0_xxx_0[i] = g_0_xxxxx_0_xxx_0[i] * pb_z + g_0_xxxxx_0_xxx_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxy_0[i] = g_0_xxxxx_0_xxy_0[i] * pb_z + g_0_xxxxx_0_xxy_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxz_0[i] = g_0_xxxxx_0_xx_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxz_0[i] * pb_z + g_0_xxxxx_0_xxz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xyy_0[i] = g_0_xxxxx_0_xyy_0[i] * pb_z + g_0_xxxxx_0_xyy_1[i] * wp_z[i];

        g_0_xxxxxz_0_xyz_0[i] = g_0_xxxxx_0_xy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyz_0[i] * pb_z + g_0_xxxxx_0_xyz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xzz_0[i] = 2.0 * g_0_xxxxx_0_xz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xzz_0[i] * pb_z + g_0_xxxxx_0_xzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_yyy_0[i] = g_0_xxxxx_0_yyy_0[i] * pb_z + g_0_xxxxx_0_yyy_1[i] * wp_z[i];

        g_0_xxxxxz_0_yyz_0[i] = g_0_xxxxx_0_yy_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyz_0[i] * pb_z + g_0_xxxxx_0_yyz_1[i] * wp_z[i];

        g_0_xxxxxz_0_yzz_0[i] = 2.0 * g_0_xxxxx_0_yz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yzz_0[i] * pb_z + g_0_xxxxx_0_yzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_zzz_0[i] = 3.0 * g_0_xxxxx_0_zz_1[i] * fi_abcd_0 + g_0_xxxxx_0_zzz_0[i] * pb_z + g_0_xxxxx_0_zzz_1[i] * wp_z[i];
    }

    /// Set up 30-40 components of targeted buffer : SISF

    auto g_0_xxxxyy_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 30);

    auto g_0_xxxxyy_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 31);

    auto g_0_xxxxyy_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 32);

    auto g_0_xxxxyy_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 33);

    auto g_0_xxxxyy_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 34);

    auto g_0_xxxxyy_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 35);

    auto g_0_xxxxyy_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 36);

    auto g_0_xxxxyy_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 37);

    auto g_0_xxxxyy_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 38);

    auto g_0_xxxxyy_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 39);

    #pragma omp simd aligned(g_0_xxxx_0_xxx_0, g_0_xxxx_0_xxx_1, g_0_xxxx_0_xxz_0, g_0_xxxx_0_xxz_1, g_0_xxxx_0_xzz_0, g_0_xxxx_0_xzz_1, g_0_xxxxy_0_xxx_0, g_0_xxxxy_0_xxx_1, g_0_xxxxy_0_xxz_0, g_0_xxxxy_0_xxz_1, g_0_xxxxy_0_xzz_0, g_0_xxxxy_0_xzz_1, g_0_xxxxyy_0_xxx_0, g_0_xxxxyy_0_xxy_0, g_0_xxxxyy_0_xxz_0, g_0_xxxxyy_0_xyy_0, g_0_xxxxyy_0_xyz_0, g_0_xxxxyy_0_xzz_0, g_0_xxxxyy_0_yyy_0, g_0_xxxxyy_0_yyz_0, g_0_xxxxyy_0_yzz_0, g_0_xxxxyy_0_zzz_0, g_0_xxxyy_0_xxy_0, g_0_xxxyy_0_xxy_1, g_0_xxxyy_0_xy_1, g_0_xxxyy_0_xyy_0, g_0_xxxyy_0_xyy_1, g_0_xxxyy_0_xyz_0, g_0_xxxyy_0_xyz_1, g_0_xxxyy_0_yy_1, g_0_xxxyy_0_yyy_0, g_0_xxxyy_0_yyy_1, g_0_xxxyy_0_yyz_0, g_0_xxxyy_0_yyz_1, g_0_xxxyy_0_yz_1, g_0_xxxyy_0_yzz_0, g_0_xxxyy_0_yzz_1, g_0_xxxyy_0_zzz_0, g_0_xxxyy_0_zzz_1, g_0_xxyy_0_xxy_0, g_0_xxyy_0_xxy_1, g_0_xxyy_0_xyy_0, g_0_xxyy_0_xyy_1, g_0_xxyy_0_xyz_0, g_0_xxyy_0_xyz_1, g_0_xxyy_0_yyy_0, g_0_xxyy_0_yyy_1, g_0_xxyy_0_yyz_0, g_0_xxyy_0_yyz_1, g_0_xxyy_0_yzz_0, g_0_xxyy_0_yzz_1, g_0_xxyy_0_zzz_0, g_0_xxyy_0_zzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxyy_0_xxx_0[i] = g_0_xxxx_0_xxx_0[i] * fi_ab_0 - g_0_xxxx_0_xxx_1[i] * fti_ab_0 + g_0_xxxxy_0_xxx_0[i] * pb_y + g_0_xxxxy_0_xxx_1[i] * wp_y[i];

        g_0_xxxxyy_0_xxy_0[i] = 3.0 * g_0_xxyy_0_xxy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxy_1[i] * fti_ab_0 + 2.0 * g_0_xxxyy_0_xy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxy_0[i] * pb_x + g_0_xxxyy_0_xxy_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxz_0[i] = g_0_xxxx_0_xxz_0[i] * fi_ab_0 - g_0_xxxx_0_xxz_1[i] * fti_ab_0 + g_0_xxxxy_0_xxz_0[i] * pb_y + g_0_xxxxy_0_xxz_1[i] * wp_y[i];

        g_0_xxxxyy_0_xyy_0[i] = 3.0 * g_0_xxyy_0_xyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xyy_1[i] * fti_ab_0 + g_0_xxxyy_0_yy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyy_0[i] * pb_x + g_0_xxxyy_0_xyy_1[i] * wp_x[i];

        g_0_xxxxyy_0_xyz_0[i] = 3.0 * g_0_xxyy_0_xyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xyz_1[i] * fti_ab_0 + g_0_xxxyy_0_yz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyz_0[i] * pb_x + g_0_xxxyy_0_xyz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xzz_0[i] = g_0_xxxx_0_xzz_0[i] * fi_ab_0 - g_0_xxxx_0_xzz_1[i] * fti_ab_0 + g_0_xxxxy_0_xzz_0[i] * pb_y + g_0_xxxxy_0_xzz_1[i] * wp_y[i];

        g_0_xxxxyy_0_yyy_0[i] = 3.0 * g_0_xxyy_0_yyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yyy_1[i] * fti_ab_0 + g_0_xxxyy_0_yyy_0[i] * pb_x + g_0_xxxyy_0_yyy_1[i] * wp_x[i];

        g_0_xxxxyy_0_yyz_0[i] = 3.0 * g_0_xxyy_0_yyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yyz_1[i] * fti_ab_0 + g_0_xxxyy_0_yyz_0[i] * pb_x + g_0_xxxyy_0_yyz_1[i] * wp_x[i];

        g_0_xxxxyy_0_yzz_0[i] = 3.0 * g_0_xxyy_0_yzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yzz_1[i] * fti_ab_0 + g_0_xxxyy_0_yzz_0[i] * pb_x + g_0_xxxyy_0_yzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_zzz_0[i] = 3.0 * g_0_xxyy_0_zzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_zzz_1[i] * fti_ab_0 + g_0_xxxyy_0_zzz_0[i] * pb_x + g_0_xxxyy_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 40-50 components of targeted buffer : SISF

    auto g_0_xxxxyz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 40);

    auto g_0_xxxxyz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 41);

    auto g_0_xxxxyz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 42);

    auto g_0_xxxxyz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 43);

    auto g_0_xxxxyz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 44);

    auto g_0_xxxxyz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 45);

    auto g_0_xxxxyz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 46);

    auto g_0_xxxxyz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 47);

    auto g_0_xxxxyz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 48);

    auto g_0_xxxxyz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 49);

    #pragma omp simd aligned(g_0_xxxxy_0_xxy_0, g_0_xxxxy_0_xxy_1, g_0_xxxxy_0_xyy_0, g_0_xxxxy_0_xyy_1, g_0_xxxxy_0_yyy_0, g_0_xxxxy_0_yyy_1, g_0_xxxxyz_0_xxx_0, g_0_xxxxyz_0_xxy_0, g_0_xxxxyz_0_xxz_0, g_0_xxxxyz_0_xyy_0, g_0_xxxxyz_0_xyz_0, g_0_xxxxyz_0_xzz_0, g_0_xxxxyz_0_yyy_0, g_0_xxxxyz_0_yyz_0, g_0_xxxxyz_0_yzz_0, g_0_xxxxyz_0_zzz_0, g_0_xxxxz_0_xxx_0, g_0_xxxxz_0_xxx_1, g_0_xxxxz_0_xxz_0, g_0_xxxxz_0_xxz_1, g_0_xxxxz_0_xyz_0, g_0_xxxxz_0_xyz_1, g_0_xxxxz_0_xz_1, g_0_xxxxz_0_xzz_0, g_0_xxxxz_0_xzz_1, g_0_xxxxz_0_yyz_0, g_0_xxxxz_0_yyz_1, g_0_xxxxz_0_yz_1, g_0_xxxxz_0_yzz_0, g_0_xxxxz_0_yzz_1, g_0_xxxxz_0_zz_1, g_0_xxxxz_0_zzz_0, g_0_xxxxz_0_zzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyz_0_xxx_0[i] = g_0_xxxxz_0_xxx_0[i] * pb_y + g_0_xxxxz_0_xxx_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxy_0[i] = g_0_xxxxy_0_xxy_0[i] * pb_z + g_0_xxxxy_0_xxy_1[i] * wp_z[i];

        g_0_xxxxyz_0_xxz_0[i] = g_0_xxxxz_0_xxz_0[i] * pb_y + g_0_xxxxz_0_xxz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xyy_0[i] = g_0_xxxxy_0_xyy_0[i] * pb_z + g_0_xxxxy_0_xyy_1[i] * wp_z[i];

        g_0_xxxxyz_0_xyz_0[i] = g_0_xxxxz_0_xz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xyz_0[i] * pb_y + g_0_xxxxz_0_xyz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xzz_0[i] = g_0_xxxxz_0_xzz_0[i] * pb_y + g_0_xxxxz_0_xzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_yyy_0[i] = g_0_xxxxy_0_yyy_0[i] * pb_z + g_0_xxxxy_0_yyy_1[i] * wp_z[i];

        g_0_xxxxyz_0_yyz_0[i] = 2.0 * g_0_xxxxz_0_yz_1[i] * fi_abcd_0 + g_0_xxxxz_0_yyz_0[i] * pb_y + g_0_xxxxz_0_yyz_1[i] * wp_y[i];

        g_0_xxxxyz_0_yzz_0[i] = g_0_xxxxz_0_zz_1[i] * fi_abcd_0 + g_0_xxxxz_0_yzz_0[i] * pb_y + g_0_xxxxz_0_yzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_zzz_0[i] = g_0_xxxxz_0_zzz_0[i] * pb_y + g_0_xxxxz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 50-60 components of targeted buffer : SISF

    auto g_0_xxxxzz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 50);

    auto g_0_xxxxzz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 51);

    auto g_0_xxxxzz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 52);

    auto g_0_xxxxzz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 53);

    auto g_0_xxxxzz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 54);

    auto g_0_xxxxzz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 55);

    auto g_0_xxxxzz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 56);

    auto g_0_xxxxzz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 57);

    auto g_0_xxxxzz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 58);

    auto g_0_xxxxzz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 59);

    #pragma omp simd aligned(g_0_xxxx_0_xxx_0, g_0_xxxx_0_xxx_1, g_0_xxxx_0_xxy_0, g_0_xxxx_0_xxy_1, g_0_xxxx_0_xyy_0, g_0_xxxx_0_xyy_1, g_0_xxxxz_0_xxx_0, g_0_xxxxz_0_xxx_1, g_0_xxxxz_0_xxy_0, g_0_xxxxz_0_xxy_1, g_0_xxxxz_0_xyy_0, g_0_xxxxz_0_xyy_1, g_0_xxxxzz_0_xxx_0, g_0_xxxxzz_0_xxy_0, g_0_xxxxzz_0_xxz_0, g_0_xxxxzz_0_xyy_0, g_0_xxxxzz_0_xyz_0, g_0_xxxxzz_0_xzz_0, g_0_xxxxzz_0_yyy_0, g_0_xxxxzz_0_yyz_0, g_0_xxxxzz_0_yzz_0, g_0_xxxxzz_0_zzz_0, g_0_xxxzz_0_xxz_0, g_0_xxxzz_0_xxz_1, g_0_xxxzz_0_xyz_0, g_0_xxxzz_0_xyz_1, g_0_xxxzz_0_xz_1, g_0_xxxzz_0_xzz_0, g_0_xxxzz_0_xzz_1, g_0_xxxzz_0_yyy_0, g_0_xxxzz_0_yyy_1, g_0_xxxzz_0_yyz_0, g_0_xxxzz_0_yyz_1, g_0_xxxzz_0_yz_1, g_0_xxxzz_0_yzz_0, g_0_xxxzz_0_yzz_1, g_0_xxxzz_0_zz_1, g_0_xxxzz_0_zzz_0, g_0_xxxzz_0_zzz_1, g_0_xxzz_0_xxz_0, g_0_xxzz_0_xxz_1, g_0_xxzz_0_xyz_0, g_0_xxzz_0_xyz_1, g_0_xxzz_0_xzz_0, g_0_xxzz_0_xzz_1, g_0_xxzz_0_yyy_0, g_0_xxzz_0_yyy_1, g_0_xxzz_0_yyz_0, g_0_xxzz_0_yyz_1, g_0_xxzz_0_yzz_0, g_0_xxzz_0_yzz_1, g_0_xxzz_0_zzz_0, g_0_xxzz_0_zzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxzz_0_xxx_0[i] = g_0_xxxx_0_xxx_0[i] * fi_ab_0 - g_0_xxxx_0_xxx_1[i] * fti_ab_0 + g_0_xxxxz_0_xxx_0[i] * pb_z + g_0_xxxxz_0_xxx_1[i] * wp_z[i];

        g_0_xxxxzz_0_xxy_0[i] = g_0_xxxx_0_xxy_0[i] * fi_ab_0 - g_0_xxxx_0_xxy_1[i] * fti_ab_0 + g_0_xxxxz_0_xxy_0[i] * pb_z + g_0_xxxxz_0_xxy_1[i] * wp_z[i];

        g_0_xxxxzz_0_xxz_0[i] = 3.0 * g_0_xxzz_0_xxz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxz_1[i] * fti_ab_0 + 2.0 * g_0_xxxzz_0_xz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxz_0[i] * pb_x + g_0_xxxzz_0_xxz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xyy_0[i] = g_0_xxxx_0_xyy_0[i] * fi_ab_0 - g_0_xxxx_0_xyy_1[i] * fti_ab_0 + g_0_xxxxz_0_xyy_0[i] * pb_z + g_0_xxxxz_0_xyy_1[i] * wp_z[i];

        g_0_xxxxzz_0_xyz_0[i] = 3.0 * g_0_xxzz_0_xyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xyz_1[i] * fti_ab_0 + g_0_xxxzz_0_yz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyz_0[i] * pb_x + g_0_xxxzz_0_xyz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xzz_0[i] = 3.0 * g_0_xxzz_0_xzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xzz_1[i] * fti_ab_0 + g_0_xxxzz_0_zz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xzz_0[i] * pb_x + g_0_xxxzz_0_xzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_yyy_0[i] = 3.0 * g_0_xxzz_0_yyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yyy_1[i] * fti_ab_0 + g_0_xxxzz_0_yyy_0[i] * pb_x + g_0_xxxzz_0_yyy_1[i] * wp_x[i];

        g_0_xxxxzz_0_yyz_0[i] = 3.0 * g_0_xxzz_0_yyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yyz_1[i] * fti_ab_0 + g_0_xxxzz_0_yyz_0[i] * pb_x + g_0_xxxzz_0_yyz_1[i] * wp_x[i];

        g_0_xxxxzz_0_yzz_0[i] = 3.0 * g_0_xxzz_0_yzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yzz_1[i] * fti_ab_0 + g_0_xxxzz_0_yzz_0[i] * pb_x + g_0_xxxzz_0_yzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_zzz_0[i] = 3.0 * g_0_xxzz_0_zzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_zzz_1[i] * fti_ab_0 + g_0_xxxzz_0_zzz_0[i] * pb_x + g_0_xxxzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 60-70 components of targeted buffer : SISF

    auto g_0_xxxyyy_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 60);

    auto g_0_xxxyyy_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 61);

    auto g_0_xxxyyy_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 62);

    auto g_0_xxxyyy_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 63);

    auto g_0_xxxyyy_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 64);

    auto g_0_xxxyyy_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 65);

    auto g_0_xxxyyy_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 66);

    auto g_0_xxxyyy_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 67);

    auto g_0_xxxyyy_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 68);

    auto g_0_xxxyyy_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 69);

    #pragma omp simd aligned(g_0_xxxy_0_xxx_0, g_0_xxxy_0_xxx_1, g_0_xxxy_0_xxz_0, g_0_xxxy_0_xxz_1, g_0_xxxy_0_xzz_0, g_0_xxxy_0_xzz_1, g_0_xxxyy_0_xxx_0, g_0_xxxyy_0_xxx_1, g_0_xxxyy_0_xxz_0, g_0_xxxyy_0_xxz_1, g_0_xxxyy_0_xzz_0, g_0_xxxyy_0_xzz_1, g_0_xxxyyy_0_xxx_0, g_0_xxxyyy_0_xxy_0, g_0_xxxyyy_0_xxz_0, g_0_xxxyyy_0_xyy_0, g_0_xxxyyy_0_xyz_0, g_0_xxxyyy_0_xzz_0, g_0_xxxyyy_0_yyy_0, g_0_xxxyyy_0_yyz_0, g_0_xxxyyy_0_yzz_0, g_0_xxxyyy_0_zzz_0, g_0_xxyyy_0_xxy_0, g_0_xxyyy_0_xxy_1, g_0_xxyyy_0_xy_1, g_0_xxyyy_0_xyy_0, g_0_xxyyy_0_xyy_1, g_0_xxyyy_0_xyz_0, g_0_xxyyy_0_xyz_1, g_0_xxyyy_0_yy_1, g_0_xxyyy_0_yyy_0, g_0_xxyyy_0_yyy_1, g_0_xxyyy_0_yyz_0, g_0_xxyyy_0_yyz_1, g_0_xxyyy_0_yz_1, g_0_xxyyy_0_yzz_0, g_0_xxyyy_0_yzz_1, g_0_xxyyy_0_zzz_0, g_0_xxyyy_0_zzz_1, g_0_xyyy_0_xxy_0, g_0_xyyy_0_xxy_1, g_0_xyyy_0_xyy_0, g_0_xyyy_0_xyy_1, g_0_xyyy_0_xyz_0, g_0_xyyy_0_xyz_1, g_0_xyyy_0_yyy_0, g_0_xyyy_0_yyy_1, g_0_xyyy_0_yyz_0, g_0_xyyy_0_yyz_1, g_0_xyyy_0_yzz_0, g_0_xyyy_0_yzz_1, g_0_xyyy_0_zzz_0, g_0_xyyy_0_zzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyy_0_xxx_0[i] = 2.0 * g_0_xxxy_0_xxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xxx_1[i] * fti_ab_0 + g_0_xxxyy_0_xxx_0[i] * pb_y + g_0_xxxyy_0_xxx_1[i] * wp_y[i];

        g_0_xxxyyy_0_xxy_0[i] = 2.0 * g_0_xyyy_0_xxy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxy_1[i] * fti_ab_0 + 2.0 * g_0_xxyyy_0_xy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxy_0[i] * pb_x + g_0_xxyyy_0_xxy_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxz_0[i] = 2.0 * g_0_xxxy_0_xxz_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xxz_1[i] * fti_ab_0 + g_0_xxxyy_0_xxz_0[i] * pb_y + g_0_xxxyy_0_xxz_1[i] * wp_y[i];

        g_0_xxxyyy_0_xyy_0[i] = 2.0 * g_0_xyyy_0_xyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xyy_1[i] * fti_ab_0 + g_0_xxyyy_0_yy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyy_0[i] * pb_x + g_0_xxyyy_0_xyy_1[i] * wp_x[i];

        g_0_xxxyyy_0_xyz_0[i] = 2.0 * g_0_xyyy_0_xyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xyz_1[i] * fti_ab_0 + g_0_xxyyy_0_yz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyz_0[i] * pb_x + g_0_xxyyy_0_xyz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xzz_0[i] = 2.0 * g_0_xxxy_0_xzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xzz_1[i] * fti_ab_0 + g_0_xxxyy_0_xzz_0[i] * pb_y + g_0_xxxyy_0_xzz_1[i] * wp_y[i];

        g_0_xxxyyy_0_yyy_0[i] = 2.0 * g_0_xyyy_0_yyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yyy_1[i] * fti_ab_0 + g_0_xxyyy_0_yyy_0[i] * pb_x + g_0_xxyyy_0_yyy_1[i] * wp_x[i];

        g_0_xxxyyy_0_yyz_0[i] = 2.0 * g_0_xyyy_0_yyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yyz_1[i] * fti_ab_0 + g_0_xxyyy_0_yyz_0[i] * pb_x + g_0_xxyyy_0_yyz_1[i] * wp_x[i];

        g_0_xxxyyy_0_yzz_0[i] = 2.0 * g_0_xyyy_0_yzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yzz_1[i] * fti_ab_0 + g_0_xxyyy_0_yzz_0[i] * pb_x + g_0_xxyyy_0_yzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_zzz_0[i] = 2.0 * g_0_xyyy_0_zzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_zzz_1[i] * fti_ab_0 + g_0_xxyyy_0_zzz_0[i] * pb_x + g_0_xxyyy_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 70-80 components of targeted buffer : SISF

    auto g_0_xxxyyz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 70);

    auto g_0_xxxyyz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 71);

    auto g_0_xxxyyz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 72);

    auto g_0_xxxyyz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 73);

    auto g_0_xxxyyz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 74);

    auto g_0_xxxyyz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 75);

    auto g_0_xxxyyz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 76);

    auto g_0_xxxyyz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 77);

    auto g_0_xxxyyz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 78);

    auto g_0_xxxyyz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 79);

    #pragma omp simd aligned(g_0_xxxyy_0_xx_1, g_0_xxxyy_0_xxx_0, g_0_xxxyy_0_xxx_1, g_0_xxxyy_0_xxy_0, g_0_xxxyy_0_xxy_1, g_0_xxxyy_0_xxz_0, g_0_xxxyy_0_xxz_1, g_0_xxxyy_0_xy_1, g_0_xxxyy_0_xyy_0, g_0_xxxyy_0_xyy_1, g_0_xxxyy_0_xyz_0, g_0_xxxyy_0_xyz_1, g_0_xxxyy_0_xz_1, g_0_xxxyy_0_xzz_0, g_0_xxxyy_0_xzz_1, g_0_xxxyy_0_yy_1, g_0_xxxyy_0_yyy_0, g_0_xxxyy_0_yyy_1, g_0_xxxyy_0_yyz_0, g_0_xxxyy_0_yyz_1, g_0_xxxyy_0_yz_1, g_0_xxxyy_0_yzz_0, g_0_xxxyy_0_yzz_1, g_0_xxxyy_0_zz_1, g_0_xxxyy_0_zzz_0, g_0_xxxyy_0_zzz_1, g_0_xxxyyz_0_xxx_0, g_0_xxxyyz_0_xxy_0, g_0_xxxyyz_0_xxz_0, g_0_xxxyyz_0_xyy_0, g_0_xxxyyz_0_xyz_0, g_0_xxxyyz_0_xzz_0, g_0_xxxyyz_0_yyy_0, g_0_xxxyyz_0_yyz_0, g_0_xxxyyz_0_yzz_0, g_0_xxxyyz_0_zzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyz_0_xxx_0[i] = g_0_xxxyy_0_xxx_0[i] * pb_z + g_0_xxxyy_0_xxx_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxy_0[i] = g_0_xxxyy_0_xxy_0[i] * pb_z + g_0_xxxyy_0_xxy_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxz_0[i] = g_0_xxxyy_0_xx_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxz_0[i] * pb_z + g_0_xxxyy_0_xxz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xyy_0[i] = g_0_xxxyy_0_xyy_0[i] * pb_z + g_0_xxxyy_0_xyy_1[i] * wp_z[i];

        g_0_xxxyyz_0_xyz_0[i] = g_0_xxxyy_0_xy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyz_0[i] * pb_z + g_0_xxxyy_0_xyz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xzz_0[i] = 2.0 * g_0_xxxyy_0_xz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xzz_0[i] * pb_z + g_0_xxxyy_0_xzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_yyy_0[i] = g_0_xxxyy_0_yyy_0[i] * pb_z + g_0_xxxyy_0_yyy_1[i] * wp_z[i];

        g_0_xxxyyz_0_yyz_0[i] = g_0_xxxyy_0_yy_1[i] * fi_abcd_0 + g_0_xxxyy_0_yyz_0[i] * pb_z + g_0_xxxyy_0_yyz_1[i] * wp_z[i];

        g_0_xxxyyz_0_yzz_0[i] = 2.0 * g_0_xxxyy_0_yz_1[i] * fi_abcd_0 + g_0_xxxyy_0_yzz_0[i] * pb_z + g_0_xxxyy_0_yzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_zzz_0[i] = 3.0 * g_0_xxxyy_0_zz_1[i] * fi_abcd_0 + g_0_xxxyy_0_zzz_0[i] * pb_z + g_0_xxxyy_0_zzz_1[i] * wp_z[i];
    }

    /// Set up 80-90 components of targeted buffer : SISF

    auto g_0_xxxyzz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 80);

    auto g_0_xxxyzz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 81);

    auto g_0_xxxyzz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 82);

    auto g_0_xxxyzz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 83);

    auto g_0_xxxyzz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 84);

    auto g_0_xxxyzz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 85);

    auto g_0_xxxyzz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 86);

    auto g_0_xxxyzz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 87);

    auto g_0_xxxyzz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 88);

    auto g_0_xxxyzz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 89);

    #pragma omp simd aligned(g_0_xxxyzz_0_xxx_0, g_0_xxxyzz_0_xxy_0, g_0_xxxyzz_0_xxz_0, g_0_xxxyzz_0_xyy_0, g_0_xxxyzz_0_xyz_0, g_0_xxxyzz_0_xzz_0, g_0_xxxyzz_0_yyy_0, g_0_xxxyzz_0_yyz_0, g_0_xxxyzz_0_yzz_0, g_0_xxxyzz_0_zzz_0, g_0_xxxzz_0_xx_1, g_0_xxxzz_0_xxx_0, g_0_xxxzz_0_xxx_1, g_0_xxxzz_0_xxy_0, g_0_xxxzz_0_xxy_1, g_0_xxxzz_0_xxz_0, g_0_xxxzz_0_xxz_1, g_0_xxxzz_0_xy_1, g_0_xxxzz_0_xyy_0, g_0_xxxzz_0_xyy_1, g_0_xxxzz_0_xyz_0, g_0_xxxzz_0_xyz_1, g_0_xxxzz_0_xz_1, g_0_xxxzz_0_xzz_0, g_0_xxxzz_0_xzz_1, g_0_xxxzz_0_yy_1, g_0_xxxzz_0_yyy_0, g_0_xxxzz_0_yyy_1, g_0_xxxzz_0_yyz_0, g_0_xxxzz_0_yyz_1, g_0_xxxzz_0_yz_1, g_0_xxxzz_0_yzz_0, g_0_xxxzz_0_yzz_1, g_0_xxxzz_0_zz_1, g_0_xxxzz_0_zzz_0, g_0_xxxzz_0_zzz_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyzz_0_xxx_0[i] = g_0_xxxzz_0_xxx_0[i] * pb_y + g_0_xxxzz_0_xxx_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxy_0[i] = g_0_xxxzz_0_xx_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxy_0[i] * pb_y + g_0_xxxzz_0_xxy_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxz_0[i] = g_0_xxxzz_0_xxz_0[i] * pb_y + g_0_xxxzz_0_xxz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xyy_0[i] = 2.0 * g_0_xxxzz_0_xy_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyy_0[i] * pb_y + g_0_xxxzz_0_xyy_1[i] * wp_y[i];

        g_0_xxxyzz_0_xyz_0[i] = g_0_xxxzz_0_xz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyz_0[i] * pb_y + g_0_xxxzz_0_xyz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xzz_0[i] = g_0_xxxzz_0_xzz_0[i] * pb_y + g_0_xxxzz_0_xzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_yyy_0[i] = 3.0 * g_0_xxxzz_0_yy_1[i] * fi_abcd_0 + g_0_xxxzz_0_yyy_0[i] * pb_y + g_0_xxxzz_0_yyy_1[i] * wp_y[i];

        g_0_xxxyzz_0_yyz_0[i] = 2.0 * g_0_xxxzz_0_yz_1[i] * fi_abcd_0 + g_0_xxxzz_0_yyz_0[i] * pb_y + g_0_xxxzz_0_yyz_1[i] * wp_y[i];

        g_0_xxxyzz_0_yzz_0[i] = g_0_xxxzz_0_zz_1[i] * fi_abcd_0 + g_0_xxxzz_0_yzz_0[i] * pb_y + g_0_xxxzz_0_yzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_zzz_0[i] = g_0_xxxzz_0_zzz_0[i] * pb_y + g_0_xxxzz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 90-100 components of targeted buffer : SISF

    auto g_0_xxxzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 90);

    auto g_0_xxxzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 91);

    auto g_0_xxxzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 92);

    auto g_0_xxxzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 93);

    auto g_0_xxxzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 94);

    auto g_0_xxxzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 95);

    auto g_0_xxxzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 96);

    auto g_0_xxxzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 97);

    auto g_0_xxxzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 98);

    auto g_0_xxxzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 99);

    #pragma omp simd aligned(g_0_xxxz_0_xxx_0, g_0_xxxz_0_xxx_1, g_0_xxxz_0_xxy_0, g_0_xxxz_0_xxy_1, g_0_xxxz_0_xyy_0, g_0_xxxz_0_xyy_1, g_0_xxxzz_0_xxx_0, g_0_xxxzz_0_xxx_1, g_0_xxxzz_0_xxy_0, g_0_xxxzz_0_xxy_1, g_0_xxxzz_0_xyy_0, g_0_xxxzz_0_xyy_1, g_0_xxxzzz_0_xxx_0, g_0_xxxzzz_0_xxy_0, g_0_xxxzzz_0_xxz_0, g_0_xxxzzz_0_xyy_0, g_0_xxxzzz_0_xyz_0, g_0_xxxzzz_0_xzz_0, g_0_xxxzzz_0_yyy_0, g_0_xxxzzz_0_yyz_0, g_0_xxxzzz_0_yzz_0, g_0_xxxzzz_0_zzz_0, g_0_xxzzz_0_xxz_0, g_0_xxzzz_0_xxz_1, g_0_xxzzz_0_xyz_0, g_0_xxzzz_0_xyz_1, g_0_xxzzz_0_xz_1, g_0_xxzzz_0_xzz_0, g_0_xxzzz_0_xzz_1, g_0_xxzzz_0_yyy_0, g_0_xxzzz_0_yyy_1, g_0_xxzzz_0_yyz_0, g_0_xxzzz_0_yyz_1, g_0_xxzzz_0_yz_1, g_0_xxzzz_0_yzz_0, g_0_xxzzz_0_yzz_1, g_0_xxzzz_0_zz_1, g_0_xxzzz_0_zzz_0, g_0_xxzzz_0_zzz_1, g_0_xzzz_0_xxz_0, g_0_xzzz_0_xxz_1, g_0_xzzz_0_xyz_0, g_0_xzzz_0_xyz_1, g_0_xzzz_0_xzz_0, g_0_xzzz_0_xzz_1, g_0_xzzz_0_yyy_0, g_0_xzzz_0_yyy_1, g_0_xzzz_0_yyz_0, g_0_xzzz_0_yyz_1, g_0_xzzz_0_yzz_0, g_0_xzzz_0_yzz_1, g_0_xzzz_0_zzz_0, g_0_xzzz_0_zzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxzzz_0_xxx_0[i] = 2.0 * g_0_xxxz_0_xxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xxx_1[i] * fti_ab_0 + g_0_xxxzz_0_xxx_0[i] * pb_z + g_0_xxxzz_0_xxx_1[i] * wp_z[i];

        g_0_xxxzzz_0_xxy_0[i] = 2.0 * g_0_xxxz_0_xxy_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xxy_1[i] * fti_ab_0 + g_0_xxxzz_0_xxy_0[i] * pb_z + g_0_xxxzz_0_xxy_1[i] * wp_z[i];

        g_0_xxxzzz_0_xxz_0[i] = 2.0 * g_0_xzzz_0_xxz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxz_1[i] * fti_ab_0 + 2.0 * g_0_xxzzz_0_xz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxz_0[i] * pb_x + g_0_xxzzz_0_xxz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xyy_0[i] = 2.0 * g_0_xxxz_0_xyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xyy_1[i] * fti_ab_0 + g_0_xxxzz_0_xyy_0[i] * pb_z + g_0_xxxzz_0_xyy_1[i] * wp_z[i];

        g_0_xxxzzz_0_xyz_0[i] = 2.0 * g_0_xzzz_0_xyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xyz_1[i] * fti_ab_0 + g_0_xxzzz_0_yz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyz_0[i] * pb_x + g_0_xxzzz_0_xyz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xzz_0[i] = 2.0 * g_0_xzzz_0_xzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xzz_1[i] * fti_ab_0 + g_0_xxzzz_0_zz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xzz_0[i] * pb_x + g_0_xxzzz_0_xzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_yyy_0[i] = 2.0 * g_0_xzzz_0_yyy_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yyy_1[i] * fti_ab_0 + g_0_xxzzz_0_yyy_0[i] * pb_x + g_0_xxzzz_0_yyy_1[i] * wp_x[i];

        g_0_xxxzzz_0_yyz_0[i] = 2.0 * g_0_xzzz_0_yyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yyz_1[i] * fti_ab_0 + g_0_xxzzz_0_yyz_0[i] * pb_x + g_0_xxzzz_0_yyz_1[i] * wp_x[i];

        g_0_xxxzzz_0_yzz_0[i] = 2.0 * g_0_xzzz_0_yzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yzz_1[i] * fti_ab_0 + g_0_xxzzz_0_yzz_0[i] * pb_x + g_0_xxzzz_0_yzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_zzz_0[i] = 2.0 * g_0_xzzz_0_zzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_zzz_1[i] * fti_ab_0 + g_0_xxzzz_0_zzz_0[i] * pb_x + g_0_xxzzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 100-110 components of targeted buffer : SISF

    auto g_0_xxyyyy_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 100);

    auto g_0_xxyyyy_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 101);

    auto g_0_xxyyyy_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 102);

    auto g_0_xxyyyy_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 103);

    auto g_0_xxyyyy_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 104);

    auto g_0_xxyyyy_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 105);

    auto g_0_xxyyyy_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 106);

    auto g_0_xxyyyy_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 107);

    auto g_0_xxyyyy_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 108);

    auto g_0_xxyyyy_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 109);

    #pragma omp simd aligned(g_0_xxyy_0_xxx_0, g_0_xxyy_0_xxx_1, g_0_xxyy_0_xxz_0, g_0_xxyy_0_xxz_1, g_0_xxyy_0_xzz_0, g_0_xxyy_0_xzz_1, g_0_xxyyy_0_xxx_0, g_0_xxyyy_0_xxx_1, g_0_xxyyy_0_xxz_0, g_0_xxyyy_0_xxz_1, g_0_xxyyy_0_xzz_0, g_0_xxyyy_0_xzz_1, g_0_xxyyyy_0_xxx_0, g_0_xxyyyy_0_xxy_0, g_0_xxyyyy_0_xxz_0, g_0_xxyyyy_0_xyy_0, g_0_xxyyyy_0_xyz_0, g_0_xxyyyy_0_xzz_0, g_0_xxyyyy_0_yyy_0, g_0_xxyyyy_0_yyz_0, g_0_xxyyyy_0_yzz_0, g_0_xxyyyy_0_zzz_0, g_0_xyyyy_0_xxy_0, g_0_xyyyy_0_xxy_1, g_0_xyyyy_0_xy_1, g_0_xyyyy_0_xyy_0, g_0_xyyyy_0_xyy_1, g_0_xyyyy_0_xyz_0, g_0_xyyyy_0_xyz_1, g_0_xyyyy_0_yy_1, g_0_xyyyy_0_yyy_0, g_0_xyyyy_0_yyy_1, g_0_xyyyy_0_yyz_0, g_0_xyyyy_0_yyz_1, g_0_xyyyy_0_yz_1, g_0_xyyyy_0_yzz_0, g_0_xyyyy_0_yzz_1, g_0_xyyyy_0_zzz_0, g_0_xyyyy_0_zzz_1, g_0_yyyy_0_xxy_0, g_0_yyyy_0_xxy_1, g_0_yyyy_0_xyy_0, g_0_yyyy_0_xyy_1, g_0_yyyy_0_xyz_0, g_0_yyyy_0_xyz_1, g_0_yyyy_0_yyy_0, g_0_yyyy_0_yyy_1, g_0_yyyy_0_yyz_0, g_0_yyyy_0_yyz_1, g_0_yyyy_0_yzz_0, g_0_yyyy_0_yzz_1, g_0_yyyy_0_zzz_0, g_0_yyyy_0_zzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyy_0_xxx_0[i] = 3.0 * g_0_xxyy_0_xxx_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxx_1[i] * fti_ab_0 + g_0_xxyyy_0_xxx_0[i] * pb_y + g_0_xxyyy_0_xxx_1[i] * wp_y[i];

        g_0_xxyyyy_0_xxy_0[i] = g_0_yyyy_0_xxy_0[i] * fi_ab_0 - g_0_yyyy_0_xxy_1[i] * fti_ab_0 + 2.0 * g_0_xyyyy_0_xy_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxy_0[i] * pb_x + g_0_xyyyy_0_xxy_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxz_0[i] = 3.0 * g_0_xxyy_0_xxz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxz_1[i] * fti_ab_0 + g_0_xxyyy_0_xxz_0[i] * pb_y + g_0_xxyyy_0_xxz_1[i] * wp_y[i];

        g_0_xxyyyy_0_xyy_0[i] = g_0_yyyy_0_xyy_0[i] * fi_ab_0 - g_0_yyyy_0_xyy_1[i] * fti_ab_0 + g_0_xyyyy_0_yy_1[i] * fi_abcd_0 + g_0_xyyyy_0_xyy_0[i] * pb_x + g_0_xyyyy_0_xyy_1[i] * wp_x[i];

        g_0_xxyyyy_0_xyz_0[i] = g_0_yyyy_0_xyz_0[i] * fi_ab_0 - g_0_yyyy_0_xyz_1[i] * fti_ab_0 + g_0_xyyyy_0_yz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xyz_0[i] * pb_x + g_0_xyyyy_0_xyz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xzz_0[i] = 3.0 * g_0_xxyy_0_xzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xzz_1[i] * fti_ab_0 + g_0_xxyyy_0_xzz_0[i] * pb_y + g_0_xxyyy_0_xzz_1[i] * wp_y[i];

        g_0_xxyyyy_0_yyy_0[i] = g_0_yyyy_0_yyy_0[i] * fi_ab_0 - g_0_yyyy_0_yyy_1[i] * fti_ab_0 + g_0_xyyyy_0_yyy_0[i] * pb_x + g_0_xyyyy_0_yyy_1[i] * wp_x[i];

        g_0_xxyyyy_0_yyz_0[i] = g_0_yyyy_0_yyz_0[i] * fi_ab_0 - g_0_yyyy_0_yyz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyz_0[i] * pb_x + g_0_xyyyy_0_yyz_1[i] * wp_x[i];

        g_0_xxyyyy_0_yzz_0[i] = g_0_yyyy_0_yzz_0[i] * fi_ab_0 - g_0_yyyy_0_yzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yzz_0[i] * pb_x + g_0_xyyyy_0_yzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_zzz_0[i] = g_0_yyyy_0_zzz_0[i] * fi_ab_0 - g_0_yyyy_0_zzz_1[i] * fti_ab_0 + g_0_xyyyy_0_zzz_0[i] * pb_x + g_0_xyyyy_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 110-120 components of targeted buffer : SISF

    auto g_0_xxyyyz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 110);

    auto g_0_xxyyyz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 111);

    auto g_0_xxyyyz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 112);

    auto g_0_xxyyyz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 113);

    auto g_0_xxyyyz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 114);

    auto g_0_xxyyyz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 115);

    auto g_0_xxyyyz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 116);

    auto g_0_xxyyyz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 117);

    auto g_0_xxyyyz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 118);

    auto g_0_xxyyyz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 119);

    #pragma omp simd aligned(g_0_xxyyy_0_xx_1, g_0_xxyyy_0_xxx_0, g_0_xxyyy_0_xxx_1, g_0_xxyyy_0_xxy_0, g_0_xxyyy_0_xxy_1, g_0_xxyyy_0_xxz_0, g_0_xxyyy_0_xxz_1, g_0_xxyyy_0_xy_1, g_0_xxyyy_0_xyy_0, g_0_xxyyy_0_xyy_1, g_0_xxyyy_0_xyz_0, g_0_xxyyy_0_xyz_1, g_0_xxyyy_0_xz_1, g_0_xxyyy_0_xzz_0, g_0_xxyyy_0_xzz_1, g_0_xxyyy_0_yy_1, g_0_xxyyy_0_yyy_0, g_0_xxyyy_0_yyy_1, g_0_xxyyy_0_yyz_0, g_0_xxyyy_0_yyz_1, g_0_xxyyy_0_yz_1, g_0_xxyyy_0_yzz_0, g_0_xxyyy_0_yzz_1, g_0_xxyyy_0_zz_1, g_0_xxyyy_0_zzz_0, g_0_xxyyy_0_zzz_1, g_0_xxyyyz_0_xxx_0, g_0_xxyyyz_0_xxy_0, g_0_xxyyyz_0_xxz_0, g_0_xxyyyz_0_xyy_0, g_0_xxyyyz_0_xyz_0, g_0_xxyyyz_0_xzz_0, g_0_xxyyyz_0_yyy_0, g_0_xxyyyz_0_yyz_0, g_0_xxyyyz_0_yzz_0, g_0_xxyyyz_0_zzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyz_0_xxx_0[i] = g_0_xxyyy_0_xxx_0[i] * pb_z + g_0_xxyyy_0_xxx_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxy_0[i] = g_0_xxyyy_0_xxy_0[i] * pb_z + g_0_xxyyy_0_xxy_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxz_0[i] = g_0_xxyyy_0_xx_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxz_0[i] * pb_z + g_0_xxyyy_0_xxz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xyy_0[i] = g_0_xxyyy_0_xyy_0[i] * pb_z + g_0_xxyyy_0_xyy_1[i] * wp_z[i];

        g_0_xxyyyz_0_xyz_0[i] = g_0_xxyyy_0_xy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyz_0[i] * pb_z + g_0_xxyyy_0_xyz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xzz_0[i] = 2.0 * g_0_xxyyy_0_xz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xzz_0[i] * pb_z + g_0_xxyyy_0_xzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_yyy_0[i] = g_0_xxyyy_0_yyy_0[i] * pb_z + g_0_xxyyy_0_yyy_1[i] * wp_z[i];

        g_0_xxyyyz_0_yyz_0[i] = g_0_xxyyy_0_yy_1[i] * fi_abcd_0 + g_0_xxyyy_0_yyz_0[i] * pb_z + g_0_xxyyy_0_yyz_1[i] * wp_z[i];

        g_0_xxyyyz_0_yzz_0[i] = 2.0 * g_0_xxyyy_0_yz_1[i] * fi_abcd_0 + g_0_xxyyy_0_yzz_0[i] * pb_z + g_0_xxyyy_0_yzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_zzz_0[i] = 3.0 * g_0_xxyyy_0_zz_1[i] * fi_abcd_0 + g_0_xxyyy_0_zzz_0[i] * pb_z + g_0_xxyyy_0_zzz_1[i] * wp_z[i];
    }

    /// Set up 120-130 components of targeted buffer : SISF

    auto g_0_xxyyzz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 120);

    auto g_0_xxyyzz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 121);

    auto g_0_xxyyzz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 122);

    auto g_0_xxyyzz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 123);

    auto g_0_xxyyzz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 124);

    auto g_0_xxyyzz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 125);

    auto g_0_xxyyzz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 126);

    auto g_0_xxyyzz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 127);

    auto g_0_xxyyzz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 128);

    auto g_0_xxyyzz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 129);

    #pragma omp simd aligned(g_0_xxyy_0_xxy_0, g_0_xxyy_0_xxy_1, g_0_xxyy_0_xyy_0, g_0_xxyy_0_xyy_1, g_0_xxyyz_0_xxy_0, g_0_xxyyz_0_xxy_1, g_0_xxyyz_0_xyy_0, g_0_xxyyz_0_xyy_1, g_0_xxyyzz_0_xxx_0, g_0_xxyyzz_0_xxy_0, g_0_xxyyzz_0_xxz_0, g_0_xxyyzz_0_xyy_0, g_0_xxyyzz_0_xyz_0, g_0_xxyyzz_0_xzz_0, g_0_xxyyzz_0_yyy_0, g_0_xxyyzz_0_yyz_0, g_0_xxyyzz_0_yzz_0, g_0_xxyyzz_0_zzz_0, g_0_xxyzz_0_xxx_0, g_0_xxyzz_0_xxx_1, g_0_xxyzz_0_xxz_0, g_0_xxyzz_0_xxz_1, g_0_xxyzz_0_xzz_0, g_0_xxyzz_0_xzz_1, g_0_xxzz_0_xxx_0, g_0_xxzz_0_xxx_1, g_0_xxzz_0_xxz_0, g_0_xxzz_0_xxz_1, g_0_xxzz_0_xzz_0, g_0_xxzz_0_xzz_1, g_0_xyyzz_0_xyz_0, g_0_xyyzz_0_xyz_1, g_0_xyyzz_0_yyy_0, g_0_xyyzz_0_yyy_1, g_0_xyyzz_0_yyz_0, g_0_xyyzz_0_yyz_1, g_0_xyyzz_0_yz_1, g_0_xyyzz_0_yzz_0, g_0_xyyzz_0_yzz_1, g_0_xyyzz_0_zzz_0, g_0_xyyzz_0_zzz_1, g_0_yyzz_0_xyz_0, g_0_yyzz_0_xyz_1, g_0_yyzz_0_yyy_0, g_0_yyzz_0_yyy_1, g_0_yyzz_0_yyz_0, g_0_yyzz_0_yyz_1, g_0_yyzz_0_yzz_0, g_0_yyzz_0_yzz_1, g_0_yyzz_0_zzz_0, g_0_yyzz_0_zzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyzz_0_xxx_0[i] = g_0_xxzz_0_xxx_0[i] * fi_ab_0 - g_0_xxzz_0_xxx_1[i] * fti_ab_0 + g_0_xxyzz_0_xxx_0[i] * pb_y + g_0_xxyzz_0_xxx_1[i] * wp_y[i];

        g_0_xxyyzz_0_xxy_0[i] = g_0_xxyy_0_xxy_0[i] * fi_ab_0 - g_0_xxyy_0_xxy_1[i] * fti_ab_0 + g_0_xxyyz_0_xxy_0[i] * pb_z + g_0_xxyyz_0_xxy_1[i] * wp_z[i];

        g_0_xxyyzz_0_xxz_0[i] = g_0_xxzz_0_xxz_0[i] * fi_ab_0 - g_0_xxzz_0_xxz_1[i] * fti_ab_0 + g_0_xxyzz_0_xxz_0[i] * pb_y + g_0_xxyzz_0_xxz_1[i] * wp_y[i];

        g_0_xxyyzz_0_xyy_0[i] = g_0_xxyy_0_xyy_0[i] * fi_ab_0 - g_0_xxyy_0_xyy_1[i] * fti_ab_0 + g_0_xxyyz_0_xyy_0[i] * pb_z + g_0_xxyyz_0_xyy_1[i] * wp_z[i];

        g_0_xxyyzz_0_xyz_0[i] = g_0_yyzz_0_xyz_0[i] * fi_ab_0 - g_0_yyzz_0_xyz_1[i] * fti_ab_0 + g_0_xyyzz_0_yz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xyz_0[i] * pb_x + g_0_xyyzz_0_xyz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xzz_0[i] = g_0_xxzz_0_xzz_0[i] * fi_ab_0 - g_0_xxzz_0_xzz_1[i] * fti_ab_0 + g_0_xxyzz_0_xzz_0[i] * pb_y + g_0_xxyzz_0_xzz_1[i] * wp_y[i];

        g_0_xxyyzz_0_yyy_0[i] = g_0_yyzz_0_yyy_0[i] * fi_ab_0 - g_0_yyzz_0_yyy_1[i] * fti_ab_0 + g_0_xyyzz_0_yyy_0[i] * pb_x + g_0_xyyzz_0_yyy_1[i] * wp_x[i];

        g_0_xxyyzz_0_yyz_0[i] = g_0_yyzz_0_yyz_0[i] * fi_ab_0 - g_0_yyzz_0_yyz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyz_0[i] * pb_x + g_0_xyyzz_0_yyz_1[i] * wp_x[i];

        g_0_xxyyzz_0_yzz_0[i] = g_0_yyzz_0_yzz_0[i] * fi_ab_0 - g_0_yyzz_0_yzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yzz_0[i] * pb_x + g_0_xyyzz_0_yzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_zzz_0[i] = g_0_yyzz_0_zzz_0[i] * fi_ab_0 - g_0_yyzz_0_zzz_1[i] * fti_ab_0 + g_0_xyyzz_0_zzz_0[i] * pb_x + g_0_xyyzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 130-140 components of targeted buffer : SISF

    auto g_0_xxyzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 130);

    auto g_0_xxyzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 131);

    auto g_0_xxyzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 132);

    auto g_0_xxyzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 133);

    auto g_0_xxyzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 134);

    auto g_0_xxyzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 135);

    auto g_0_xxyzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 136);

    auto g_0_xxyzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 137);

    auto g_0_xxyzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 138);

    auto g_0_xxyzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 139);

    #pragma omp simd aligned(g_0_xxyzzz_0_xxx_0, g_0_xxyzzz_0_xxy_0, g_0_xxyzzz_0_xxz_0, g_0_xxyzzz_0_xyy_0, g_0_xxyzzz_0_xyz_0, g_0_xxyzzz_0_xzz_0, g_0_xxyzzz_0_yyy_0, g_0_xxyzzz_0_yyz_0, g_0_xxyzzz_0_yzz_0, g_0_xxyzzz_0_zzz_0, g_0_xxzzz_0_xx_1, g_0_xxzzz_0_xxx_0, g_0_xxzzz_0_xxx_1, g_0_xxzzz_0_xxy_0, g_0_xxzzz_0_xxy_1, g_0_xxzzz_0_xxz_0, g_0_xxzzz_0_xxz_1, g_0_xxzzz_0_xy_1, g_0_xxzzz_0_xyy_0, g_0_xxzzz_0_xyy_1, g_0_xxzzz_0_xyz_0, g_0_xxzzz_0_xyz_1, g_0_xxzzz_0_xz_1, g_0_xxzzz_0_xzz_0, g_0_xxzzz_0_xzz_1, g_0_xxzzz_0_yy_1, g_0_xxzzz_0_yyy_0, g_0_xxzzz_0_yyy_1, g_0_xxzzz_0_yyz_0, g_0_xxzzz_0_yyz_1, g_0_xxzzz_0_yz_1, g_0_xxzzz_0_yzz_0, g_0_xxzzz_0_yzz_1, g_0_xxzzz_0_zz_1, g_0_xxzzz_0_zzz_0, g_0_xxzzz_0_zzz_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyzzz_0_xxx_0[i] = g_0_xxzzz_0_xxx_0[i] * pb_y + g_0_xxzzz_0_xxx_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxy_0[i] = g_0_xxzzz_0_xx_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxy_0[i] * pb_y + g_0_xxzzz_0_xxy_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxz_0[i] = g_0_xxzzz_0_xxz_0[i] * pb_y + g_0_xxzzz_0_xxz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xyy_0[i] = 2.0 * g_0_xxzzz_0_xy_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyy_0[i] * pb_y + g_0_xxzzz_0_xyy_1[i] * wp_y[i];

        g_0_xxyzzz_0_xyz_0[i] = g_0_xxzzz_0_xz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyz_0[i] * pb_y + g_0_xxzzz_0_xyz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xzz_0[i] = g_0_xxzzz_0_xzz_0[i] * pb_y + g_0_xxzzz_0_xzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_yyy_0[i] = 3.0 * g_0_xxzzz_0_yy_1[i] * fi_abcd_0 + g_0_xxzzz_0_yyy_0[i] * pb_y + g_0_xxzzz_0_yyy_1[i] * wp_y[i];

        g_0_xxyzzz_0_yyz_0[i] = 2.0 * g_0_xxzzz_0_yz_1[i] * fi_abcd_0 + g_0_xxzzz_0_yyz_0[i] * pb_y + g_0_xxzzz_0_yyz_1[i] * wp_y[i];

        g_0_xxyzzz_0_yzz_0[i] = g_0_xxzzz_0_zz_1[i] * fi_abcd_0 + g_0_xxzzz_0_yzz_0[i] * pb_y + g_0_xxzzz_0_yzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_zzz_0[i] = g_0_xxzzz_0_zzz_0[i] * pb_y + g_0_xxzzz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 140-150 components of targeted buffer : SISF

    auto g_0_xxzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 140);

    auto g_0_xxzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 141);

    auto g_0_xxzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 142);

    auto g_0_xxzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 143);

    auto g_0_xxzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 144);

    auto g_0_xxzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 145);

    auto g_0_xxzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 146);

    auto g_0_xxzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 147);

    auto g_0_xxzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 148);

    auto g_0_xxzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 149);

    #pragma omp simd aligned(g_0_xxzz_0_xxx_0, g_0_xxzz_0_xxx_1, g_0_xxzz_0_xxy_0, g_0_xxzz_0_xxy_1, g_0_xxzz_0_xyy_0, g_0_xxzz_0_xyy_1, g_0_xxzzz_0_xxx_0, g_0_xxzzz_0_xxx_1, g_0_xxzzz_0_xxy_0, g_0_xxzzz_0_xxy_1, g_0_xxzzz_0_xyy_0, g_0_xxzzz_0_xyy_1, g_0_xxzzzz_0_xxx_0, g_0_xxzzzz_0_xxy_0, g_0_xxzzzz_0_xxz_0, g_0_xxzzzz_0_xyy_0, g_0_xxzzzz_0_xyz_0, g_0_xxzzzz_0_xzz_0, g_0_xxzzzz_0_yyy_0, g_0_xxzzzz_0_yyz_0, g_0_xxzzzz_0_yzz_0, g_0_xxzzzz_0_zzz_0, g_0_xzzzz_0_xxz_0, g_0_xzzzz_0_xxz_1, g_0_xzzzz_0_xyz_0, g_0_xzzzz_0_xyz_1, g_0_xzzzz_0_xz_1, g_0_xzzzz_0_xzz_0, g_0_xzzzz_0_xzz_1, g_0_xzzzz_0_yyy_0, g_0_xzzzz_0_yyy_1, g_0_xzzzz_0_yyz_0, g_0_xzzzz_0_yyz_1, g_0_xzzzz_0_yz_1, g_0_xzzzz_0_yzz_0, g_0_xzzzz_0_yzz_1, g_0_xzzzz_0_zz_1, g_0_xzzzz_0_zzz_0, g_0_xzzzz_0_zzz_1, g_0_zzzz_0_xxz_0, g_0_zzzz_0_xxz_1, g_0_zzzz_0_xyz_0, g_0_zzzz_0_xyz_1, g_0_zzzz_0_xzz_0, g_0_zzzz_0_xzz_1, g_0_zzzz_0_yyy_0, g_0_zzzz_0_yyy_1, g_0_zzzz_0_yyz_0, g_0_zzzz_0_yyz_1, g_0_zzzz_0_yzz_0, g_0_zzzz_0_yzz_1, g_0_zzzz_0_zzz_0, g_0_zzzz_0_zzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzzzz_0_xxx_0[i] = 3.0 * g_0_xxzz_0_xxx_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxx_1[i] * fti_ab_0 + g_0_xxzzz_0_xxx_0[i] * pb_z + g_0_xxzzz_0_xxx_1[i] * wp_z[i];

        g_0_xxzzzz_0_xxy_0[i] = 3.0 * g_0_xxzz_0_xxy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxy_1[i] * fti_ab_0 + g_0_xxzzz_0_xxy_0[i] * pb_z + g_0_xxzzz_0_xxy_1[i] * wp_z[i];

        g_0_xxzzzz_0_xxz_0[i] = g_0_zzzz_0_xxz_0[i] * fi_ab_0 - g_0_zzzz_0_xxz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzz_0_xz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxz_0[i] * pb_x + g_0_xzzzz_0_xxz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xyy_0[i] = 3.0 * g_0_xxzz_0_xyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xyy_1[i] * fti_ab_0 + g_0_xxzzz_0_xyy_0[i] * pb_z + g_0_xxzzz_0_xyy_1[i] * wp_z[i];

        g_0_xxzzzz_0_xyz_0[i] = g_0_zzzz_0_xyz_0[i] * fi_ab_0 - g_0_zzzz_0_xyz_1[i] * fti_ab_0 + g_0_xzzzz_0_yz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xyz_0[i] * pb_x + g_0_xzzzz_0_xyz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xzz_0[i] = g_0_zzzz_0_xzz_0[i] * fi_ab_0 - g_0_zzzz_0_xzz_1[i] * fti_ab_0 + g_0_xzzzz_0_zz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xzz_0[i] * pb_x + g_0_xzzzz_0_xzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_yyy_0[i] = g_0_zzzz_0_yyy_0[i] * fi_ab_0 - g_0_zzzz_0_yyy_1[i] * fti_ab_0 + g_0_xzzzz_0_yyy_0[i] * pb_x + g_0_xzzzz_0_yyy_1[i] * wp_x[i];

        g_0_xxzzzz_0_yyz_0[i] = g_0_zzzz_0_yyz_0[i] * fi_ab_0 - g_0_zzzz_0_yyz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyz_0[i] * pb_x + g_0_xzzzz_0_yyz_1[i] * wp_x[i];

        g_0_xxzzzz_0_yzz_0[i] = g_0_zzzz_0_yzz_0[i] * fi_ab_0 - g_0_zzzz_0_yzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yzz_0[i] * pb_x + g_0_xzzzz_0_yzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_zzz_0[i] = g_0_zzzz_0_zzz_0[i] * fi_ab_0 - g_0_zzzz_0_zzz_1[i] * fti_ab_0 + g_0_xzzzz_0_zzz_0[i] * pb_x + g_0_xzzzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 150-160 components of targeted buffer : SISF

    auto g_0_xyyyyy_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 150);

    auto g_0_xyyyyy_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 151);

    auto g_0_xyyyyy_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 152);

    auto g_0_xyyyyy_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 153);

    auto g_0_xyyyyy_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 154);

    auto g_0_xyyyyy_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 155);

    auto g_0_xyyyyy_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 156);

    auto g_0_xyyyyy_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 157);

    auto g_0_xyyyyy_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 158);

    auto g_0_xyyyyy_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 159);

    #pragma omp simd aligned(g_0_xyyyyy_0_xxx_0, g_0_xyyyyy_0_xxy_0, g_0_xyyyyy_0_xxz_0, g_0_xyyyyy_0_xyy_0, g_0_xyyyyy_0_xyz_0, g_0_xyyyyy_0_xzz_0, g_0_xyyyyy_0_yyy_0, g_0_xyyyyy_0_yyz_0, g_0_xyyyyy_0_yzz_0, g_0_xyyyyy_0_zzz_0, g_0_yyyyy_0_xx_1, g_0_yyyyy_0_xxx_0, g_0_yyyyy_0_xxx_1, g_0_yyyyy_0_xxy_0, g_0_yyyyy_0_xxy_1, g_0_yyyyy_0_xxz_0, g_0_yyyyy_0_xxz_1, g_0_yyyyy_0_xy_1, g_0_yyyyy_0_xyy_0, g_0_yyyyy_0_xyy_1, g_0_yyyyy_0_xyz_0, g_0_yyyyy_0_xyz_1, g_0_yyyyy_0_xz_1, g_0_yyyyy_0_xzz_0, g_0_yyyyy_0_xzz_1, g_0_yyyyy_0_yy_1, g_0_yyyyy_0_yyy_0, g_0_yyyyy_0_yyy_1, g_0_yyyyy_0_yyz_0, g_0_yyyyy_0_yyz_1, g_0_yyyyy_0_yz_1, g_0_yyyyy_0_yzz_0, g_0_yyyyy_0_yzz_1, g_0_yyyyy_0_zz_1, g_0_yyyyy_0_zzz_0, g_0_yyyyy_0_zzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyy_0_xxx_0[i] = 3.0 * g_0_yyyyy_0_xx_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxx_0[i] * pb_x + g_0_yyyyy_0_xxx_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxy_0[i] = 2.0 * g_0_yyyyy_0_xy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxy_0[i] * pb_x + g_0_yyyyy_0_xxy_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxz_0[i] = 2.0 * g_0_yyyyy_0_xz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxz_0[i] * pb_x + g_0_yyyyy_0_xxz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xyy_0[i] = g_0_yyyyy_0_yy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyy_0[i] * pb_x + g_0_yyyyy_0_xyy_1[i] * wp_x[i];

        g_0_xyyyyy_0_xyz_0[i] = g_0_yyyyy_0_yz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyz_0[i] * pb_x + g_0_yyyyy_0_xyz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xzz_0[i] = g_0_yyyyy_0_zz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xzz_0[i] * pb_x + g_0_yyyyy_0_xzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_yyy_0[i] = g_0_yyyyy_0_yyy_0[i] * pb_x + g_0_yyyyy_0_yyy_1[i] * wp_x[i];

        g_0_xyyyyy_0_yyz_0[i] = g_0_yyyyy_0_yyz_0[i] * pb_x + g_0_yyyyy_0_yyz_1[i] * wp_x[i];

        g_0_xyyyyy_0_yzz_0[i] = g_0_yyyyy_0_yzz_0[i] * pb_x + g_0_yyyyy_0_yzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_zzz_0[i] = g_0_yyyyy_0_zzz_0[i] * pb_x + g_0_yyyyy_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 160-170 components of targeted buffer : SISF

    auto g_0_xyyyyz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 160);

    auto g_0_xyyyyz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 161);

    auto g_0_xyyyyz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 162);

    auto g_0_xyyyyz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 163);

    auto g_0_xyyyyz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 164);

    auto g_0_xyyyyz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 165);

    auto g_0_xyyyyz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 166);

    auto g_0_xyyyyz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 167);

    auto g_0_xyyyyz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 168);

    auto g_0_xyyyyz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 169);

    #pragma omp simd aligned(g_0_xyyyy_0_xxx_0, g_0_xyyyy_0_xxx_1, g_0_xyyyy_0_xxy_0, g_0_xyyyy_0_xxy_1, g_0_xyyyy_0_xyy_0, g_0_xyyyy_0_xyy_1, g_0_xyyyyz_0_xxx_0, g_0_xyyyyz_0_xxy_0, g_0_xyyyyz_0_xxz_0, g_0_xyyyyz_0_xyy_0, g_0_xyyyyz_0_xyz_0, g_0_xyyyyz_0_xzz_0, g_0_xyyyyz_0_yyy_0, g_0_xyyyyz_0_yyz_0, g_0_xyyyyz_0_yzz_0, g_0_xyyyyz_0_zzz_0, g_0_yyyyz_0_xxz_0, g_0_yyyyz_0_xxz_1, g_0_yyyyz_0_xyz_0, g_0_yyyyz_0_xyz_1, g_0_yyyyz_0_xz_1, g_0_yyyyz_0_xzz_0, g_0_yyyyz_0_xzz_1, g_0_yyyyz_0_yyy_0, g_0_yyyyz_0_yyy_1, g_0_yyyyz_0_yyz_0, g_0_yyyyz_0_yyz_1, g_0_yyyyz_0_yz_1, g_0_yyyyz_0_yzz_0, g_0_yyyyz_0_yzz_1, g_0_yyyyz_0_zz_1, g_0_yyyyz_0_zzz_0, g_0_yyyyz_0_zzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyz_0_xxx_0[i] = g_0_xyyyy_0_xxx_0[i] * pb_z + g_0_xyyyy_0_xxx_1[i] * wp_z[i];

        g_0_xyyyyz_0_xxy_0[i] = g_0_xyyyy_0_xxy_0[i] * pb_z + g_0_xyyyy_0_xxy_1[i] * wp_z[i];

        g_0_xyyyyz_0_xxz_0[i] = 2.0 * g_0_yyyyz_0_xz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxz_0[i] * pb_x + g_0_yyyyz_0_xxz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xyy_0[i] = g_0_xyyyy_0_xyy_0[i] * pb_z + g_0_xyyyy_0_xyy_1[i] * wp_z[i];

        g_0_xyyyyz_0_xyz_0[i] = g_0_yyyyz_0_yz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xyz_0[i] * pb_x + g_0_yyyyz_0_xyz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xzz_0[i] = g_0_yyyyz_0_zz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xzz_0[i] * pb_x + g_0_yyyyz_0_xzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_yyy_0[i] = g_0_yyyyz_0_yyy_0[i] * pb_x + g_0_yyyyz_0_yyy_1[i] * wp_x[i];

        g_0_xyyyyz_0_yyz_0[i] = g_0_yyyyz_0_yyz_0[i] * pb_x + g_0_yyyyz_0_yyz_1[i] * wp_x[i];

        g_0_xyyyyz_0_yzz_0[i] = g_0_yyyyz_0_yzz_0[i] * pb_x + g_0_yyyyz_0_yzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_zzz_0[i] = g_0_yyyyz_0_zzz_0[i] * pb_x + g_0_yyyyz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 170-180 components of targeted buffer : SISF

    auto g_0_xyyyzz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 170);

    auto g_0_xyyyzz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 171);

    auto g_0_xyyyzz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 172);

    auto g_0_xyyyzz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 173);

    auto g_0_xyyyzz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 174);

    auto g_0_xyyyzz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 175);

    auto g_0_xyyyzz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 176);

    auto g_0_xyyyzz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 177);

    auto g_0_xyyyzz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 178);

    auto g_0_xyyyzz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 179);

    #pragma omp simd aligned(g_0_xyyyzz_0_xxx_0, g_0_xyyyzz_0_xxy_0, g_0_xyyyzz_0_xxz_0, g_0_xyyyzz_0_xyy_0, g_0_xyyyzz_0_xyz_0, g_0_xyyyzz_0_xzz_0, g_0_xyyyzz_0_yyy_0, g_0_xyyyzz_0_yyz_0, g_0_xyyyzz_0_yzz_0, g_0_xyyyzz_0_zzz_0, g_0_yyyzz_0_xx_1, g_0_yyyzz_0_xxx_0, g_0_yyyzz_0_xxx_1, g_0_yyyzz_0_xxy_0, g_0_yyyzz_0_xxy_1, g_0_yyyzz_0_xxz_0, g_0_yyyzz_0_xxz_1, g_0_yyyzz_0_xy_1, g_0_yyyzz_0_xyy_0, g_0_yyyzz_0_xyy_1, g_0_yyyzz_0_xyz_0, g_0_yyyzz_0_xyz_1, g_0_yyyzz_0_xz_1, g_0_yyyzz_0_xzz_0, g_0_yyyzz_0_xzz_1, g_0_yyyzz_0_yy_1, g_0_yyyzz_0_yyy_0, g_0_yyyzz_0_yyy_1, g_0_yyyzz_0_yyz_0, g_0_yyyzz_0_yyz_1, g_0_yyyzz_0_yz_1, g_0_yyyzz_0_yzz_0, g_0_yyyzz_0_yzz_1, g_0_yyyzz_0_zz_1, g_0_yyyzz_0_zzz_0, g_0_yyyzz_0_zzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyzz_0_xxx_0[i] = 3.0 * g_0_yyyzz_0_xx_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxx_0[i] * pb_x + g_0_yyyzz_0_xxx_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxy_0[i] = 2.0 * g_0_yyyzz_0_xy_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxy_0[i] * pb_x + g_0_yyyzz_0_xxy_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxz_0[i] = 2.0 * g_0_yyyzz_0_xz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxz_0[i] * pb_x + g_0_yyyzz_0_xxz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xyy_0[i] = g_0_yyyzz_0_yy_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyy_0[i] * pb_x + g_0_yyyzz_0_xyy_1[i] * wp_x[i];

        g_0_xyyyzz_0_xyz_0[i] = g_0_yyyzz_0_yz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyz_0[i] * pb_x + g_0_yyyzz_0_xyz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xzz_0[i] = g_0_yyyzz_0_zz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xzz_0[i] * pb_x + g_0_yyyzz_0_xzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_yyy_0[i] = g_0_yyyzz_0_yyy_0[i] * pb_x + g_0_yyyzz_0_yyy_1[i] * wp_x[i];

        g_0_xyyyzz_0_yyz_0[i] = g_0_yyyzz_0_yyz_0[i] * pb_x + g_0_yyyzz_0_yyz_1[i] * wp_x[i];

        g_0_xyyyzz_0_yzz_0[i] = g_0_yyyzz_0_yzz_0[i] * pb_x + g_0_yyyzz_0_yzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_zzz_0[i] = g_0_yyyzz_0_zzz_0[i] * pb_x + g_0_yyyzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 180-190 components of targeted buffer : SISF

    auto g_0_xyyzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 180);

    auto g_0_xyyzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 181);

    auto g_0_xyyzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 182);

    auto g_0_xyyzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 183);

    auto g_0_xyyzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 184);

    auto g_0_xyyzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 185);

    auto g_0_xyyzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 186);

    auto g_0_xyyzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 187);

    auto g_0_xyyzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 188);

    auto g_0_xyyzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 189);

    #pragma omp simd aligned(g_0_xyyzzz_0_xxx_0, g_0_xyyzzz_0_xxy_0, g_0_xyyzzz_0_xxz_0, g_0_xyyzzz_0_xyy_0, g_0_xyyzzz_0_xyz_0, g_0_xyyzzz_0_xzz_0, g_0_xyyzzz_0_yyy_0, g_0_xyyzzz_0_yyz_0, g_0_xyyzzz_0_yzz_0, g_0_xyyzzz_0_zzz_0, g_0_yyzzz_0_xx_1, g_0_yyzzz_0_xxx_0, g_0_yyzzz_0_xxx_1, g_0_yyzzz_0_xxy_0, g_0_yyzzz_0_xxy_1, g_0_yyzzz_0_xxz_0, g_0_yyzzz_0_xxz_1, g_0_yyzzz_0_xy_1, g_0_yyzzz_0_xyy_0, g_0_yyzzz_0_xyy_1, g_0_yyzzz_0_xyz_0, g_0_yyzzz_0_xyz_1, g_0_yyzzz_0_xz_1, g_0_yyzzz_0_xzz_0, g_0_yyzzz_0_xzz_1, g_0_yyzzz_0_yy_1, g_0_yyzzz_0_yyy_0, g_0_yyzzz_0_yyy_1, g_0_yyzzz_0_yyz_0, g_0_yyzzz_0_yyz_1, g_0_yyzzz_0_yz_1, g_0_yyzzz_0_yzz_0, g_0_yyzzz_0_yzz_1, g_0_yyzzz_0_zz_1, g_0_yyzzz_0_zzz_0, g_0_yyzzz_0_zzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyzzz_0_xxx_0[i] = 3.0 * g_0_yyzzz_0_xx_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxx_0[i] * pb_x + g_0_yyzzz_0_xxx_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxy_0[i] = 2.0 * g_0_yyzzz_0_xy_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxy_0[i] * pb_x + g_0_yyzzz_0_xxy_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxz_0[i] = 2.0 * g_0_yyzzz_0_xz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxz_0[i] * pb_x + g_0_yyzzz_0_xxz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xyy_0[i] = g_0_yyzzz_0_yy_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyy_0[i] * pb_x + g_0_yyzzz_0_xyy_1[i] * wp_x[i];

        g_0_xyyzzz_0_xyz_0[i] = g_0_yyzzz_0_yz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyz_0[i] * pb_x + g_0_yyzzz_0_xyz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xzz_0[i] = g_0_yyzzz_0_zz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xzz_0[i] * pb_x + g_0_yyzzz_0_xzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_yyy_0[i] = g_0_yyzzz_0_yyy_0[i] * pb_x + g_0_yyzzz_0_yyy_1[i] * wp_x[i];

        g_0_xyyzzz_0_yyz_0[i] = g_0_yyzzz_0_yyz_0[i] * pb_x + g_0_yyzzz_0_yyz_1[i] * wp_x[i];

        g_0_xyyzzz_0_yzz_0[i] = g_0_yyzzz_0_yzz_0[i] * pb_x + g_0_yyzzz_0_yzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_zzz_0[i] = g_0_yyzzz_0_zzz_0[i] * pb_x + g_0_yyzzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 190-200 components of targeted buffer : SISF

    auto g_0_xyzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 190);

    auto g_0_xyzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 191);

    auto g_0_xyzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 192);

    auto g_0_xyzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 193);

    auto g_0_xyzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 194);

    auto g_0_xyzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 195);

    auto g_0_xyzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 196);

    auto g_0_xyzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 197);

    auto g_0_xyzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 198);

    auto g_0_xyzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 199);

    #pragma omp simd aligned(g_0_xyzzzz_0_xxx_0, g_0_xyzzzz_0_xxy_0, g_0_xyzzzz_0_xxz_0, g_0_xyzzzz_0_xyy_0, g_0_xyzzzz_0_xyz_0, g_0_xyzzzz_0_xzz_0, g_0_xyzzzz_0_yyy_0, g_0_xyzzzz_0_yyz_0, g_0_xyzzzz_0_yzz_0, g_0_xyzzzz_0_zzz_0, g_0_xzzzz_0_xxx_0, g_0_xzzzz_0_xxx_1, g_0_xzzzz_0_xxz_0, g_0_xzzzz_0_xxz_1, g_0_xzzzz_0_xzz_0, g_0_xzzzz_0_xzz_1, g_0_yzzzz_0_xxy_0, g_0_yzzzz_0_xxy_1, g_0_yzzzz_0_xy_1, g_0_yzzzz_0_xyy_0, g_0_yzzzz_0_xyy_1, g_0_yzzzz_0_xyz_0, g_0_yzzzz_0_xyz_1, g_0_yzzzz_0_yy_1, g_0_yzzzz_0_yyy_0, g_0_yzzzz_0_yyy_1, g_0_yzzzz_0_yyz_0, g_0_yzzzz_0_yyz_1, g_0_yzzzz_0_yz_1, g_0_yzzzz_0_yzz_0, g_0_yzzzz_0_yzz_1, g_0_yzzzz_0_zzz_0, g_0_yzzzz_0_zzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzzzz_0_xxx_0[i] = g_0_xzzzz_0_xxx_0[i] * pb_y + g_0_xzzzz_0_xxx_1[i] * wp_y[i];

        g_0_xyzzzz_0_xxy_0[i] = 2.0 * g_0_yzzzz_0_xy_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxy_0[i] * pb_x + g_0_yzzzz_0_xxy_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxz_0[i] = g_0_xzzzz_0_xxz_0[i] * pb_y + g_0_xzzzz_0_xxz_1[i] * wp_y[i];

        g_0_xyzzzz_0_xyy_0[i] = g_0_yzzzz_0_yy_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyy_0[i] * pb_x + g_0_yzzzz_0_xyy_1[i] * wp_x[i];

        g_0_xyzzzz_0_xyz_0[i] = g_0_yzzzz_0_yz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyz_0[i] * pb_x + g_0_yzzzz_0_xyz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xzz_0[i] = g_0_xzzzz_0_xzz_0[i] * pb_y + g_0_xzzzz_0_xzz_1[i] * wp_y[i];

        g_0_xyzzzz_0_yyy_0[i] = g_0_yzzzz_0_yyy_0[i] * pb_x + g_0_yzzzz_0_yyy_1[i] * wp_x[i];

        g_0_xyzzzz_0_yyz_0[i] = g_0_yzzzz_0_yyz_0[i] * pb_x + g_0_yzzzz_0_yyz_1[i] * wp_x[i];

        g_0_xyzzzz_0_yzz_0[i] = g_0_yzzzz_0_yzz_0[i] * pb_x + g_0_yzzzz_0_yzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_zzz_0[i] = g_0_yzzzz_0_zzz_0[i] * pb_x + g_0_yzzzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 200-210 components of targeted buffer : SISF

    auto g_0_xzzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 200);

    auto g_0_xzzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 201);

    auto g_0_xzzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 202);

    auto g_0_xzzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 203);

    auto g_0_xzzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 204);

    auto g_0_xzzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 205);

    auto g_0_xzzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 206);

    auto g_0_xzzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 207);

    auto g_0_xzzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 208);

    auto g_0_xzzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 209);

    #pragma omp simd aligned(g_0_xzzzzz_0_xxx_0, g_0_xzzzzz_0_xxy_0, g_0_xzzzzz_0_xxz_0, g_0_xzzzzz_0_xyy_0, g_0_xzzzzz_0_xyz_0, g_0_xzzzzz_0_xzz_0, g_0_xzzzzz_0_yyy_0, g_0_xzzzzz_0_yyz_0, g_0_xzzzzz_0_yzz_0, g_0_xzzzzz_0_zzz_0, g_0_zzzzz_0_xx_1, g_0_zzzzz_0_xxx_0, g_0_zzzzz_0_xxx_1, g_0_zzzzz_0_xxy_0, g_0_zzzzz_0_xxy_1, g_0_zzzzz_0_xxz_0, g_0_zzzzz_0_xxz_1, g_0_zzzzz_0_xy_1, g_0_zzzzz_0_xyy_0, g_0_zzzzz_0_xyy_1, g_0_zzzzz_0_xyz_0, g_0_zzzzz_0_xyz_1, g_0_zzzzz_0_xz_1, g_0_zzzzz_0_xzz_0, g_0_zzzzz_0_xzz_1, g_0_zzzzz_0_yy_1, g_0_zzzzz_0_yyy_0, g_0_zzzzz_0_yyy_1, g_0_zzzzz_0_yyz_0, g_0_zzzzz_0_yyz_1, g_0_zzzzz_0_yz_1, g_0_zzzzz_0_yzz_0, g_0_zzzzz_0_yzz_1, g_0_zzzzz_0_zz_1, g_0_zzzzz_0_zzz_0, g_0_zzzzz_0_zzz_1, wp_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzzzz_0_xxx_0[i] = 3.0 * g_0_zzzzz_0_xx_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxx_0[i] * pb_x + g_0_zzzzz_0_xxx_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxy_0[i] = 2.0 * g_0_zzzzz_0_xy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxy_0[i] * pb_x + g_0_zzzzz_0_xxy_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxz_0[i] = 2.0 * g_0_zzzzz_0_xz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxz_0[i] * pb_x + g_0_zzzzz_0_xxz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xyy_0[i] = g_0_zzzzz_0_yy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyy_0[i] * pb_x + g_0_zzzzz_0_xyy_1[i] * wp_x[i];

        g_0_xzzzzz_0_xyz_0[i] = g_0_zzzzz_0_yz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyz_0[i] * pb_x + g_0_zzzzz_0_xyz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xzz_0[i] = g_0_zzzzz_0_zz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xzz_0[i] * pb_x + g_0_zzzzz_0_xzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_yyy_0[i] = g_0_zzzzz_0_yyy_0[i] * pb_x + g_0_zzzzz_0_yyy_1[i] * wp_x[i];

        g_0_xzzzzz_0_yyz_0[i] = g_0_zzzzz_0_yyz_0[i] * pb_x + g_0_zzzzz_0_yyz_1[i] * wp_x[i];

        g_0_xzzzzz_0_yzz_0[i] = g_0_zzzzz_0_yzz_0[i] * pb_x + g_0_zzzzz_0_yzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_zzz_0[i] = g_0_zzzzz_0_zzz_0[i] * pb_x + g_0_zzzzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 210-220 components of targeted buffer : SISF

    auto g_0_yyyyyy_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 210);

    auto g_0_yyyyyy_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 211);

    auto g_0_yyyyyy_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 212);

    auto g_0_yyyyyy_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 213);

    auto g_0_yyyyyy_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 214);

    auto g_0_yyyyyy_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 215);

    auto g_0_yyyyyy_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 216);

    auto g_0_yyyyyy_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 217);

    auto g_0_yyyyyy_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 218);

    auto g_0_yyyyyy_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 219);

    #pragma omp simd aligned(g_0_yyyy_0_xxx_0, g_0_yyyy_0_xxx_1, g_0_yyyy_0_xxy_0, g_0_yyyy_0_xxy_1, g_0_yyyy_0_xxz_0, g_0_yyyy_0_xxz_1, g_0_yyyy_0_xyy_0, g_0_yyyy_0_xyy_1, g_0_yyyy_0_xyz_0, g_0_yyyy_0_xyz_1, g_0_yyyy_0_xzz_0, g_0_yyyy_0_xzz_1, g_0_yyyy_0_yyy_0, g_0_yyyy_0_yyy_1, g_0_yyyy_0_yyz_0, g_0_yyyy_0_yyz_1, g_0_yyyy_0_yzz_0, g_0_yyyy_0_yzz_1, g_0_yyyy_0_zzz_0, g_0_yyyy_0_zzz_1, g_0_yyyyy_0_xx_1, g_0_yyyyy_0_xxx_0, g_0_yyyyy_0_xxx_1, g_0_yyyyy_0_xxy_0, g_0_yyyyy_0_xxy_1, g_0_yyyyy_0_xxz_0, g_0_yyyyy_0_xxz_1, g_0_yyyyy_0_xy_1, g_0_yyyyy_0_xyy_0, g_0_yyyyy_0_xyy_1, g_0_yyyyy_0_xyz_0, g_0_yyyyy_0_xyz_1, g_0_yyyyy_0_xz_1, g_0_yyyyy_0_xzz_0, g_0_yyyyy_0_xzz_1, g_0_yyyyy_0_yy_1, g_0_yyyyy_0_yyy_0, g_0_yyyyy_0_yyy_1, g_0_yyyyy_0_yyz_0, g_0_yyyyy_0_yyz_1, g_0_yyyyy_0_yz_1, g_0_yyyyy_0_yzz_0, g_0_yyyyy_0_yzz_1, g_0_yyyyy_0_zz_1, g_0_yyyyy_0_zzz_0, g_0_yyyyy_0_zzz_1, g_0_yyyyyy_0_xxx_0, g_0_yyyyyy_0_xxy_0, g_0_yyyyyy_0_xxz_0, g_0_yyyyyy_0_xyy_0, g_0_yyyyyy_0_xyz_0, g_0_yyyyyy_0_xzz_0, g_0_yyyyyy_0_yyy_0, g_0_yyyyyy_0_yyz_0, g_0_yyyyyy_0_yzz_0, g_0_yyyyyy_0_zzz_0, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyy_0_xxx_0[i] = 5.0 * g_0_yyyy_0_xxx_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxx_1[i] * fti_ab_0 + g_0_yyyyy_0_xxx_0[i] * pb_y + g_0_yyyyy_0_xxx_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxy_0[i] = 5.0 * g_0_yyyy_0_xxy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxy_1[i] * fti_ab_0 + g_0_yyyyy_0_xx_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxy_0[i] * pb_y + g_0_yyyyy_0_xxy_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxz_0[i] = 5.0 * g_0_yyyy_0_xxz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxz_1[i] * fti_ab_0 + g_0_yyyyy_0_xxz_0[i] * pb_y + g_0_yyyyy_0_xxz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xyy_0[i] = 5.0 * g_0_yyyy_0_xyy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xyy_1[i] * fti_ab_0 + 2.0 * g_0_yyyyy_0_xy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyy_0[i] * pb_y + g_0_yyyyy_0_xyy_1[i] * wp_y[i];

        g_0_yyyyyy_0_xyz_0[i] = 5.0 * g_0_yyyy_0_xyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xyz_1[i] * fti_ab_0 + g_0_yyyyy_0_xz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyz_0[i] * pb_y + g_0_yyyyy_0_xyz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xzz_0[i] = 5.0 * g_0_yyyy_0_xzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xzz_1[i] * fti_ab_0 + g_0_yyyyy_0_xzz_0[i] * pb_y + g_0_yyyyy_0_xzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_yyy_0[i] = 5.0 * g_0_yyyy_0_yyy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yyy_1[i] * fti_ab_0 + 3.0 * g_0_yyyyy_0_yy_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyy_0[i] * pb_y + g_0_yyyyy_0_yyy_1[i] * wp_y[i];

        g_0_yyyyyy_0_yyz_0[i] = 5.0 * g_0_yyyy_0_yyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyy_0_yz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyz_0[i] * pb_y + g_0_yyyyy_0_yyz_1[i] * wp_y[i];

        g_0_yyyyyy_0_yzz_0[i] = 5.0 * g_0_yyyy_0_yzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yzz_1[i] * fti_ab_0 + g_0_yyyyy_0_zz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yzz_0[i] * pb_y + g_0_yyyyy_0_yzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_zzz_0[i] = 5.0 * g_0_yyyy_0_zzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_zzz_1[i] * fti_ab_0 + g_0_yyyyy_0_zzz_0[i] * pb_y + g_0_yyyyy_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 220-230 components of targeted buffer : SISF

    auto g_0_yyyyyz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 220);

    auto g_0_yyyyyz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 221);

    auto g_0_yyyyyz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 222);

    auto g_0_yyyyyz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 223);

    auto g_0_yyyyyz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 224);

    auto g_0_yyyyyz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 225);

    auto g_0_yyyyyz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 226);

    auto g_0_yyyyyz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 227);

    auto g_0_yyyyyz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 228);

    auto g_0_yyyyyz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 229);

    #pragma omp simd aligned(g_0_yyyyy_0_xx_1, g_0_yyyyy_0_xxx_0, g_0_yyyyy_0_xxx_1, g_0_yyyyy_0_xxy_0, g_0_yyyyy_0_xxy_1, g_0_yyyyy_0_xxz_0, g_0_yyyyy_0_xxz_1, g_0_yyyyy_0_xy_1, g_0_yyyyy_0_xyy_0, g_0_yyyyy_0_xyy_1, g_0_yyyyy_0_xyz_0, g_0_yyyyy_0_xyz_1, g_0_yyyyy_0_xz_1, g_0_yyyyy_0_xzz_0, g_0_yyyyy_0_xzz_1, g_0_yyyyy_0_yy_1, g_0_yyyyy_0_yyy_0, g_0_yyyyy_0_yyy_1, g_0_yyyyy_0_yyz_0, g_0_yyyyy_0_yyz_1, g_0_yyyyy_0_yz_1, g_0_yyyyy_0_yzz_0, g_0_yyyyy_0_yzz_1, g_0_yyyyy_0_zz_1, g_0_yyyyy_0_zzz_0, g_0_yyyyy_0_zzz_1, g_0_yyyyyz_0_xxx_0, g_0_yyyyyz_0_xxy_0, g_0_yyyyyz_0_xxz_0, g_0_yyyyyz_0_xyy_0, g_0_yyyyyz_0_xyz_0, g_0_yyyyyz_0_xzz_0, g_0_yyyyyz_0_yyy_0, g_0_yyyyyz_0_yyz_0, g_0_yyyyyz_0_yzz_0, g_0_yyyyyz_0_zzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyyz_0_xxx_0[i] = g_0_yyyyy_0_xxx_0[i] * pb_z + g_0_yyyyy_0_xxx_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxy_0[i] = g_0_yyyyy_0_xxy_0[i] * pb_z + g_0_yyyyy_0_xxy_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxz_0[i] = g_0_yyyyy_0_xx_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxz_0[i] * pb_z + g_0_yyyyy_0_xxz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xyy_0[i] = g_0_yyyyy_0_xyy_0[i] * pb_z + g_0_yyyyy_0_xyy_1[i] * wp_z[i];

        g_0_yyyyyz_0_xyz_0[i] = g_0_yyyyy_0_xy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyz_0[i] * pb_z + g_0_yyyyy_0_xyz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xzz_0[i] = 2.0 * g_0_yyyyy_0_xz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xzz_0[i] * pb_z + g_0_yyyyy_0_xzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_yyy_0[i] = g_0_yyyyy_0_yyy_0[i] * pb_z + g_0_yyyyy_0_yyy_1[i] * wp_z[i];

        g_0_yyyyyz_0_yyz_0[i] = g_0_yyyyy_0_yy_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyz_0[i] * pb_z + g_0_yyyyy_0_yyz_1[i] * wp_z[i];

        g_0_yyyyyz_0_yzz_0[i] = 2.0 * g_0_yyyyy_0_yz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yzz_0[i] * pb_z + g_0_yyyyy_0_yzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_zzz_0[i] = 3.0 * g_0_yyyyy_0_zz_1[i] * fi_abcd_0 + g_0_yyyyy_0_zzz_0[i] * pb_z + g_0_yyyyy_0_zzz_1[i] * wp_z[i];
    }

    /// Set up 230-240 components of targeted buffer : SISF

    auto g_0_yyyyzz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 230);

    auto g_0_yyyyzz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 231);

    auto g_0_yyyyzz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 232);

    auto g_0_yyyyzz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 233);

    auto g_0_yyyyzz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 234);

    auto g_0_yyyyzz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 235);

    auto g_0_yyyyzz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 236);

    auto g_0_yyyyzz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 237);

    auto g_0_yyyyzz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 238);

    auto g_0_yyyyzz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 239);

    #pragma omp simd aligned(g_0_yyyy_0_xxy_0, g_0_yyyy_0_xxy_1, g_0_yyyy_0_xyy_0, g_0_yyyy_0_xyy_1, g_0_yyyy_0_yyy_0, g_0_yyyy_0_yyy_1, g_0_yyyyz_0_xxy_0, g_0_yyyyz_0_xxy_1, g_0_yyyyz_0_xyy_0, g_0_yyyyz_0_xyy_1, g_0_yyyyz_0_yyy_0, g_0_yyyyz_0_yyy_1, g_0_yyyyzz_0_xxx_0, g_0_yyyyzz_0_xxy_0, g_0_yyyyzz_0_xxz_0, g_0_yyyyzz_0_xyy_0, g_0_yyyyzz_0_xyz_0, g_0_yyyyzz_0_xzz_0, g_0_yyyyzz_0_yyy_0, g_0_yyyyzz_0_yyz_0, g_0_yyyyzz_0_yzz_0, g_0_yyyyzz_0_zzz_0, g_0_yyyzz_0_xxx_0, g_0_yyyzz_0_xxx_1, g_0_yyyzz_0_xxz_0, g_0_yyyzz_0_xxz_1, g_0_yyyzz_0_xyz_0, g_0_yyyzz_0_xyz_1, g_0_yyyzz_0_xz_1, g_0_yyyzz_0_xzz_0, g_0_yyyzz_0_xzz_1, g_0_yyyzz_0_yyz_0, g_0_yyyzz_0_yyz_1, g_0_yyyzz_0_yz_1, g_0_yyyzz_0_yzz_0, g_0_yyyzz_0_yzz_1, g_0_yyyzz_0_zz_1, g_0_yyyzz_0_zzz_0, g_0_yyyzz_0_zzz_1, g_0_yyzz_0_xxx_0, g_0_yyzz_0_xxx_1, g_0_yyzz_0_xxz_0, g_0_yyzz_0_xxz_1, g_0_yyzz_0_xyz_0, g_0_yyzz_0_xyz_1, g_0_yyzz_0_xzz_0, g_0_yyzz_0_xzz_1, g_0_yyzz_0_yyz_0, g_0_yyzz_0_yyz_1, g_0_yyzz_0_yzz_0, g_0_yyzz_0_yzz_1, g_0_yyzz_0_zzz_0, g_0_yyzz_0_zzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyzz_0_xxx_0[i] = 3.0 * g_0_yyzz_0_xxx_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxx_1[i] * fti_ab_0 + g_0_yyyzz_0_xxx_0[i] * pb_y + g_0_yyyzz_0_xxx_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxy_0[i] = g_0_yyyy_0_xxy_0[i] * fi_ab_0 - g_0_yyyy_0_xxy_1[i] * fti_ab_0 + g_0_yyyyz_0_xxy_0[i] * pb_z + g_0_yyyyz_0_xxy_1[i] * wp_z[i];

        g_0_yyyyzz_0_xxz_0[i] = 3.0 * g_0_yyzz_0_xxz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxz_1[i] * fti_ab_0 + g_0_yyyzz_0_xxz_0[i] * pb_y + g_0_yyyzz_0_xxz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xyy_0[i] = g_0_yyyy_0_xyy_0[i] * fi_ab_0 - g_0_yyyy_0_xyy_1[i] * fti_ab_0 + g_0_yyyyz_0_xyy_0[i] * pb_z + g_0_yyyyz_0_xyy_1[i] * wp_z[i];

        g_0_yyyyzz_0_xyz_0[i] = 3.0 * g_0_yyzz_0_xyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xyz_1[i] * fti_ab_0 + g_0_yyyzz_0_xz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyz_0[i] * pb_y + g_0_yyyzz_0_xyz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xzz_0[i] = 3.0 * g_0_yyzz_0_xzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xzz_1[i] * fti_ab_0 + g_0_yyyzz_0_xzz_0[i] * pb_y + g_0_yyyzz_0_xzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_yyy_0[i] = g_0_yyyy_0_yyy_0[i] * fi_ab_0 - g_0_yyyy_0_yyy_1[i] * fti_ab_0 + g_0_yyyyz_0_yyy_0[i] * pb_z + g_0_yyyyz_0_yyy_1[i] * wp_z[i];

        g_0_yyyyzz_0_yyz_0[i] = 3.0 * g_0_yyzz_0_yyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyzz_0_yz_1[i] * fi_abcd_0 + g_0_yyyzz_0_yyz_0[i] * pb_y + g_0_yyyzz_0_yyz_1[i] * wp_y[i];

        g_0_yyyyzz_0_yzz_0[i] = 3.0 * g_0_yyzz_0_yzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yzz_1[i] * fti_ab_0 + g_0_yyyzz_0_zz_1[i] * fi_abcd_0 + g_0_yyyzz_0_yzz_0[i] * pb_y + g_0_yyyzz_0_yzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_zzz_0[i] = 3.0 * g_0_yyzz_0_zzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_zzz_1[i] * fti_ab_0 + g_0_yyyzz_0_zzz_0[i] * pb_y + g_0_yyyzz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 240-250 components of targeted buffer : SISF

    auto g_0_yyyzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 240);

    auto g_0_yyyzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 241);

    auto g_0_yyyzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 242);

    auto g_0_yyyzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 243);

    auto g_0_yyyzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 244);

    auto g_0_yyyzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 245);

    auto g_0_yyyzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 246);

    auto g_0_yyyzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 247);

    auto g_0_yyyzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 248);

    auto g_0_yyyzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 249);

    #pragma omp simd aligned(g_0_yyyz_0_xxy_0, g_0_yyyz_0_xxy_1, g_0_yyyz_0_xyy_0, g_0_yyyz_0_xyy_1, g_0_yyyz_0_yyy_0, g_0_yyyz_0_yyy_1, g_0_yyyzz_0_xxy_0, g_0_yyyzz_0_xxy_1, g_0_yyyzz_0_xyy_0, g_0_yyyzz_0_xyy_1, g_0_yyyzz_0_yyy_0, g_0_yyyzz_0_yyy_1, g_0_yyyzzz_0_xxx_0, g_0_yyyzzz_0_xxy_0, g_0_yyyzzz_0_xxz_0, g_0_yyyzzz_0_xyy_0, g_0_yyyzzz_0_xyz_0, g_0_yyyzzz_0_xzz_0, g_0_yyyzzz_0_yyy_0, g_0_yyyzzz_0_yyz_0, g_0_yyyzzz_0_yzz_0, g_0_yyyzzz_0_zzz_0, g_0_yyzzz_0_xxx_0, g_0_yyzzz_0_xxx_1, g_0_yyzzz_0_xxz_0, g_0_yyzzz_0_xxz_1, g_0_yyzzz_0_xyz_0, g_0_yyzzz_0_xyz_1, g_0_yyzzz_0_xz_1, g_0_yyzzz_0_xzz_0, g_0_yyzzz_0_xzz_1, g_0_yyzzz_0_yyz_0, g_0_yyzzz_0_yyz_1, g_0_yyzzz_0_yz_1, g_0_yyzzz_0_yzz_0, g_0_yyzzz_0_yzz_1, g_0_yyzzz_0_zz_1, g_0_yyzzz_0_zzz_0, g_0_yyzzz_0_zzz_1, g_0_yzzz_0_xxx_0, g_0_yzzz_0_xxx_1, g_0_yzzz_0_xxz_0, g_0_yzzz_0_xxz_1, g_0_yzzz_0_xyz_0, g_0_yzzz_0_xyz_1, g_0_yzzz_0_xzz_0, g_0_yzzz_0_xzz_1, g_0_yzzz_0_yyz_0, g_0_yzzz_0_yyz_1, g_0_yzzz_0_yzz_0, g_0_yzzz_0_yzz_1, g_0_yzzz_0_zzz_0, g_0_yzzz_0_zzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyzzz_0_xxx_0[i] = 2.0 * g_0_yzzz_0_xxx_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxx_1[i] * fti_ab_0 + g_0_yyzzz_0_xxx_0[i] * pb_y + g_0_yyzzz_0_xxx_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxy_0[i] = 2.0 * g_0_yyyz_0_xxy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_xxy_1[i] * fti_ab_0 + g_0_yyyzz_0_xxy_0[i] * pb_z + g_0_yyyzz_0_xxy_1[i] * wp_z[i];

        g_0_yyyzzz_0_xxz_0[i] = 2.0 * g_0_yzzz_0_xxz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxz_1[i] * fti_ab_0 + g_0_yyzzz_0_xxz_0[i] * pb_y + g_0_yyzzz_0_xxz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xyy_0[i] = 2.0 * g_0_yyyz_0_xyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_xyy_1[i] * fti_ab_0 + g_0_yyyzz_0_xyy_0[i] * pb_z + g_0_yyyzz_0_xyy_1[i] * wp_z[i];

        g_0_yyyzzz_0_xyz_0[i] = 2.0 * g_0_yzzz_0_xyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xyz_1[i] * fti_ab_0 + g_0_yyzzz_0_xz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyz_0[i] * pb_y + g_0_yyzzz_0_xyz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xzz_0[i] = 2.0 * g_0_yzzz_0_xzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xzz_1[i] * fti_ab_0 + g_0_yyzzz_0_xzz_0[i] * pb_y + g_0_yyzzz_0_xzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_yyy_0[i] = 2.0 * g_0_yyyz_0_yyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_yyy_1[i] * fti_ab_0 + g_0_yyyzz_0_yyy_0[i] * pb_z + g_0_yyyzz_0_yyy_1[i] * wp_z[i];

        g_0_yyyzzz_0_yyz_0[i] = 2.0 * g_0_yzzz_0_yyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_yyz_1[i] * fti_ab_0 + 2.0 * g_0_yyzzz_0_yz_1[i] * fi_abcd_0 + g_0_yyzzz_0_yyz_0[i] * pb_y + g_0_yyzzz_0_yyz_1[i] * wp_y[i];

        g_0_yyyzzz_0_yzz_0[i] = 2.0 * g_0_yzzz_0_yzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_yzz_1[i] * fti_ab_0 + g_0_yyzzz_0_zz_1[i] * fi_abcd_0 + g_0_yyzzz_0_yzz_0[i] * pb_y + g_0_yyzzz_0_yzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_zzz_0[i] = 2.0 * g_0_yzzz_0_zzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_zzz_1[i] * fti_ab_0 + g_0_yyzzz_0_zzz_0[i] * pb_y + g_0_yyzzz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 250-260 components of targeted buffer : SISF

    auto g_0_yyzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 250);

    auto g_0_yyzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 251);

    auto g_0_yyzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 252);

    auto g_0_yyzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 253);

    auto g_0_yyzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 254);

    auto g_0_yyzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 255);

    auto g_0_yyzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 256);

    auto g_0_yyzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 257);

    auto g_0_yyzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 258);

    auto g_0_yyzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 259);

    #pragma omp simd aligned(g_0_yyzz_0_xxy_0, g_0_yyzz_0_xxy_1, g_0_yyzz_0_xyy_0, g_0_yyzz_0_xyy_1, g_0_yyzz_0_yyy_0, g_0_yyzz_0_yyy_1, g_0_yyzzz_0_xxy_0, g_0_yyzzz_0_xxy_1, g_0_yyzzz_0_xyy_0, g_0_yyzzz_0_xyy_1, g_0_yyzzz_0_yyy_0, g_0_yyzzz_0_yyy_1, g_0_yyzzzz_0_xxx_0, g_0_yyzzzz_0_xxy_0, g_0_yyzzzz_0_xxz_0, g_0_yyzzzz_0_xyy_0, g_0_yyzzzz_0_xyz_0, g_0_yyzzzz_0_xzz_0, g_0_yyzzzz_0_yyy_0, g_0_yyzzzz_0_yyz_0, g_0_yyzzzz_0_yzz_0, g_0_yyzzzz_0_zzz_0, g_0_yzzzz_0_xxx_0, g_0_yzzzz_0_xxx_1, g_0_yzzzz_0_xxz_0, g_0_yzzzz_0_xxz_1, g_0_yzzzz_0_xyz_0, g_0_yzzzz_0_xyz_1, g_0_yzzzz_0_xz_1, g_0_yzzzz_0_xzz_0, g_0_yzzzz_0_xzz_1, g_0_yzzzz_0_yyz_0, g_0_yzzzz_0_yyz_1, g_0_yzzzz_0_yz_1, g_0_yzzzz_0_yzz_0, g_0_yzzzz_0_yzz_1, g_0_yzzzz_0_zz_1, g_0_yzzzz_0_zzz_0, g_0_yzzzz_0_zzz_1, g_0_zzzz_0_xxx_0, g_0_zzzz_0_xxx_1, g_0_zzzz_0_xxz_0, g_0_zzzz_0_xxz_1, g_0_zzzz_0_xyz_0, g_0_zzzz_0_xyz_1, g_0_zzzz_0_xzz_0, g_0_zzzz_0_xzz_1, g_0_zzzz_0_yyz_0, g_0_zzzz_0_yyz_1, g_0_zzzz_0_yzz_0, g_0_zzzz_0_yzz_1, g_0_zzzz_0_zzz_0, g_0_zzzz_0_zzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzzzz_0_xxx_0[i] = g_0_zzzz_0_xxx_0[i] * fi_ab_0 - g_0_zzzz_0_xxx_1[i] * fti_ab_0 + g_0_yzzzz_0_xxx_0[i] * pb_y + g_0_yzzzz_0_xxx_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxy_0[i] = 3.0 * g_0_yyzz_0_xxy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxy_1[i] * fti_ab_0 + g_0_yyzzz_0_xxy_0[i] * pb_z + g_0_yyzzz_0_xxy_1[i] * wp_z[i];

        g_0_yyzzzz_0_xxz_0[i] = g_0_zzzz_0_xxz_0[i] * fi_ab_0 - g_0_zzzz_0_xxz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxz_0[i] * pb_y + g_0_yzzzz_0_xxz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xyy_0[i] = 3.0 * g_0_yyzz_0_xyy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xyy_1[i] * fti_ab_0 + g_0_yyzzz_0_xyy_0[i] * pb_z + g_0_yyzzz_0_xyy_1[i] * wp_z[i];

        g_0_yyzzzz_0_xyz_0[i] = g_0_zzzz_0_xyz_0[i] * fi_ab_0 - g_0_zzzz_0_xyz_1[i] * fti_ab_0 + g_0_yzzzz_0_xz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyz_0[i] * pb_y + g_0_yzzzz_0_xyz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xzz_0[i] = g_0_zzzz_0_xzz_0[i] * fi_ab_0 - g_0_zzzz_0_xzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xzz_0[i] * pb_y + g_0_yzzzz_0_xzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_yyy_0[i] = 3.0 * g_0_yyzz_0_yyy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yyy_1[i] * fti_ab_0 + g_0_yyzzz_0_yyy_0[i] * pb_z + g_0_yyzzz_0_yyy_1[i] * wp_z[i];

        g_0_yyzzzz_0_yyz_0[i] = g_0_zzzz_0_yyz_0[i] * fi_ab_0 - g_0_zzzz_0_yyz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzz_0_yz_1[i] * fi_abcd_0 + g_0_yzzzz_0_yyz_0[i] * pb_y + g_0_yzzzz_0_yyz_1[i] * wp_y[i];

        g_0_yyzzzz_0_yzz_0[i] = g_0_zzzz_0_yzz_0[i] * fi_ab_0 - g_0_zzzz_0_yzz_1[i] * fti_ab_0 + g_0_yzzzz_0_zz_1[i] * fi_abcd_0 + g_0_yzzzz_0_yzz_0[i] * pb_y + g_0_yzzzz_0_yzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_zzz_0[i] = g_0_zzzz_0_zzz_0[i] * fi_ab_0 - g_0_zzzz_0_zzz_1[i] * fti_ab_0 + g_0_yzzzz_0_zzz_0[i] * pb_y + g_0_yzzzz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 260-270 components of targeted buffer : SISF

    auto g_0_yzzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 260);

    auto g_0_yzzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 261);

    auto g_0_yzzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 262);

    auto g_0_yzzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 263);

    auto g_0_yzzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 264);

    auto g_0_yzzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 265);

    auto g_0_yzzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 266);

    auto g_0_yzzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 267);

    auto g_0_yzzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 268);

    auto g_0_yzzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 269);

    #pragma omp simd aligned(g_0_yzzzzz_0_xxx_0, g_0_yzzzzz_0_xxy_0, g_0_yzzzzz_0_xxz_0, g_0_yzzzzz_0_xyy_0, g_0_yzzzzz_0_xyz_0, g_0_yzzzzz_0_xzz_0, g_0_yzzzzz_0_yyy_0, g_0_yzzzzz_0_yyz_0, g_0_yzzzzz_0_yzz_0, g_0_yzzzzz_0_zzz_0, g_0_zzzzz_0_xx_1, g_0_zzzzz_0_xxx_0, g_0_zzzzz_0_xxx_1, g_0_zzzzz_0_xxy_0, g_0_zzzzz_0_xxy_1, g_0_zzzzz_0_xxz_0, g_0_zzzzz_0_xxz_1, g_0_zzzzz_0_xy_1, g_0_zzzzz_0_xyy_0, g_0_zzzzz_0_xyy_1, g_0_zzzzz_0_xyz_0, g_0_zzzzz_0_xyz_1, g_0_zzzzz_0_xz_1, g_0_zzzzz_0_xzz_0, g_0_zzzzz_0_xzz_1, g_0_zzzzz_0_yy_1, g_0_zzzzz_0_yyy_0, g_0_zzzzz_0_yyy_1, g_0_zzzzz_0_yyz_0, g_0_zzzzz_0_yyz_1, g_0_zzzzz_0_yz_1, g_0_zzzzz_0_yzz_0, g_0_zzzzz_0_yzz_1, g_0_zzzzz_0_zz_1, g_0_zzzzz_0_zzz_0, g_0_zzzzz_0_zzz_1, wp_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzzzz_0_xxx_0[i] = g_0_zzzzz_0_xxx_0[i] * pb_y + g_0_zzzzz_0_xxx_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxy_0[i] = g_0_zzzzz_0_xx_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxy_0[i] * pb_y + g_0_zzzzz_0_xxy_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxz_0[i] = g_0_zzzzz_0_xxz_0[i] * pb_y + g_0_zzzzz_0_xxz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xyy_0[i] = 2.0 * g_0_zzzzz_0_xy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyy_0[i] * pb_y + g_0_zzzzz_0_xyy_1[i] * wp_y[i];

        g_0_yzzzzz_0_xyz_0[i] = g_0_zzzzz_0_xz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyz_0[i] * pb_y + g_0_zzzzz_0_xyz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xzz_0[i] = g_0_zzzzz_0_xzz_0[i] * pb_y + g_0_zzzzz_0_xzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_yyy_0[i] = 3.0 * g_0_zzzzz_0_yy_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyy_0[i] * pb_y + g_0_zzzzz_0_yyy_1[i] * wp_y[i];

        g_0_yzzzzz_0_yyz_0[i] = 2.0 * g_0_zzzzz_0_yz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyz_0[i] * pb_y + g_0_zzzzz_0_yyz_1[i] * wp_y[i];

        g_0_yzzzzz_0_yzz_0[i] = g_0_zzzzz_0_zz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yzz_0[i] * pb_y + g_0_zzzzz_0_yzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_zzz_0[i] = g_0_zzzzz_0_zzz_0[i] * pb_y + g_0_zzzzz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 270-280 components of targeted buffer : SISF

    auto g_0_zzzzzz_0_xxx_0 = pbuffer.data(idx_eri_0_sisf + 270);

    auto g_0_zzzzzz_0_xxy_0 = pbuffer.data(idx_eri_0_sisf + 271);

    auto g_0_zzzzzz_0_xxz_0 = pbuffer.data(idx_eri_0_sisf + 272);

    auto g_0_zzzzzz_0_xyy_0 = pbuffer.data(idx_eri_0_sisf + 273);

    auto g_0_zzzzzz_0_xyz_0 = pbuffer.data(idx_eri_0_sisf + 274);

    auto g_0_zzzzzz_0_xzz_0 = pbuffer.data(idx_eri_0_sisf + 275);

    auto g_0_zzzzzz_0_yyy_0 = pbuffer.data(idx_eri_0_sisf + 276);

    auto g_0_zzzzzz_0_yyz_0 = pbuffer.data(idx_eri_0_sisf + 277);

    auto g_0_zzzzzz_0_yzz_0 = pbuffer.data(idx_eri_0_sisf + 278);

    auto g_0_zzzzzz_0_zzz_0 = pbuffer.data(idx_eri_0_sisf + 279);

    #pragma omp simd aligned(g_0_zzzz_0_xxx_0, g_0_zzzz_0_xxx_1, g_0_zzzz_0_xxy_0, g_0_zzzz_0_xxy_1, g_0_zzzz_0_xxz_0, g_0_zzzz_0_xxz_1, g_0_zzzz_0_xyy_0, g_0_zzzz_0_xyy_1, g_0_zzzz_0_xyz_0, g_0_zzzz_0_xyz_1, g_0_zzzz_0_xzz_0, g_0_zzzz_0_xzz_1, g_0_zzzz_0_yyy_0, g_0_zzzz_0_yyy_1, g_0_zzzz_0_yyz_0, g_0_zzzz_0_yyz_1, g_0_zzzz_0_yzz_0, g_0_zzzz_0_yzz_1, g_0_zzzz_0_zzz_0, g_0_zzzz_0_zzz_1, g_0_zzzzz_0_xx_1, g_0_zzzzz_0_xxx_0, g_0_zzzzz_0_xxx_1, g_0_zzzzz_0_xxy_0, g_0_zzzzz_0_xxy_1, g_0_zzzzz_0_xxz_0, g_0_zzzzz_0_xxz_1, g_0_zzzzz_0_xy_1, g_0_zzzzz_0_xyy_0, g_0_zzzzz_0_xyy_1, g_0_zzzzz_0_xyz_0, g_0_zzzzz_0_xyz_1, g_0_zzzzz_0_xz_1, g_0_zzzzz_0_xzz_0, g_0_zzzzz_0_xzz_1, g_0_zzzzz_0_yy_1, g_0_zzzzz_0_yyy_0, g_0_zzzzz_0_yyy_1, g_0_zzzzz_0_yyz_0, g_0_zzzzz_0_yyz_1, g_0_zzzzz_0_yz_1, g_0_zzzzz_0_yzz_0, g_0_zzzzz_0_yzz_1, g_0_zzzzz_0_zz_1, g_0_zzzzz_0_zzz_0, g_0_zzzzz_0_zzz_1, g_0_zzzzzz_0_xxx_0, g_0_zzzzzz_0_xxy_0, g_0_zzzzzz_0_xxz_0, g_0_zzzzzz_0_xyy_0, g_0_zzzzzz_0_xyz_0, g_0_zzzzzz_0_xzz_0, g_0_zzzzzz_0_yyy_0, g_0_zzzzzz_0_yyz_0, g_0_zzzzzz_0_yzz_0, g_0_zzzzzz_0_zzz_0, wp_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzzzz_0_xxx_0[i] = 5.0 * g_0_zzzz_0_xxx_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxx_1[i] * fti_ab_0 + g_0_zzzzz_0_xxx_0[i] * pb_z + g_0_zzzzz_0_xxx_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxy_0[i] = 5.0 * g_0_zzzz_0_xxy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxy_1[i] * fti_ab_0 + g_0_zzzzz_0_xxy_0[i] * pb_z + g_0_zzzzz_0_xxy_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxz_0[i] = 5.0 * g_0_zzzz_0_xxz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxz_1[i] * fti_ab_0 + g_0_zzzzz_0_xx_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxz_0[i] * pb_z + g_0_zzzzz_0_xxz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xyy_0[i] = 5.0 * g_0_zzzz_0_xyy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xyy_1[i] * fti_ab_0 + g_0_zzzzz_0_xyy_0[i] * pb_z + g_0_zzzzz_0_xyy_1[i] * wp_z[i];

        g_0_zzzzzz_0_xyz_0[i] = 5.0 * g_0_zzzz_0_xyz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xyz_1[i] * fti_ab_0 + g_0_zzzzz_0_xy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyz_0[i] * pb_z + g_0_zzzzz_0_xyz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xzz_0[i] = 5.0 * g_0_zzzz_0_xzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzz_0_xz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xzz_0[i] * pb_z + g_0_zzzzz_0_xzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_yyy_0[i] = 5.0 * g_0_zzzz_0_yyy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yyy_1[i] * fti_ab_0 + g_0_zzzzz_0_yyy_0[i] * pb_z + g_0_zzzzz_0_yyy_1[i] * wp_z[i];

        g_0_zzzzzz_0_yyz_0[i] = 5.0 * g_0_zzzz_0_yyz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yyz_1[i] * fti_ab_0 + g_0_zzzzz_0_yy_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyz_0[i] * pb_z + g_0_zzzzz_0_yyz_1[i] * wp_z[i];

        g_0_zzzzzz_0_yzz_0[i] = 5.0 * g_0_zzzz_0_yzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzz_0_yz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yzz_0[i] * pb_z + g_0_zzzzz_0_yzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_zzz_0[i] = 5.0 * g_0_zzzz_0_zzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_zzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzzz_0_zz_1[i] * fi_abcd_0 + g_0_zzzzz_0_zzz_0[i] * pb_z + g_0_zzzzz_0_zzz_1[i] * wp_z[i];
    }
}

} // erirec namespace

