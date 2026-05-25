#include "LocalCorePotentialPrimRecHD.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_hd(CSimdArray<double>& pbuffer, 
                                  const size_t idx_hd,
                                  const size_t idx_fd,
                                  const size_t idx_gp,
                                  const size_t idx_gd,
                                  const CSimdArray<double>& factors,
                                  const size_t idx_ra,
                                  const size_t idx_zeta) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(RA) distances

    auto ra_x = factors.data(idx_ra);

    auto ra_y = factors.data(idx_ra + 1);

    auto ra_z = factors.data(idx_ra + 2);

    // Set up inverted 1/2xi

    auto fxi = factors.data(idx_zeta);

    // Set up components of auxiliary buffer : FD

    auto tg_xxx_xx = pbuffer.data(idx_fd);

    auto tg_xxx_xy = pbuffer.data(idx_fd + 1);

    auto tg_xxx_xz = pbuffer.data(idx_fd + 2);

    auto tg_xxx_yy = pbuffer.data(idx_fd + 3);

    auto tg_xxx_yz = pbuffer.data(idx_fd + 4);

    auto tg_xxx_zz = pbuffer.data(idx_fd + 5);

    auto tg_xxy_xx = pbuffer.data(idx_fd + 6);

    auto tg_xxy_xz = pbuffer.data(idx_fd + 8);

    auto tg_xxy_yy = pbuffer.data(idx_fd + 9);

    auto tg_xxy_yz = pbuffer.data(idx_fd + 10);

    auto tg_xxz_xx = pbuffer.data(idx_fd + 12);

    auto tg_xxz_xy = pbuffer.data(idx_fd + 13);

    auto tg_xxz_xz = pbuffer.data(idx_fd + 14);

    auto tg_xxz_yz = pbuffer.data(idx_fd + 16);

    auto tg_xxz_zz = pbuffer.data(idx_fd + 17);

    auto tg_xyy_xy = pbuffer.data(idx_fd + 19);

    auto tg_xyy_yy = pbuffer.data(idx_fd + 21);

    auto tg_xyy_yz = pbuffer.data(idx_fd + 22);

    auto tg_xyy_zz = pbuffer.data(idx_fd + 23);

    auto tg_xyz_yz = pbuffer.data(idx_fd + 28);

    auto tg_xzz_xz = pbuffer.data(idx_fd + 32);

    auto tg_xzz_yy = pbuffer.data(idx_fd + 33);

    auto tg_xzz_yz = pbuffer.data(idx_fd + 34);

    auto tg_xzz_zz = pbuffer.data(idx_fd + 35);

    auto tg_yyy_xx = pbuffer.data(idx_fd + 36);

    auto tg_yyy_xy = pbuffer.data(idx_fd + 37);

    auto tg_yyy_xz = pbuffer.data(idx_fd + 38);

    auto tg_yyy_yy = pbuffer.data(idx_fd + 39);

    auto tg_yyy_yz = pbuffer.data(idx_fd + 40);

    auto tg_yyy_zz = pbuffer.data(idx_fd + 41);

    auto tg_yyz_xy = pbuffer.data(idx_fd + 43);

    auto tg_yyz_xz = pbuffer.data(idx_fd + 44);

    auto tg_yyz_yy = pbuffer.data(idx_fd + 45);

    auto tg_yyz_yz = pbuffer.data(idx_fd + 46);

    auto tg_yyz_zz = pbuffer.data(idx_fd + 47);

    auto tg_yzz_xx = pbuffer.data(idx_fd + 48);

    auto tg_yzz_xz = pbuffer.data(idx_fd + 50);

    auto tg_yzz_yy = pbuffer.data(idx_fd + 51);

    auto tg_yzz_yz = pbuffer.data(idx_fd + 52);

    auto tg_yzz_zz = pbuffer.data(idx_fd + 53);

    auto tg_zzz_xx = pbuffer.data(idx_fd + 54);

    auto tg_zzz_xy = pbuffer.data(idx_fd + 55);

    auto tg_zzz_xz = pbuffer.data(idx_fd + 56);

    auto tg_zzz_yy = pbuffer.data(idx_fd + 57);

    auto tg_zzz_yz = pbuffer.data(idx_fd + 58);

    auto tg_zzz_zz = pbuffer.data(idx_fd + 59);

    // Set up components of auxiliary buffer : GP

    auto tg_xxxx_x = pbuffer.data(idx_gp);

    auto tg_xxxx_y = pbuffer.data(idx_gp + 1);

    auto tg_xxxx_z = pbuffer.data(idx_gp + 2);

    auto tg_xxyy_y = pbuffer.data(idx_gp + 10);

    auto tg_xxzz_x = pbuffer.data(idx_gp + 15);

    auto tg_xxzz_z = pbuffer.data(idx_gp + 17);

    auto tg_xyyy_y = pbuffer.data(idx_gp + 19);

    auto tg_xzzz_z = pbuffer.data(idx_gp + 29);

    auto tg_yyyy_x = pbuffer.data(idx_gp + 30);

    auto tg_yyyy_y = pbuffer.data(idx_gp + 31);

    auto tg_yyyy_z = pbuffer.data(idx_gp + 32);

    auto tg_yyyz_z = pbuffer.data(idx_gp + 35);

    auto tg_yyzz_x = pbuffer.data(idx_gp + 36);

    auto tg_yyzz_y = pbuffer.data(idx_gp + 37);

    auto tg_yyzz_z = pbuffer.data(idx_gp + 38);

    auto tg_yzzz_y = pbuffer.data(idx_gp + 40);

    auto tg_yzzz_z = pbuffer.data(idx_gp + 41);

    auto tg_zzzz_x = pbuffer.data(idx_gp + 42);

    auto tg_zzzz_y = pbuffer.data(idx_gp + 43);

    auto tg_zzzz_z = pbuffer.data(idx_gp + 44);

    // Set up components of auxiliary buffer : GD

    auto tg_xxxx_xx = pbuffer.data(idx_gd);

    auto tg_xxxx_xy = pbuffer.data(idx_gd + 1);

    auto tg_xxxx_xz = pbuffer.data(idx_gd + 2);

    auto tg_xxxx_yy = pbuffer.data(idx_gd + 3);

    auto tg_xxxx_yz = pbuffer.data(idx_gd + 4);

    auto tg_xxxx_zz = pbuffer.data(idx_gd + 5);

    auto tg_xxxy_xx = pbuffer.data(idx_gd + 6);

    auto tg_xxxy_xy = pbuffer.data(idx_gd + 7);

    auto tg_xxxy_xz = pbuffer.data(idx_gd + 8);

    auto tg_xxxy_yy = pbuffer.data(idx_gd + 9);

    auto tg_xxxy_yz = pbuffer.data(idx_gd + 10);

    auto tg_xxxz_xx = pbuffer.data(idx_gd + 12);

    auto tg_xxxz_xy = pbuffer.data(idx_gd + 13);

    auto tg_xxxz_xz = pbuffer.data(idx_gd + 14);

    auto tg_xxxz_yz = pbuffer.data(idx_gd + 16);

    auto tg_xxxz_zz = pbuffer.data(idx_gd + 17);

    auto tg_xxyy_xx = pbuffer.data(idx_gd + 18);

    auto tg_xxyy_xy = pbuffer.data(idx_gd + 19);

    auto tg_xxyy_xz = pbuffer.data(idx_gd + 20);

    auto tg_xxyy_yy = pbuffer.data(idx_gd + 21);

    auto tg_xxyy_yz = pbuffer.data(idx_gd + 22);

    auto tg_xxyy_zz = pbuffer.data(idx_gd + 23);

    auto tg_xxyz_xz = pbuffer.data(idx_gd + 26);

    auto tg_xxyz_yz = pbuffer.data(idx_gd + 28);

    auto tg_xxzz_xx = pbuffer.data(idx_gd + 30);

    auto tg_xxzz_xy = pbuffer.data(idx_gd + 31);

    auto tg_xxzz_xz = pbuffer.data(idx_gd + 32);

    auto tg_xxzz_yy = pbuffer.data(idx_gd + 33);

    auto tg_xxzz_yz = pbuffer.data(idx_gd + 34);

    auto tg_xxzz_zz = pbuffer.data(idx_gd + 35);

    auto tg_xyyy_xx = pbuffer.data(idx_gd + 36);

    auto tg_xyyy_xy = pbuffer.data(idx_gd + 37);

    auto tg_xyyy_yy = pbuffer.data(idx_gd + 39);

    auto tg_xyyy_yz = pbuffer.data(idx_gd + 40);

    auto tg_xyyy_zz = pbuffer.data(idx_gd + 41);

    auto tg_xyyz_yz = pbuffer.data(idx_gd + 46);

    auto tg_xyyz_zz = pbuffer.data(idx_gd + 47);

    auto tg_xyzz_yy = pbuffer.data(idx_gd + 51);

    auto tg_xyzz_yz = pbuffer.data(idx_gd + 52);

    auto tg_xzzz_xx = pbuffer.data(idx_gd + 54);

    auto tg_xzzz_xz = pbuffer.data(idx_gd + 56);

    auto tg_xzzz_yy = pbuffer.data(idx_gd + 57);

    auto tg_xzzz_yz = pbuffer.data(idx_gd + 58);

    auto tg_xzzz_zz = pbuffer.data(idx_gd + 59);

    auto tg_yyyy_xx = pbuffer.data(idx_gd + 60);

    auto tg_yyyy_xy = pbuffer.data(idx_gd + 61);

    auto tg_yyyy_xz = pbuffer.data(idx_gd + 62);

    auto tg_yyyy_yy = pbuffer.data(idx_gd + 63);

    auto tg_yyyy_yz = pbuffer.data(idx_gd + 64);

    auto tg_yyyy_zz = pbuffer.data(idx_gd + 65);

    auto tg_yyyz_xy = pbuffer.data(idx_gd + 67);

    auto tg_yyyz_xz = pbuffer.data(idx_gd + 68);

    auto tg_yyyz_yy = pbuffer.data(idx_gd + 69);

    auto tg_yyyz_yz = pbuffer.data(idx_gd + 70);

    auto tg_yyyz_zz = pbuffer.data(idx_gd + 71);

    auto tg_yyzz_xx = pbuffer.data(idx_gd + 72);

    auto tg_yyzz_xy = pbuffer.data(idx_gd + 73);

    auto tg_yyzz_xz = pbuffer.data(idx_gd + 74);

    auto tg_yyzz_yy = pbuffer.data(idx_gd + 75);

    auto tg_yyzz_yz = pbuffer.data(idx_gd + 76);

    auto tg_yyzz_zz = pbuffer.data(idx_gd + 77);

    auto tg_yzzz_xx = pbuffer.data(idx_gd + 78);

    auto tg_yzzz_xy = pbuffer.data(idx_gd + 79);

    auto tg_yzzz_xz = pbuffer.data(idx_gd + 80);

    auto tg_yzzz_yy = pbuffer.data(idx_gd + 81);

    auto tg_yzzz_yz = pbuffer.data(idx_gd + 82);

    auto tg_yzzz_zz = pbuffer.data(idx_gd + 83);

    auto tg_zzzz_xx = pbuffer.data(idx_gd + 84);

    auto tg_zzzz_xy = pbuffer.data(idx_gd + 85);

    auto tg_zzzz_xz = pbuffer.data(idx_gd + 86);

    auto tg_zzzz_yy = pbuffer.data(idx_gd + 87);

    auto tg_zzzz_yz = pbuffer.data(idx_gd + 88);

    auto tg_zzzz_zz = pbuffer.data(idx_gd + 89);

    // Set up components of targeted buffer : HD

    auto tg_xxxxx_xx = pbuffer.data(idx_hd);

    auto tg_xxxxx_xy = pbuffer.data(idx_hd + 1);

    auto tg_xxxxx_xz = pbuffer.data(idx_hd + 2);

    auto tg_xxxxx_yy = pbuffer.data(idx_hd + 3);

    auto tg_xxxxx_yz = pbuffer.data(idx_hd + 4);

    auto tg_xxxxx_zz = pbuffer.data(idx_hd + 5);

    auto tg_xxxxy_xx = pbuffer.data(idx_hd + 6);

    auto tg_xxxxy_xy = pbuffer.data(idx_hd + 7);

    auto tg_xxxxy_xz = pbuffer.data(idx_hd + 8);

    auto tg_xxxxy_yy = pbuffer.data(idx_hd + 9);

    auto tg_xxxxy_yz = pbuffer.data(idx_hd + 10);

    auto tg_xxxxy_zz = pbuffer.data(idx_hd + 11);

    auto tg_xxxxz_xx = pbuffer.data(idx_hd + 12);

    auto tg_xxxxz_xy = pbuffer.data(idx_hd + 13);

    auto tg_xxxxz_xz = pbuffer.data(idx_hd + 14);

    auto tg_xxxxz_yy = pbuffer.data(idx_hd + 15);

    auto tg_xxxxz_yz = pbuffer.data(idx_hd + 16);

    auto tg_xxxxz_zz = pbuffer.data(idx_hd + 17);

    auto tg_xxxyy_xx = pbuffer.data(idx_hd + 18);

    auto tg_xxxyy_xy = pbuffer.data(idx_hd + 19);

    auto tg_xxxyy_xz = pbuffer.data(idx_hd + 20);

    auto tg_xxxyy_yy = pbuffer.data(idx_hd + 21);

    auto tg_xxxyy_yz = pbuffer.data(idx_hd + 22);

    auto tg_xxxyy_zz = pbuffer.data(idx_hd + 23);

    auto tg_xxxyz_xx = pbuffer.data(idx_hd + 24);

    auto tg_xxxyz_xy = pbuffer.data(idx_hd + 25);

    auto tg_xxxyz_xz = pbuffer.data(idx_hd + 26);

    auto tg_xxxyz_yy = pbuffer.data(idx_hd + 27);

    auto tg_xxxyz_yz = pbuffer.data(idx_hd + 28);

    auto tg_xxxyz_zz = pbuffer.data(idx_hd + 29);

    auto tg_xxxzz_xx = pbuffer.data(idx_hd + 30);

    auto tg_xxxzz_xy = pbuffer.data(idx_hd + 31);

    auto tg_xxxzz_xz = pbuffer.data(idx_hd + 32);

    auto tg_xxxzz_yy = pbuffer.data(idx_hd + 33);

    auto tg_xxxzz_yz = pbuffer.data(idx_hd + 34);

    auto tg_xxxzz_zz = pbuffer.data(idx_hd + 35);

    auto tg_xxyyy_xx = pbuffer.data(idx_hd + 36);

    auto tg_xxyyy_xy = pbuffer.data(idx_hd + 37);

    auto tg_xxyyy_xz = pbuffer.data(idx_hd + 38);

    auto tg_xxyyy_yy = pbuffer.data(idx_hd + 39);

    auto tg_xxyyy_yz = pbuffer.data(idx_hd + 40);

    auto tg_xxyyy_zz = pbuffer.data(idx_hd + 41);

    auto tg_xxyyz_xx = pbuffer.data(idx_hd + 42);

    auto tg_xxyyz_xy = pbuffer.data(idx_hd + 43);

    auto tg_xxyyz_xz = pbuffer.data(idx_hd + 44);

    auto tg_xxyyz_yy = pbuffer.data(idx_hd + 45);

    auto tg_xxyyz_yz = pbuffer.data(idx_hd + 46);

    auto tg_xxyyz_zz = pbuffer.data(idx_hd + 47);

    auto tg_xxyzz_xx = pbuffer.data(idx_hd + 48);

    auto tg_xxyzz_xy = pbuffer.data(idx_hd + 49);

    auto tg_xxyzz_xz = pbuffer.data(idx_hd + 50);

    auto tg_xxyzz_yy = pbuffer.data(idx_hd + 51);

    auto tg_xxyzz_yz = pbuffer.data(idx_hd + 52);

    auto tg_xxyzz_zz = pbuffer.data(idx_hd + 53);

    auto tg_xxzzz_xx = pbuffer.data(idx_hd + 54);

    auto tg_xxzzz_xy = pbuffer.data(idx_hd + 55);

    auto tg_xxzzz_xz = pbuffer.data(idx_hd + 56);

    auto tg_xxzzz_yy = pbuffer.data(idx_hd + 57);

    auto tg_xxzzz_yz = pbuffer.data(idx_hd + 58);

    auto tg_xxzzz_zz = pbuffer.data(idx_hd + 59);

    auto tg_xyyyy_xx = pbuffer.data(idx_hd + 60);

    auto tg_xyyyy_xy = pbuffer.data(idx_hd + 61);

    auto tg_xyyyy_xz = pbuffer.data(idx_hd + 62);

    auto tg_xyyyy_yy = pbuffer.data(idx_hd + 63);

    auto tg_xyyyy_yz = pbuffer.data(idx_hd + 64);

    auto tg_xyyyy_zz = pbuffer.data(idx_hd + 65);

    auto tg_xyyyz_xx = pbuffer.data(idx_hd + 66);

    auto tg_xyyyz_xy = pbuffer.data(idx_hd + 67);

    auto tg_xyyyz_xz = pbuffer.data(idx_hd + 68);

    auto tg_xyyyz_yy = pbuffer.data(idx_hd + 69);

    auto tg_xyyyz_yz = pbuffer.data(idx_hd + 70);

    auto tg_xyyyz_zz = pbuffer.data(idx_hd + 71);

    auto tg_xyyzz_xx = pbuffer.data(idx_hd + 72);

    auto tg_xyyzz_xy = pbuffer.data(idx_hd + 73);

    auto tg_xyyzz_xz = pbuffer.data(idx_hd + 74);

    auto tg_xyyzz_yy = pbuffer.data(idx_hd + 75);

    auto tg_xyyzz_yz = pbuffer.data(idx_hd + 76);

    auto tg_xyyzz_zz = pbuffer.data(idx_hd + 77);

    auto tg_xyzzz_xx = pbuffer.data(idx_hd + 78);

    auto tg_xyzzz_xy = pbuffer.data(idx_hd + 79);

    auto tg_xyzzz_xz = pbuffer.data(idx_hd + 80);

    auto tg_xyzzz_yy = pbuffer.data(idx_hd + 81);

    auto tg_xyzzz_yz = pbuffer.data(idx_hd + 82);

    auto tg_xyzzz_zz = pbuffer.data(idx_hd + 83);

    auto tg_xzzzz_xx = pbuffer.data(idx_hd + 84);

    auto tg_xzzzz_xy = pbuffer.data(idx_hd + 85);

    auto tg_xzzzz_xz = pbuffer.data(idx_hd + 86);

    auto tg_xzzzz_yy = pbuffer.data(idx_hd + 87);

    auto tg_xzzzz_yz = pbuffer.data(idx_hd + 88);

    auto tg_xzzzz_zz = pbuffer.data(idx_hd + 89);

    auto tg_yyyyy_xx = pbuffer.data(idx_hd + 90);

    auto tg_yyyyy_xy = pbuffer.data(idx_hd + 91);

    auto tg_yyyyy_xz = pbuffer.data(idx_hd + 92);

    auto tg_yyyyy_yy = pbuffer.data(idx_hd + 93);

    auto tg_yyyyy_yz = pbuffer.data(idx_hd + 94);

    auto tg_yyyyy_zz = pbuffer.data(idx_hd + 95);

    auto tg_yyyyz_xx = pbuffer.data(idx_hd + 96);

    auto tg_yyyyz_xy = pbuffer.data(idx_hd + 97);

    auto tg_yyyyz_xz = pbuffer.data(idx_hd + 98);

    auto tg_yyyyz_yy = pbuffer.data(idx_hd + 99);

    auto tg_yyyyz_yz = pbuffer.data(idx_hd + 100);

    auto tg_yyyyz_zz = pbuffer.data(idx_hd + 101);

    auto tg_yyyzz_xx = pbuffer.data(idx_hd + 102);

    auto tg_yyyzz_xy = pbuffer.data(idx_hd + 103);

    auto tg_yyyzz_xz = pbuffer.data(idx_hd + 104);

    auto tg_yyyzz_yy = pbuffer.data(idx_hd + 105);

    auto tg_yyyzz_yz = pbuffer.data(idx_hd + 106);

    auto tg_yyyzz_zz = pbuffer.data(idx_hd + 107);

    auto tg_yyzzz_xx = pbuffer.data(idx_hd + 108);

    auto tg_yyzzz_xy = pbuffer.data(idx_hd + 109);

    auto tg_yyzzz_xz = pbuffer.data(idx_hd + 110);

    auto tg_yyzzz_yy = pbuffer.data(idx_hd + 111);

    auto tg_yyzzz_yz = pbuffer.data(idx_hd + 112);

    auto tg_yyzzz_zz = pbuffer.data(idx_hd + 113);

    auto tg_yzzzz_xx = pbuffer.data(idx_hd + 114);

    auto tg_yzzzz_xy = pbuffer.data(idx_hd + 115);

    auto tg_yzzzz_xz = pbuffer.data(idx_hd + 116);

    auto tg_yzzzz_yy = pbuffer.data(idx_hd + 117);

    auto tg_yzzzz_yz = pbuffer.data(idx_hd + 118);

    auto tg_yzzzz_zz = pbuffer.data(idx_hd + 119);

    auto tg_zzzzz_xx = pbuffer.data(idx_hd + 120);

    auto tg_zzzzz_xy = pbuffer.data(idx_hd + 121);

    auto tg_zzzzz_xz = pbuffer.data(idx_hd + 122);

    auto tg_zzzzz_yy = pbuffer.data(idx_hd + 123);

    auto tg_zzzzz_yz = pbuffer.data(idx_hd + 124);

    auto tg_zzzzz_zz = pbuffer.data(idx_hd + 125);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_xxx_xx, tg_xxx_xy, tg_xxx_xz, tg_xxx_yy, tg_xxx_yz, tg_xxx_zz, tg_xxxx_x, tg_xxxx_xx, tg_xxxx_xy, tg_xxxx_xz, tg_xxxx_y, tg_xxxx_yy, tg_xxxx_yz, tg_xxxx_z, tg_xxxx_zz, tg_xxxxx_xx, tg_xxxxx_xy, tg_xxxxx_xz, tg_xxxxx_yy, tg_xxxxx_yz, tg_xxxxx_zz, tg_xxxxy_xx, tg_xxxxy_xy, tg_xxxxy_xz, tg_xxxxy_yy, tg_xxxxy_yz, tg_xxxxy_zz, tg_xxxxz_xx, tg_xxxxz_xy, tg_xxxxz_xz, tg_xxxxz_yy, tg_xxxxz_yz, tg_xxxxz_zz, tg_xxxy_xx, tg_xxxy_xy, tg_xxxy_xz, tg_xxxy_yy, tg_xxxy_yz, tg_xxxyy_xx, tg_xxxyy_xy, tg_xxxyy_xz, tg_xxxyy_yy, tg_xxxyy_yz, tg_xxxyy_zz, tg_xxxyz_xx, tg_xxxyz_xy, tg_xxxyz_xz, tg_xxxyz_yy, tg_xxxyz_yz, tg_xxxyz_zz, tg_xxxz_xx, tg_xxxz_xy, tg_xxxz_xz, tg_xxxz_yz, tg_xxxz_zz, tg_xxxzz_xx, tg_xxxzz_xy, tg_xxxzz_xz, tg_xxxzz_yy, tg_xxxzz_yz, tg_xxxzz_zz, tg_xxy_xx, tg_xxy_xz, tg_xxy_yy, tg_xxy_yz, tg_xxyy_xx, tg_xxyy_xy, tg_xxyy_xz, tg_xxyy_y, tg_xxyy_yy, tg_xxyy_yz, tg_xxyy_zz, tg_xxyyy_xx, tg_xxyyy_xy, tg_xxyyy_xz, tg_xxyyy_yy, tg_xxyyy_yz, tg_xxyyy_zz, tg_xxyyz_xx, tg_xxyyz_xy, tg_xxyyz_xz, tg_xxyyz_yy, tg_xxyyz_yz, tg_xxyyz_zz, tg_xxyz_xz, tg_xxyz_yz, tg_xxyzz_xx, tg_xxyzz_xy, tg_xxyzz_xz, tg_xxyzz_yy, tg_xxyzz_yz, tg_xxyzz_zz, tg_xxz_xx, tg_xxz_xy, tg_xxz_xz, tg_xxz_yz, tg_xxz_zz, tg_xxzz_x, tg_xxzz_xx, tg_xxzz_xy, tg_xxzz_xz, tg_xxzz_yy, tg_xxzz_yz, tg_xxzz_z, tg_xxzz_zz, tg_xxzzz_xx, tg_xxzzz_xy, tg_xxzzz_xz, tg_xxzzz_yy, tg_xxzzz_yz, tg_xxzzz_zz, tg_xyy_xy, tg_xyy_yy, tg_xyy_yz, tg_xyy_zz, tg_xyyy_xx, tg_xyyy_xy, tg_xyyy_y, tg_xyyy_yy, tg_xyyy_yz, tg_xyyy_zz, tg_xyyyy_xx, tg_xyyyy_xy, tg_xyyyy_xz, tg_xyyyy_yy, tg_xyyyy_yz, tg_xyyyy_zz, tg_xyyyz_xx, tg_xyyyz_xy, tg_xyyyz_xz, tg_xyyyz_yy, tg_xyyyz_yz, tg_xyyyz_zz, tg_xyyz_yz, tg_xyyz_zz, tg_xyyzz_xx, tg_xyyzz_xy, tg_xyyzz_xz, tg_xyyzz_yy, tg_xyyzz_yz, tg_xyyzz_zz, tg_xyz_yz, tg_xyzz_yy, tg_xyzz_yz, tg_xyzzz_xx, tg_xyzzz_xy, tg_xyzzz_xz, tg_xyzzz_yy, tg_xyzzz_yz, tg_xyzzz_zz, tg_xzz_xz, tg_xzz_yy, tg_xzz_yz, tg_xzz_zz, tg_xzzz_xx, tg_xzzz_xz, tg_xzzz_yy, tg_xzzz_yz, tg_xzzz_z, tg_xzzz_zz, tg_xzzzz_xx, tg_xzzzz_xy, tg_xzzzz_xz, tg_xzzzz_yy, tg_xzzzz_yz, tg_xzzzz_zz, tg_yyy_xx, tg_yyy_xy, tg_yyy_xz, tg_yyy_yy, tg_yyy_yz, tg_yyy_zz, tg_yyyy_x, tg_yyyy_xx, tg_yyyy_xy, tg_yyyy_xz, tg_yyyy_y, tg_yyyy_yy, tg_yyyy_yz, tg_yyyy_z, tg_yyyy_zz, tg_yyyyy_xx, tg_yyyyy_xy, tg_yyyyy_xz, tg_yyyyy_yy, tg_yyyyy_yz, tg_yyyyy_zz, tg_yyyyz_xx, tg_yyyyz_xy, tg_yyyyz_xz, tg_yyyyz_yy, tg_yyyyz_yz, tg_yyyyz_zz, tg_yyyz_xy, tg_yyyz_xz, tg_yyyz_yy, tg_yyyz_yz, tg_yyyz_z, tg_yyyz_zz, tg_yyyzz_xx, tg_yyyzz_xy, tg_yyyzz_xz, tg_yyyzz_yy, tg_yyyzz_yz, tg_yyyzz_zz, tg_yyz_xy, tg_yyz_xz, tg_yyz_yy, tg_yyz_yz, tg_yyz_zz, tg_yyzz_x, tg_yyzz_xx, tg_yyzz_xy, tg_yyzz_xz, tg_yyzz_y, tg_yyzz_yy, tg_yyzz_yz, tg_yyzz_z, tg_yyzz_zz, tg_yyzzz_xx, tg_yyzzz_xy, tg_yyzzz_xz, tg_yyzzz_yy, tg_yyzzz_yz, tg_yyzzz_zz, tg_yzz_xx, tg_yzz_xz, tg_yzz_yy, tg_yzz_yz, tg_yzz_zz, tg_yzzz_xx, tg_yzzz_xy, tg_yzzz_xz, tg_yzzz_y, tg_yzzz_yy, tg_yzzz_yz, tg_yzzz_z, tg_yzzz_zz, tg_yzzzz_xx, tg_yzzzz_xy, tg_yzzzz_xz, tg_yzzzz_yy, tg_yzzzz_yz, tg_yzzzz_zz, tg_zzz_xx, tg_zzz_xy, tg_zzz_xz, tg_zzz_yy, tg_zzz_yz, tg_zzz_zz, tg_zzzz_x, tg_zzzz_xx, tg_zzzz_xy, tg_zzzz_xz, tg_zzzz_y, tg_zzzz_yy, tg_zzzz_yz, tg_zzzz_z, tg_zzzz_zz, tg_zzzzz_xx, tg_zzzzz_xy, tg_zzzzz_xz, tg_zzzzz_yy, tg_zzzzz_yz, tg_zzzzz_zz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xxxxx_xx[i] = 4.0 * tg_xxx_xx[i] * fxi[i] + 2.0 * tg_xxxx_x[i] * fxi[i] + tg_xxxx_xx[i] * ra_x[i];

        tg_xxxxx_xy[i] = 4.0 * tg_xxx_xy[i] * fxi[i] + tg_xxxx_y[i] * fxi[i] + tg_xxxx_xy[i] * ra_x[i];

        tg_xxxxx_xz[i] = 4.0 * tg_xxx_xz[i] * fxi[i] + tg_xxxx_z[i] * fxi[i] + tg_xxxx_xz[i] * ra_x[i];

        tg_xxxxx_yy[i] = 4.0 * tg_xxx_yy[i] * fxi[i] + tg_xxxx_yy[i] * ra_x[i];

        tg_xxxxx_yz[i] = 4.0 * tg_xxx_yz[i] * fxi[i] + tg_xxxx_yz[i] * ra_x[i];

        tg_xxxxx_zz[i] = 4.0 * tg_xxx_zz[i] * fxi[i] + tg_xxxx_zz[i] * ra_x[i];

        tg_xxxxy_xx[i] = tg_xxxx_xx[i] * ra_y[i];

        tg_xxxxy_xy[i] = tg_xxxx_x[i] * fxi[i] + tg_xxxx_xy[i] * ra_y[i];

        tg_xxxxy_xz[i] = tg_xxxx_xz[i] * ra_y[i];

        tg_xxxxy_yy[i] = 3.0 * tg_xxy_yy[i] * fxi[i] + tg_xxxy_yy[i] * ra_x[i];

        tg_xxxxy_yz[i] = 3.0 * tg_xxy_yz[i] * fxi[i] + tg_xxxy_yz[i] * ra_x[i];

        tg_xxxxy_zz[i] = tg_xxxx_zz[i] * ra_y[i];

        tg_xxxxz_xx[i] = tg_xxxx_xx[i] * ra_z[i];

        tg_xxxxz_xy[i] = tg_xxxx_xy[i] * ra_z[i];

        tg_xxxxz_xz[i] = tg_xxxx_x[i] * fxi[i] + tg_xxxx_xz[i] * ra_z[i];

        tg_xxxxz_yy[i] = tg_xxxx_yy[i] * ra_z[i];

        tg_xxxxz_yz[i] = 3.0 * tg_xxz_yz[i] * fxi[i] + tg_xxxz_yz[i] * ra_x[i];

        tg_xxxxz_zz[i] = 3.0 * tg_xxz_zz[i] * fxi[i] + tg_xxxz_zz[i] * ra_x[i];

        tg_xxxyy_xx[i] = tg_xxx_xx[i] * fxi[i] + tg_xxxy_xx[i] * ra_y[i];

        tg_xxxyy_xy[i] = 2.0 * tg_xyy_xy[i] * fxi[i] + tg_xxyy_y[i] * fxi[i] + tg_xxyy_xy[i] * ra_x[i];

        tg_xxxyy_xz[i] = tg_xxx_xz[i] * fxi[i] + tg_xxxy_xz[i] * ra_y[i];

        tg_xxxyy_yy[i] = 2.0 * tg_xyy_yy[i] * fxi[i] + tg_xxyy_yy[i] * ra_x[i];

        tg_xxxyy_yz[i] = 2.0 * tg_xyy_yz[i] * fxi[i] + tg_xxyy_yz[i] * ra_x[i];

        tg_xxxyy_zz[i] = 2.0 * tg_xyy_zz[i] * fxi[i] + tg_xxyy_zz[i] * ra_x[i];

        tg_xxxyz_xx[i] = tg_xxxz_xx[i] * ra_y[i];

        tg_xxxyz_xy[i] = tg_xxxy_xy[i] * ra_z[i];

        tg_xxxyz_xz[i] = tg_xxxz_xz[i] * ra_y[i];

        tg_xxxyz_yy[i] = tg_xxxy_yy[i] * ra_z[i];

        tg_xxxyz_yz[i] = 2.0 * tg_xyz_yz[i] * fxi[i] + tg_xxyz_yz[i] * ra_x[i];

        tg_xxxyz_zz[i] = tg_xxxz_zz[i] * ra_y[i];

        tg_xxxzz_xx[i] = tg_xxx_xx[i] * fxi[i] + tg_xxxz_xx[i] * ra_z[i];

        tg_xxxzz_xy[i] = tg_xxx_xy[i] * fxi[i] + tg_xxxz_xy[i] * ra_z[i];

        tg_xxxzz_xz[i] = 2.0 * tg_xzz_xz[i] * fxi[i] + tg_xxzz_z[i] * fxi[i] + tg_xxzz_xz[i] * ra_x[i];

        tg_xxxzz_yy[i] = 2.0 * tg_xzz_yy[i] * fxi[i] + tg_xxzz_yy[i] * ra_x[i];

        tg_xxxzz_yz[i] = 2.0 * tg_xzz_yz[i] * fxi[i] + tg_xxzz_yz[i] * ra_x[i];

        tg_xxxzz_zz[i] = 2.0 * tg_xzz_zz[i] * fxi[i] + tg_xxzz_zz[i] * ra_x[i];

        tg_xxyyy_xx[i] = 2.0 * tg_xxy_xx[i] * fxi[i] + tg_xxyy_xx[i] * ra_y[i];

        tg_xxyyy_xy[i] = tg_yyy_xy[i] * fxi[i] + tg_xyyy_y[i] * fxi[i] + tg_xyyy_xy[i] * ra_x[i];

        tg_xxyyy_xz[i] = 2.0 * tg_xxy_xz[i] * fxi[i] + tg_xxyy_xz[i] * ra_y[i];

        tg_xxyyy_yy[i] = tg_yyy_yy[i] * fxi[i] + tg_xyyy_yy[i] * ra_x[i];

        tg_xxyyy_yz[i] = tg_yyy_yz[i] * fxi[i] + tg_xyyy_yz[i] * ra_x[i];

        tg_xxyyy_zz[i] = tg_yyy_zz[i] * fxi[i] + tg_xyyy_zz[i] * ra_x[i];

        tg_xxyyz_xx[i] = tg_xxyy_xx[i] * ra_z[i];

        tg_xxyyz_xy[i] = tg_xxyy_xy[i] * ra_z[i];

        tg_xxyyz_xz[i] = tg_xxz_xz[i] * fxi[i] + tg_xxyz_xz[i] * ra_y[i];

        tg_xxyyz_yy[i] = tg_xxyy_yy[i] * ra_z[i];

        tg_xxyyz_yz[i] = tg_yyz_yz[i] * fxi[i] + tg_xyyz_yz[i] * ra_x[i];

        tg_xxyyz_zz[i] = tg_yyz_zz[i] * fxi[i] + tg_xyyz_zz[i] * ra_x[i];

        tg_xxyzz_xx[i] = tg_xxzz_xx[i] * ra_y[i];

        tg_xxyzz_xy[i] = tg_xxzz_x[i] * fxi[i] + tg_xxzz_xy[i] * ra_y[i];

        tg_xxyzz_xz[i] = tg_xxzz_xz[i] * ra_y[i];

        tg_xxyzz_yy[i] = tg_yzz_yy[i] * fxi[i] + tg_xyzz_yy[i] * ra_x[i];

        tg_xxyzz_yz[i] = tg_yzz_yz[i] * fxi[i] + tg_xyzz_yz[i] * ra_x[i];

        tg_xxyzz_zz[i] = tg_xxzz_zz[i] * ra_y[i];

        tg_xxzzz_xx[i] = 2.0 * tg_xxz_xx[i] * fxi[i] + tg_xxzz_xx[i] * ra_z[i];

        tg_xxzzz_xy[i] = 2.0 * tg_xxz_xy[i] * fxi[i] + tg_xxzz_xy[i] * ra_z[i];

        tg_xxzzz_xz[i] = tg_zzz_xz[i] * fxi[i] + tg_xzzz_z[i] * fxi[i] + tg_xzzz_xz[i] * ra_x[i];

        tg_xxzzz_yy[i] = tg_zzz_yy[i] * fxi[i] + tg_xzzz_yy[i] * ra_x[i];

        tg_xxzzz_yz[i] = tg_zzz_yz[i] * fxi[i] + tg_xzzz_yz[i] * ra_x[i];

        tg_xxzzz_zz[i] = tg_zzz_zz[i] * fxi[i] + tg_xzzz_zz[i] * ra_x[i];

        tg_xyyyy_xx[i] = 2.0 * tg_yyyy_x[i] * fxi[i] + tg_yyyy_xx[i] * ra_x[i];

        tg_xyyyy_xy[i] = tg_yyyy_y[i] * fxi[i] + tg_yyyy_xy[i] * ra_x[i];

        tg_xyyyy_xz[i] = tg_yyyy_z[i] * fxi[i] + tg_yyyy_xz[i] * ra_x[i];

        tg_xyyyy_yy[i] = tg_yyyy_yy[i] * ra_x[i];

        tg_xyyyy_yz[i] = tg_yyyy_yz[i] * ra_x[i];

        tg_xyyyy_zz[i] = tg_yyyy_zz[i] * ra_x[i];

        tg_xyyyz_xx[i] = tg_xyyy_xx[i] * ra_z[i];

        tg_xyyyz_xy[i] = tg_xyyy_xy[i] * ra_z[i];

        tg_xyyyz_xz[i] = tg_yyyz_z[i] * fxi[i] + tg_yyyz_xz[i] * ra_x[i];

        tg_xyyyz_yy[i] = tg_yyyz_yy[i] * ra_x[i];

        tg_xyyyz_yz[i] = tg_yyyz_yz[i] * ra_x[i];

        tg_xyyyz_zz[i] = tg_yyyz_zz[i] * ra_x[i];

        tg_xyyzz_xx[i] = 2.0 * tg_yyzz_x[i] * fxi[i] + tg_yyzz_xx[i] * ra_x[i];

        tg_xyyzz_xy[i] = tg_yyzz_y[i] * fxi[i] + tg_yyzz_xy[i] * ra_x[i];

        tg_xyyzz_xz[i] = tg_yyzz_z[i] * fxi[i] + tg_yyzz_xz[i] * ra_x[i];

        tg_xyyzz_yy[i] = tg_yyzz_yy[i] * ra_x[i];

        tg_xyyzz_yz[i] = tg_yyzz_yz[i] * ra_x[i];

        tg_xyyzz_zz[i] = tg_yyzz_zz[i] * ra_x[i];

        tg_xyzzz_xx[i] = tg_xzzz_xx[i] * ra_y[i];

        tg_xyzzz_xy[i] = tg_yzzz_y[i] * fxi[i] + tg_yzzz_xy[i] * ra_x[i];

        tg_xyzzz_xz[i] = tg_xzzz_xz[i] * ra_y[i];

        tg_xyzzz_yy[i] = tg_yzzz_yy[i] * ra_x[i];

        tg_xyzzz_yz[i] = tg_yzzz_yz[i] * ra_x[i];

        tg_xyzzz_zz[i] = tg_yzzz_zz[i] * ra_x[i];

        tg_xzzzz_xx[i] = 2.0 * tg_zzzz_x[i] * fxi[i] + tg_zzzz_xx[i] * ra_x[i];

        tg_xzzzz_xy[i] = tg_zzzz_y[i] * fxi[i] + tg_zzzz_xy[i] * ra_x[i];

        tg_xzzzz_xz[i] = tg_zzzz_z[i] * fxi[i] + tg_zzzz_xz[i] * ra_x[i];

        tg_xzzzz_yy[i] = tg_zzzz_yy[i] * ra_x[i];

        tg_xzzzz_yz[i] = tg_zzzz_yz[i] * ra_x[i];

        tg_xzzzz_zz[i] = tg_zzzz_zz[i] * ra_x[i];

        tg_yyyyy_xx[i] = 4.0 * tg_yyy_xx[i] * fxi[i] + tg_yyyy_xx[i] * ra_y[i];

        tg_yyyyy_xy[i] = 4.0 * tg_yyy_xy[i] * fxi[i] + tg_yyyy_x[i] * fxi[i] + tg_yyyy_xy[i] * ra_y[i];

        tg_yyyyy_xz[i] = 4.0 * tg_yyy_xz[i] * fxi[i] + tg_yyyy_xz[i] * ra_y[i];

        tg_yyyyy_yy[i] = 4.0 * tg_yyy_yy[i] * fxi[i] + 2.0 * tg_yyyy_y[i] * fxi[i] + tg_yyyy_yy[i] * ra_y[i];

        tg_yyyyy_yz[i] = 4.0 * tg_yyy_yz[i] * fxi[i] + tg_yyyy_z[i] * fxi[i] + tg_yyyy_yz[i] * ra_y[i];

        tg_yyyyy_zz[i] = 4.0 * tg_yyy_zz[i] * fxi[i] + tg_yyyy_zz[i] * ra_y[i];

        tg_yyyyz_xx[i] = tg_yyyy_xx[i] * ra_z[i];

        tg_yyyyz_xy[i] = tg_yyyy_xy[i] * ra_z[i];

        tg_yyyyz_xz[i] = 3.0 * tg_yyz_xz[i] * fxi[i] + tg_yyyz_xz[i] * ra_y[i];

        tg_yyyyz_yy[i] = tg_yyyy_yy[i] * ra_z[i];

        tg_yyyyz_yz[i] = tg_yyyy_y[i] * fxi[i] + tg_yyyy_yz[i] * ra_z[i];

        tg_yyyyz_zz[i] = 3.0 * tg_yyz_zz[i] * fxi[i] + tg_yyyz_zz[i] * ra_y[i];

        tg_yyyzz_xx[i] = 2.0 * tg_yzz_xx[i] * fxi[i] + tg_yyzz_xx[i] * ra_y[i];

        tg_yyyzz_xy[i] = tg_yyy_xy[i] * fxi[i] + tg_yyyz_xy[i] * ra_z[i];

        tg_yyyzz_xz[i] = 2.0 * tg_yzz_xz[i] * fxi[i] + tg_yyzz_xz[i] * ra_y[i];

        tg_yyyzz_yy[i] = tg_yyy_yy[i] * fxi[i] + tg_yyyz_yy[i] * ra_z[i];

        tg_yyyzz_yz[i] = 2.0 * tg_yzz_yz[i] * fxi[i] + tg_yyzz_z[i] * fxi[i] + tg_yyzz_yz[i] * ra_y[i];

        tg_yyyzz_zz[i] = 2.0 * tg_yzz_zz[i] * fxi[i] + tg_yyzz_zz[i] * ra_y[i];

        tg_yyzzz_xx[i] = tg_zzz_xx[i] * fxi[i] + tg_yzzz_xx[i] * ra_y[i];

        tg_yyzzz_xy[i] = 2.0 * tg_yyz_xy[i] * fxi[i] + tg_yyzz_xy[i] * ra_z[i];

        tg_yyzzz_xz[i] = tg_zzz_xz[i] * fxi[i] + tg_yzzz_xz[i] * ra_y[i];

        tg_yyzzz_yy[i] = 2.0 * tg_yyz_yy[i] * fxi[i] + tg_yyzz_yy[i] * ra_z[i];

        tg_yyzzz_yz[i] = tg_zzz_yz[i] * fxi[i] + tg_yzzz_z[i] * fxi[i] + tg_yzzz_yz[i] * ra_y[i];

        tg_yyzzz_zz[i] = tg_zzz_zz[i] * fxi[i] + tg_yzzz_zz[i] * ra_y[i];

        tg_yzzzz_xx[i] = tg_zzzz_xx[i] * ra_y[i];

        tg_yzzzz_xy[i] = tg_zzzz_x[i] * fxi[i] + tg_zzzz_xy[i] * ra_y[i];

        tg_yzzzz_xz[i] = tg_zzzz_xz[i] * ra_y[i];

        tg_yzzzz_yy[i] = 2.0 * tg_zzzz_y[i] * fxi[i] + tg_zzzz_yy[i] * ra_y[i];

        tg_yzzzz_yz[i] = tg_zzzz_z[i] * fxi[i] + tg_zzzz_yz[i] * ra_y[i];

        tg_yzzzz_zz[i] = tg_zzzz_zz[i] * ra_y[i];

        tg_zzzzz_xx[i] = 4.0 * tg_zzz_xx[i] * fxi[i] + tg_zzzz_xx[i] * ra_z[i];

        tg_zzzzz_xy[i] = 4.0 * tg_zzz_xy[i] * fxi[i] + tg_zzzz_xy[i] * ra_z[i];

        tg_zzzzz_xz[i] = 4.0 * tg_zzz_xz[i] * fxi[i] + tg_zzzz_x[i] * fxi[i] + tg_zzzz_xz[i] * ra_z[i];

        tg_zzzzz_yy[i] = 4.0 * tg_zzz_yy[i] * fxi[i] + tg_zzzz_yy[i] * ra_z[i];

        tg_zzzzz_yz[i] = 4.0 * tg_zzz_yz[i] * fxi[i] + tg_zzzz_y[i] * fxi[i] + tg_zzzz_yz[i] * ra_z[i];

        tg_zzzzz_zz[i] = 4.0 * tg_zzz_zz[i] * fxi[i] + 2.0 * tg_zzzz_z[i] * fxi[i] + tg_zzzz_zz[i] * ra_z[i];
    }
}

} // t2lecp namespace

