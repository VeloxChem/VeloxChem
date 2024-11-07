#include "ElectricDipoleMomentumPrimRecDG.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_dg(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_dg,
                                      const size_t              idx_dip_sg,
                                      const size_t              idx_dip_pf,
                                      const size_t              idx_ovl_pg,
                                      const size_t              idx_dip_pg,
                                      const CSimdArray<double>& factors,
                                      const size_t              idx_rpa,
                                      const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : SG

    auto tr_x_0_xxxx = pbuffer.data(idx_dip_sg);

    auto tr_x_0_xxxy = pbuffer.data(idx_dip_sg + 1);

    auto tr_x_0_xxxz = pbuffer.data(idx_dip_sg + 2);

    auto tr_x_0_xxyy = pbuffer.data(idx_dip_sg + 3);

    auto tr_x_0_xxyz = pbuffer.data(idx_dip_sg + 4);

    auto tr_x_0_xxzz = pbuffer.data(idx_dip_sg + 5);

    auto tr_x_0_xyyy = pbuffer.data(idx_dip_sg + 6);

    auto tr_x_0_xyyz = pbuffer.data(idx_dip_sg + 7);

    auto tr_x_0_xyzz = pbuffer.data(idx_dip_sg + 8);

    auto tr_x_0_xzzz = pbuffer.data(idx_dip_sg + 9);

    auto tr_x_0_yyyy = pbuffer.data(idx_dip_sg + 10);

    auto tr_x_0_yyyz = pbuffer.data(idx_dip_sg + 11);

    auto tr_x_0_yyzz = pbuffer.data(idx_dip_sg + 12);

    auto tr_x_0_yzzz = pbuffer.data(idx_dip_sg + 13);

    auto tr_x_0_zzzz = pbuffer.data(idx_dip_sg + 14);

    auto tr_y_0_xxxx = pbuffer.data(idx_dip_sg + 15);

    auto tr_y_0_xxxy = pbuffer.data(idx_dip_sg + 16);

    auto tr_y_0_xxxz = pbuffer.data(idx_dip_sg + 17);

    auto tr_y_0_xxyy = pbuffer.data(idx_dip_sg + 18);

    auto tr_y_0_xxyz = pbuffer.data(idx_dip_sg + 19);

    auto tr_y_0_xxzz = pbuffer.data(idx_dip_sg + 20);

    auto tr_y_0_xyyy = pbuffer.data(idx_dip_sg + 21);

    auto tr_y_0_xyyz = pbuffer.data(idx_dip_sg + 22);

    auto tr_y_0_xyzz = pbuffer.data(idx_dip_sg + 23);

    auto tr_y_0_xzzz = pbuffer.data(idx_dip_sg + 24);

    auto tr_y_0_yyyy = pbuffer.data(idx_dip_sg + 25);

    auto tr_y_0_yyyz = pbuffer.data(idx_dip_sg + 26);

    auto tr_y_0_yyzz = pbuffer.data(idx_dip_sg + 27);

    auto tr_y_0_yzzz = pbuffer.data(idx_dip_sg + 28);

    auto tr_y_0_zzzz = pbuffer.data(idx_dip_sg + 29);

    auto tr_z_0_xxxx = pbuffer.data(idx_dip_sg + 30);

    auto tr_z_0_xxxy = pbuffer.data(idx_dip_sg + 31);

    auto tr_z_0_xxxz = pbuffer.data(idx_dip_sg + 32);

    auto tr_z_0_xxyy = pbuffer.data(idx_dip_sg + 33);

    auto tr_z_0_xxyz = pbuffer.data(idx_dip_sg + 34);

    auto tr_z_0_xxzz = pbuffer.data(idx_dip_sg + 35);

    auto tr_z_0_xyyy = pbuffer.data(idx_dip_sg + 36);

    auto tr_z_0_xyyz = pbuffer.data(idx_dip_sg + 37);

    auto tr_z_0_xyzz = pbuffer.data(idx_dip_sg + 38);

    auto tr_z_0_xzzz = pbuffer.data(idx_dip_sg + 39);

    auto tr_z_0_yyyy = pbuffer.data(idx_dip_sg + 40);

    auto tr_z_0_yyyz = pbuffer.data(idx_dip_sg + 41);

    auto tr_z_0_yyzz = pbuffer.data(idx_dip_sg + 42);

    auto tr_z_0_yzzz = pbuffer.data(idx_dip_sg + 43);

    auto tr_z_0_zzzz = pbuffer.data(idx_dip_sg + 44);

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

    // Set up components of auxiliary buffer : PG

    auto ts_x_xxxx = pbuffer.data(idx_ovl_pg);

    auto ts_x_xxxy = pbuffer.data(idx_ovl_pg + 1);

    auto ts_x_xxxz = pbuffer.data(idx_ovl_pg + 2);

    auto ts_x_xxyy = pbuffer.data(idx_ovl_pg + 3);

    auto ts_x_xxyz = pbuffer.data(idx_ovl_pg + 4);

    auto ts_x_xxzz = pbuffer.data(idx_ovl_pg + 5);

    auto ts_x_xyyy = pbuffer.data(idx_ovl_pg + 6);

    auto ts_x_xyyz = pbuffer.data(idx_ovl_pg + 7);

    auto ts_x_xyzz = pbuffer.data(idx_ovl_pg + 8);

    auto ts_x_xzzz = pbuffer.data(idx_ovl_pg + 9);

    auto ts_x_yyyy = pbuffer.data(idx_ovl_pg + 10);

    auto ts_x_yyyz = pbuffer.data(idx_ovl_pg + 11);

    auto ts_x_yyzz = pbuffer.data(idx_ovl_pg + 12);

    auto ts_x_yzzz = pbuffer.data(idx_ovl_pg + 13);

    auto ts_x_zzzz = pbuffer.data(idx_ovl_pg + 14);

    auto ts_y_xxxx = pbuffer.data(idx_ovl_pg + 15);

    auto ts_y_xxxy = pbuffer.data(idx_ovl_pg + 16);

    auto ts_y_xxxz = pbuffer.data(idx_ovl_pg + 17);

    auto ts_y_xxyy = pbuffer.data(idx_ovl_pg + 18);

    auto ts_y_xxyz = pbuffer.data(idx_ovl_pg + 19);

    auto ts_y_xxzz = pbuffer.data(idx_ovl_pg + 20);

    auto ts_y_xyyy = pbuffer.data(idx_ovl_pg + 21);

    auto ts_y_xyyz = pbuffer.data(idx_ovl_pg + 22);

    auto ts_y_xyzz = pbuffer.data(idx_ovl_pg + 23);

    auto ts_y_xzzz = pbuffer.data(idx_ovl_pg + 24);

    auto ts_y_yyyy = pbuffer.data(idx_ovl_pg + 25);

    auto ts_y_yyyz = pbuffer.data(idx_ovl_pg + 26);

    auto ts_y_yyzz = pbuffer.data(idx_ovl_pg + 27);

    auto ts_y_yzzz = pbuffer.data(idx_ovl_pg + 28);

    auto ts_y_zzzz = pbuffer.data(idx_ovl_pg + 29);

    auto ts_z_xxxx = pbuffer.data(idx_ovl_pg + 30);

    auto ts_z_xxxy = pbuffer.data(idx_ovl_pg + 31);

    auto ts_z_xxxz = pbuffer.data(idx_ovl_pg + 32);

    auto ts_z_xxyy = pbuffer.data(idx_ovl_pg + 33);

    auto ts_z_xxyz = pbuffer.data(idx_ovl_pg + 34);

    auto ts_z_xxzz = pbuffer.data(idx_ovl_pg + 35);

    auto ts_z_xyyy = pbuffer.data(idx_ovl_pg + 36);

    auto ts_z_xyyz = pbuffer.data(idx_ovl_pg + 37);

    auto ts_z_xyzz = pbuffer.data(idx_ovl_pg + 38);

    auto ts_z_xzzz = pbuffer.data(idx_ovl_pg + 39);

    auto ts_z_yyyy = pbuffer.data(idx_ovl_pg + 40);

    auto ts_z_yyyz = pbuffer.data(idx_ovl_pg + 41);

    auto ts_z_yyzz = pbuffer.data(idx_ovl_pg + 42);

    auto ts_z_yzzz = pbuffer.data(idx_ovl_pg + 43);

    auto ts_z_zzzz = pbuffer.data(idx_ovl_pg + 44);

    // Set up components of auxiliary buffer : PG

    auto tr_x_x_xxxx = pbuffer.data(idx_dip_pg);

    auto tr_x_x_xxxy = pbuffer.data(idx_dip_pg + 1);

    auto tr_x_x_xxxz = pbuffer.data(idx_dip_pg + 2);

    auto tr_x_x_xxyy = pbuffer.data(idx_dip_pg + 3);

    auto tr_x_x_xxyz = pbuffer.data(idx_dip_pg + 4);

    auto tr_x_x_xxzz = pbuffer.data(idx_dip_pg + 5);

    auto tr_x_x_xyyy = pbuffer.data(idx_dip_pg + 6);

    auto tr_x_x_xyyz = pbuffer.data(idx_dip_pg + 7);

    auto tr_x_x_xyzz = pbuffer.data(idx_dip_pg + 8);

    auto tr_x_x_xzzz = pbuffer.data(idx_dip_pg + 9);

    auto tr_x_x_yyyy = pbuffer.data(idx_dip_pg + 10);

    auto tr_x_x_yyyz = pbuffer.data(idx_dip_pg + 11);

    auto tr_x_x_yyzz = pbuffer.data(idx_dip_pg + 12);

    auto tr_x_x_yzzz = pbuffer.data(idx_dip_pg + 13);

    auto tr_x_x_zzzz = pbuffer.data(idx_dip_pg + 14);

    auto tr_x_y_xxxx = pbuffer.data(idx_dip_pg + 15);

    auto tr_x_y_xxxy = pbuffer.data(idx_dip_pg + 16);

    auto tr_x_y_xxxz = pbuffer.data(idx_dip_pg + 17);

    auto tr_x_y_xxyy = pbuffer.data(idx_dip_pg + 18);

    auto tr_x_y_xxyz = pbuffer.data(idx_dip_pg + 19);

    auto tr_x_y_xxzz = pbuffer.data(idx_dip_pg + 20);

    auto tr_x_y_xyyy = pbuffer.data(idx_dip_pg + 21);

    auto tr_x_y_xyyz = pbuffer.data(idx_dip_pg + 22);

    auto tr_x_y_xyzz = pbuffer.data(idx_dip_pg + 23);

    auto tr_x_y_xzzz = pbuffer.data(idx_dip_pg + 24);

    auto tr_x_y_yyyy = pbuffer.data(idx_dip_pg + 25);

    auto tr_x_y_yyyz = pbuffer.data(idx_dip_pg + 26);

    auto tr_x_y_yyzz = pbuffer.data(idx_dip_pg + 27);

    auto tr_x_y_yzzz = pbuffer.data(idx_dip_pg + 28);

    auto tr_x_y_zzzz = pbuffer.data(idx_dip_pg + 29);

    auto tr_x_z_xxxx = pbuffer.data(idx_dip_pg + 30);

    auto tr_x_z_xxxy = pbuffer.data(idx_dip_pg + 31);

    auto tr_x_z_xxxz = pbuffer.data(idx_dip_pg + 32);

    auto tr_x_z_xxyy = pbuffer.data(idx_dip_pg + 33);

    auto tr_x_z_xxyz = pbuffer.data(idx_dip_pg + 34);

    auto tr_x_z_xxzz = pbuffer.data(idx_dip_pg + 35);

    auto tr_x_z_xyyy = pbuffer.data(idx_dip_pg + 36);

    auto tr_x_z_xyyz = pbuffer.data(idx_dip_pg + 37);

    auto tr_x_z_xyzz = pbuffer.data(idx_dip_pg + 38);

    auto tr_x_z_xzzz = pbuffer.data(idx_dip_pg + 39);

    auto tr_x_z_yyyy = pbuffer.data(idx_dip_pg + 40);

    auto tr_x_z_yyyz = pbuffer.data(idx_dip_pg + 41);

    auto tr_x_z_yyzz = pbuffer.data(idx_dip_pg + 42);

    auto tr_x_z_yzzz = pbuffer.data(idx_dip_pg + 43);

    auto tr_x_z_zzzz = pbuffer.data(idx_dip_pg + 44);

    auto tr_y_x_xxxx = pbuffer.data(idx_dip_pg + 45);

    auto tr_y_x_xxxy = pbuffer.data(idx_dip_pg + 46);

    auto tr_y_x_xxxz = pbuffer.data(idx_dip_pg + 47);

    auto tr_y_x_xxyy = pbuffer.data(idx_dip_pg + 48);

    auto tr_y_x_xxyz = pbuffer.data(idx_dip_pg + 49);

    auto tr_y_x_xxzz = pbuffer.data(idx_dip_pg + 50);

    auto tr_y_x_xyyy = pbuffer.data(idx_dip_pg + 51);

    auto tr_y_x_xyyz = pbuffer.data(idx_dip_pg + 52);

    auto tr_y_x_xyzz = pbuffer.data(idx_dip_pg + 53);

    auto tr_y_x_xzzz = pbuffer.data(idx_dip_pg + 54);

    auto tr_y_x_yyyy = pbuffer.data(idx_dip_pg + 55);

    auto tr_y_x_yyyz = pbuffer.data(idx_dip_pg + 56);

    auto tr_y_x_yyzz = pbuffer.data(idx_dip_pg + 57);

    auto tr_y_x_yzzz = pbuffer.data(idx_dip_pg + 58);

    auto tr_y_x_zzzz = pbuffer.data(idx_dip_pg + 59);

    auto tr_y_y_xxxx = pbuffer.data(idx_dip_pg + 60);

    auto tr_y_y_xxxy = pbuffer.data(idx_dip_pg + 61);

    auto tr_y_y_xxxz = pbuffer.data(idx_dip_pg + 62);

    auto tr_y_y_xxyy = pbuffer.data(idx_dip_pg + 63);

    auto tr_y_y_xxyz = pbuffer.data(idx_dip_pg + 64);

    auto tr_y_y_xxzz = pbuffer.data(idx_dip_pg + 65);

    auto tr_y_y_xyyy = pbuffer.data(idx_dip_pg + 66);

    auto tr_y_y_xyyz = pbuffer.data(idx_dip_pg + 67);

    auto tr_y_y_xyzz = pbuffer.data(idx_dip_pg + 68);

    auto tr_y_y_xzzz = pbuffer.data(idx_dip_pg + 69);

    auto tr_y_y_yyyy = pbuffer.data(idx_dip_pg + 70);

    auto tr_y_y_yyyz = pbuffer.data(idx_dip_pg + 71);

    auto tr_y_y_yyzz = pbuffer.data(idx_dip_pg + 72);

    auto tr_y_y_yzzz = pbuffer.data(idx_dip_pg + 73);

    auto tr_y_y_zzzz = pbuffer.data(idx_dip_pg + 74);

    auto tr_y_z_xxxx = pbuffer.data(idx_dip_pg + 75);

    auto tr_y_z_xxxy = pbuffer.data(idx_dip_pg + 76);

    auto tr_y_z_xxxz = pbuffer.data(idx_dip_pg + 77);

    auto tr_y_z_xxyy = pbuffer.data(idx_dip_pg + 78);

    auto tr_y_z_xxyz = pbuffer.data(idx_dip_pg + 79);

    auto tr_y_z_xxzz = pbuffer.data(idx_dip_pg + 80);

    auto tr_y_z_xyyy = pbuffer.data(idx_dip_pg + 81);

    auto tr_y_z_xyyz = pbuffer.data(idx_dip_pg + 82);

    auto tr_y_z_xyzz = pbuffer.data(idx_dip_pg + 83);

    auto tr_y_z_xzzz = pbuffer.data(idx_dip_pg + 84);

    auto tr_y_z_yyyy = pbuffer.data(idx_dip_pg + 85);

    auto tr_y_z_yyyz = pbuffer.data(idx_dip_pg + 86);

    auto tr_y_z_yyzz = pbuffer.data(idx_dip_pg + 87);

    auto tr_y_z_yzzz = pbuffer.data(idx_dip_pg + 88);

    auto tr_y_z_zzzz = pbuffer.data(idx_dip_pg + 89);

    auto tr_z_x_xxxx = pbuffer.data(idx_dip_pg + 90);

    auto tr_z_x_xxxy = pbuffer.data(idx_dip_pg + 91);

    auto tr_z_x_xxxz = pbuffer.data(idx_dip_pg + 92);

    auto tr_z_x_xxyy = pbuffer.data(idx_dip_pg + 93);

    auto tr_z_x_xxyz = pbuffer.data(idx_dip_pg + 94);

    auto tr_z_x_xxzz = pbuffer.data(idx_dip_pg + 95);

    auto tr_z_x_xyyy = pbuffer.data(idx_dip_pg + 96);

    auto tr_z_x_xyyz = pbuffer.data(idx_dip_pg + 97);

    auto tr_z_x_xyzz = pbuffer.data(idx_dip_pg + 98);

    auto tr_z_x_xzzz = pbuffer.data(idx_dip_pg + 99);

    auto tr_z_x_yyyy = pbuffer.data(idx_dip_pg + 100);

    auto tr_z_x_yyyz = pbuffer.data(idx_dip_pg + 101);

    auto tr_z_x_yyzz = pbuffer.data(idx_dip_pg + 102);

    auto tr_z_x_yzzz = pbuffer.data(idx_dip_pg + 103);

    auto tr_z_x_zzzz = pbuffer.data(idx_dip_pg + 104);

    auto tr_z_y_xxxx = pbuffer.data(idx_dip_pg + 105);

    auto tr_z_y_xxxy = pbuffer.data(idx_dip_pg + 106);

    auto tr_z_y_xxxz = pbuffer.data(idx_dip_pg + 107);

    auto tr_z_y_xxyy = pbuffer.data(idx_dip_pg + 108);

    auto tr_z_y_xxyz = pbuffer.data(idx_dip_pg + 109);

    auto tr_z_y_xxzz = pbuffer.data(idx_dip_pg + 110);

    auto tr_z_y_xyyy = pbuffer.data(idx_dip_pg + 111);

    auto tr_z_y_xyyz = pbuffer.data(idx_dip_pg + 112);

    auto tr_z_y_xyzz = pbuffer.data(idx_dip_pg + 113);

    auto tr_z_y_xzzz = pbuffer.data(idx_dip_pg + 114);

    auto tr_z_y_yyyy = pbuffer.data(idx_dip_pg + 115);

    auto tr_z_y_yyyz = pbuffer.data(idx_dip_pg + 116);

    auto tr_z_y_yyzz = pbuffer.data(idx_dip_pg + 117);

    auto tr_z_y_yzzz = pbuffer.data(idx_dip_pg + 118);

    auto tr_z_y_zzzz = pbuffer.data(idx_dip_pg + 119);

    auto tr_z_z_xxxx = pbuffer.data(idx_dip_pg + 120);

    auto tr_z_z_xxxy = pbuffer.data(idx_dip_pg + 121);

    auto tr_z_z_xxxz = pbuffer.data(idx_dip_pg + 122);

    auto tr_z_z_xxyy = pbuffer.data(idx_dip_pg + 123);

    auto tr_z_z_xxyz = pbuffer.data(idx_dip_pg + 124);

    auto tr_z_z_xxzz = pbuffer.data(idx_dip_pg + 125);

    auto tr_z_z_xyyy = pbuffer.data(idx_dip_pg + 126);

    auto tr_z_z_xyyz = pbuffer.data(idx_dip_pg + 127);

    auto tr_z_z_xyzz = pbuffer.data(idx_dip_pg + 128);

    auto tr_z_z_xzzz = pbuffer.data(idx_dip_pg + 129);

    auto tr_z_z_yyyy = pbuffer.data(idx_dip_pg + 130);

    auto tr_z_z_yyyz = pbuffer.data(idx_dip_pg + 131);

    auto tr_z_z_yyzz = pbuffer.data(idx_dip_pg + 132);

    auto tr_z_z_yzzz = pbuffer.data(idx_dip_pg + 133);

    auto tr_z_z_zzzz = pbuffer.data(idx_dip_pg + 134);

    // Set up 0-15 components of targeted buffer : DG

    auto tr_x_xx_xxxx = pbuffer.data(idx_dip_dg);

    auto tr_x_xx_xxxy = pbuffer.data(idx_dip_dg + 1);

    auto tr_x_xx_xxxz = pbuffer.data(idx_dip_dg + 2);

    auto tr_x_xx_xxyy = pbuffer.data(idx_dip_dg + 3);

    auto tr_x_xx_xxyz = pbuffer.data(idx_dip_dg + 4);

    auto tr_x_xx_xxzz = pbuffer.data(idx_dip_dg + 5);

    auto tr_x_xx_xyyy = pbuffer.data(idx_dip_dg + 6);

    auto tr_x_xx_xyyz = pbuffer.data(idx_dip_dg + 7);

    auto tr_x_xx_xyzz = pbuffer.data(idx_dip_dg + 8);

    auto tr_x_xx_xzzz = pbuffer.data(idx_dip_dg + 9);

    auto tr_x_xx_yyyy = pbuffer.data(idx_dip_dg + 10);

    auto tr_x_xx_yyyz = pbuffer.data(idx_dip_dg + 11);

    auto tr_x_xx_yyzz = pbuffer.data(idx_dip_dg + 12);

    auto tr_x_xx_yzzz = pbuffer.data(idx_dip_dg + 13);

    auto tr_x_xx_zzzz = pbuffer.data(idx_dip_dg + 14);

#pragma omp simd aligned(pa_x,             \
                             tr_x_0_xxxx,  \
                             tr_x_0_xxxy,  \
                             tr_x_0_xxxz,  \
                             tr_x_0_xxyy,  \
                             tr_x_0_xxyz,  \
                             tr_x_0_xxzz,  \
                             tr_x_0_xyyy,  \
                             tr_x_0_xyyz,  \
                             tr_x_0_xyzz,  \
                             tr_x_0_xzzz,  \
                             tr_x_0_yyyy,  \
                             tr_x_0_yyyz,  \
                             tr_x_0_yyzz,  \
                             tr_x_0_yzzz,  \
                             tr_x_0_zzzz,  \
                             tr_x_x_xxx,   \
                             tr_x_x_xxxx,  \
                             tr_x_x_xxxy,  \
                             tr_x_x_xxxz,  \
                             tr_x_x_xxy,   \
                             tr_x_x_xxyy,  \
                             tr_x_x_xxyz,  \
                             tr_x_x_xxz,   \
                             tr_x_x_xxzz,  \
                             tr_x_x_xyy,   \
                             tr_x_x_xyyy,  \
                             tr_x_x_xyyz,  \
                             tr_x_x_xyz,   \
                             tr_x_x_xyzz,  \
                             tr_x_x_xzz,   \
                             tr_x_x_xzzz,  \
                             tr_x_x_yyy,   \
                             tr_x_x_yyyy,  \
                             tr_x_x_yyyz,  \
                             tr_x_x_yyz,   \
                             tr_x_x_yyzz,  \
                             tr_x_x_yzz,   \
                             tr_x_x_yzzz,  \
                             tr_x_x_zzz,   \
                             tr_x_x_zzzz,  \
                             tr_x_xx_xxxx, \
                             tr_x_xx_xxxy, \
                             tr_x_xx_xxxz, \
                             tr_x_xx_xxyy, \
                             tr_x_xx_xxyz, \
                             tr_x_xx_xxzz, \
                             tr_x_xx_xyyy, \
                             tr_x_xx_xyyz, \
                             tr_x_xx_xyzz, \
                             tr_x_xx_xzzz, \
                             tr_x_xx_yyyy, \
                             tr_x_xx_yyyz, \
                             tr_x_xx_yyzz, \
                             tr_x_xx_yzzz, \
                             tr_x_xx_zzzz, \
                             ts_x_xxxx,    \
                             ts_x_xxxy,    \
                             ts_x_xxxz,    \
                             ts_x_xxyy,    \
                             ts_x_xxyz,    \
                             ts_x_xxzz,    \
                             ts_x_xyyy,    \
                             ts_x_xyyz,    \
                             ts_x_xyzz,    \
                             ts_x_xzzz,    \
                             ts_x_yyyy,    \
                             ts_x_yyyz,    \
                             ts_x_yyzz,    \
                             ts_x_yzzz,    \
                             ts_x_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xx_xxxx[i] = tr_x_0_xxxx[i] * fe_0 + 4.0 * tr_x_x_xxx[i] * fe_0 + ts_x_xxxx[i] * fe_0 + tr_x_x_xxxx[i] * pa_x[i];

        tr_x_xx_xxxy[i] = tr_x_0_xxxy[i] * fe_0 + 3.0 * tr_x_x_xxy[i] * fe_0 + ts_x_xxxy[i] * fe_0 + tr_x_x_xxxy[i] * pa_x[i];

        tr_x_xx_xxxz[i] = tr_x_0_xxxz[i] * fe_0 + 3.0 * tr_x_x_xxz[i] * fe_0 + ts_x_xxxz[i] * fe_0 + tr_x_x_xxxz[i] * pa_x[i];

        tr_x_xx_xxyy[i] = tr_x_0_xxyy[i] * fe_0 + 2.0 * tr_x_x_xyy[i] * fe_0 + ts_x_xxyy[i] * fe_0 + tr_x_x_xxyy[i] * pa_x[i];

        tr_x_xx_xxyz[i] = tr_x_0_xxyz[i] * fe_0 + 2.0 * tr_x_x_xyz[i] * fe_0 + ts_x_xxyz[i] * fe_0 + tr_x_x_xxyz[i] * pa_x[i];

        tr_x_xx_xxzz[i] = tr_x_0_xxzz[i] * fe_0 + 2.0 * tr_x_x_xzz[i] * fe_0 + ts_x_xxzz[i] * fe_0 + tr_x_x_xxzz[i] * pa_x[i];

        tr_x_xx_xyyy[i] = tr_x_0_xyyy[i] * fe_0 + tr_x_x_yyy[i] * fe_0 + ts_x_xyyy[i] * fe_0 + tr_x_x_xyyy[i] * pa_x[i];

        tr_x_xx_xyyz[i] = tr_x_0_xyyz[i] * fe_0 + tr_x_x_yyz[i] * fe_0 + ts_x_xyyz[i] * fe_0 + tr_x_x_xyyz[i] * pa_x[i];

        tr_x_xx_xyzz[i] = tr_x_0_xyzz[i] * fe_0 + tr_x_x_yzz[i] * fe_0 + ts_x_xyzz[i] * fe_0 + tr_x_x_xyzz[i] * pa_x[i];

        tr_x_xx_xzzz[i] = tr_x_0_xzzz[i] * fe_0 + tr_x_x_zzz[i] * fe_0 + ts_x_xzzz[i] * fe_0 + tr_x_x_xzzz[i] * pa_x[i];

        tr_x_xx_yyyy[i] = tr_x_0_yyyy[i] * fe_0 + ts_x_yyyy[i] * fe_0 + tr_x_x_yyyy[i] * pa_x[i];

        tr_x_xx_yyyz[i] = tr_x_0_yyyz[i] * fe_0 + ts_x_yyyz[i] * fe_0 + tr_x_x_yyyz[i] * pa_x[i];

        tr_x_xx_yyzz[i] = tr_x_0_yyzz[i] * fe_0 + ts_x_yyzz[i] * fe_0 + tr_x_x_yyzz[i] * pa_x[i];

        tr_x_xx_yzzz[i] = tr_x_0_yzzz[i] * fe_0 + ts_x_yzzz[i] * fe_0 + tr_x_x_yzzz[i] * pa_x[i];

        tr_x_xx_zzzz[i] = tr_x_0_zzzz[i] * fe_0 + ts_x_zzzz[i] * fe_0 + tr_x_x_zzzz[i] * pa_x[i];
    }

    // Set up 15-30 components of targeted buffer : DG

    auto tr_x_xy_xxxx = pbuffer.data(idx_dip_dg + 15);

    auto tr_x_xy_xxxy = pbuffer.data(idx_dip_dg + 16);

    auto tr_x_xy_xxxz = pbuffer.data(idx_dip_dg + 17);

    auto tr_x_xy_xxyy = pbuffer.data(idx_dip_dg + 18);

    auto tr_x_xy_xxyz = pbuffer.data(idx_dip_dg + 19);

    auto tr_x_xy_xxzz = pbuffer.data(idx_dip_dg + 20);

    auto tr_x_xy_xyyy = pbuffer.data(idx_dip_dg + 21);

    auto tr_x_xy_xyyz = pbuffer.data(idx_dip_dg + 22);

    auto tr_x_xy_xyzz = pbuffer.data(idx_dip_dg + 23);

    auto tr_x_xy_xzzz = pbuffer.data(idx_dip_dg + 24);

    auto tr_x_xy_yyyy = pbuffer.data(idx_dip_dg + 25);

    auto tr_x_xy_yyyz = pbuffer.data(idx_dip_dg + 26);

    auto tr_x_xy_yyzz = pbuffer.data(idx_dip_dg + 27);

    auto tr_x_xy_yzzz = pbuffer.data(idx_dip_dg + 28);

    auto tr_x_xy_zzzz = pbuffer.data(idx_dip_dg + 29);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             tr_x_x_xxx,   \
                             tr_x_x_xxxx,  \
                             tr_x_x_xxxy,  \
                             tr_x_x_xxxz,  \
                             tr_x_x_xxy,   \
                             tr_x_x_xxyy,  \
                             tr_x_x_xxyz,  \
                             tr_x_x_xxz,   \
                             tr_x_x_xxzz,  \
                             tr_x_x_xyy,   \
                             tr_x_x_xyyy,  \
                             tr_x_x_xyyz,  \
                             tr_x_x_xyz,   \
                             tr_x_x_xyzz,  \
                             tr_x_x_xzz,   \
                             tr_x_x_xzzz,  \
                             tr_x_x_zzzz,  \
                             tr_x_xy_xxxx, \
                             tr_x_xy_xxxy, \
                             tr_x_xy_xxxz, \
                             tr_x_xy_xxyy, \
                             tr_x_xy_xxyz, \
                             tr_x_xy_xxzz, \
                             tr_x_xy_xyyy, \
                             tr_x_xy_xyyz, \
                             tr_x_xy_xyzz, \
                             tr_x_xy_xzzz, \
                             tr_x_xy_yyyy, \
                             tr_x_xy_yyyz, \
                             tr_x_xy_yyzz, \
                             tr_x_xy_yzzz, \
                             tr_x_xy_zzzz, \
                             tr_x_y_yyyy,  \
                             tr_x_y_yyyz,  \
                             tr_x_y_yyzz,  \
                             tr_x_y_yzzz,  \
                             ts_y_yyyy,    \
                             ts_y_yyyz,    \
                             ts_y_yyzz,    \
                             ts_y_yzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xy_xxxx[i] = tr_x_x_xxxx[i] * pa_y[i];

        tr_x_xy_xxxy[i] = tr_x_x_xxx[i] * fe_0 + tr_x_x_xxxy[i] * pa_y[i];

        tr_x_xy_xxxz[i] = tr_x_x_xxxz[i] * pa_y[i];

        tr_x_xy_xxyy[i] = 2.0 * tr_x_x_xxy[i] * fe_0 + tr_x_x_xxyy[i] * pa_y[i];

        tr_x_xy_xxyz[i] = tr_x_x_xxz[i] * fe_0 + tr_x_x_xxyz[i] * pa_y[i];

        tr_x_xy_xxzz[i] = tr_x_x_xxzz[i] * pa_y[i];

        tr_x_xy_xyyy[i] = 3.0 * tr_x_x_xyy[i] * fe_0 + tr_x_x_xyyy[i] * pa_y[i];

        tr_x_xy_xyyz[i] = 2.0 * tr_x_x_xyz[i] * fe_0 + tr_x_x_xyyz[i] * pa_y[i];

        tr_x_xy_xyzz[i] = tr_x_x_xzz[i] * fe_0 + tr_x_x_xyzz[i] * pa_y[i];

        tr_x_xy_xzzz[i] = tr_x_x_xzzz[i] * pa_y[i];

        tr_x_xy_yyyy[i] = ts_y_yyyy[i] * fe_0 + tr_x_y_yyyy[i] * pa_x[i];

        tr_x_xy_yyyz[i] = ts_y_yyyz[i] * fe_0 + tr_x_y_yyyz[i] * pa_x[i];

        tr_x_xy_yyzz[i] = ts_y_yyzz[i] * fe_0 + tr_x_y_yyzz[i] * pa_x[i];

        tr_x_xy_yzzz[i] = ts_y_yzzz[i] * fe_0 + tr_x_y_yzzz[i] * pa_x[i];

        tr_x_xy_zzzz[i] = tr_x_x_zzzz[i] * pa_y[i];
    }

    // Set up 30-45 components of targeted buffer : DG

    auto tr_x_xz_xxxx = pbuffer.data(idx_dip_dg + 30);

    auto tr_x_xz_xxxy = pbuffer.data(idx_dip_dg + 31);

    auto tr_x_xz_xxxz = pbuffer.data(idx_dip_dg + 32);

    auto tr_x_xz_xxyy = pbuffer.data(idx_dip_dg + 33);

    auto tr_x_xz_xxyz = pbuffer.data(idx_dip_dg + 34);

    auto tr_x_xz_xxzz = pbuffer.data(idx_dip_dg + 35);

    auto tr_x_xz_xyyy = pbuffer.data(idx_dip_dg + 36);

    auto tr_x_xz_xyyz = pbuffer.data(idx_dip_dg + 37);

    auto tr_x_xz_xyzz = pbuffer.data(idx_dip_dg + 38);

    auto tr_x_xz_xzzz = pbuffer.data(idx_dip_dg + 39);

    auto tr_x_xz_yyyy = pbuffer.data(idx_dip_dg + 40);

    auto tr_x_xz_yyyz = pbuffer.data(idx_dip_dg + 41);

    auto tr_x_xz_yyzz = pbuffer.data(idx_dip_dg + 42);

    auto tr_x_xz_yzzz = pbuffer.data(idx_dip_dg + 43);

    auto tr_x_xz_zzzz = pbuffer.data(idx_dip_dg + 44);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             tr_x_x_xxx,   \
                             tr_x_x_xxxx,  \
                             tr_x_x_xxxy,  \
                             tr_x_x_xxxz,  \
                             tr_x_x_xxy,   \
                             tr_x_x_xxyy,  \
                             tr_x_x_xxyz,  \
                             tr_x_x_xxz,   \
                             tr_x_x_xxzz,  \
                             tr_x_x_xyy,   \
                             tr_x_x_xyyy,  \
                             tr_x_x_xyyz,  \
                             tr_x_x_xyz,   \
                             tr_x_x_xyzz,  \
                             tr_x_x_xzz,   \
                             tr_x_x_xzzz,  \
                             tr_x_x_yyyy,  \
                             tr_x_xz_xxxx, \
                             tr_x_xz_xxxy, \
                             tr_x_xz_xxxz, \
                             tr_x_xz_xxyy, \
                             tr_x_xz_xxyz, \
                             tr_x_xz_xxzz, \
                             tr_x_xz_xyyy, \
                             tr_x_xz_xyyz, \
                             tr_x_xz_xyzz, \
                             tr_x_xz_xzzz, \
                             tr_x_xz_yyyy, \
                             tr_x_xz_yyyz, \
                             tr_x_xz_yyzz, \
                             tr_x_xz_yzzz, \
                             tr_x_xz_zzzz, \
                             tr_x_z_yyyz,  \
                             tr_x_z_yyzz,  \
                             tr_x_z_yzzz,  \
                             tr_x_z_zzzz,  \
                             ts_z_yyyz,    \
                             ts_z_yyzz,    \
                             ts_z_yzzz,    \
                             ts_z_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xz_xxxx[i] = tr_x_x_xxxx[i] * pa_z[i];

        tr_x_xz_xxxy[i] = tr_x_x_xxxy[i] * pa_z[i];

        tr_x_xz_xxxz[i] = tr_x_x_xxx[i] * fe_0 + tr_x_x_xxxz[i] * pa_z[i];

        tr_x_xz_xxyy[i] = tr_x_x_xxyy[i] * pa_z[i];

        tr_x_xz_xxyz[i] = tr_x_x_xxy[i] * fe_0 + tr_x_x_xxyz[i] * pa_z[i];

        tr_x_xz_xxzz[i] = 2.0 * tr_x_x_xxz[i] * fe_0 + tr_x_x_xxzz[i] * pa_z[i];

        tr_x_xz_xyyy[i] = tr_x_x_xyyy[i] * pa_z[i];

        tr_x_xz_xyyz[i] = tr_x_x_xyy[i] * fe_0 + tr_x_x_xyyz[i] * pa_z[i];

        tr_x_xz_xyzz[i] = 2.0 * tr_x_x_xyz[i] * fe_0 + tr_x_x_xyzz[i] * pa_z[i];

        tr_x_xz_xzzz[i] = 3.0 * tr_x_x_xzz[i] * fe_0 + tr_x_x_xzzz[i] * pa_z[i];

        tr_x_xz_yyyy[i] = tr_x_x_yyyy[i] * pa_z[i];

        tr_x_xz_yyyz[i] = ts_z_yyyz[i] * fe_0 + tr_x_z_yyyz[i] * pa_x[i];

        tr_x_xz_yyzz[i] = ts_z_yyzz[i] * fe_0 + tr_x_z_yyzz[i] * pa_x[i];

        tr_x_xz_yzzz[i] = ts_z_yzzz[i] * fe_0 + tr_x_z_yzzz[i] * pa_x[i];

        tr_x_xz_zzzz[i] = ts_z_zzzz[i] * fe_0 + tr_x_z_zzzz[i] * pa_x[i];
    }

    // Set up 45-60 components of targeted buffer : DG

    auto tr_x_yy_xxxx = pbuffer.data(idx_dip_dg + 45);

    auto tr_x_yy_xxxy = pbuffer.data(idx_dip_dg + 46);

    auto tr_x_yy_xxxz = pbuffer.data(idx_dip_dg + 47);

    auto tr_x_yy_xxyy = pbuffer.data(idx_dip_dg + 48);

    auto tr_x_yy_xxyz = pbuffer.data(idx_dip_dg + 49);

    auto tr_x_yy_xxzz = pbuffer.data(idx_dip_dg + 50);

    auto tr_x_yy_xyyy = pbuffer.data(idx_dip_dg + 51);

    auto tr_x_yy_xyyz = pbuffer.data(idx_dip_dg + 52);

    auto tr_x_yy_xyzz = pbuffer.data(idx_dip_dg + 53);

    auto tr_x_yy_xzzz = pbuffer.data(idx_dip_dg + 54);

    auto tr_x_yy_yyyy = pbuffer.data(idx_dip_dg + 55);

    auto tr_x_yy_yyyz = pbuffer.data(idx_dip_dg + 56);

    auto tr_x_yy_yyzz = pbuffer.data(idx_dip_dg + 57);

    auto tr_x_yy_yzzz = pbuffer.data(idx_dip_dg + 58);

    auto tr_x_yy_zzzz = pbuffer.data(idx_dip_dg + 59);

#pragma omp simd aligned(pa_y,             \
                             tr_x_0_xxxx,  \
                             tr_x_0_xxxy,  \
                             tr_x_0_xxxz,  \
                             tr_x_0_xxyy,  \
                             tr_x_0_xxyz,  \
                             tr_x_0_xxzz,  \
                             tr_x_0_xyyy,  \
                             tr_x_0_xyyz,  \
                             tr_x_0_xyzz,  \
                             tr_x_0_xzzz,  \
                             tr_x_0_yyyy,  \
                             tr_x_0_yyyz,  \
                             tr_x_0_yyzz,  \
                             tr_x_0_yzzz,  \
                             tr_x_0_zzzz,  \
                             tr_x_y_xxx,   \
                             tr_x_y_xxxx,  \
                             tr_x_y_xxxy,  \
                             tr_x_y_xxxz,  \
                             tr_x_y_xxy,   \
                             tr_x_y_xxyy,  \
                             tr_x_y_xxyz,  \
                             tr_x_y_xxz,   \
                             tr_x_y_xxzz,  \
                             tr_x_y_xyy,   \
                             tr_x_y_xyyy,  \
                             tr_x_y_xyyz,  \
                             tr_x_y_xyz,   \
                             tr_x_y_xyzz,  \
                             tr_x_y_xzz,   \
                             tr_x_y_xzzz,  \
                             tr_x_y_yyy,   \
                             tr_x_y_yyyy,  \
                             tr_x_y_yyyz,  \
                             tr_x_y_yyz,   \
                             tr_x_y_yyzz,  \
                             tr_x_y_yzz,   \
                             tr_x_y_yzzz,  \
                             tr_x_y_zzz,   \
                             tr_x_y_zzzz,  \
                             tr_x_yy_xxxx, \
                             tr_x_yy_xxxy, \
                             tr_x_yy_xxxz, \
                             tr_x_yy_xxyy, \
                             tr_x_yy_xxyz, \
                             tr_x_yy_xxzz, \
                             tr_x_yy_xyyy, \
                             tr_x_yy_xyyz, \
                             tr_x_yy_xyzz, \
                             tr_x_yy_xzzz, \
                             tr_x_yy_yyyy, \
                             tr_x_yy_yyyz, \
                             tr_x_yy_yyzz, \
                             tr_x_yy_yzzz, \
                             tr_x_yy_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yy_xxxx[i] = tr_x_0_xxxx[i] * fe_0 + tr_x_y_xxxx[i] * pa_y[i];

        tr_x_yy_xxxy[i] = tr_x_0_xxxy[i] * fe_0 + tr_x_y_xxx[i] * fe_0 + tr_x_y_xxxy[i] * pa_y[i];

        tr_x_yy_xxxz[i] = tr_x_0_xxxz[i] * fe_0 + tr_x_y_xxxz[i] * pa_y[i];

        tr_x_yy_xxyy[i] = tr_x_0_xxyy[i] * fe_0 + 2.0 * tr_x_y_xxy[i] * fe_0 + tr_x_y_xxyy[i] * pa_y[i];

        tr_x_yy_xxyz[i] = tr_x_0_xxyz[i] * fe_0 + tr_x_y_xxz[i] * fe_0 + tr_x_y_xxyz[i] * pa_y[i];

        tr_x_yy_xxzz[i] = tr_x_0_xxzz[i] * fe_0 + tr_x_y_xxzz[i] * pa_y[i];

        tr_x_yy_xyyy[i] = tr_x_0_xyyy[i] * fe_0 + 3.0 * tr_x_y_xyy[i] * fe_0 + tr_x_y_xyyy[i] * pa_y[i];

        tr_x_yy_xyyz[i] = tr_x_0_xyyz[i] * fe_0 + 2.0 * tr_x_y_xyz[i] * fe_0 + tr_x_y_xyyz[i] * pa_y[i];

        tr_x_yy_xyzz[i] = tr_x_0_xyzz[i] * fe_0 + tr_x_y_xzz[i] * fe_0 + tr_x_y_xyzz[i] * pa_y[i];

        tr_x_yy_xzzz[i] = tr_x_0_xzzz[i] * fe_0 + tr_x_y_xzzz[i] * pa_y[i];

        tr_x_yy_yyyy[i] = tr_x_0_yyyy[i] * fe_0 + 4.0 * tr_x_y_yyy[i] * fe_0 + tr_x_y_yyyy[i] * pa_y[i];

        tr_x_yy_yyyz[i] = tr_x_0_yyyz[i] * fe_0 + 3.0 * tr_x_y_yyz[i] * fe_0 + tr_x_y_yyyz[i] * pa_y[i];

        tr_x_yy_yyzz[i] = tr_x_0_yyzz[i] * fe_0 + 2.0 * tr_x_y_yzz[i] * fe_0 + tr_x_y_yyzz[i] * pa_y[i];

        tr_x_yy_yzzz[i] = tr_x_0_yzzz[i] * fe_0 + tr_x_y_zzz[i] * fe_0 + tr_x_y_yzzz[i] * pa_y[i];

        tr_x_yy_zzzz[i] = tr_x_0_zzzz[i] * fe_0 + tr_x_y_zzzz[i] * pa_y[i];
    }

    // Set up 60-75 components of targeted buffer : DG

    auto tr_x_yz_xxxx = pbuffer.data(idx_dip_dg + 60);

    auto tr_x_yz_xxxy = pbuffer.data(idx_dip_dg + 61);

    auto tr_x_yz_xxxz = pbuffer.data(idx_dip_dg + 62);

    auto tr_x_yz_xxyy = pbuffer.data(idx_dip_dg + 63);

    auto tr_x_yz_xxyz = pbuffer.data(idx_dip_dg + 64);

    auto tr_x_yz_xxzz = pbuffer.data(idx_dip_dg + 65);

    auto tr_x_yz_xyyy = pbuffer.data(idx_dip_dg + 66);

    auto tr_x_yz_xyyz = pbuffer.data(idx_dip_dg + 67);

    auto tr_x_yz_xyzz = pbuffer.data(idx_dip_dg + 68);

    auto tr_x_yz_xzzz = pbuffer.data(idx_dip_dg + 69);

    auto tr_x_yz_yyyy = pbuffer.data(idx_dip_dg + 70);

    auto tr_x_yz_yyyz = pbuffer.data(idx_dip_dg + 71);

    auto tr_x_yz_yyzz = pbuffer.data(idx_dip_dg + 72);

    auto tr_x_yz_yzzz = pbuffer.data(idx_dip_dg + 73);

    auto tr_x_yz_zzzz = pbuffer.data(idx_dip_dg + 74);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             tr_x_y_xxxy,  \
                             tr_x_y_xxyy,  \
                             tr_x_y_xyyy,  \
                             tr_x_y_yyyy,  \
                             tr_x_yz_xxxx, \
                             tr_x_yz_xxxy, \
                             tr_x_yz_xxxz, \
                             tr_x_yz_xxyy, \
                             tr_x_yz_xxyz, \
                             tr_x_yz_xxzz, \
                             tr_x_yz_xyyy, \
                             tr_x_yz_xyyz, \
                             tr_x_yz_xyzz, \
                             tr_x_yz_xzzz, \
                             tr_x_yz_yyyy, \
                             tr_x_yz_yyyz, \
                             tr_x_yz_yyzz, \
                             tr_x_yz_yzzz, \
                             tr_x_yz_zzzz, \
                             tr_x_z_xxxx,  \
                             tr_x_z_xxxz,  \
                             tr_x_z_xxyz,  \
                             tr_x_z_xxz,   \
                             tr_x_z_xxzz,  \
                             tr_x_z_xyyz,  \
                             tr_x_z_xyz,   \
                             tr_x_z_xyzz,  \
                             tr_x_z_xzz,   \
                             tr_x_z_xzzz,  \
                             tr_x_z_yyyz,  \
                             tr_x_z_yyz,   \
                             tr_x_z_yyzz,  \
                             tr_x_z_yzz,   \
                             tr_x_z_yzzz,  \
                             tr_x_z_zzz,   \
                             tr_x_z_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yz_xxxx[i] = tr_x_z_xxxx[i] * pa_y[i];

        tr_x_yz_xxxy[i] = tr_x_y_xxxy[i] * pa_z[i];

        tr_x_yz_xxxz[i] = tr_x_z_xxxz[i] * pa_y[i];

        tr_x_yz_xxyy[i] = tr_x_y_xxyy[i] * pa_z[i];

        tr_x_yz_xxyz[i] = tr_x_z_xxz[i] * fe_0 + tr_x_z_xxyz[i] * pa_y[i];

        tr_x_yz_xxzz[i] = tr_x_z_xxzz[i] * pa_y[i];

        tr_x_yz_xyyy[i] = tr_x_y_xyyy[i] * pa_z[i];

        tr_x_yz_xyyz[i] = 2.0 * tr_x_z_xyz[i] * fe_0 + tr_x_z_xyyz[i] * pa_y[i];

        tr_x_yz_xyzz[i] = tr_x_z_xzz[i] * fe_0 + tr_x_z_xyzz[i] * pa_y[i];

        tr_x_yz_xzzz[i] = tr_x_z_xzzz[i] * pa_y[i];

        tr_x_yz_yyyy[i] = tr_x_y_yyyy[i] * pa_z[i];

        tr_x_yz_yyyz[i] = 3.0 * tr_x_z_yyz[i] * fe_0 + tr_x_z_yyyz[i] * pa_y[i];

        tr_x_yz_yyzz[i] = 2.0 * tr_x_z_yzz[i] * fe_0 + tr_x_z_yyzz[i] * pa_y[i];

        tr_x_yz_yzzz[i] = tr_x_z_zzz[i] * fe_0 + tr_x_z_yzzz[i] * pa_y[i];

        tr_x_yz_zzzz[i] = tr_x_z_zzzz[i] * pa_y[i];
    }

    // Set up 75-90 components of targeted buffer : DG

    auto tr_x_zz_xxxx = pbuffer.data(idx_dip_dg + 75);

    auto tr_x_zz_xxxy = pbuffer.data(idx_dip_dg + 76);

    auto tr_x_zz_xxxz = pbuffer.data(idx_dip_dg + 77);

    auto tr_x_zz_xxyy = pbuffer.data(idx_dip_dg + 78);

    auto tr_x_zz_xxyz = pbuffer.data(idx_dip_dg + 79);

    auto tr_x_zz_xxzz = pbuffer.data(idx_dip_dg + 80);

    auto tr_x_zz_xyyy = pbuffer.data(idx_dip_dg + 81);

    auto tr_x_zz_xyyz = pbuffer.data(idx_dip_dg + 82);

    auto tr_x_zz_xyzz = pbuffer.data(idx_dip_dg + 83);

    auto tr_x_zz_xzzz = pbuffer.data(idx_dip_dg + 84);

    auto tr_x_zz_yyyy = pbuffer.data(idx_dip_dg + 85);

    auto tr_x_zz_yyyz = pbuffer.data(idx_dip_dg + 86);

    auto tr_x_zz_yyzz = pbuffer.data(idx_dip_dg + 87);

    auto tr_x_zz_yzzz = pbuffer.data(idx_dip_dg + 88);

    auto tr_x_zz_zzzz = pbuffer.data(idx_dip_dg + 89);

#pragma omp simd aligned(pa_z,             \
                             tr_x_0_xxxx,  \
                             tr_x_0_xxxy,  \
                             tr_x_0_xxxz,  \
                             tr_x_0_xxyy,  \
                             tr_x_0_xxyz,  \
                             tr_x_0_xxzz,  \
                             tr_x_0_xyyy,  \
                             tr_x_0_xyyz,  \
                             tr_x_0_xyzz,  \
                             tr_x_0_xzzz,  \
                             tr_x_0_yyyy,  \
                             tr_x_0_yyyz,  \
                             tr_x_0_yyzz,  \
                             tr_x_0_yzzz,  \
                             tr_x_0_zzzz,  \
                             tr_x_z_xxx,   \
                             tr_x_z_xxxx,  \
                             tr_x_z_xxxy,  \
                             tr_x_z_xxxz,  \
                             tr_x_z_xxy,   \
                             tr_x_z_xxyy,  \
                             tr_x_z_xxyz,  \
                             tr_x_z_xxz,   \
                             tr_x_z_xxzz,  \
                             tr_x_z_xyy,   \
                             tr_x_z_xyyy,  \
                             tr_x_z_xyyz,  \
                             tr_x_z_xyz,   \
                             tr_x_z_xyzz,  \
                             tr_x_z_xzz,   \
                             tr_x_z_xzzz,  \
                             tr_x_z_yyy,   \
                             tr_x_z_yyyy,  \
                             tr_x_z_yyyz,  \
                             tr_x_z_yyz,   \
                             tr_x_z_yyzz,  \
                             tr_x_z_yzz,   \
                             tr_x_z_yzzz,  \
                             tr_x_z_zzz,   \
                             tr_x_z_zzzz,  \
                             tr_x_zz_xxxx, \
                             tr_x_zz_xxxy, \
                             tr_x_zz_xxxz, \
                             tr_x_zz_xxyy, \
                             tr_x_zz_xxyz, \
                             tr_x_zz_xxzz, \
                             tr_x_zz_xyyy, \
                             tr_x_zz_xyyz, \
                             tr_x_zz_xyzz, \
                             tr_x_zz_xzzz, \
                             tr_x_zz_yyyy, \
                             tr_x_zz_yyyz, \
                             tr_x_zz_yyzz, \
                             tr_x_zz_yzzz, \
                             tr_x_zz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_zz_xxxx[i] = tr_x_0_xxxx[i] * fe_0 + tr_x_z_xxxx[i] * pa_z[i];

        tr_x_zz_xxxy[i] = tr_x_0_xxxy[i] * fe_0 + tr_x_z_xxxy[i] * pa_z[i];

        tr_x_zz_xxxz[i] = tr_x_0_xxxz[i] * fe_0 + tr_x_z_xxx[i] * fe_0 + tr_x_z_xxxz[i] * pa_z[i];

        tr_x_zz_xxyy[i] = tr_x_0_xxyy[i] * fe_0 + tr_x_z_xxyy[i] * pa_z[i];

        tr_x_zz_xxyz[i] = tr_x_0_xxyz[i] * fe_0 + tr_x_z_xxy[i] * fe_0 + tr_x_z_xxyz[i] * pa_z[i];

        tr_x_zz_xxzz[i] = tr_x_0_xxzz[i] * fe_0 + 2.0 * tr_x_z_xxz[i] * fe_0 + tr_x_z_xxzz[i] * pa_z[i];

        tr_x_zz_xyyy[i] = tr_x_0_xyyy[i] * fe_0 + tr_x_z_xyyy[i] * pa_z[i];

        tr_x_zz_xyyz[i] = tr_x_0_xyyz[i] * fe_0 + tr_x_z_xyy[i] * fe_0 + tr_x_z_xyyz[i] * pa_z[i];

        tr_x_zz_xyzz[i] = tr_x_0_xyzz[i] * fe_0 + 2.0 * tr_x_z_xyz[i] * fe_0 + tr_x_z_xyzz[i] * pa_z[i];

        tr_x_zz_xzzz[i] = tr_x_0_xzzz[i] * fe_0 + 3.0 * tr_x_z_xzz[i] * fe_0 + tr_x_z_xzzz[i] * pa_z[i];

        tr_x_zz_yyyy[i] = tr_x_0_yyyy[i] * fe_0 + tr_x_z_yyyy[i] * pa_z[i];

        tr_x_zz_yyyz[i] = tr_x_0_yyyz[i] * fe_0 + tr_x_z_yyy[i] * fe_0 + tr_x_z_yyyz[i] * pa_z[i];

        tr_x_zz_yyzz[i] = tr_x_0_yyzz[i] * fe_0 + 2.0 * tr_x_z_yyz[i] * fe_0 + tr_x_z_yyzz[i] * pa_z[i];

        tr_x_zz_yzzz[i] = tr_x_0_yzzz[i] * fe_0 + 3.0 * tr_x_z_yzz[i] * fe_0 + tr_x_z_yzzz[i] * pa_z[i];

        tr_x_zz_zzzz[i] = tr_x_0_zzzz[i] * fe_0 + 4.0 * tr_x_z_zzz[i] * fe_0 + tr_x_z_zzzz[i] * pa_z[i];
    }

    // Set up 90-105 components of targeted buffer : DG

    auto tr_y_xx_xxxx = pbuffer.data(idx_dip_dg + 90);

    auto tr_y_xx_xxxy = pbuffer.data(idx_dip_dg + 91);

    auto tr_y_xx_xxxz = pbuffer.data(idx_dip_dg + 92);

    auto tr_y_xx_xxyy = pbuffer.data(idx_dip_dg + 93);

    auto tr_y_xx_xxyz = pbuffer.data(idx_dip_dg + 94);

    auto tr_y_xx_xxzz = pbuffer.data(idx_dip_dg + 95);

    auto tr_y_xx_xyyy = pbuffer.data(idx_dip_dg + 96);

    auto tr_y_xx_xyyz = pbuffer.data(idx_dip_dg + 97);

    auto tr_y_xx_xyzz = pbuffer.data(idx_dip_dg + 98);

    auto tr_y_xx_xzzz = pbuffer.data(idx_dip_dg + 99);

    auto tr_y_xx_yyyy = pbuffer.data(idx_dip_dg + 100);

    auto tr_y_xx_yyyz = pbuffer.data(idx_dip_dg + 101);

    auto tr_y_xx_yyzz = pbuffer.data(idx_dip_dg + 102);

    auto tr_y_xx_yzzz = pbuffer.data(idx_dip_dg + 103);

    auto tr_y_xx_zzzz = pbuffer.data(idx_dip_dg + 104);

#pragma omp simd aligned(pa_x,             \
                             tr_y_0_xxxx,  \
                             tr_y_0_xxxy,  \
                             tr_y_0_xxxz,  \
                             tr_y_0_xxyy,  \
                             tr_y_0_xxyz,  \
                             tr_y_0_xxzz,  \
                             tr_y_0_xyyy,  \
                             tr_y_0_xyyz,  \
                             tr_y_0_xyzz,  \
                             tr_y_0_xzzz,  \
                             tr_y_0_yyyy,  \
                             tr_y_0_yyyz,  \
                             tr_y_0_yyzz,  \
                             tr_y_0_yzzz,  \
                             tr_y_0_zzzz,  \
                             tr_y_x_xxx,   \
                             tr_y_x_xxxx,  \
                             tr_y_x_xxxy,  \
                             tr_y_x_xxxz,  \
                             tr_y_x_xxy,   \
                             tr_y_x_xxyy,  \
                             tr_y_x_xxyz,  \
                             tr_y_x_xxz,   \
                             tr_y_x_xxzz,  \
                             tr_y_x_xyy,   \
                             tr_y_x_xyyy,  \
                             tr_y_x_xyyz,  \
                             tr_y_x_xyz,   \
                             tr_y_x_xyzz,  \
                             tr_y_x_xzz,   \
                             tr_y_x_xzzz,  \
                             tr_y_x_yyy,   \
                             tr_y_x_yyyy,  \
                             tr_y_x_yyyz,  \
                             tr_y_x_yyz,   \
                             tr_y_x_yyzz,  \
                             tr_y_x_yzz,   \
                             tr_y_x_yzzz,  \
                             tr_y_x_zzz,   \
                             tr_y_x_zzzz,  \
                             tr_y_xx_xxxx, \
                             tr_y_xx_xxxy, \
                             tr_y_xx_xxxz, \
                             tr_y_xx_xxyy, \
                             tr_y_xx_xxyz, \
                             tr_y_xx_xxzz, \
                             tr_y_xx_xyyy, \
                             tr_y_xx_xyyz, \
                             tr_y_xx_xyzz, \
                             tr_y_xx_xzzz, \
                             tr_y_xx_yyyy, \
                             tr_y_xx_yyyz, \
                             tr_y_xx_yyzz, \
                             tr_y_xx_yzzz, \
                             tr_y_xx_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xx_xxxx[i] = tr_y_0_xxxx[i] * fe_0 + 4.0 * tr_y_x_xxx[i] * fe_0 + tr_y_x_xxxx[i] * pa_x[i];

        tr_y_xx_xxxy[i] = tr_y_0_xxxy[i] * fe_0 + 3.0 * tr_y_x_xxy[i] * fe_0 + tr_y_x_xxxy[i] * pa_x[i];

        tr_y_xx_xxxz[i] = tr_y_0_xxxz[i] * fe_0 + 3.0 * tr_y_x_xxz[i] * fe_0 + tr_y_x_xxxz[i] * pa_x[i];

        tr_y_xx_xxyy[i] = tr_y_0_xxyy[i] * fe_0 + 2.0 * tr_y_x_xyy[i] * fe_0 + tr_y_x_xxyy[i] * pa_x[i];

        tr_y_xx_xxyz[i] = tr_y_0_xxyz[i] * fe_0 + 2.0 * tr_y_x_xyz[i] * fe_0 + tr_y_x_xxyz[i] * pa_x[i];

        tr_y_xx_xxzz[i] = tr_y_0_xxzz[i] * fe_0 + 2.0 * tr_y_x_xzz[i] * fe_0 + tr_y_x_xxzz[i] * pa_x[i];

        tr_y_xx_xyyy[i] = tr_y_0_xyyy[i] * fe_0 + tr_y_x_yyy[i] * fe_0 + tr_y_x_xyyy[i] * pa_x[i];

        tr_y_xx_xyyz[i] = tr_y_0_xyyz[i] * fe_0 + tr_y_x_yyz[i] * fe_0 + tr_y_x_xyyz[i] * pa_x[i];

        tr_y_xx_xyzz[i] = tr_y_0_xyzz[i] * fe_0 + tr_y_x_yzz[i] * fe_0 + tr_y_x_xyzz[i] * pa_x[i];

        tr_y_xx_xzzz[i] = tr_y_0_xzzz[i] * fe_0 + tr_y_x_zzz[i] * fe_0 + tr_y_x_xzzz[i] * pa_x[i];

        tr_y_xx_yyyy[i] = tr_y_0_yyyy[i] * fe_0 + tr_y_x_yyyy[i] * pa_x[i];

        tr_y_xx_yyyz[i] = tr_y_0_yyyz[i] * fe_0 + tr_y_x_yyyz[i] * pa_x[i];

        tr_y_xx_yyzz[i] = tr_y_0_yyzz[i] * fe_0 + tr_y_x_yyzz[i] * pa_x[i];

        tr_y_xx_yzzz[i] = tr_y_0_yzzz[i] * fe_0 + tr_y_x_yzzz[i] * pa_x[i];

        tr_y_xx_zzzz[i] = tr_y_0_zzzz[i] * fe_0 + tr_y_x_zzzz[i] * pa_x[i];
    }

    // Set up 105-120 components of targeted buffer : DG

    auto tr_y_xy_xxxx = pbuffer.data(idx_dip_dg + 105);

    auto tr_y_xy_xxxy = pbuffer.data(idx_dip_dg + 106);

    auto tr_y_xy_xxxz = pbuffer.data(idx_dip_dg + 107);

    auto tr_y_xy_xxyy = pbuffer.data(idx_dip_dg + 108);

    auto tr_y_xy_xxyz = pbuffer.data(idx_dip_dg + 109);

    auto tr_y_xy_xxzz = pbuffer.data(idx_dip_dg + 110);

    auto tr_y_xy_xyyy = pbuffer.data(idx_dip_dg + 111);

    auto tr_y_xy_xyyz = pbuffer.data(idx_dip_dg + 112);

    auto tr_y_xy_xyzz = pbuffer.data(idx_dip_dg + 113);

    auto tr_y_xy_xzzz = pbuffer.data(idx_dip_dg + 114);

    auto tr_y_xy_yyyy = pbuffer.data(idx_dip_dg + 115);

    auto tr_y_xy_yyyz = pbuffer.data(idx_dip_dg + 116);

    auto tr_y_xy_yyzz = pbuffer.data(idx_dip_dg + 117);

    auto tr_y_xy_yzzz = pbuffer.data(idx_dip_dg + 118);

    auto tr_y_xy_zzzz = pbuffer.data(idx_dip_dg + 119);

#pragma omp simd aligned(pa_x,             \
                             tr_y_xy_xxxx, \
                             tr_y_xy_xxxy, \
                             tr_y_xy_xxxz, \
                             tr_y_xy_xxyy, \
                             tr_y_xy_xxyz, \
                             tr_y_xy_xxzz, \
                             tr_y_xy_xyyy, \
                             tr_y_xy_xyyz, \
                             tr_y_xy_xyzz, \
                             tr_y_xy_xzzz, \
                             tr_y_xy_yyyy, \
                             tr_y_xy_yyyz, \
                             tr_y_xy_yyzz, \
                             tr_y_xy_yzzz, \
                             tr_y_xy_zzzz, \
                             tr_y_y_xxx,   \
                             tr_y_y_xxxx,  \
                             tr_y_y_xxxy,  \
                             tr_y_y_xxxz,  \
                             tr_y_y_xxy,   \
                             tr_y_y_xxyy,  \
                             tr_y_y_xxyz,  \
                             tr_y_y_xxz,   \
                             tr_y_y_xxzz,  \
                             tr_y_y_xyy,   \
                             tr_y_y_xyyy,  \
                             tr_y_y_xyyz,  \
                             tr_y_y_xyz,   \
                             tr_y_y_xyzz,  \
                             tr_y_y_xzz,   \
                             tr_y_y_xzzz,  \
                             tr_y_y_yyy,   \
                             tr_y_y_yyyy,  \
                             tr_y_y_yyyz,  \
                             tr_y_y_yyz,   \
                             tr_y_y_yyzz,  \
                             tr_y_y_yzz,   \
                             tr_y_y_yzzz,  \
                             tr_y_y_zzz,   \
                             tr_y_y_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xy_xxxx[i] = 4.0 * tr_y_y_xxx[i] * fe_0 + tr_y_y_xxxx[i] * pa_x[i];

        tr_y_xy_xxxy[i] = 3.0 * tr_y_y_xxy[i] * fe_0 + tr_y_y_xxxy[i] * pa_x[i];

        tr_y_xy_xxxz[i] = 3.0 * tr_y_y_xxz[i] * fe_0 + tr_y_y_xxxz[i] * pa_x[i];

        tr_y_xy_xxyy[i] = 2.0 * tr_y_y_xyy[i] * fe_0 + tr_y_y_xxyy[i] * pa_x[i];

        tr_y_xy_xxyz[i] = 2.0 * tr_y_y_xyz[i] * fe_0 + tr_y_y_xxyz[i] * pa_x[i];

        tr_y_xy_xxzz[i] = 2.0 * tr_y_y_xzz[i] * fe_0 + tr_y_y_xxzz[i] * pa_x[i];

        tr_y_xy_xyyy[i] = tr_y_y_yyy[i] * fe_0 + tr_y_y_xyyy[i] * pa_x[i];

        tr_y_xy_xyyz[i] = tr_y_y_yyz[i] * fe_0 + tr_y_y_xyyz[i] * pa_x[i];

        tr_y_xy_xyzz[i] = tr_y_y_yzz[i] * fe_0 + tr_y_y_xyzz[i] * pa_x[i];

        tr_y_xy_xzzz[i] = tr_y_y_zzz[i] * fe_0 + tr_y_y_xzzz[i] * pa_x[i];

        tr_y_xy_yyyy[i] = tr_y_y_yyyy[i] * pa_x[i];

        tr_y_xy_yyyz[i] = tr_y_y_yyyz[i] * pa_x[i];

        tr_y_xy_yyzz[i] = tr_y_y_yyzz[i] * pa_x[i];

        tr_y_xy_yzzz[i] = tr_y_y_yzzz[i] * pa_x[i];

        tr_y_xy_zzzz[i] = tr_y_y_zzzz[i] * pa_x[i];
    }

    // Set up 120-135 components of targeted buffer : DG

    auto tr_y_xz_xxxx = pbuffer.data(idx_dip_dg + 120);

    auto tr_y_xz_xxxy = pbuffer.data(idx_dip_dg + 121);

    auto tr_y_xz_xxxz = pbuffer.data(idx_dip_dg + 122);

    auto tr_y_xz_xxyy = pbuffer.data(idx_dip_dg + 123);

    auto tr_y_xz_xxyz = pbuffer.data(idx_dip_dg + 124);

    auto tr_y_xz_xxzz = pbuffer.data(idx_dip_dg + 125);

    auto tr_y_xz_xyyy = pbuffer.data(idx_dip_dg + 126);

    auto tr_y_xz_xyyz = pbuffer.data(idx_dip_dg + 127);

    auto tr_y_xz_xyzz = pbuffer.data(idx_dip_dg + 128);

    auto tr_y_xz_xzzz = pbuffer.data(idx_dip_dg + 129);

    auto tr_y_xz_yyyy = pbuffer.data(idx_dip_dg + 130);

    auto tr_y_xz_yyyz = pbuffer.data(idx_dip_dg + 131);

    auto tr_y_xz_yyzz = pbuffer.data(idx_dip_dg + 132);

    auto tr_y_xz_yzzz = pbuffer.data(idx_dip_dg + 133);

    auto tr_y_xz_zzzz = pbuffer.data(idx_dip_dg + 134);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             tr_y_x_xxxx,  \
                             tr_y_x_xxxy,  \
                             tr_y_x_xxyy,  \
                             tr_y_x_xyyy,  \
                             tr_y_xz_xxxx, \
                             tr_y_xz_xxxy, \
                             tr_y_xz_xxxz, \
                             tr_y_xz_xxyy, \
                             tr_y_xz_xxyz, \
                             tr_y_xz_xxzz, \
                             tr_y_xz_xyyy, \
                             tr_y_xz_xyyz, \
                             tr_y_xz_xyzz, \
                             tr_y_xz_xzzz, \
                             tr_y_xz_yyyy, \
                             tr_y_xz_yyyz, \
                             tr_y_xz_yyzz, \
                             tr_y_xz_yzzz, \
                             tr_y_xz_zzzz, \
                             tr_y_z_xxxz,  \
                             tr_y_z_xxyz,  \
                             tr_y_z_xxz,   \
                             tr_y_z_xxzz,  \
                             tr_y_z_xyyz,  \
                             tr_y_z_xyz,   \
                             tr_y_z_xyzz,  \
                             tr_y_z_xzz,   \
                             tr_y_z_xzzz,  \
                             tr_y_z_yyyy,  \
                             tr_y_z_yyyz,  \
                             tr_y_z_yyz,   \
                             tr_y_z_yyzz,  \
                             tr_y_z_yzz,   \
                             tr_y_z_yzzz,  \
                             tr_y_z_zzz,   \
                             tr_y_z_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xz_xxxx[i] = tr_y_x_xxxx[i] * pa_z[i];

        tr_y_xz_xxxy[i] = tr_y_x_xxxy[i] * pa_z[i];

        tr_y_xz_xxxz[i] = 3.0 * tr_y_z_xxz[i] * fe_0 + tr_y_z_xxxz[i] * pa_x[i];

        tr_y_xz_xxyy[i] = tr_y_x_xxyy[i] * pa_z[i];

        tr_y_xz_xxyz[i] = 2.0 * tr_y_z_xyz[i] * fe_0 + tr_y_z_xxyz[i] * pa_x[i];

        tr_y_xz_xxzz[i] = 2.0 * tr_y_z_xzz[i] * fe_0 + tr_y_z_xxzz[i] * pa_x[i];

        tr_y_xz_xyyy[i] = tr_y_x_xyyy[i] * pa_z[i];

        tr_y_xz_xyyz[i] = tr_y_z_yyz[i] * fe_0 + tr_y_z_xyyz[i] * pa_x[i];

        tr_y_xz_xyzz[i] = tr_y_z_yzz[i] * fe_0 + tr_y_z_xyzz[i] * pa_x[i];

        tr_y_xz_xzzz[i] = tr_y_z_zzz[i] * fe_0 + tr_y_z_xzzz[i] * pa_x[i];

        tr_y_xz_yyyy[i] = tr_y_z_yyyy[i] * pa_x[i];

        tr_y_xz_yyyz[i] = tr_y_z_yyyz[i] * pa_x[i];

        tr_y_xz_yyzz[i] = tr_y_z_yyzz[i] * pa_x[i];

        tr_y_xz_yzzz[i] = tr_y_z_yzzz[i] * pa_x[i];

        tr_y_xz_zzzz[i] = tr_y_z_zzzz[i] * pa_x[i];
    }

    // Set up 135-150 components of targeted buffer : DG

    auto tr_y_yy_xxxx = pbuffer.data(idx_dip_dg + 135);

    auto tr_y_yy_xxxy = pbuffer.data(idx_dip_dg + 136);

    auto tr_y_yy_xxxz = pbuffer.data(idx_dip_dg + 137);

    auto tr_y_yy_xxyy = pbuffer.data(idx_dip_dg + 138);

    auto tr_y_yy_xxyz = pbuffer.data(idx_dip_dg + 139);

    auto tr_y_yy_xxzz = pbuffer.data(idx_dip_dg + 140);

    auto tr_y_yy_xyyy = pbuffer.data(idx_dip_dg + 141);

    auto tr_y_yy_xyyz = pbuffer.data(idx_dip_dg + 142);

    auto tr_y_yy_xyzz = pbuffer.data(idx_dip_dg + 143);

    auto tr_y_yy_xzzz = pbuffer.data(idx_dip_dg + 144);

    auto tr_y_yy_yyyy = pbuffer.data(idx_dip_dg + 145);

    auto tr_y_yy_yyyz = pbuffer.data(idx_dip_dg + 146);

    auto tr_y_yy_yyzz = pbuffer.data(idx_dip_dg + 147);

    auto tr_y_yy_yzzz = pbuffer.data(idx_dip_dg + 148);

    auto tr_y_yy_zzzz = pbuffer.data(idx_dip_dg + 149);

#pragma omp simd aligned(pa_y,             \
                             tr_y_0_xxxx,  \
                             tr_y_0_xxxy,  \
                             tr_y_0_xxxz,  \
                             tr_y_0_xxyy,  \
                             tr_y_0_xxyz,  \
                             tr_y_0_xxzz,  \
                             tr_y_0_xyyy,  \
                             tr_y_0_xyyz,  \
                             tr_y_0_xyzz,  \
                             tr_y_0_xzzz,  \
                             tr_y_0_yyyy,  \
                             tr_y_0_yyyz,  \
                             tr_y_0_yyzz,  \
                             tr_y_0_yzzz,  \
                             tr_y_0_zzzz,  \
                             tr_y_y_xxx,   \
                             tr_y_y_xxxx,  \
                             tr_y_y_xxxy,  \
                             tr_y_y_xxxz,  \
                             tr_y_y_xxy,   \
                             tr_y_y_xxyy,  \
                             tr_y_y_xxyz,  \
                             tr_y_y_xxz,   \
                             tr_y_y_xxzz,  \
                             tr_y_y_xyy,   \
                             tr_y_y_xyyy,  \
                             tr_y_y_xyyz,  \
                             tr_y_y_xyz,   \
                             tr_y_y_xyzz,  \
                             tr_y_y_xzz,   \
                             tr_y_y_xzzz,  \
                             tr_y_y_yyy,   \
                             tr_y_y_yyyy,  \
                             tr_y_y_yyyz,  \
                             tr_y_y_yyz,   \
                             tr_y_y_yyzz,  \
                             tr_y_y_yzz,   \
                             tr_y_y_yzzz,  \
                             tr_y_y_zzz,   \
                             tr_y_y_zzzz,  \
                             tr_y_yy_xxxx, \
                             tr_y_yy_xxxy, \
                             tr_y_yy_xxxz, \
                             tr_y_yy_xxyy, \
                             tr_y_yy_xxyz, \
                             tr_y_yy_xxzz, \
                             tr_y_yy_xyyy, \
                             tr_y_yy_xyyz, \
                             tr_y_yy_xyzz, \
                             tr_y_yy_xzzz, \
                             tr_y_yy_yyyy, \
                             tr_y_yy_yyyz, \
                             tr_y_yy_yyzz, \
                             tr_y_yy_yzzz, \
                             tr_y_yy_zzzz, \
                             ts_y_xxxx,    \
                             ts_y_xxxy,    \
                             ts_y_xxxz,    \
                             ts_y_xxyy,    \
                             ts_y_xxyz,    \
                             ts_y_xxzz,    \
                             ts_y_xyyy,    \
                             ts_y_xyyz,    \
                             ts_y_xyzz,    \
                             ts_y_xzzz,    \
                             ts_y_yyyy,    \
                             ts_y_yyyz,    \
                             ts_y_yyzz,    \
                             ts_y_yzzz,    \
                             ts_y_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yy_xxxx[i] = tr_y_0_xxxx[i] * fe_0 + ts_y_xxxx[i] * fe_0 + tr_y_y_xxxx[i] * pa_y[i];

        tr_y_yy_xxxy[i] = tr_y_0_xxxy[i] * fe_0 + tr_y_y_xxx[i] * fe_0 + ts_y_xxxy[i] * fe_0 + tr_y_y_xxxy[i] * pa_y[i];

        tr_y_yy_xxxz[i] = tr_y_0_xxxz[i] * fe_0 + ts_y_xxxz[i] * fe_0 + tr_y_y_xxxz[i] * pa_y[i];

        tr_y_yy_xxyy[i] = tr_y_0_xxyy[i] * fe_0 + 2.0 * tr_y_y_xxy[i] * fe_0 + ts_y_xxyy[i] * fe_0 + tr_y_y_xxyy[i] * pa_y[i];

        tr_y_yy_xxyz[i] = tr_y_0_xxyz[i] * fe_0 + tr_y_y_xxz[i] * fe_0 + ts_y_xxyz[i] * fe_0 + tr_y_y_xxyz[i] * pa_y[i];

        tr_y_yy_xxzz[i] = tr_y_0_xxzz[i] * fe_0 + ts_y_xxzz[i] * fe_0 + tr_y_y_xxzz[i] * pa_y[i];

        tr_y_yy_xyyy[i] = tr_y_0_xyyy[i] * fe_0 + 3.0 * tr_y_y_xyy[i] * fe_0 + ts_y_xyyy[i] * fe_0 + tr_y_y_xyyy[i] * pa_y[i];

        tr_y_yy_xyyz[i] = tr_y_0_xyyz[i] * fe_0 + 2.0 * tr_y_y_xyz[i] * fe_0 + ts_y_xyyz[i] * fe_0 + tr_y_y_xyyz[i] * pa_y[i];

        tr_y_yy_xyzz[i] = tr_y_0_xyzz[i] * fe_0 + tr_y_y_xzz[i] * fe_0 + ts_y_xyzz[i] * fe_0 + tr_y_y_xyzz[i] * pa_y[i];

        tr_y_yy_xzzz[i] = tr_y_0_xzzz[i] * fe_0 + ts_y_xzzz[i] * fe_0 + tr_y_y_xzzz[i] * pa_y[i];

        tr_y_yy_yyyy[i] = tr_y_0_yyyy[i] * fe_0 + 4.0 * tr_y_y_yyy[i] * fe_0 + ts_y_yyyy[i] * fe_0 + tr_y_y_yyyy[i] * pa_y[i];

        tr_y_yy_yyyz[i] = tr_y_0_yyyz[i] * fe_0 + 3.0 * tr_y_y_yyz[i] * fe_0 + ts_y_yyyz[i] * fe_0 + tr_y_y_yyyz[i] * pa_y[i];

        tr_y_yy_yyzz[i] = tr_y_0_yyzz[i] * fe_0 + 2.0 * tr_y_y_yzz[i] * fe_0 + ts_y_yyzz[i] * fe_0 + tr_y_y_yyzz[i] * pa_y[i];

        tr_y_yy_yzzz[i] = tr_y_0_yzzz[i] * fe_0 + tr_y_y_zzz[i] * fe_0 + ts_y_yzzz[i] * fe_0 + tr_y_y_yzzz[i] * pa_y[i];

        tr_y_yy_zzzz[i] = tr_y_0_zzzz[i] * fe_0 + ts_y_zzzz[i] * fe_0 + tr_y_y_zzzz[i] * pa_y[i];
    }

    // Set up 150-165 components of targeted buffer : DG

    auto tr_y_yz_xxxx = pbuffer.data(idx_dip_dg + 150);

    auto tr_y_yz_xxxy = pbuffer.data(idx_dip_dg + 151);

    auto tr_y_yz_xxxz = pbuffer.data(idx_dip_dg + 152);

    auto tr_y_yz_xxyy = pbuffer.data(idx_dip_dg + 153);

    auto tr_y_yz_xxyz = pbuffer.data(idx_dip_dg + 154);

    auto tr_y_yz_xxzz = pbuffer.data(idx_dip_dg + 155);

    auto tr_y_yz_xyyy = pbuffer.data(idx_dip_dg + 156);

    auto tr_y_yz_xyyz = pbuffer.data(idx_dip_dg + 157);

    auto tr_y_yz_xyzz = pbuffer.data(idx_dip_dg + 158);

    auto tr_y_yz_xzzz = pbuffer.data(idx_dip_dg + 159);

    auto tr_y_yz_yyyy = pbuffer.data(idx_dip_dg + 160);

    auto tr_y_yz_yyyz = pbuffer.data(idx_dip_dg + 161);

    auto tr_y_yz_yyzz = pbuffer.data(idx_dip_dg + 162);

    auto tr_y_yz_yzzz = pbuffer.data(idx_dip_dg + 163);

    auto tr_y_yz_zzzz = pbuffer.data(idx_dip_dg + 164);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             tr_y_y_xxxx,  \
                             tr_y_y_xxxy,  \
                             tr_y_y_xxy,   \
                             tr_y_y_xxyy,  \
                             tr_y_y_xxyz,  \
                             tr_y_y_xyy,   \
                             tr_y_y_xyyy,  \
                             tr_y_y_xyyz,  \
                             tr_y_y_xyz,   \
                             tr_y_y_xyzz,  \
                             tr_y_y_yyy,   \
                             tr_y_y_yyyy,  \
                             tr_y_y_yyyz,  \
                             tr_y_y_yyz,   \
                             tr_y_y_yyzz,  \
                             tr_y_y_yzz,   \
                             tr_y_y_yzzz,  \
                             tr_y_yz_xxxx, \
                             tr_y_yz_xxxy, \
                             tr_y_yz_xxxz, \
                             tr_y_yz_xxyy, \
                             tr_y_yz_xxyz, \
                             tr_y_yz_xxzz, \
                             tr_y_yz_xyyy, \
                             tr_y_yz_xyyz, \
                             tr_y_yz_xyzz, \
                             tr_y_yz_xzzz, \
                             tr_y_yz_yyyy, \
                             tr_y_yz_yyyz, \
                             tr_y_yz_yyzz, \
                             tr_y_yz_yzzz, \
                             tr_y_yz_zzzz, \
                             tr_y_z_xxxz,  \
                             tr_y_z_xxzz,  \
                             tr_y_z_xzzz,  \
                             tr_y_z_zzzz,  \
                             ts_z_xxxz,    \
                             ts_z_xxzz,    \
                             ts_z_xzzz,    \
                             ts_z_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yz_xxxx[i] = tr_y_y_xxxx[i] * pa_z[i];

        tr_y_yz_xxxy[i] = tr_y_y_xxxy[i] * pa_z[i];

        tr_y_yz_xxxz[i] = ts_z_xxxz[i] * fe_0 + tr_y_z_xxxz[i] * pa_y[i];

        tr_y_yz_xxyy[i] = tr_y_y_xxyy[i] * pa_z[i];

        tr_y_yz_xxyz[i] = tr_y_y_xxy[i] * fe_0 + tr_y_y_xxyz[i] * pa_z[i];

        tr_y_yz_xxzz[i] = ts_z_xxzz[i] * fe_0 + tr_y_z_xxzz[i] * pa_y[i];

        tr_y_yz_xyyy[i] = tr_y_y_xyyy[i] * pa_z[i];

        tr_y_yz_xyyz[i] = tr_y_y_xyy[i] * fe_0 + tr_y_y_xyyz[i] * pa_z[i];

        tr_y_yz_xyzz[i] = 2.0 * tr_y_y_xyz[i] * fe_0 + tr_y_y_xyzz[i] * pa_z[i];

        tr_y_yz_xzzz[i] = ts_z_xzzz[i] * fe_0 + tr_y_z_xzzz[i] * pa_y[i];

        tr_y_yz_yyyy[i] = tr_y_y_yyyy[i] * pa_z[i];

        tr_y_yz_yyyz[i] = tr_y_y_yyy[i] * fe_0 + tr_y_y_yyyz[i] * pa_z[i];

        tr_y_yz_yyzz[i] = 2.0 * tr_y_y_yyz[i] * fe_0 + tr_y_y_yyzz[i] * pa_z[i];

        tr_y_yz_yzzz[i] = 3.0 * tr_y_y_yzz[i] * fe_0 + tr_y_y_yzzz[i] * pa_z[i];

        tr_y_yz_zzzz[i] = ts_z_zzzz[i] * fe_0 + tr_y_z_zzzz[i] * pa_y[i];
    }

    // Set up 165-180 components of targeted buffer : DG

    auto tr_y_zz_xxxx = pbuffer.data(idx_dip_dg + 165);

    auto tr_y_zz_xxxy = pbuffer.data(idx_dip_dg + 166);

    auto tr_y_zz_xxxz = pbuffer.data(idx_dip_dg + 167);

    auto tr_y_zz_xxyy = pbuffer.data(idx_dip_dg + 168);

    auto tr_y_zz_xxyz = pbuffer.data(idx_dip_dg + 169);

    auto tr_y_zz_xxzz = pbuffer.data(idx_dip_dg + 170);

    auto tr_y_zz_xyyy = pbuffer.data(idx_dip_dg + 171);

    auto tr_y_zz_xyyz = pbuffer.data(idx_dip_dg + 172);

    auto tr_y_zz_xyzz = pbuffer.data(idx_dip_dg + 173);

    auto tr_y_zz_xzzz = pbuffer.data(idx_dip_dg + 174);

    auto tr_y_zz_yyyy = pbuffer.data(idx_dip_dg + 175);

    auto tr_y_zz_yyyz = pbuffer.data(idx_dip_dg + 176);

    auto tr_y_zz_yyzz = pbuffer.data(idx_dip_dg + 177);

    auto tr_y_zz_yzzz = pbuffer.data(idx_dip_dg + 178);

    auto tr_y_zz_zzzz = pbuffer.data(idx_dip_dg + 179);

#pragma omp simd aligned(pa_z,             \
                             tr_y_0_xxxx,  \
                             tr_y_0_xxxy,  \
                             tr_y_0_xxxz,  \
                             tr_y_0_xxyy,  \
                             tr_y_0_xxyz,  \
                             tr_y_0_xxzz,  \
                             tr_y_0_xyyy,  \
                             tr_y_0_xyyz,  \
                             tr_y_0_xyzz,  \
                             tr_y_0_xzzz,  \
                             tr_y_0_yyyy,  \
                             tr_y_0_yyyz,  \
                             tr_y_0_yyzz,  \
                             tr_y_0_yzzz,  \
                             tr_y_0_zzzz,  \
                             tr_y_z_xxx,   \
                             tr_y_z_xxxx,  \
                             tr_y_z_xxxy,  \
                             tr_y_z_xxxz,  \
                             tr_y_z_xxy,   \
                             tr_y_z_xxyy,  \
                             tr_y_z_xxyz,  \
                             tr_y_z_xxz,   \
                             tr_y_z_xxzz,  \
                             tr_y_z_xyy,   \
                             tr_y_z_xyyy,  \
                             tr_y_z_xyyz,  \
                             tr_y_z_xyz,   \
                             tr_y_z_xyzz,  \
                             tr_y_z_xzz,   \
                             tr_y_z_xzzz,  \
                             tr_y_z_yyy,   \
                             tr_y_z_yyyy,  \
                             tr_y_z_yyyz,  \
                             tr_y_z_yyz,   \
                             tr_y_z_yyzz,  \
                             tr_y_z_yzz,   \
                             tr_y_z_yzzz,  \
                             tr_y_z_zzz,   \
                             tr_y_z_zzzz,  \
                             tr_y_zz_xxxx, \
                             tr_y_zz_xxxy, \
                             tr_y_zz_xxxz, \
                             tr_y_zz_xxyy, \
                             tr_y_zz_xxyz, \
                             tr_y_zz_xxzz, \
                             tr_y_zz_xyyy, \
                             tr_y_zz_xyyz, \
                             tr_y_zz_xyzz, \
                             tr_y_zz_xzzz, \
                             tr_y_zz_yyyy, \
                             tr_y_zz_yyyz, \
                             tr_y_zz_yyzz, \
                             tr_y_zz_yzzz, \
                             tr_y_zz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_zz_xxxx[i] = tr_y_0_xxxx[i] * fe_0 + tr_y_z_xxxx[i] * pa_z[i];

        tr_y_zz_xxxy[i] = tr_y_0_xxxy[i] * fe_0 + tr_y_z_xxxy[i] * pa_z[i];

        tr_y_zz_xxxz[i] = tr_y_0_xxxz[i] * fe_0 + tr_y_z_xxx[i] * fe_0 + tr_y_z_xxxz[i] * pa_z[i];

        tr_y_zz_xxyy[i] = tr_y_0_xxyy[i] * fe_0 + tr_y_z_xxyy[i] * pa_z[i];

        tr_y_zz_xxyz[i] = tr_y_0_xxyz[i] * fe_0 + tr_y_z_xxy[i] * fe_0 + tr_y_z_xxyz[i] * pa_z[i];

        tr_y_zz_xxzz[i] = tr_y_0_xxzz[i] * fe_0 + 2.0 * tr_y_z_xxz[i] * fe_0 + tr_y_z_xxzz[i] * pa_z[i];

        tr_y_zz_xyyy[i] = tr_y_0_xyyy[i] * fe_0 + tr_y_z_xyyy[i] * pa_z[i];

        tr_y_zz_xyyz[i] = tr_y_0_xyyz[i] * fe_0 + tr_y_z_xyy[i] * fe_0 + tr_y_z_xyyz[i] * pa_z[i];

        tr_y_zz_xyzz[i] = tr_y_0_xyzz[i] * fe_0 + 2.0 * tr_y_z_xyz[i] * fe_0 + tr_y_z_xyzz[i] * pa_z[i];

        tr_y_zz_xzzz[i] = tr_y_0_xzzz[i] * fe_0 + 3.0 * tr_y_z_xzz[i] * fe_0 + tr_y_z_xzzz[i] * pa_z[i];

        tr_y_zz_yyyy[i] = tr_y_0_yyyy[i] * fe_0 + tr_y_z_yyyy[i] * pa_z[i];

        tr_y_zz_yyyz[i] = tr_y_0_yyyz[i] * fe_0 + tr_y_z_yyy[i] * fe_0 + tr_y_z_yyyz[i] * pa_z[i];

        tr_y_zz_yyzz[i] = tr_y_0_yyzz[i] * fe_0 + 2.0 * tr_y_z_yyz[i] * fe_0 + tr_y_z_yyzz[i] * pa_z[i];

        tr_y_zz_yzzz[i] = tr_y_0_yzzz[i] * fe_0 + 3.0 * tr_y_z_yzz[i] * fe_0 + tr_y_z_yzzz[i] * pa_z[i];

        tr_y_zz_zzzz[i] = tr_y_0_zzzz[i] * fe_0 + 4.0 * tr_y_z_zzz[i] * fe_0 + tr_y_z_zzzz[i] * pa_z[i];
    }

    // Set up 180-195 components of targeted buffer : DG

    auto tr_z_xx_xxxx = pbuffer.data(idx_dip_dg + 180);

    auto tr_z_xx_xxxy = pbuffer.data(idx_dip_dg + 181);

    auto tr_z_xx_xxxz = pbuffer.data(idx_dip_dg + 182);

    auto tr_z_xx_xxyy = pbuffer.data(idx_dip_dg + 183);

    auto tr_z_xx_xxyz = pbuffer.data(idx_dip_dg + 184);

    auto tr_z_xx_xxzz = pbuffer.data(idx_dip_dg + 185);

    auto tr_z_xx_xyyy = pbuffer.data(idx_dip_dg + 186);

    auto tr_z_xx_xyyz = pbuffer.data(idx_dip_dg + 187);

    auto tr_z_xx_xyzz = pbuffer.data(idx_dip_dg + 188);

    auto tr_z_xx_xzzz = pbuffer.data(idx_dip_dg + 189);

    auto tr_z_xx_yyyy = pbuffer.data(idx_dip_dg + 190);

    auto tr_z_xx_yyyz = pbuffer.data(idx_dip_dg + 191);

    auto tr_z_xx_yyzz = pbuffer.data(idx_dip_dg + 192);

    auto tr_z_xx_yzzz = pbuffer.data(idx_dip_dg + 193);

    auto tr_z_xx_zzzz = pbuffer.data(idx_dip_dg + 194);

#pragma omp simd aligned(pa_x,             \
                             tr_z_0_xxxx,  \
                             tr_z_0_xxxy,  \
                             tr_z_0_xxxz,  \
                             tr_z_0_xxyy,  \
                             tr_z_0_xxyz,  \
                             tr_z_0_xxzz,  \
                             tr_z_0_xyyy,  \
                             tr_z_0_xyyz,  \
                             tr_z_0_xyzz,  \
                             tr_z_0_xzzz,  \
                             tr_z_0_yyyy,  \
                             tr_z_0_yyyz,  \
                             tr_z_0_yyzz,  \
                             tr_z_0_yzzz,  \
                             tr_z_0_zzzz,  \
                             tr_z_x_xxx,   \
                             tr_z_x_xxxx,  \
                             tr_z_x_xxxy,  \
                             tr_z_x_xxxz,  \
                             tr_z_x_xxy,   \
                             tr_z_x_xxyy,  \
                             tr_z_x_xxyz,  \
                             tr_z_x_xxz,   \
                             tr_z_x_xxzz,  \
                             tr_z_x_xyy,   \
                             tr_z_x_xyyy,  \
                             tr_z_x_xyyz,  \
                             tr_z_x_xyz,   \
                             tr_z_x_xyzz,  \
                             tr_z_x_xzz,   \
                             tr_z_x_xzzz,  \
                             tr_z_x_yyy,   \
                             tr_z_x_yyyy,  \
                             tr_z_x_yyyz,  \
                             tr_z_x_yyz,   \
                             tr_z_x_yyzz,  \
                             tr_z_x_yzz,   \
                             tr_z_x_yzzz,  \
                             tr_z_x_zzz,   \
                             tr_z_x_zzzz,  \
                             tr_z_xx_xxxx, \
                             tr_z_xx_xxxy, \
                             tr_z_xx_xxxz, \
                             tr_z_xx_xxyy, \
                             tr_z_xx_xxyz, \
                             tr_z_xx_xxzz, \
                             tr_z_xx_xyyy, \
                             tr_z_xx_xyyz, \
                             tr_z_xx_xyzz, \
                             tr_z_xx_xzzz, \
                             tr_z_xx_yyyy, \
                             tr_z_xx_yyyz, \
                             tr_z_xx_yyzz, \
                             tr_z_xx_yzzz, \
                             tr_z_xx_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xx_xxxx[i] = tr_z_0_xxxx[i] * fe_0 + 4.0 * tr_z_x_xxx[i] * fe_0 + tr_z_x_xxxx[i] * pa_x[i];

        tr_z_xx_xxxy[i] = tr_z_0_xxxy[i] * fe_0 + 3.0 * tr_z_x_xxy[i] * fe_0 + tr_z_x_xxxy[i] * pa_x[i];

        tr_z_xx_xxxz[i] = tr_z_0_xxxz[i] * fe_0 + 3.0 * tr_z_x_xxz[i] * fe_0 + tr_z_x_xxxz[i] * pa_x[i];

        tr_z_xx_xxyy[i] = tr_z_0_xxyy[i] * fe_0 + 2.0 * tr_z_x_xyy[i] * fe_0 + tr_z_x_xxyy[i] * pa_x[i];

        tr_z_xx_xxyz[i] = tr_z_0_xxyz[i] * fe_0 + 2.0 * tr_z_x_xyz[i] * fe_0 + tr_z_x_xxyz[i] * pa_x[i];

        tr_z_xx_xxzz[i] = tr_z_0_xxzz[i] * fe_0 + 2.0 * tr_z_x_xzz[i] * fe_0 + tr_z_x_xxzz[i] * pa_x[i];

        tr_z_xx_xyyy[i] = tr_z_0_xyyy[i] * fe_0 + tr_z_x_yyy[i] * fe_0 + tr_z_x_xyyy[i] * pa_x[i];

        tr_z_xx_xyyz[i] = tr_z_0_xyyz[i] * fe_0 + tr_z_x_yyz[i] * fe_0 + tr_z_x_xyyz[i] * pa_x[i];

        tr_z_xx_xyzz[i] = tr_z_0_xyzz[i] * fe_0 + tr_z_x_yzz[i] * fe_0 + tr_z_x_xyzz[i] * pa_x[i];

        tr_z_xx_xzzz[i] = tr_z_0_xzzz[i] * fe_0 + tr_z_x_zzz[i] * fe_0 + tr_z_x_xzzz[i] * pa_x[i];

        tr_z_xx_yyyy[i] = tr_z_0_yyyy[i] * fe_0 + tr_z_x_yyyy[i] * pa_x[i];

        tr_z_xx_yyyz[i] = tr_z_0_yyyz[i] * fe_0 + tr_z_x_yyyz[i] * pa_x[i];

        tr_z_xx_yyzz[i] = tr_z_0_yyzz[i] * fe_0 + tr_z_x_yyzz[i] * pa_x[i];

        tr_z_xx_yzzz[i] = tr_z_0_yzzz[i] * fe_0 + tr_z_x_yzzz[i] * pa_x[i];

        tr_z_xx_zzzz[i] = tr_z_0_zzzz[i] * fe_0 + tr_z_x_zzzz[i] * pa_x[i];
    }

    // Set up 195-210 components of targeted buffer : DG

    auto tr_z_xy_xxxx = pbuffer.data(idx_dip_dg + 195);

    auto tr_z_xy_xxxy = pbuffer.data(idx_dip_dg + 196);

    auto tr_z_xy_xxxz = pbuffer.data(idx_dip_dg + 197);

    auto tr_z_xy_xxyy = pbuffer.data(idx_dip_dg + 198);

    auto tr_z_xy_xxyz = pbuffer.data(idx_dip_dg + 199);

    auto tr_z_xy_xxzz = pbuffer.data(idx_dip_dg + 200);

    auto tr_z_xy_xyyy = pbuffer.data(idx_dip_dg + 201);

    auto tr_z_xy_xyyz = pbuffer.data(idx_dip_dg + 202);

    auto tr_z_xy_xyzz = pbuffer.data(idx_dip_dg + 203);

    auto tr_z_xy_xzzz = pbuffer.data(idx_dip_dg + 204);

    auto tr_z_xy_yyyy = pbuffer.data(idx_dip_dg + 205);

    auto tr_z_xy_yyyz = pbuffer.data(idx_dip_dg + 206);

    auto tr_z_xy_yyzz = pbuffer.data(idx_dip_dg + 207);

    auto tr_z_xy_yzzz = pbuffer.data(idx_dip_dg + 208);

    auto tr_z_xy_zzzz = pbuffer.data(idx_dip_dg + 209);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             tr_z_x_xxxx,  \
                             tr_z_x_xxxz,  \
                             tr_z_x_xxzz,  \
                             tr_z_x_xzzz,  \
                             tr_z_xy_xxxx, \
                             tr_z_xy_xxxy, \
                             tr_z_xy_xxxz, \
                             tr_z_xy_xxyy, \
                             tr_z_xy_xxyz, \
                             tr_z_xy_xxzz, \
                             tr_z_xy_xyyy, \
                             tr_z_xy_xyyz, \
                             tr_z_xy_xyzz, \
                             tr_z_xy_xzzz, \
                             tr_z_xy_yyyy, \
                             tr_z_xy_yyyz, \
                             tr_z_xy_yyzz, \
                             tr_z_xy_yzzz, \
                             tr_z_xy_zzzz, \
                             tr_z_y_xxxy,  \
                             tr_z_y_xxy,   \
                             tr_z_y_xxyy,  \
                             tr_z_y_xxyz,  \
                             tr_z_y_xyy,   \
                             tr_z_y_xyyy,  \
                             tr_z_y_xyyz,  \
                             tr_z_y_xyz,   \
                             tr_z_y_xyzz,  \
                             tr_z_y_yyy,   \
                             tr_z_y_yyyy,  \
                             tr_z_y_yyyz,  \
                             tr_z_y_yyz,   \
                             tr_z_y_yyzz,  \
                             tr_z_y_yzz,   \
                             tr_z_y_yzzz,  \
                             tr_z_y_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xy_xxxx[i] = tr_z_x_xxxx[i] * pa_y[i];

        tr_z_xy_xxxy[i] = 3.0 * tr_z_y_xxy[i] * fe_0 + tr_z_y_xxxy[i] * pa_x[i];

        tr_z_xy_xxxz[i] = tr_z_x_xxxz[i] * pa_y[i];

        tr_z_xy_xxyy[i] = 2.0 * tr_z_y_xyy[i] * fe_0 + tr_z_y_xxyy[i] * pa_x[i];

        tr_z_xy_xxyz[i] = 2.0 * tr_z_y_xyz[i] * fe_0 + tr_z_y_xxyz[i] * pa_x[i];

        tr_z_xy_xxzz[i] = tr_z_x_xxzz[i] * pa_y[i];

        tr_z_xy_xyyy[i] = tr_z_y_yyy[i] * fe_0 + tr_z_y_xyyy[i] * pa_x[i];

        tr_z_xy_xyyz[i] = tr_z_y_yyz[i] * fe_0 + tr_z_y_xyyz[i] * pa_x[i];

        tr_z_xy_xyzz[i] = tr_z_y_yzz[i] * fe_0 + tr_z_y_xyzz[i] * pa_x[i];

        tr_z_xy_xzzz[i] = tr_z_x_xzzz[i] * pa_y[i];

        tr_z_xy_yyyy[i] = tr_z_y_yyyy[i] * pa_x[i];

        tr_z_xy_yyyz[i] = tr_z_y_yyyz[i] * pa_x[i];

        tr_z_xy_yyzz[i] = tr_z_y_yyzz[i] * pa_x[i];

        tr_z_xy_yzzz[i] = tr_z_y_yzzz[i] * pa_x[i];

        tr_z_xy_zzzz[i] = tr_z_y_zzzz[i] * pa_x[i];
    }

    // Set up 210-225 components of targeted buffer : DG

    auto tr_z_xz_xxxx = pbuffer.data(idx_dip_dg + 210);

    auto tr_z_xz_xxxy = pbuffer.data(idx_dip_dg + 211);

    auto tr_z_xz_xxxz = pbuffer.data(idx_dip_dg + 212);

    auto tr_z_xz_xxyy = pbuffer.data(idx_dip_dg + 213);

    auto tr_z_xz_xxyz = pbuffer.data(idx_dip_dg + 214);

    auto tr_z_xz_xxzz = pbuffer.data(idx_dip_dg + 215);

    auto tr_z_xz_xyyy = pbuffer.data(idx_dip_dg + 216);

    auto tr_z_xz_xyyz = pbuffer.data(idx_dip_dg + 217);

    auto tr_z_xz_xyzz = pbuffer.data(idx_dip_dg + 218);

    auto tr_z_xz_xzzz = pbuffer.data(idx_dip_dg + 219);

    auto tr_z_xz_yyyy = pbuffer.data(idx_dip_dg + 220);

    auto tr_z_xz_yyyz = pbuffer.data(idx_dip_dg + 221);

    auto tr_z_xz_yyzz = pbuffer.data(idx_dip_dg + 222);

    auto tr_z_xz_yzzz = pbuffer.data(idx_dip_dg + 223);

    auto tr_z_xz_zzzz = pbuffer.data(idx_dip_dg + 224);

#pragma omp simd aligned(pa_x,             \
                             tr_z_xz_xxxx, \
                             tr_z_xz_xxxy, \
                             tr_z_xz_xxxz, \
                             tr_z_xz_xxyy, \
                             tr_z_xz_xxyz, \
                             tr_z_xz_xxzz, \
                             tr_z_xz_xyyy, \
                             tr_z_xz_xyyz, \
                             tr_z_xz_xyzz, \
                             tr_z_xz_xzzz, \
                             tr_z_xz_yyyy, \
                             tr_z_xz_yyyz, \
                             tr_z_xz_yyzz, \
                             tr_z_xz_yzzz, \
                             tr_z_xz_zzzz, \
                             tr_z_z_xxx,   \
                             tr_z_z_xxxx,  \
                             tr_z_z_xxxy,  \
                             tr_z_z_xxxz,  \
                             tr_z_z_xxy,   \
                             tr_z_z_xxyy,  \
                             tr_z_z_xxyz,  \
                             tr_z_z_xxz,   \
                             tr_z_z_xxzz,  \
                             tr_z_z_xyy,   \
                             tr_z_z_xyyy,  \
                             tr_z_z_xyyz,  \
                             tr_z_z_xyz,   \
                             tr_z_z_xyzz,  \
                             tr_z_z_xzz,   \
                             tr_z_z_xzzz,  \
                             tr_z_z_yyy,   \
                             tr_z_z_yyyy,  \
                             tr_z_z_yyyz,  \
                             tr_z_z_yyz,   \
                             tr_z_z_yyzz,  \
                             tr_z_z_yzz,   \
                             tr_z_z_yzzz,  \
                             tr_z_z_zzz,   \
                             tr_z_z_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xz_xxxx[i] = 4.0 * tr_z_z_xxx[i] * fe_0 + tr_z_z_xxxx[i] * pa_x[i];

        tr_z_xz_xxxy[i] = 3.0 * tr_z_z_xxy[i] * fe_0 + tr_z_z_xxxy[i] * pa_x[i];

        tr_z_xz_xxxz[i] = 3.0 * tr_z_z_xxz[i] * fe_0 + tr_z_z_xxxz[i] * pa_x[i];

        tr_z_xz_xxyy[i] = 2.0 * tr_z_z_xyy[i] * fe_0 + tr_z_z_xxyy[i] * pa_x[i];

        tr_z_xz_xxyz[i] = 2.0 * tr_z_z_xyz[i] * fe_0 + tr_z_z_xxyz[i] * pa_x[i];

        tr_z_xz_xxzz[i] = 2.0 * tr_z_z_xzz[i] * fe_0 + tr_z_z_xxzz[i] * pa_x[i];

        tr_z_xz_xyyy[i] = tr_z_z_yyy[i] * fe_0 + tr_z_z_xyyy[i] * pa_x[i];

        tr_z_xz_xyyz[i] = tr_z_z_yyz[i] * fe_0 + tr_z_z_xyyz[i] * pa_x[i];

        tr_z_xz_xyzz[i] = tr_z_z_yzz[i] * fe_0 + tr_z_z_xyzz[i] * pa_x[i];

        tr_z_xz_xzzz[i] = tr_z_z_zzz[i] * fe_0 + tr_z_z_xzzz[i] * pa_x[i];

        tr_z_xz_yyyy[i] = tr_z_z_yyyy[i] * pa_x[i];

        tr_z_xz_yyyz[i] = tr_z_z_yyyz[i] * pa_x[i];

        tr_z_xz_yyzz[i] = tr_z_z_yyzz[i] * pa_x[i];

        tr_z_xz_yzzz[i] = tr_z_z_yzzz[i] * pa_x[i];

        tr_z_xz_zzzz[i] = tr_z_z_zzzz[i] * pa_x[i];
    }

    // Set up 225-240 components of targeted buffer : DG

    auto tr_z_yy_xxxx = pbuffer.data(idx_dip_dg + 225);

    auto tr_z_yy_xxxy = pbuffer.data(idx_dip_dg + 226);

    auto tr_z_yy_xxxz = pbuffer.data(idx_dip_dg + 227);

    auto tr_z_yy_xxyy = pbuffer.data(idx_dip_dg + 228);

    auto tr_z_yy_xxyz = pbuffer.data(idx_dip_dg + 229);

    auto tr_z_yy_xxzz = pbuffer.data(idx_dip_dg + 230);

    auto tr_z_yy_xyyy = pbuffer.data(idx_dip_dg + 231);

    auto tr_z_yy_xyyz = pbuffer.data(idx_dip_dg + 232);

    auto tr_z_yy_xyzz = pbuffer.data(idx_dip_dg + 233);

    auto tr_z_yy_xzzz = pbuffer.data(idx_dip_dg + 234);

    auto tr_z_yy_yyyy = pbuffer.data(idx_dip_dg + 235);

    auto tr_z_yy_yyyz = pbuffer.data(idx_dip_dg + 236);

    auto tr_z_yy_yyzz = pbuffer.data(idx_dip_dg + 237);

    auto tr_z_yy_yzzz = pbuffer.data(idx_dip_dg + 238);

    auto tr_z_yy_zzzz = pbuffer.data(idx_dip_dg + 239);

#pragma omp simd aligned(pa_y,             \
                             tr_z_0_xxxx,  \
                             tr_z_0_xxxy,  \
                             tr_z_0_xxxz,  \
                             tr_z_0_xxyy,  \
                             tr_z_0_xxyz,  \
                             tr_z_0_xxzz,  \
                             tr_z_0_xyyy,  \
                             tr_z_0_xyyz,  \
                             tr_z_0_xyzz,  \
                             tr_z_0_xzzz,  \
                             tr_z_0_yyyy,  \
                             tr_z_0_yyyz,  \
                             tr_z_0_yyzz,  \
                             tr_z_0_yzzz,  \
                             tr_z_0_zzzz,  \
                             tr_z_y_xxx,   \
                             tr_z_y_xxxx,  \
                             tr_z_y_xxxy,  \
                             tr_z_y_xxxz,  \
                             tr_z_y_xxy,   \
                             tr_z_y_xxyy,  \
                             tr_z_y_xxyz,  \
                             tr_z_y_xxz,   \
                             tr_z_y_xxzz,  \
                             tr_z_y_xyy,   \
                             tr_z_y_xyyy,  \
                             tr_z_y_xyyz,  \
                             tr_z_y_xyz,   \
                             tr_z_y_xyzz,  \
                             tr_z_y_xzz,   \
                             tr_z_y_xzzz,  \
                             tr_z_y_yyy,   \
                             tr_z_y_yyyy,  \
                             tr_z_y_yyyz,  \
                             tr_z_y_yyz,   \
                             tr_z_y_yyzz,  \
                             tr_z_y_yzz,   \
                             tr_z_y_yzzz,  \
                             tr_z_y_zzz,   \
                             tr_z_y_zzzz,  \
                             tr_z_yy_xxxx, \
                             tr_z_yy_xxxy, \
                             tr_z_yy_xxxz, \
                             tr_z_yy_xxyy, \
                             tr_z_yy_xxyz, \
                             tr_z_yy_xxzz, \
                             tr_z_yy_xyyy, \
                             tr_z_yy_xyyz, \
                             tr_z_yy_xyzz, \
                             tr_z_yy_xzzz, \
                             tr_z_yy_yyyy, \
                             tr_z_yy_yyyz, \
                             tr_z_yy_yyzz, \
                             tr_z_yy_yzzz, \
                             tr_z_yy_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yy_xxxx[i] = tr_z_0_xxxx[i] * fe_0 + tr_z_y_xxxx[i] * pa_y[i];

        tr_z_yy_xxxy[i] = tr_z_0_xxxy[i] * fe_0 + tr_z_y_xxx[i] * fe_0 + tr_z_y_xxxy[i] * pa_y[i];

        tr_z_yy_xxxz[i] = tr_z_0_xxxz[i] * fe_0 + tr_z_y_xxxz[i] * pa_y[i];

        tr_z_yy_xxyy[i] = tr_z_0_xxyy[i] * fe_0 + 2.0 * tr_z_y_xxy[i] * fe_0 + tr_z_y_xxyy[i] * pa_y[i];

        tr_z_yy_xxyz[i] = tr_z_0_xxyz[i] * fe_0 + tr_z_y_xxz[i] * fe_0 + tr_z_y_xxyz[i] * pa_y[i];

        tr_z_yy_xxzz[i] = tr_z_0_xxzz[i] * fe_0 + tr_z_y_xxzz[i] * pa_y[i];

        tr_z_yy_xyyy[i] = tr_z_0_xyyy[i] * fe_0 + 3.0 * tr_z_y_xyy[i] * fe_0 + tr_z_y_xyyy[i] * pa_y[i];

        tr_z_yy_xyyz[i] = tr_z_0_xyyz[i] * fe_0 + 2.0 * tr_z_y_xyz[i] * fe_0 + tr_z_y_xyyz[i] * pa_y[i];

        tr_z_yy_xyzz[i] = tr_z_0_xyzz[i] * fe_0 + tr_z_y_xzz[i] * fe_0 + tr_z_y_xyzz[i] * pa_y[i];

        tr_z_yy_xzzz[i] = tr_z_0_xzzz[i] * fe_0 + tr_z_y_xzzz[i] * pa_y[i];

        tr_z_yy_yyyy[i] = tr_z_0_yyyy[i] * fe_0 + 4.0 * tr_z_y_yyy[i] * fe_0 + tr_z_y_yyyy[i] * pa_y[i];

        tr_z_yy_yyyz[i] = tr_z_0_yyyz[i] * fe_0 + 3.0 * tr_z_y_yyz[i] * fe_0 + tr_z_y_yyyz[i] * pa_y[i];

        tr_z_yy_yyzz[i] = tr_z_0_yyzz[i] * fe_0 + 2.0 * tr_z_y_yzz[i] * fe_0 + tr_z_y_yyzz[i] * pa_y[i];

        tr_z_yy_yzzz[i] = tr_z_0_yzzz[i] * fe_0 + tr_z_y_zzz[i] * fe_0 + tr_z_y_yzzz[i] * pa_y[i];

        tr_z_yy_zzzz[i] = tr_z_0_zzzz[i] * fe_0 + tr_z_y_zzzz[i] * pa_y[i];
    }

    // Set up 240-255 components of targeted buffer : DG

    auto tr_z_yz_xxxx = pbuffer.data(idx_dip_dg + 240);

    auto tr_z_yz_xxxy = pbuffer.data(idx_dip_dg + 241);

    auto tr_z_yz_xxxz = pbuffer.data(idx_dip_dg + 242);

    auto tr_z_yz_xxyy = pbuffer.data(idx_dip_dg + 243);

    auto tr_z_yz_xxyz = pbuffer.data(idx_dip_dg + 244);

    auto tr_z_yz_xxzz = pbuffer.data(idx_dip_dg + 245);

    auto tr_z_yz_xyyy = pbuffer.data(idx_dip_dg + 246);

    auto tr_z_yz_xyyz = pbuffer.data(idx_dip_dg + 247);

    auto tr_z_yz_xyzz = pbuffer.data(idx_dip_dg + 248);

    auto tr_z_yz_xzzz = pbuffer.data(idx_dip_dg + 249);

    auto tr_z_yz_yyyy = pbuffer.data(idx_dip_dg + 250);

    auto tr_z_yz_yyyz = pbuffer.data(idx_dip_dg + 251);

    auto tr_z_yz_yyzz = pbuffer.data(idx_dip_dg + 252);

    auto tr_z_yz_yzzz = pbuffer.data(idx_dip_dg + 253);

    auto tr_z_yz_zzzz = pbuffer.data(idx_dip_dg + 254);

#pragma omp simd aligned(pa_y,             \
                             tr_z_yz_xxxx, \
                             tr_z_yz_xxxy, \
                             tr_z_yz_xxxz, \
                             tr_z_yz_xxyy, \
                             tr_z_yz_xxyz, \
                             tr_z_yz_xxzz, \
                             tr_z_yz_xyyy, \
                             tr_z_yz_xyyz, \
                             tr_z_yz_xyzz, \
                             tr_z_yz_xzzz, \
                             tr_z_yz_yyyy, \
                             tr_z_yz_yyyz, \
                             tr_z_yz_yyzz, \
                             tr_z_yz_yzzz, \
                             tr_z_yz_zzzz, \
                             tr_z_z_xxx,   \
                             tr_z_z_xxxx,  \
                             tr_z_z_xxxy,  \
                             tr_z_z_xxxz,  \
                             tr_z_z_xxy,   \
                             tr_z_z_xxyy,  \
                             tr_z_z_xxyz,  \
                             tr_z_z_xxz,   \
                             tr_z_z_xxzz,  \
                             tr_z_z_xyy,   \
                             tr_z_z_xyyy,  \
                             tr_z_z_xyyz,  \
                             tr_z_z_xyz,   \
                             tr_z_z_xyzz,  \
                             tr_z_z_xzz,   \
                             tr_z_z_xzzz,  \
                             tr_z_z_yyy,   \
                             tr_z_z_yyyy,  \
                             tr_z_z_yyyz,  \
                             tr_z_z_yyz,   \
                             tr_z_z_yyzz,  \
                             tr_z_z_yzz,   \
                             tr_z_z_yzzz,  \
                             tr_z_z_zzz,   \
                             tr_z_z_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yz_xxxx[i] = tr_z_z_xxxx[i] * pa_y[i];

        tr_z_yz_xxxy[i] = tr_z_z_xxx[i] * fe_0 + tr_z_z_xxxy[i] * pa_y[i];

        tr_z_yz_xxxz[i] = tr_z_z_xxxz[i] * pa_y[i];

        tr_z_yz_xxyy[i] = 2.0 * tr_z_z_xxy[i] * fe_0 + tr_z_z_xxyy[i] * pa_y[i];

        tr_z_yz_xxyz[i] = tr_z_z_xxz[i] * fe_0 + tr_z_z_xxyz[i] * pa_y[i];

        tr_z_yz_xxzz[i] = tr_z_z_xxzz[i] * pa_y[i];

        tr_z_yz_xyyy[i] = 3.0 * tr_z_z_xyy[i] * fe_0 + tr_z_z_xyyy[i] * pa_y[i];

        tr_z_yz_xyyz[i] = 2.0 * tr_z_z_xyz[i] * fe_0 + tr_z_z_xyyz[i] * pa_y[i];

        tr_z_yz_xyzz[i] = tr_z_z_xzz[i] * fe_0 + tr_z_z_xyzz[i] * pa_y[i];

        tr_z_yz_xzzz[i] = tr_z_z_xzzz[i] * pa_y[i];

        tr_z_yz_yyyy[i] = 4.0 * tr_z_z_yyy[i] * fe_0 + tr_z_z_yyyy[i] * pa_y[i];

        tr_z_yz_yyyz[i] = 3.0 * tr_z_z_yyz[i] * fe_0 + tr_z_z_yyyz[i] * pa_y[i];

        tr_z_yz_yyzz[i] = 2.0 * tr_z_z_yzz[i] * fe_0 + tr_z_z_yyzz[i] * pa_y[i];

        tr_z_yz_yzzz[i] = tr_z_z_zzz[i] * fe_0 + tr_z_z_yzzz[i] * pa_y[i];

        tr_z_yz_zzzz[i] = tr_z_z_zzzz[i] * pa_y[i];
    }

    // Set up 255-270 components of targeted buffer : DG

    auto tr_z_zz_xxxx = pbuffer.data(idx_dip_dg + 255);

    auto tr_z_zz_xxxy = pbuffer.data(idx_dip_dg + 256);

    auto tr_z_zz_xxxz = pbuffer.data(idx_dip_dg + 257);

    auto tr_z_zz_xxyy = pbuffer.data(idx_dip_dg + 258);

    auto tr_z_zz_xxyz = pbuffer.data(idx_dip_dg + 259);

    auto tr_z_zz_xxzz = pbuffer.data(idx_dip_dg + 260);

    auto tr_z_zz_xyyy = pbuffer.data(idx_dip_dg + 261);

    auto tr_z_zz_xyyz = pbuffer.data(idx_dip_dg + 262);

    auto tr_z_zz_xyzz = pbuffer.data(idx_dip_dg + 263);

    auto tr_z_zz_xzzz = pbuffer.data(idx_dip_dg + 264);

    auto tr_z_zz_yyyy = pbuffer.data(idx_dip_dg + 265);

    auto tr_z_zz_yyyz = pbuffer.data(idx_dip_dg + 266);

    auto tr_z_zz_yyzz = pbuffer.data(idx_dip_dg + 267);

    auto tr_z_zz_yzzz = pbuffer.data(idx_dip_dg + 268);

    auto tr_z_zz_zzzz = pbuffer.data(idx_dip_dg + 269);

#pragma omp simd aligned(pa_z,             \
                             tr_z_0_xxxx,  \
                             tr_z_0_xxxy,  \
                             tr_z_0_xxxz,  \
                             tr_z_0_xxyy,  \
                             tr_z_0_xxyz,  \
                             tr_z_0_xxzz,  \
                             tr_z_0_xyyy,  \
                             tr_z_0_xyyz,  \
                             tr_z_0_xyzz,  \
                             tr_z_0_xzzz,  \
                             tr_z_0_yyyy,  \
                             tr_z_0_yyyz,  \
                             tr_z_0_yyzz,  \
                             tr_z_0_yzzz,  \
                             tr_z_0_zzzz,  \
                             tr_z_z_xxx,   \
                             tr_z_z_xxxx,  \
                             tr_z_z_xxxy,  \
                             tr_z_z_xxxz,  \
                             tr_z_z_xxy,   \
                             tr_z_z_xxyy,  \
                             tr_z_z_xxyz,  \
                             tr_z_z_xxz,   \
                             tr_z_z_xxzz,  \
                             tr_z_z_xyy,   \
                             tr_z_z_xyyy,  \
                             tr_z_z_xyyz,  \
                             tr_z_z_xyz,   \
                             tr_z_z_xyzz,  \
                             tr_z_z_xzz,   \
                             tr_z_z_xzzz,  \
                             tr_z_z_yyy,   \
                             tr_z_z_yyyy,  \
                             tr_z_z_yyyz,  \
                             tr_z_z_yyz,   \
                             tr_z_z_yyzz,  \
                             tr_z_z_yzz,   \
                             tr_z_z_yzzz,  \
                             tr_z_z_zzz,   \
                             tr_z_z_zzzz,  \
                             tr_z_zz_xxxx, \
                             tr_z_zz_xxxy, \
                             tr_z_zz_xxxz, \
                             tr_z_zz_xxyy, \
                             tr_z_zz_xxyz, \
                             tr_z_zz_xxzz, \
                             tr_z_zz_xyyy, \
                             tr_z_zz_xyyz, \
                             tr_z_zz_xyzz, \
                             tr_z_zz_xzzz, \
                             tr_z_zz_yyyy, \
                             tr_z_zz_yyyz, \
                             tr_z_zz_yyzz, \
                             tr_z_zz_yzzz, \
                             tr_z_zz_zzzz, \
                             ts_z_xxxx,    \
                             ts_z_xxxy,    \
                             ts_z_xxxz,    \
                             ts_z_xxyy,    \
                             ts_z_xxyz,    \
                             ts_z_xxzz,    \
                             ts_z_xyyy,    \
                             ts_z_xyyz,    \
                             ts_z_xyzz,    \
                             ts_z_xzzz,    \
                             ts_z_yyyy,    \
                             ts_z_yyyz,    \
                             ts_z_yyzz,    \
                             ts_z_yzzz,    \
                             ts_z_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_zz_xxxx[i] = tr_z_0_xxxx[i] * fe_0 + ts_z_xxxx[i] * fe_0 + tr_z_z_xxxx[i] * pa_z[i];

        tr_z_zz_xxxy[i] = tr_z_0_xxxy[i] * fe_0 + ts_z_xxxy[i] * fe_0 + tr_z_z_xxxy[i] * pa_z[i];

        tr_z_zz_xxxz[i] = tr_z_0_xxxz[i] * fe_0 + tr_z_z_xxx[i] * fe_0 + ts_z_xxxz[i] * fe_0 + tr_z_z_xxxz[i] * pa_z[i];

        tr_z_zz_xxyy[i] = tr_z_0_xxyy[i] * fe_0 + ts_z_xxyy[i] * fe_0 + tr_z_z_xxyy[i] * pa_z[i];

        tr_z_zz_xxyz[i] = tr_z_0_xxyz[i] * fe_0 + tr_z_z_xxy[i] * fe_0 + ts_z_xxyz[i] * fe_0 + tr_z_z_xxyz[i] * pa_z[i];

        tr_z_zz_xxzz[i] = tr_z_0_xxzz[i] * fe_0 + 2.0 * tr_z_z_xxz[i] * fe_0 + ts_z_xxzz[i] * fe_0 + tr_z_z_xxzz[i] * pa_z[i];

        tr_z_zz_xyyy[i] = tr_z_0_xyyy[i] * fe_0 + ts_z_xyyy[i] * fe_0 + tr_z_z_xyyy[i] * pa_z[i];

        tr_z_zz_xyyz[i] = tr_z_0_xyyz[i] * fe_0 + tr_z_z_xyy[i] * fe_0 + ts_z_xyyz[i] * fe_0 + tr_z_z_xyyz[i] * pa_z[i];

        tr_z_zz_xyzz[i] = tr_z_0_xyzz[i] * fe_0 + 2.0 * tr_z_z_xyz[i] * fe_0 + ts_z_xyzz[i] * fe_0 + tr_z_z_xyzz[i] * pa_z[i];

        tr_z_zz_xzzz[i] = tr_z_0_xzzz[i] * fe_0 + 3.0 * tr_z_z_xzz[i] * fe_0 + ts_z_xzzz[i] * fe_0 + tr_z_z_xzzz[i] * pa_z[i];

        tr_z_zz_yyyy[i] = tr_z_0_yyyy[i] * fe_0 + ts_z_yyyy[i] * fe_0 + tr_z_z_yyyy[i] * pa_z[i];

        tr_z_zz_yyyz[i] = tr_z_0_yyyz[i] * fe_0 + tr_z_z_yyy[i] * fe_0 + ts_z_yyyz[i] * fe_0 + tr_z_z_yyyz[i] * pa_z[i];

        tr_z_zz_yyzz[i] = tr_z_0_yyzz[i] * fe_0 + 2.0 * tr_z_z_yyz[i] * fe_0 + ts_z_yyzz[i] * fe_0 + tr_z_z_yyzz[i] * pa_z[i];

        tr_z_zz_yzzz[i] = tr_z_0_yzzz[i] * fe_0 + 3.0 * tr_z_z_yzz[i] * fe_0 + ts_z_yzzz[i] * fe_0 + tr_z_z_yzzz[i] * pa_z[i];

        tr_z_zz_zzzz[i] = tr_z_0_zzzz[i] * fe_0 + 4.0 * tr_z_z_zzz[i] * fe_0 + ts_z_zzzz[i] * fe_0 + tr_z_z_zzzz[i] * pa_z[i];
    }
}

}  // namespace diprec
