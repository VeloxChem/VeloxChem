#include "NuclearPotentialGeom020PrimRecDG.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_020_dg(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_020_0_dg,
                                        const size_t              idx_npot_geom_020_0_sg,
                                        const size_t              idx_npot_geom_020_1_sg,
                                        const size_t              idx_npot_geom_020_0_pf,
                                        const size_t              idx_npot_geom_020_1_pf,
                                        const size_t              idx_npot_geom_010_1_pg,
                                        const size_t              idx_npot_geom_020_0_pg,
                                        const size_t              idx_npot_geom_020_1_pg,
                                        const CSimdArray<double>& factors,
                                        const size_t              idx_rpa,
                                        const size_t              idx_rpc,
                                        const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up R(PC) distances

    auto pc_x = factors.data(idx_rpc);

    auto pc_y = factors.data(idx_rpc + 1);

    auto pc_z = factors.data(idx_rpc + 2);

    // Set up components of auxiliary buffer : SG

    auto ta2_xx_0_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_sg);

    auto ta2_xx_0_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 1);

    auto ta2_xx_0_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 2);

    auto ta2_xx_0_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 3);

    auto ta2_xx_0_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 4);

    auto ta2_xx_0_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 5);

    auto ta2_xx_0_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 6);

    auto ta2_xx_0_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 7);

    auto ta2_xx_0_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 8);

    auto ta2_xx_0_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 9);

    auto ta2_xx_0_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 10);

    auto ta2_xx_0_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 11);

    auto ta2_xx_0_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 12);

    auto ta2_xx_0_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 13);

    auto ta2_xx_0_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 14);

    auto ta2_xy_0_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_sg + 15);

    auto ta2_xy_0_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 16);

    auto ta2_xy_0_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 17);

    auto ta2_xy_0_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 18);

    auto ta2_xy_0_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 19);

    auto ta2_xy_0_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 20);

    auto ta2_xy_0_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 21);

    auto ta2_xy_0_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 22);

    auto ta2_xy_0_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 23);

    auto ta2_xy_0_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 24);

    auto ta2_xy_0_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 25);

    auto ta2_xy_0_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 26);

    auto ta2_xy_0_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 27);

    auto ta2_xy_0_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 28);

    auto ta2_xy_0_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 29);

    auto ta2_xz_0_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_sg + 30);

    auto ta2_xz_0_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 31);

    auto ta2_xz_0_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 32);

    auto ta2_xz_0_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 33);

    auto ta2_xz_0_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 34);

    auto ta2_xz_0_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 35);

    auto ta2_xz_0_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 36);

    auto ta2_xz_0_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 37);

    auto ta2_xz_0_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 38);

    auto ta2_xz_0_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 39);

    auto ta2_xz_0_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 40);

    auto ta2_xz_0_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 41);

    auto ta2_xz_0_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 42);

    auto ta2_xz_0_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 43);

    auto ta2_xz_0_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 44);

    auto ta2_yy_0_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_sg + 45);

    auto ta2_yy_0_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 46);

    auto ta2_yy_0_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 47);

    auto ta2_yy_0_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 48);

    auto ta2_yy_0_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 49);

    auto ta2_yy_0_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 50);

    auto ta2_yy_0_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 51);

    auto ta2_yy_0_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 52);

    auto ta2_yy_0_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 53);

    auto ta2_yy_0_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 54);

    auto ta2_yy_0_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 55);

    auto ta2_yy_0_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 56);

    auto ta2_yy_0_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 57);

    auto ta2_yy_0_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 58);

    auto ta2_yy_0_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 59);

    auto ta2_yz_0_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_sg + 60);

    auto ta2_yz_0_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 61);

    auto ta2_yz_0_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 62);

    auto ta2_yz_0_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 63);

    auto ta2_yz_0_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 64);

    auto ta2_yz_0_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 65);

    auto ta2_yz_0_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 66);

    auto ta2_yz_0_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 67);

    auto ta2_yz_0_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 68);

    auto ta2_yz_0_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 69);

    auto ta2_yz_0_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 70);

    auto ta2_yz_0_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 71);

    auto ta2_yz_0_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 72);

    auto ta2_yz_0_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 73);

    auto ta2_yz_0_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 74);

    auto ta2_zz_0_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_sg + 75);

    auto ta2_zz_0_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 76);

    auto ta2_zz_0_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 77);

    auto ta2_zz_0_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 78);

    auto ta2_zz_0_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 79);

    auto ta2_zz_0_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 80);

    auto ta2_zz_0_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 81);

    auto ta2_zz_0_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 82);

    auto ta2_zz_0_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 83);

    auto ta2_zz_0_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 84);

    auto ta2_zz_0_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_sg + 85);

    auto ta2_zz_0_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 86);

    auto ta2_zz_0_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 87);

    auto ta2_zz_0_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 88);

    auto ta2_zz_0_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_sg + 89);

    // Set up components of auxiliary buffer : SG

    auto ta2_xx_0_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_sg);

    auto ta2_xx_0_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 1);

    auto ta2_xx_0_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 2);

    auto ta2_xx_0_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 3);

    auto ta2_xx_0_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 4);

    auto ta2_xx_0_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 5);

    auto ta2_xx_0_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 6);

    auto ta2_xx_0_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 7);

    auto ta2_xx_0_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 8);

    auto ta2_xx_0_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 9);

    auto ta2_xx_0_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 10);

    auto ta2_xx_0_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 11);

    auto ta2_xx_0_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 12);

    auto ta2_xx_0_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 13);

    auto ta2_xx_0_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 14);

    auto ta2_xy_0_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_sg + 15);

    auto ta2_xy_0_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 16);

    auto ta2_xy_0_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 17);

    auto ta2_xy_0_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 18);

    auto ta2_xy_0_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 19);

    auto ta2_xy_0_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 20);

    auto ta2_xy_0_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 21);

    auto ta2_xy_0_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 22);

    auto ta2_xy_0_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 23);

    auto ta2_xy_0_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 24);

    auto ta2_xy_0_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 25);

    auto ta2_xy_0_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 26);

    auto ta2_xy_0_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 27);

    auto ta2_xy_0_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 28);

    auto ta2_xy_0_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 29);

    auto ta2_xz_0_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_sg + 30);

    auto ta2_xz_0_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 31);

    auto ta2_xz_0_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 32);

    auto ta2_xz_0_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 33);

    auto ta2_xz_0_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 34);

    auto ta2_xz_0_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 35);

    auto ta2_xz_0_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 36);

    auto ta2_xz_0_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 37);

    auto ta2_xz_0_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 38);

    auto ta2_xz_0_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 39);

    auto ta2_xz_0_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 40);

    auto ta2_xz_0_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 41);

    auto ta2_xz_0_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 42);

    auto ta2_xz_0_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 43);

    auto ta2_xz_0_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 44);

    auto ta2_yy_0_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_sg + 45);

    auto ta2_yy_0_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 46);

    auto ta2_yy_0_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 47);

    auto ta2_yy_0_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 48);

    auto ta2_yy_0_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 49);

    auto ta2_yy_0_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 50);

    auto ta2_yy_0_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 51);

    auto ta2_yy_0_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 52);

    auto ta2_yy_0_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 53);

    auto ta2_yy_0_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 54);

    auto ta2_yy_0_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 55);

    auto ta2_yy_0_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 56);

    auto ta2_yy_0_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 57);

    auto ta2_yy_0_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 58);

    auto ta2_yy_0_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 59);

    auto ta2_yz_0_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_sg + 60);

    auto ta2_yz_0_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 61);

    auto ta2_yz_0_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 62);

    auto ta2_yz_0_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 63);

    auto ta2_yz_0_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 64);

    auto ta2_yz_0_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 65);

    auto ta2_yz_0_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 66);

    auto ta2_yz_0_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 67);

    auto ta2_yz_0_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 68);

    auto ta2_yz_0_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 69);

    auto ta2_yz_0_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 70);

    auto ta2_yz_0_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 71);

    auto ta2_yz_0_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 72);

    auto ta2_yz_0_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 73);

    auto ta2_yz_0_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 74);

    auto ta2_zz_0_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_sg + 75);

    auto ta2_zz_0_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 76);

    auto ta2_zz_0_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 77);

    auto ta2_zz_0_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 78);

    auto ta2_zz_0_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 79);

    auto ta2_zz_0_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 80);

    auto ta2_zz_0_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 81);

    auto ta2_zz_0_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 82);

    auto ta2_zz_0_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 83);

    auto ta2_zz_0_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 84);

    auto ta2_zz_0_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_sg + 85);

    auto ta2_zz_0_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 86);

    auto ta2_zz_0_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 87);

    auto ta2_zz_0_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 88);

    auto ta2_zz_0_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_sg + 89);

    // Set up components of auxiliary buffer : PF

    auto ta2_xx_x_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf);

    auto ta2_xx_x_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 1);

    auto ta2_xx_x_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 2);

    auto ta2_xx_x_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 3);

    auto ta2_xx_x_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 4);

    auto ta2_xx_x_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 5);

    auto ta2_xx_x_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 6);

    auto ta2_xx_x_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 7);

    auto ta2_xx_x_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 8);

    auto ta2_xx_x_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 9);

    auto ta2_xx_y_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 10);

    auto ta2_xx_y_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 11);

    auto ta2_xx_y_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 12);

    auto ta2_xx_y_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 13);

    auto ta2_xx_y_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 14);

    auto ta2_xx_y_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 15);

    auto ta2_xx_y_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 16);

    auto ta2_xx_y_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 17);

    auto ta2_xx_y_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 18);

    auto ta2_xx_y_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 19);

    auto ta2_xx_z_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 20);

    auto ta2_xx_z_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 21);

    auto ta2_xx_z_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 22);

    auto ta2_xx_z_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 23);

    auto ta2_xx_z_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 24);

    auto ta2_xx_z_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 25);

    auto ta2_xx_z_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 26);

    auto ta2_xx_z_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 27);

    auto ta2_xx_z_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 28);

    auto ta2_xx_z_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 29);

    auto ta2_xy_x_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 30);

    auto ta2_xy_x_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 31);

    auto ta2_xy_x_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 32);

    auto ta2_xy_x_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 33);

    auto ta2_xy_x_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 34);

    auto ta2_xy_x_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 35);

    auto ta2_xy_x_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 36);

    auto ta2_xy_x_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 37);

    auto ta2_xy_x_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 38);

    auto ta2_xy_x_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 39);

    auto ta2_xy_y_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 40);

    auto ta2_xy_y_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 41);

    auto ta2_xy_y_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 42);

    auto ta2_xy_y_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 43);

    auto ta2_xy_y_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 44);

    auto ta2_xy_y_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 45);

    auto ta2_xy_y_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 46);

    auto ta2_xy_y_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 47);

    auto ta2_xy_y_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 48);

    auto ta2_xy_y_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 49);

    auto ta2_xy_z_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 50);

    auto ta2_xy_z_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 51);

    auto ta2_xy_z_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 52);

    auto ta2_xy_z_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 53);

    auto ta2_xy_z_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 54);

    auto ta2_xy_z_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 55);

    auto ta2_xy_z_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 56);

    auto ta2_xy_z_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 57);

    auto ta2_xy_z_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 58);

    auto ta2_xy_z_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 59);

    auto ta2_xz_x_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 60);

    auto ta2_xz_x_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 61);

    auto ta2_xz_x_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 62);

    auto ta2_xz_x_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 63);

    auto ta2_xz_x_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 64);

    auto ta2_xz_x_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 65);

    auto ta2_xz_x_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 66);

    auto ta2_xz_x_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 67);

    auto ta2_xz_x_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 68);

    auto ta2_xz_x_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 69);

    auto ta2_xz_y_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 70);

    auto ta2_xz_y_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 71);

    auto ta2_xz_y_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 72);

    auto ta2_xz_y_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 73);

    auto ta2_xz_y_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 74);

    auto ta2_xz_y_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 75);

    auto ta2_xz_y_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 76);

    auto ta2_xz_y_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 77);

    auto ta2_xz_y_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 78);

    auto ta2_xz_y_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 79);

    auto ta2_xz_z_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 80);

    auto ta2_xz_z_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 81);

    auto ta2_xz_z_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 82);

    auto ta2_xz_z_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 83);

    auto ta2_xz_z_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 84);

    auto ta2_xz_z_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 85);

    auto ta2_xz_z_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 86);

    auto ta2_xz_z_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 87);

    auto ta2_xz_z_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 88);

    auto ta2_xz_z_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 89);

    auto ta2_yy_x_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 90);

    auto ta2_yy_x_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 91);

    auto ta2_yy_x_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 92);

    auto ta2_yy_x_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 93);

    auto ta2_yy_x_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 94);

    auto ta2_yy_x_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 95);

    auto ta2_yy_x_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 96);

    auto ta2_yy_x_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 97);

    auto ta2_yy_x_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 98);

    auto ta2_yy_x_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 99);

    auto ta2_yy_y_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 100);

    auto ta2_yy_y_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 101);

    auto ta2_yy_y_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 102);

    auto ta2_yy_y_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 103);

    auto ta2_yy_y_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 104);

    auto ta2_yy_y_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 105);

    auto ta2_yy_y_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 106);

    auto ta2_yy_y_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 107);

    auto ta2_yy_y_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 108);

    auto ta2_yy_y_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 109);

    auto ta2_yy_z_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 110);

    auto ta2_yy_z_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 111);

    auto ta2_yy_z_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 112);

    auto ta2_yy_z_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 113);

    auto ta2_yy_z_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 114);

    auto ta2_yy_z_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 115);

    auto ta2_yy_z_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 116);

    auto ta2_yy_z_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 117);

    auto ta2_yy_z_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 118);

    auto ta2_yy_z_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 119);

    auto ta2_yz_x_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 120);

    auto ta2_yz_x_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 121);

    auto ta2_yz_x_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 122);

    auto ta2_yz_x_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 123);

    auto ta2_yz_x_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 124);

    auto ta2_yz_x_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 125);

    auto ta2_yz_x_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 126);

    auto ta2_yz_x_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 127);

    auto ta2_yz_x_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 128);

    auto ta2_yz_x_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 129);

    auto ta2_yz_y_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 130);

    auto ta2_yz_y_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 131);

    auto ta2_yz_y_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 132);

    auto ta2_yz_y_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 133);

    auto ta2_yz_y_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 134);

    auto ta2_yz_y_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 135);

    auto ta2_yz_y_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 136);

    auto ta2_yz_y_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 137);

    auto ta2_yz_y_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 138);

    auto ta2_yz_y_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 139);

    auto ta2_yz_z_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 140);

    auto ta2_yz_z_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 141);

    auto ta2_yz_z_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 142);

    auto ta2_yz_z_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 143);

    auto ta2_yz_z_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 144);

    auto ta2_yz_z_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 145);

    auto ta2_yz_z_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 146);

    auto ta2_yz_z_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 147);

    auto ta2_yz_z_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 148);

    auto ta2_yz_z_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 149);

    auto ta2_zz_x_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 150);

    auto ta2_zz_x_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 151);

    auto ta2_zz_x_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 152);

    auto ta2_zz_x_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 153);

    auto ta2_zz_x_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 154);

    auto ta2_zz_x_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 155);

    auto ta2_zz_x_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 156);

    auto ta2_zz_x_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 157);

    auto ta2_zz_x_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 158);

    auto ta2_zz_x_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 159);

    auto ta2_zz_y_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 160);

    auto ta2_zz_y_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 161);

    auto ta2_zz_y_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 162);

    auto ta2_zz_y_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 163);

    auto ta2_zz_y_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 164);

    auto ta2_zz_y_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 165);

    auto ta2_zz_y_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 166);

    auto ta2_zz_y_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 167);

    auto ta2_zz_y_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 168);

    auto ta2_zz_y_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 169);

    auto ta2_zz_z_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 170);

    auto ta2_zz_z_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 171);

    auto ta2_zz_z_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 172);

    auto ta2_zz_z_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 173);

    auto ta2_zz_z_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 174);

    auto ta2_zz_z_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 175);

    auto ta2_zz_z_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 176);

    auto ta2_zz_z_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 177);

    auto ta2_zz_z_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 178);

    auto ta2_zz_z_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 179);

    // Set up components of auxiliary buffer : PF

    auto ta2_xx_x_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf);

    auto ta2_xx_x_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 1);

    auto ta2_xx_x_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 2);

    auto ta2_xx_x_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 3);

    auto ta2_xx_x_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 4);

    auto ta2_xx_x_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 5);

    auto ta2_xx_x_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 6);

    auto ta2_xx_x_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 7);

    auto ta2_xx_x_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 8);

    auto ta2_xx_x_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 9);

    auto ta2_xx_y_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 10);

    auto ta2_xx_y_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 11);

    auto ta2_xx_y_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 12);

    auto ta2_xx_y_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 13);

    auto ta2_xx_y_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 14);

    auto ta2_xx_y_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 15);

    auto ta2_xx_y_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 16);

    auto ta2_xx_y_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 17);

    auto ta2_xx_y_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 18);

    auto ta2_xx_y_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 19);

    auto ta2_xx_z_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 20);

    auto ta2_xx_z_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 21);

    auto ta2_xx_z_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 22);

    auto ta2_xx_z_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 23);

    auto ta2_xx_z_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 24);

    auto ta2_xx_z_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 25);

    auto ta2_xx_z_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 26);

    auto ta2_xx_z_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 27);

    auto ta2_xx_z_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 28);

    auto ta2_xx_z_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 29);

    auto ta2_xy_x_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 30);

    auto ta2_xy_x_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 31);

    auto ta2_xy_x_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 32);

    auto ta2_xy_x_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 33);

    auto ta2_xy_x_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 34);

    auto ta2_xy_x_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 35);

    auto ta2_xy_x_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 36);

    auto ta2_xy_x_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 37);

    auto ta2_xy_x_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 38);

    auto ta2_xy_x_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 39);

    auto ta2_xy_y_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 40);

    auto ta2_xy_y_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 41);

    auto ta2_xy_y_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 42);

    auto ta2_xy_y_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 43);

    auto ta2_xy_y_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 44);

    auto ta2_xy_y_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 45);

    auto ta2_xy_y_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 46);

    auto ta2_xy_y_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 47);

    auto ta2_xy_y_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 48);

    auto ta2_xy_y_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 49);

    auto ta2_xy_z_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 50);

    auto ta2_xy_z_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 51);

    auto ta2_xy_z_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 52);

    auto ta2_xy_z_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 53);

    auto ta2_xy_z_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 54);

    auto ta2_xy_z_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 55);

    auto ta2_xy_z_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 56);

    auto ta2_xy_z_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 57);

    auto ta2_xy_z_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 58);

    auto ta2_xy_z_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 59);

    auto ta2_xz_x_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 60);

    auto ta2_xz_x_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 61);

    auto ta2_xz_x_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 62);

    auto ta2_xz_x_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 63);

    auto ta2_xz_x_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 64);

    auto ta2_xz_x_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 65);

    auto ta2_xz_x_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 66);

    auto ta2_xz_x_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 67);

    auto ta2_xz_x_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 68);

    auto ta2_xz_x_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 69);

    auto ta2_xz_y_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 70);

    auto ta2_xz_y_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 71);

    auto ta2_xz_y_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 72);

    auto ta2_xz_y_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 73);

    auto ta2_xz_y_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 74);

    auto ta2_xz_y_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 75);

    auto ta2_xz_y_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 76);

    auto ta2_xz_y_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 77);

    auto ta2_xz_y_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 78);

    auto ta2_xz_y_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 79);

    auto ta2_xz_z_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 80);

    auto ta2_xz_z_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 81);

    auto ta2_xz_z_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 82);

    auto ta2_xz_z_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 83);

    auto ta2_xz_z_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 84);

    auto ta2_xz_z_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 85);

    auto ta2_xz_z_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 86);

    auto ta2_xz_z_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 87);

    auto ta2_xz_z_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 88);

    auto ta2_xz_z_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 89);

    auto ta2_yy_x_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 90);

    auto ta2_yy_x_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 91);

    auto ta2_yy_x_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 92);

    auto ta2_yy_x_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 93);

    auto ta2_yy_x_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 94);

    auto ta2_yy_x_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 95);

    auto ta2_yy_x_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 96);

    auto ta2_yy_x_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 97);

    auto ta2_yy_x_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 98);

    auto ta2_yy_x_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 99);

    auto ta2_yy_y_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 100);

    auto ta2_yy_y_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 101);

    auto ta2_yy_y_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 102);

    auto ta2_yy_y_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 103);

    auto ta2_yy_y_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 104);

    auto ta2_yy_y_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 105);

    auto ta2_yy_y_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 106);

    auto ta2_yy_y_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 107);

    auto ta2_yy_y_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 108);

    auto ta2_yy_y_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 109);

    auto ta2_yy_z_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 110);

    auto ta2_yy_z_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 111);

    auto ta2_yy_z_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 112);

    auto ta2_yy_z_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 113);

    auto ta2_yy_z_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 114);

    auto ta2_yy_z_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 115);

    auto ta2_yy_z_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 116);

    auto ta2_yy_z_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 117);

    auto ta2_yy_z_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 118);

    auto ta2_yy_z_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 119);

    auto ta2_yz_x_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 120);

    auto ta2_yz_x_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 121);

    auto ta2_yz_x_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 122);

    auto ta2_yz_x_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 123);

    auto ta2_yz_x_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 124);

    auto ta2_yz_x_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 125);

    auto ta2_yz_x_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 126);

    auto ta2_yz_x_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 127);

    auto ta2_yz_x_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 128);

    auto ta2_yz_x_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 129);

    auto ta2_yz_y_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 130);

    auto ta2_yz_y_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 131);

    auto ta2_yz_y_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 132);

    auto ta2_yz_y_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 133);

    auto ta2_yz_y_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 134);

    auto ta2_yz_y_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 135);

    auto ta2_yz_y_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 136);

    auto ta2_yz_y_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 137);

    auto ta2_yz_y_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 138);

    auto ta2_yz_y_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 139);

    auto ta2_yz_z_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 140);

    auto ta2_yz_z_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 141);

    auto ta2_yz_z_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 142);

    auto ta2_yz_z_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 143);

    auto ta2_yz_z_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 144);

    auto ta2_yz_z_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 145);

    auto ta2_yz_z_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 146);

    auto ta2_yz_z_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 147);

    auto ta2_yz_z_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 148);

    auto ta2_yz_z_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 149);

    auto ta2_zz_x_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 150);

    auto ta2_zz_x_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 151);

    auto ta2_zz_x_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 152);

    auto ta2_zz_x_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 153);

    auto ta2_zz_x_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 154);

    auto ta2_zz_x_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 155);

    auto ta2_zz_x_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 156);

    auto ta2_zz_x_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 157);

    auto ta2_zz_x_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 158);

    auto ta2_zz_x_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 159);

    auto ta2_zz_y_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 160);

    auto ta2_zz_y_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 161);

    auto ta2_zz_y_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 162);

    auto ta2_zz_y_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 163);

    auto ta2_zz_y_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 164);

    auto ta2_zz_y_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 165);

    auto ta2_zz_y_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 166);

    auto ta2_zz_y_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 167);

    auto ta2_zz_y_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 168);

    auto ta2_zz_y_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 169);

    auto ta2_zz_z_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 170);

    auto ta2_zz_z_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 171);

    auto ta2_zz_z_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 172);

    auto ta2_zz_z_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 173);

    auto ta2_zz_z_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 174);

    auto ta2_zz_z_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 175);

    auto ta2_zz_z_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 176);

    auto ta2_zz_z_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 177);

    auto ta2_zz_z_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 178);

    auto ta2_zz_z_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 179);

    // Set up components of auxiliary buffer : PG

    auto ta1_x_x_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_pg);

    auto ta1_x_x_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 1);

    auto ta1_x_x_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 2);

    auto ta1_x_x_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 3);

    auto ta1_x_x_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 4);

    auto ta1_x_x_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 5);

    auto ta1_x_x_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 6);

    auto ta1_x_x_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 7);

    auto ta1_x_x_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 8);

    auto ta1_x_x_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 9);

    auto ta1_x_x_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 10);

    auto ta1_x_x_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 11);

    auto ta1_x_x_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 12);

    auto ta1_x_x_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 13);

    auto ta1_x_x_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 14);

    auto ta1_x_y_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_pg + 15);

    auto ta1_x_y_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 16);

    auto ta1_x_y_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 17);

    auto ta1_x_y_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 18);

    auto ta1_x_y_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 19);

    auto ta1_x_y_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 20);

    auto ta1_x_y_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 21);

    auto ta1_x_y_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 22);

    auto ta1_x_y_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 23);

    auto ta1_x_y_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 24);

    auto ta1_x_y_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 25);

    auto ta1_x_y_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 26);

    auto ta1_x_y_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 27);

    auto ta1_x_y_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 28);

    auto ta1_x_y_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 29);

    auto ta1_x_z_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_pg + 30);

    auto ta1_x_z_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 31);

    auto ta1_x_z_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 32);

    auto ta1_x_z_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 33);

    auto ta1_x_z_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 34);

    auto ta1_x_z_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 35);

    auto ta1_x_z_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 36);

    auto ta1_x_z_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 37);

    auto ta1_x_z_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 38);

    auto ta1_x_z_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 39);

    auto ta1_x_z_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 40);

    auto ta1_x_z_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 41);

    auto ta1_x_z_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 42);

    auto ta1_x_z_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 43);

    auto ta1_x_z_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 44);

    auto ta1_y_x_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_pg + 45);

    auto ta1_y_x_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 46);

    auto ta1_y_x_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 47);

    auto ta1_y_x_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 48);

    auto ta1_y_x_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 49);

    auto ta1_y_x_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 50);

    auto ta1_y_x_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 51);

    auto ta1_y_x_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 52);

    auto ta1_y_x_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 53);

    auto ta1_y_x_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 54);

    auto ta1_y_x_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 55);

    auto ta1_y_x_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 56);

    auto ta1_y_x_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 57);

    auto ta1_y_x_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 58);

    auto ta1_y_x_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 59);

    auto ta1_y_y_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_pg + 60);

    auto ta1_y_y_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 61);

    auto ta1_y_y_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 62);

    auto ta1_y_y_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 63);

    auto ta1_y_y_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 64);

    auto ta1_y_y_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 65);

    auto ta1_y_y_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 66);

    auto ta1_y_y_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 67);

    auto ta1_y_y_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 68);

    auto ta1_y_y_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 69);

    auto ta1_y_y_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 70);

    auto ta1_y_y_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 71);

    auto ta1_y_y_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 72);

    auto ta1_y_y_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 73);

    auto ta1_y_y_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 74);

    auto ta1_y_z_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_pg + 75);

    auto ta1_y_z_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 76);

    auto ta1_y_z_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 77);

    auto ta1_y_z_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 78);

    auto ta1_y_z_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 79);

    auto ta1_y_z_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 80);

    auto ta1_y_z_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 81);

    auto ta1_y_z_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 82);

    auto ta1_y_z_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 83);

    auto ta1_y_z_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 84);

    auto ta1_y_z_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 85);

    auto ta1_y_z_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 86);

    auto ta1_y_z_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 87);

    auto ta1_y_z_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 88);

    auto ta1_y_z_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 89);

    auto ta1_z_x_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_pg + 90);

    auto ta1_z_x_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 91);

    auto ta1_z_x_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 92);

    auto ta1_z_x_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 93);

    auto ta1_z_x_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 94);

    auto ta1_z_x_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 95);

    auto ta1_z_x_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 96);

    auto ta1_z_x_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 97);

    auto ta1_z_x_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 98);

    auto ta1_z_x_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 99);

    auto ta1_z_x_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 100);

    auto ta1_z_x_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 101);

    auto ta1_z_x_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 102);

    auto ta1_z_x_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 103);

    auto ta1_z_x_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 104);

    auto ta1_z_y_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_pg + 105);

    auto ta1_z_y_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 106);

    auto ta1_z_y_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 107);

    auto ta1_z_y_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 108);

    auto ta1_z_y_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 109);

    auto ta1_z_y_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 110);

    auto ta1_z_y_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 111);

    auto ta1_z_y_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 112);

    auto ta1_z_y_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 113);

    auto ta1_z_y_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 114);

    auto ta1_z_y_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 115);

    auto ta1_z_y_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 116);

    auto ta1_z_y_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 117);

    auto ta1_z_y_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 118);

    auto ta1_z_y_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 119);

    auto ta1_z_z_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_pg + 120);

    auto ta1_z_z_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 121);

    auto ta1_z_z_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 122);

    auto ta1_z_z_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 123);

    auto ta1_z_z_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 124);

    auto ta1_z_z_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 125);

    auto ta1_z_z_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 126);

    auto ta1_z_z_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 127);

    auto ta1_z_z_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 128);

    auto ta1_z_z_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 129);

    auto ta1_z_z_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 130);

    auto ta1_z_z_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 131);

    auto ta1_z_z_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 132);

    auto ta1_z_z_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 133);

    auto ta1_z_z_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 134);

    // Set up components of auxiliary buffer : PG

    auto ta2_xx_x_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg);

    auto ta2_xx_x_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 1);

    auto ta2_xx_x_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 2);

    auto ta2_xx_x_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 3);

    auto ta2_xx_x_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 4);

    auto ta2_xx_x_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 5);

    auto ta2_xx_x_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 6);

    auto ta2_xx_x_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 7);

    auto ta2_xx_x_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 8);

    auto ta2_xx_x_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 9);

    auto ta2_xx_x_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 10);

    auto ta2_xx_x_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 11);

    auto ta2_xx_x_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 12);

    auto ta2_xx_x_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 13);

    auto ta2_xx_x_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 14);

    auto ta2_xx_y_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 15);

    auto ta2_xx_y_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 16);

    auto ta2_xx_y_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 17);

    auto ta2_xx_y_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 18);

    auto ta2_xx_y_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 19);

    auto ta2_xx_y_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 20);

    auto ta2_xx_y_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 21);

    auto ta2_xx_y_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 22);

    auto ta2_xx_y_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 23);

    auto ta2_xx_y_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 24);

    auto ta2_xx_y_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 25);

    auto ta2_xx_y_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 26);

    auto ta2_xx_y_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 27);

    auto ta2_xx_y_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 28);

    auto ta2_xx_y_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 29);

    auto ta2_xx_z_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 30);

    auto ta2_xx_z_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 31);

    auto ta2_xx_z_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 32);

    auto ta2_xx_z_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 33);

    auto ta2_xx_z_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 34);

    auto ta2_xx_z_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 35);

    auto ta2_xx_z_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 36);

    auto ta2_xx_z_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 37);

    auto ta2_xx_z_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 38);

    auto ta2_xx_z_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 39);

    auto ta2_xx_z_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 40);

    auto ta2_xx_z_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 41);

    auto ta2_xx_z_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 42);

    auto ta2_xx_z_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 43);

    auto ta2_xx_z_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 44);

    auto ta2_xy_x_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 45);

    auto ta2_xy_x_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 46);

    auto ta2_xy_x_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 47);

    auto ta2_xy_x_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 48);

    auto ta2_xy_x_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 49);

    auto ta2_xy_x_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 50);

    auto ta2_xy_x_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 51);

    auto ta2_xy_x_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 52);

    auto ta2_xy_x_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 53);

    auto ta2_xy_x_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 54);

    auto ta2_xy_x_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 55);

    auto ta2_xy_x_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 56);

    auto ta2_xy_x_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 57);

    auto ta2_xy_x_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 58);

    auto ta2_xy_x_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 59);

    auto ta2_xy_y_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 60);

    auto ta2_xy_y_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 61);

    auto ta2_xy_y_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 62);

    auto ta2_xy_y_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 63);

    auto ta2_xy_y_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 64);

    auto ta2_xy_y_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 65);

    auto ta2_xy_y_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 66);

    auto ta2_xy_y_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 67);

    auto ta2_xy_y_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 68);

    auto ta2_xy_y_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 69);

    auto ta2_xy_y_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 70);

    auto ta2_xy_y_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 71);

    auto ta2_xy_y_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 72);

    auto ta2_xy_y_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 73);

    auto ta2_xy_y_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 74);

    auto ta2_xy_z_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 75);

    auto ta2_xy_z_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 76);

    auto ta2_xy_z_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 77);

    auto ta2_xy_z_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 78);

    auto ta2_xy_z_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 79);

    auto ta2_xy_z_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 80);

    auto ta2_xy_z_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 81);

    auto ta2_xy_z_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 82);

    auto ta2_xy_z_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 83);

    auto ta2_xy_z_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 84);

    auto ta2_xy_z_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 85);

    auto ta2_xy_z_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 86);

    auto ta2_xy_z_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 87);

    auto ta2_xy_z_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 88);

    auto ta2_xy_z_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 89);

    auto ta2_xz_x_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 90);

    auto ta2_xz_x_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 91);

    auto ta2_xz_x_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 92);

    auto ta2_xz_x_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 93);

    auto ta2_xz_x_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 94);

    auto ta2_xz_x_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 95);

    auto ta2_xz_x_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 96);

    auto ta2_xz_x_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 97);

    auto ta2_xz_x_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 98);

    auto ta2_xz_x_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 99);

    auto ta2_xz_x_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 100);

    auto ta2_xz_x_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 101);

    auto ta2_xz_x_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 102);

    auto ta2_xz_x_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 103);

    auto ta2_xz_x_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 104);

    auto ta2_xz_y_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 105);

    auto ta2_xz_y_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 106);

    auto ta2_xz_y_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 107);

    auto ta2_xz_y_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 108);

    auto ta2_xz_y_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 109);

    auto ta2_xz_y_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 110);

    auto ta2_xz_y_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 111);

    auto ta2_xz_y_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 112);

    auto ta2_xz_y_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 113);

    auto ta2_xz_y_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 114);

    auto ta2_xz_y_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 115);

    auto ta2_xz_y_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 116);

    auto ta2_xz_y_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 117);

    auto ta2_xz_y_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 118);

    auto ta2_xz_y_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 119);

    auto ta2_xz_z_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 120);

    auto ta2_xz_z_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 121);

    auto ta2_xz_z_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 122);

    auto ta2_xz_z_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 123);

    auto ta2_xz_z_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 124);

    auto ta2_xz_z_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 125);

    auto ta2_xz_z_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 126);

    auto ta2_xz_z_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 127);

    auto ta2_xz_z_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 128);

    auto ta2_xz_z_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 129);

    auto ta2_xz_z_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 130);

    auto ta2_xz_z_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 131);

    auto ta2_xz_z_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 132);

    auto ta2_xz_z_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 133);

    auto ta2_xz_z_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 134);

    auto ta2_yy_x_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 135);

    auto ta2_yy_x_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 136);

    auto ta2_yy_x_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 137);

    auto ta2_yy_x_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 138);

    auto ta2_yy_x_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 139);

    auto ta2_yy_x_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 140);

    auto ta2_yy_x_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 141);

    auto ta2_yy_x_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 142);

    auto ta2_yy_x_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 143);

    auto ta2_yy_x_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 144);

    auto ta2_yy_x_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 145);

    auto ta2_yy_x_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 146);

    auto ta2_yy_x_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 147);

    auto ta2_yy_x_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 148);

    auto ta2_yy_x_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 149);

    auto ta2_yy_y_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 150);

    auto ta2_yy_y_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 151);

    auto ta2_yy_y_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 152);

    auto ta2_yy_y_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 153);

    auto ta2_yy_y_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 154);

    auto ta2_yy_y_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 155);

    auto ta2_yy_y_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 156);

    auto ta2_yy_y_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 157);

    auto ta2_yy_y_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 158);

    auto ta2_yy_y_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 159);

    auto ta2_yy_y_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 160);

    auto ta2_yy_y_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 161);

    auto ta2_yy_y_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 162);

    auto ta2_yy_y_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 163);

    auto ta2_yy_y_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 164);

    auto ta2_yy_z_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 165);

    auto ta2_yy_z_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 166);

    auto ta2_yy_z_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 167);

    auto ta2_yy_z_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 168);

    auto ta2_yy_z_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 169);

    auto ta2_yy_z_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 170);

    auto ta2_yy_z_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 171);

    auto ta2_yy_z_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 172);

    auto ta2_yy_z_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 173);

    auto ta2_yy_z_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 174);

    auto ta2_yy_z_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 175);

    auto ta2_yy_z_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 176);

    auto ta2_yy_z_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 177);

    auto ta2_yy_z_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 178);

    auto ta2_yy_z_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 179);

    auto ta2_yz_x_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 180);

    auto ta2_yz_x_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 181);

    auto ta2_yz_x_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 182);

    auto ta2_yz_x_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 183);

    auto ta2_yz_x_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 184);

    auto ta2_yz_x_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 185);

    auto ta2_yz_x_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 186);

    auto ta2_yz_x_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 187);

    auto ta2_yz_x_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 188);

    auto ta2_yz_x_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 189);

    auto ta2_yz_x_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 190);

    auto ta2_yz_x_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 191);

    auto ta2_yz_x_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 192);

    auto ta2_yz_x_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 193);

    auto ta2_yz_x_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 194);

    auto ta2_yz_y_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 195);

    auto ta2_yz_y_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 196);

    auto ta2_yz_y_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 197);

    auto ta2_yz_y_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 198);

    auto ta2_yz_y_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 199);

    auto ta2_yz_y_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 200);

    auto ta2_yz_y_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 201);

    auto ta2_yz_y_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 202);

    auto ta2_yz_y_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 203);

    auto ta2_yz_y_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 204);

    auto ta2_yz_y_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 205);

    auto ta2_yz_y_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 206);

    auto ta2_yz_y_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 207);

    auto ta2_yz_y_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 208);

    auto ta2_yz_y_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 209);

    auto ta2_yz_z_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 210);

    auto ta2_yz_z_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 211);

    auto ta2_yz_z_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 212);

    auto ta2_yz_z_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 213);

    auto ta2_yz_z_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 214);

    auto ta2_yz_z_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 215);

    auto ta2_yz_z_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 216);

    auto ta2_yz_z_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 217);

    auto ta2_yz_z_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 218);

    auto ta2_yz_z_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 219);

    auto ta2_yz_z_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 220);

    auto ta2_yz_z_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 221);

    auto ta2_yz_z_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 222);

    auto ta2_yz_z_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 223);

    auto ta2_yz_z_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 224);

    auto ta2_zz_x_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 225);

    auto ta2_zz_x_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 226);

    auto ta2_zz_x_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 227);

    auto ta2_zz_x_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 228);

    auto ta2_zz_x_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 229);

    auto ta2_zz_x_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 230);

    auto ta2_zz_x_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 231);

    auto ta2_zz_x_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 232);

    auto ta2_zz_x_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 233);

    auto ta2_zz_x_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 234);

    auto ta2_zz_x_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 235);

    auto ta2_zz_x_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 236);

    auto ta2_zz_x_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 237);

    auto ta2_zz_x_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 238);

    auto ta2_zz_x_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 239);

    auto ta2_zz_y_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 240);

    auto ta2_zz_y_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 241);

    auto ta2_zz_y_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 242);

    auto ta2_zz_y_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 243);

    auto ta2_zz_y_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 244);

    auto ta2_zz_y_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 245);

    auto ta2_zz_y_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 246);

    auto ta2_zz_y_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 247);

    auto ta2_zz_y_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 248);

    auto ta2_zz_y_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 249);

    auto ta2_zz_y_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 250);

    auto ta2_zz_y_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 251);

    auto ta2_zz_y_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 252);

    auto ta2_zz_y_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 253);

    auto ta2_zz_y_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 254);

    auto ta2_zz_z_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_pg + 255);

    auto ta2_zz_z_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 256);

    auto ta2_zz_z_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 257);

    auto ta2_zz_z_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 258);

    auto ta2_zz_z_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 259);

    auto ta2_zz_z_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 260);

    auto ta2_zz_z_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 261);

    auto ta2_zz_z_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 262);

    auto ta2_zz_z_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 263);

    auto ta2_zz_z_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 264);

    auto ta2_zz_z_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_pg + 265);

    auto ta2_zz_z_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 266);

    auto ta2_zz_z_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 267);

    auto ta2_zz_z_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 268);

    auto ta2_zz_z_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_pg + 269);

    // Set up components of auxiliary buffer : PG

    auto ta2_xx_x_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg);

    auto ta2_xx_x_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 1);

    auto ta2_xx_x_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 2);

    auto ta2_xx_x_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 3);

    auto ta2_xx_x_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 4);

    auto ta2_xx_x_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 5);

    auto ta2_xx_x_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 6);

    auto ta2_xx_x_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 7);

    auto ta2_xx_x_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 8);

    auto ta2_xx_x_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 9);

    auto ta2_xx_x_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 10);

    auto ta2_xx_x_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 11);

    auto ta2_xx_x_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 12);

    auto ta2_xx_x_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 13);

    auto ta2_xx_x_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 14);

    auto ta2_xx_y_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 15);

    auto ta2_xx_y_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 16);

    auto ta2_xx_y_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 17);

    auto ta2_xx_y_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 18);

    auto ta2_xx_y_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 19);

    auto ta2_xx_y_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 20);

    auto ta2_xx_y_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 21);

    auto ta2_xx_y_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 22);

    auto ta2_xx_y_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 23);

    auto ta2_xx_y_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 24);

    auto ta2_xx_y_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 25);

    auto ta2_xx_y_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 26);

    auto ta2_xx_y_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 27);

    auto ta2_xx_y_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 28);

    auto ta2_xx_y_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 29);

    auto ta2_xx_z_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 30);

    auto ta2_xx_z_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 31);

    auto ta2_xx_z_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 32);

    auto ta2_xx_z_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 33);

    auto ta2_xx_z_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 34);

    auto ta2_xx_z_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 35);

    auto ta2_xx_z_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 36);

    auto ta2_xx_z_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 37);

    auto ta2_xx_z_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 38);

    auto ta2_xx_z_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 39);

    auto ta2_xx_z_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 40);

    auto ta2_xx_z_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 41);

    auto ta2_xx_z_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 42);

    auto ta2_xx_z_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 43);

    auto ta2_xx_z_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 44);

    auto ta2_xy_x_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 45);

    auto ta2_xy_x_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 46);

    auto ta2_xy_x_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 47);

    auto ta2_xy_x_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 48);

    auto ta2_xy_x_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 49);

    auto ta2_xy_x_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 50);

    auto ta2_xy_x_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 51);

    auto ta2_xy_x_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 52);

    auto ta2_xy_x_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 53);

    auto ta2_xy_x_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 54);

    auto ta2_xy_x_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 55);

    auto ta2_xy_x_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 56);

    auto ta2_xy_x_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 57);

    auto ta2_xy_x_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 58);

    auto ta2_xy_x_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 59);

    auto ta2_xy_y_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 60);

    auto ta2_xy_y_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 61);

    auto ta2_xy_y_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 62);

    auto ta2_xy_y_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 63);

    auto ta2_xy_y_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 64);

    auto ta2_xy_y_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 65);

    auto ta2_xy_y_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 66);

    auto ta2_xy_y_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 67);

    auto ta2_xy_y_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 68);

    auto ta2_xy_y_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 69);

    auto ta2_xy_y_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 70);

    auto ta2_xy_y_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 71);

    auto ta2_xy_y_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 72);

    auto ta2_xy_y_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 73);

    auto ta2_xy_y_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 74);

    auto ta2_xy_z_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 75);

    auto ta2_xy_z_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 76);

    auto ta2_xy_z_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 77);

    auto ta2_xy_z_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 78);

    auto ta2_xy_z_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 79);

    auto ta2_xy_z_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 80);

    auto ta2_xy_z_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 81);

    auto ta2_xy_z_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 82);

    auto ta2_xy_z_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 83);

    auto ta2_xy_z_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 84);

    auto ta2_xy_z_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 85);

    auto ta2_xy_z_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 86);

    auto ta2_xy_z_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 87);

    auto ta2_xy_z_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 88);

    auto ta2_xy_z_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 89);

    auto ta2_xz_x_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 90);

    auto ta2_xz_x_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 91);

    auto ta2_xz_x_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 92);

    auto ta2_xz_x_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 93);

    auto ta2_xz_x_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 94);

    auto ta2_xz_x_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 95);

    auto ta2_xz_x_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 96);

    auto ta2_xz_x_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 97);

    auto ta2_xz_x_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 98);

    auto ta2_xz_x_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 99);

    auto ta2_xz_x_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 100);

    auto ta2_xz_x_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 101);

    auto ta2_xz_x_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 102);

    auto ta2_xz_x_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 103);

    auto ta2_xz_x_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 104);

    auto ta2_xz_y_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 105);

    auto ta2_xz_y_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 106);

    auto ta2_xz_y_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 107);

    auto ta2_xz_y_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 108);

    auto ta2_xz_y_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 109);

    auto ta2_xz_y_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 110);

    auto ta2_xz_y_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 111);

    auto ta2_xz_y_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 112);

    auto ta2_xz_y_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 113);

    auto ta2_xz_y_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 114);

    auto ta2_xz_y_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 115);

    auto ta2_xz_y_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 116);

    auto ta2_xz_y_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 117);

    auto ta2_xz_y_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 118);

    auto ta2_xz_y_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 119);

    auto ta2_xz_z_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 120);

    auto ta2_xz_z_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 121);

    auto ta2_xz_z_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 122);

    auto ta2_xz_z_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 123);

    auto ta2_xz_z_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 124);

    auto ta2_xz_z_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 125);

    auto ta2_xz_z_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 126);

    auto ta2_xz_z_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 127);

    auto ta2_xz_z_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 128);

    auto ta2_xz_z_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 129);

    auto ta2_xz_z_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 130);

    auto ta2_xz_z_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 131);

    auto ta2_xz_z_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 132);

    auto ta2_xz_z_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 133);

    auto ta2_xz_z_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 134);

    auto ta2_yy_x_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 135);

    auto ta2_yy_x_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 136);

    auto ta2_yy_x_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 137);

    auto ta2_yy_x_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 138);

    auto ta2_yy_x_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 139);

    auto ta2_yy_x_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 140);

    auto ta2_yy_x_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 141);

    auto ta2_yy_x_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 142);

    auto ta2_yy_x_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 143);

    auto ta2_yy_x_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 144);

    auto ta2_yy_x_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 145);

    auto ta2_yy_x_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 146);

    auto ta2_yy_x_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 147);

    auto ta2_yy_x_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 148);

    auto ta2_yy_x_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 149);

    auto ta2_yy_y_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 150);

    auto ta2_yy_y_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 151);

    auto ta2_yy_y_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 152);

    auto ta2_yy_y_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 153);

    auto ta2_yy_y_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 154);

    auto ta2_yy_y_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 155);

    auto ta2_yy_y_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 156);

    auto ta2_yy_y_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 157);

    auto ta2_yy_y_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 158);

    auto ta2_yy_y_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 159);

    auto ta2_yy_y_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 160);

    auto ta2_yy_y_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 161);

    auto ta2_yy_y_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 162);

    auto ta2_yy_y_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 163);

    auto ta2_yy_y_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 164);

    auto ta2_yy_z_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 165);

    auto ta2_yy_z_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 166);

    auto ta2_yy_z_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 167);

    auto ta2_yy_z_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 168);

    auto ta2_yy_z_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 169);

    auto ta2_yy_z_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 170);

    auto ta2_yy_z_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 171);

    auto ta2_yy_z_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 172);

    auto ta2_yy_z_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 173);

    auto ta2_yy_z_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 174);

    auto ta2_yy_z_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 175);

    auto ta2_yy_z_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 176);

    auto ta2_yy_z_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 177);

    auto ta2_yy_z_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 178);

    auto ta2_yy_z_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 179);

    auto ta2_yz_x_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 180);

    auto ta2_yz_x_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 181);

    auto ta2_yz_x_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 182);

    auto ta2_yz_x_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 183);

    auto ta2_yz_x_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 184);

    auto ta2_yz_x_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 185);

    auto ta2_yz_x_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 186);

    auto ta2_yz_x_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 187);

    auto ta2_yz_x_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 188);

    auto ta2_yz_x_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 189);

    auto ta2_yz_x_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 190);

    auto ta2_yz_x_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 191);

    auto ta2_yz_x_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 192);

    auto ta2_yz_x_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 193);

    auto ta2_yz_x_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 194);

    auto ta2_yz_y_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 195);

    auto ta2_yz_y_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 196);

    auto ta2_yz_y_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 197);

    auto ta2_yz_y_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 198);

    auto ta2_yz_y_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 199);

    auto ta2_yz_y_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 200);

    auto ta2_yz_y_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 201);

    auto ta2_yz_y_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 202);

    auto ta2_yz_y_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 203);

    auto ta2_yz_y_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 204);

    auto ta2_yz_y_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 205);

    auto ta2_yz_y_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 206);

    auto ta2_yz_y_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 207);

    auto ta2_yz_y_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 208);

    auto ta2_yz_y_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 209);

    auto ta2_yz_z_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 210);

    auto ta2_yz_z_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 211);

    auto ta2_yz_z_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 212);

    auto ta2_yz_z_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 213);

    auto ta2_yz_z_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 214);

    auto ta2_yz_z_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 215);

    auto ta2_yz_z_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 216);

    auto ta2_yz_z_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 217);

    auto ta2_yz_z_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 218);

    auto ta2_yz_z_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 219);

    auto ta2_yz_z_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 220);

    auto ta2_yz_z_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 221);

    auto ta2_yz_z_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 222);

    auto ta2_yz_z_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 223);

    auto ta2_yz_z_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 224);

    auto ta2_zz_x_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 225);

    auto ta2_zz_x_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 226);

    auto ta2_zz_x_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 227);

    auto ta2_zz_x_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 228);

    auto ta2_zz_x_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 229);

    auto ta2_zz_x_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 230);

    auto ta2_zz_x_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 231);

    auto ta2_zz_x_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 232);

    auto ta2_zz_x_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 233);

    auto ta2_zz_x_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 234);

    auto ta2_zz_x_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 235);

    auto ta2_zz_x_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 236);

    auto ta2_zz_x_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 237);

    auto ta2_zz_x_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 238);

    auto ta2_zz_x_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 239);

    auto ta2_zz_y_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 240);

    auto ta2_zz_y_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 241);

    auto ta2_zz_y_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 242);

    auto ta2_zz_y_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 243);

    auto ta2_zz_y_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 244);

    auto ta2_zz_y_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 245);

    auto ta2_zz_y_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 246);

    auto ta2_zz_y_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 247);

    auto ta2_zz_y_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 248);

    auto ta2_zz_y_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 249);

    auto ta2_zz_y_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 250);

    auto ta2_zz_y_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 251);

    auto ta2_zz_y_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 252);

    auto ta2_zz_y_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 253);

    auto ta2_zz_y_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 254);

    auto ta2_zz_z_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_pg + 255);

    auto ta2_zz_z_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 256);

    auto ta2_zz_z_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 257);

    auto ta2_zz_z_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 258);

    auto ta2_zz_z_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 259);

    auto ta2_zz_z_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 260);

    auto ta2_zz_z_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 261);

    auto ta2_zz_z_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 262);

    auto ta2_zz_z_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 263);

    auto ta2_zz_z_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 264);

    auto ta2_zz_z_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_pg + 265);

    auto ta2_zz_z_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 266);

    auto ta2_zz_z_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 267);

    auto ta2_zz_z_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 268);

    auto ta2_zz_z_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_pg + 269);

    // Set up 0-15 components of targeted buffer : DG

    auto ta2_xx_xx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg);

    auto ta2_xx_xx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 1);

    auto ta2_xx_xx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 2);

    auto ta2_xx_xx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 3);

    auto ta2_xx_xx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 4);

    auto ta2_xx_xx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 5);

    auto ta2_xx_xx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 6);

    auto ta2_xx_xx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 7);

    auto ta2_xx_xx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 8);

    auto ta2_xx_xx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 9);

    auto ta2_xx_xx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 10);

    auto ta2_xx_xx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 11);

    auto ta2_xx_xx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 12);

    auto ta2_xx_xx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 13);

    auto ta2_xx_xx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 14);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_x_x_xxxx_1,   \
                             ta1_x_x_xxxy_1,   \
                             ta1_x_x_xxxz_1,   \
                             ta1_x_x_xxyy_1,   \
                             ta1_x_x_xxyz_1,   \
                             ta1_x_x_xxzz_1,   \
                             ta1_x_x_xyyy_1,   \
                             ta1_x_x_xyyz_1,   \
                             ta1_x_x_xyzz_1,   \
                             ta1_x_x_xzzz_1,   \
                             ta1_x_x_yyyy_1,   \
                             ta1_x_x_yyyz_1,   \
                             ta1_x_x_yyzz_1,   \
                             ta1_x_x_yzzz_1,   \
                             ta1_x_x_zzzz_1,   \
                             ta2_xx_0_xxxx_0,  \
                             ta2_xx_0_xxxx_1,  \
                             ta2_xx_0_xxxy_0,  \
                             ta2_xx_0_xxxy_1,  \
                             ta2_xx_0_xxxz_0,  \
                             ta2_xx_0_xxxz_1,  \
                             ta2_xx_0_xxyy_0,  \
                             ta2_xx_0_xxyy_1,  \
                             ta2_xx_0_xxyz_0,  \
                             ta2_xx_0_xxyz_1,  \
                             ta2_xx_0_xxzz_0,  \
                             ta2_xx_0_xxzz_1,  \
                             ta2_xx_0_xyyy_0,  \
                             ta2_xx_0_xyyy_1,  \
                             ta2_xx_0_xyyz_0,  \
                             ta2_xx_0_xyyz_1,  \
                             ta2_xx_0_xyzz_0,  \
                             ta2_xx_0_xyzz_1,  \
                             ta2_xx_0_xzzz_0,  \
                             ta2_xx_0_xzzz_1,  \
                             ta2_xx_0_yyyy_0,  \
                             ta2_xx_0_yyyy_1,  \
                             ta2_xx_0_yyyz_0,  \
                             ta2_xx_0_yyyz_1,  \
                             ta2_xx_0_yyzz_0,  \
                             ta2_xx_0_yyzz_1,  \
                             ta2_xx_0_yzzz_0,  \
                             ta2_xx_0_yzzz_1,  \
                             ta2_xx_0_zzzz_0,  \
                             ta2_xx_0_zzzz_1,  \
                             ta2_xx_x_xxx_0,   \
                             ta2_xx_x_xxx_1,   \
                             ta2_xx_x_xxxx_0,  \
                             ta2_xx_x_xxxx_1,  \
                             ta2_xx_x_xxxy_0,  \
                             ta2_xx_x_xxxy_1,  \
                             ta2_xx_x_xxxz_0,  \
                             ta2_xx_x_xxxz_1,  \
                             ta2_xx_x_xxy_0,   \
                             ta2_xx_x_xxy_1,   \
                             ta2_xx_x_xxyy_0,  \
                             ta2_xx_x_xxyy_1,  \
                             ta2_xx_x_xxyz_0,  \
                             ta2_xx_x_xxyz_1,  \
                             ta2_xx_x_xxz_0,   \
                             ta2_xx_x_xxz_1,   \
                             ta2_xx_x_xxzz_0,  \
                             ta2_xx_x_xxzz_1,  \
                             ta2_xx_x_xyy_0,   \
                             ta2_xx_x_xyy_1,   \
                             ta2_xx_x_xyyy_0,  \
                             ta2_xx_x_xyyy_1,  \
                             ta2_xx_x_xyyz_0,  \
                             ta2_xx_x_xyyz_1,  \
                             ta2_xx_x_xyz_0,   \
                             ta2_xx_x_xyz_1,   \
                             ta2_xx_x_xyzz_0,  \
                             ta2_xx_x_xyzz_1,  \
                             ta2_xx_x_xzz_0,   \
                             ta2_xx_x_xzz_1,   \
                             ta2_xx_x_xzzz_0,  \
                             ta2_xx_x_xzzz_1,  \
                             ta2_xx_x_yyy_0,   \
                             ta2_xx_x_yyy_1,   \
                             ta2_xx_x_yyyy_0,  \
                             ta2_xx_x_yyyy_1,  \
                             ta2_xx_x_yyyz_0,  \
                             ta2_xx_x_yyyz_1,  \
                             ta2_xx_x_yyz_0,   \
                             ta2_xx_x_yyz_1,   \
                             ta2_xx_x_yyzz_0,  \
                             ta2_xx_x_yyzz_1,  \
                             ta2_xx_x_yzz_0,   \
                             ta2_xx_x_yzz_1,   \
                             ta2_xx_x_yzzz_0,  \
                             ta2_xx_x_yzzz_1,  \
                             ta2_xx_x_zzz_0,   \
                             ta2_xx_x_zzz_1,   \
                             ta2_xx_x_zzzz_0,  \
                             ta2_xx_x_zzzz_1,  \
                             ta2_xx_xx_xxxx_0, \
                             ta2_xx_xx_xxxy_0, \
                             ta2_xx_xx_xxxz_0, \
                             ta2_xx_xx_xxyy_0, \
                             ta2_xx_xx_xxyz_0, \
                             ta2_xx_xx_xxzz_0, \
                             ta2_xx_xx_xyyy_0, \
                             ta2_xx_xx_xyyz_0, \
                             ta2_xx_xx_xyzz_0, \
                             ta2_xx_xx_xzzz_0, \
                             ta2_xx_xx_yyyy_0, \
                             ta2_xx_xx_yyyz_0, \
                             ta2_xx_xx_yyzz_0, \
                             ta2_xx_xx_yzzz_0, \
                             ta2_xx_xx_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xx_xxxx_0[i] = ta2_xx_0_xxxx_0[i] * fe_0 - ta2_xx_0_xxxx_1[i] * fe_0 + 4.0 * ta2_xx_x_xxx_0[i] * fe_0 -
                              4.0 * ta2_xx_x_xxx_1[i] * fe_0 + 2.0 * ta1_x_x_xxxx_1[i] + ta2_xx_x_xxxx_0[i] * pa_x[i] - ta2_xx_x_xxxx_1[i] * pc_x[i];

        ta2_xx_xx_xxxy_0[i] = ta2_xx_0_xxxy_0[i] * fe_0 - ta2_xx_0_xxxy_1[i] * fe_0 + 3.0 * ta2_xx_x_xxy_0[i] * fe_0 -
                              3.0 * ta2_xx_x_xxy_1[i] * fe_0 + 2.0 * ta1_x_x_xxxy_1[i] + ta2_xx_x_xxxy_0[i] * pa_x[i] - ta2_xx_x_xxxy_1[i] * pc_x[i];

        ta2_xx_xx_xxxz_0[i] = ta2_xx_0_xxxz_0[i] * fe_0 - ta2_xx_0_xxxz_1[i] * fe_0 + 3.0 * ta2_xx_x_xxz_0[i] * fe_0 -
                              3.0 * ta2_xx_x_xxz_1[i] * fe_0 + 2.0 * ta1_x_x_xxxz_1[i] + ta2_xx_x_xxxz_0[i] * pa_x[i] - ta2_xx_x_xxxz_1[i] * pc_x[i];

        ta2_xx_xx_xxyy_0[i] = ta2_xx_0_xxyy_0[i] * fe_0 - ta2_xx_0_xxyy_1[i] * fe_0 + 2.0 * ta2_xx_x_xyy_0[i] * fe_0 -
                              2.0 * ta2_xx_x_xyy_1[i] * fe_0 + 2.0 * ta1_x_x_xxyy_1[i] + ta2_xx_x_xxyy_0[i] * pa_x[i] - ta2_xx_x_xxyy_1[i] * pc_x[i];

        ta2_xx_xx_xxyz_0[i] = ta2_xx_0_xxyz_0[i] * fe_0 - ta2_xx_0_xxyz_1[i] * fe_0 + 2.0 * ta2_xx_x_xyz_0[i] * fe_0 -
                              2.0 * ta2_xx_x_xyz_1[i] * fe_0 + 2.0 * ta1_x_x_xxyz_1[i] + ta2_xx_x_xxyz_0[i] * pa_x[i] - ta2_xx_x_xxyz_1[i] * pc_x[i];

        ta2_xx_xx_xxzz_0[i] = ta2_xx_0_xxzz_0[i] * fe_0 - ta2_xx_0_xxzz_1[i] * fe_0 + 2.0 * ta2_xx_x_xzz_0[i] * fe_0 -
                              2.0 * ta2_xx_x_xzz_1[i] * fe_0 + 2.0 * ta1_x_x_xxzz_1[i] + ta2_xx_x_xxzz_0[i] * pa_x[i] - ta2_xx_x_xxzz_1[i] * pc_x[i];

        ta2_xx_xx_xyyy_0[i] = ta2_xx_0_xyyy_0[i] * fe_0 - ta2_xx_0_xyyy_1[i] * fe_0 + ta2_xx_x_yyy_0[i] * fe_0 - ta2_xx_x_yyy_1[i] * fe_0 +
                              2.0 * ta1_x_x_xyyy_1[i] + ta2_xx_x_xyyy_0[i] * pa_x[i] - ta2_xx_x_xyyy_1[i] * pc_x[i];

        ta2_xx_xx_xyyz_0[i] = ta2_xx_0_xyyz_0[i] * fe_0 - ta2_xx_0_xyyz_1[i] * fe_0 + ta2_xx_x_yyz_0[i] * fe_0 - ta2_xx_x_yyz_1[i] * fe_0 +
                              2.0 * ta1_x_x_xyyz_1[i] + ta2_xx_x_xyyz_0[i] * pa_x[i] - ta2_xx_x_xyyz_1[i] * pc_x[i];

        ta2_xx_xx_xyzz_0[i] = ta2_xx_0_xyzz_0[i] * fe_0 - ta2_xx_0_xyzz_1[i] * fe_0 + ta2_xx_x_yzz_0[i] * fe_0 - ta2_xx_x_yzz_1[i] * fe_0 +
                              2.0 * ta1_x_x_xyzz_1[i] + ta2_xx_x_xyzz_0[i] * pa_x[i] - ta2_xx_x_xyzz_1[i] * pc_x[i];

        ta2_xx_xx_xzzz_0[i] = ta2_xx_0_xzzz_0[i] * fe_0 - ta2_xx_0_xzzz_1[i] * fe_0 + ta2_xx_x_zzz_0[i] * fe_0 - ta2_xx_x_zzz_1[i] * fe_0 +
                              2.0 * ta1_x_x_xzzz_1[i] + ta2_xx_x_xzzz_0[i] * pa_x[i] - ta2_xx_x_xzzz_1[i] * pc_x[i];

        ta2_xx_xx_yyyy_0[i] = ta2_xx_0_yyyy_0[i] * fe_0 - ta2_xx_0_yyyy_1[i] * fe_0 + 2.0 * ta1_x_x_yyyy_1[i] + ta2_xx_x_yyyy_0[i] * pa_x[i] -
                              ta2_xx_x_yyyy_1[i] * pc_x[i];

        ta2_xx_xx_yyyz_0[i] = ta2_xx_0_yyyz_0[i] * fe_0 - ta2_xx_0_yyyz_1[i] * fe_0 + 2.0 * ta1_x_x_yyyz_1[i] + ta2_xx_x_yyyz_0[i] * pa_x[i] -
                              ta2_xx_x_yyyz_1[i] * pc_x[i];

        ta2_xx_xx_yyzz_0[i] = ta2_xx_0_yyzz_0[i] * fe_0 - ta2_xx_0_yyzz_1[i] * fe_0 + 2.0 * ta1_x_x_yyzz_1[i] + ta2_xx_x_yyzz_0[i] * pa_x[i] -
                              ta2_xx_x_yyzz_1[i] * pc_x[i];

        ta2_xx_xx_yzzz_0[i] = ta2_xx_0_yzzz_0[i] * fe_0 - ta2_xx_0_yzzz_1[i] * fe_0 + 2.0 * ta1_x_x_yzzz_1[i] + ta2_xx_x_yzzz_0[i] * pa_x[i] -
                              ta2_xx_x_yzzz_1[i] * pc_x[i];

        ta2_xx_xx_zzzz_0[i] = ta2_xx_0_zzzz_0[i] * fe_0 - ta2_xx_0_zzzz_1[i] * fe_0 + 2.0 * ta1_x_x_zzzz_1[i] + ta2_xx_x_zzzz_0[i] * pa_x[i] -
                              ta2_xx_x_zzzz_1[i] * pc_x[i];
    }

    // Set up 15-30 components of targeted buffer : DG

    auto ta2_xx_xy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 15);

    auto ta2_xx_xy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 16);

    auto ta2_xx_xy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 17);

    auto ta2_xx_xy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 18);

    auto ta2_xx_xy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 19);

    auto ta2_xx_xy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 20);

    auto ta2_xx_xy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 21);

    auto ta2_xx_xy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 22);

    auto ta2_xx_xy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 23);

    auto ta2_xx_xy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 24);

    auto ta2_xx_xy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 25);

    auto ta2_xx_xy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 26);

    auto ta2_xx_xy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 27);

    auto ta2_xx_xy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 28);

    auto ta2_xx_xy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 29);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_x_y_yyyy_1,   \
                             ta1_x_y_yyyz_1,   \
                             ta1_x_y_yyzz_1,   \
                             ta1_x_y_yzzz_1,   \
                             ta2_xx_x_xxx_0,   \
                             ta2_xx_x_xxx_1,   \
                             ta2_xx_x_xxxx_0,  \
                             ta2_xx_x_xxxx_1,  \
                             ta2_xx_x_xxxy_0,  \
                             ta2_xx_x_xxxy_1,  \
                             ta2_xx_x_xxxz_0,  \
                             ta2_xx_x_xxxz_1,  \
                             ta2_xx_x_xxy_0,   \
                             ta2_xx_x_xxy_1,   \
                             ta2_xx_x_xxyy_0,  \
                             ta2_xx_x_xxyy_1,  \
                             ta2_xx_x_xxyz_0,  \
                             ta2_xx_x_xxyz_1,  \
                             ta2_xx_x_xxz_0,   \
                             ta2_xx_x_xxz_1,   \
                             ta2_xx_x_xxzz_0,  \
                             ta2_xx_x_xxzz_1,  \
                             ta2_xx_x_xyy_0,   \
                             ta2_xx_x_xyy_1,   \
                             ta2_xx_x_xyyy_0,  \
                             ta2_xx_x_xyyy_1,  \
                             ta2_xx_x_xyyz_0,  \
                             ta2_xx_x_xyyz_1,  \
                             ta2_xx_x_xyz_0,   \
                             ta2_xx_x_xyz_1,   \
                             ta2_xx_x_xyzz_0,  \
                             ta2_xx_x_xyzz_1,  \
                             ta2_xx_x_xzz_0,   \
                             ta2_xx_x_xzz_1,   \
                             ta2_xx_x_xzzz_0,  \
                             ta2_xx_x_xzzz_1,  \
                             ta2_xx_x_zzzz_0,  \
                             ta2_xx_x_zzzz_1,  \
                             ta2_xx_xy_xxxx_0, \
                             ta2_xx_xy_xxxy_0, \
                             ta2_xx_xy_xxxz_0, \
                             ta2_xx_xy_xxyy_0, \
                             ta2_xx_xy_xxyz_0, \
                             ta2_xx_xy_xxzz_0, \
                             ta2_xx_xy_xyyy_0, \
                             ta2_xx_xy_xyyz_0, \
                             ta2_xx_xy_xyzz_0, \
                             ta2_xx_xy_xzzz_0, \
                             ta2_xx_xy_yyyy_0, \
                             ta2_xx_xy_yyyz_0, \
                             ta2_xx_xy_yyzz_0, \
                             ta2_xx_xy_yzzz_0, \
                             ta2_xx_xy_zzzz_0, \
                             ta2_xx_y_yyyy_0,  \
                             ta2_xx_y_yyyy_1,  \
                             ta2_xx_y_yyyz_0,  \
                             ta2_xx_y_yyyz_1,  \
                             ta2_xx_y_yyzz_0,  \
                             ta2_xx_y_yyzz_1,  \
                             ta2_xx_y_yzzz_0,  \
                             ta2_xx_y_yzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xy_xxxx_0[i] = ta2_xx_x_xxxx_0[i] * pa_y[i] - ta2_xx_x_xxxx_1[i] * pc_y[i];

        ta2_xx_xy_xxxy_0[i] = ta2_xx_x_xxx_0[i] * fe_0 - ta2_xx_x_xxx_1[i] * fe_0 + ta2_xx_x_xxxy_0[i] * pa_y[i] - ta2_xx_x_xxxy_1[i] * pc_y[i];

        ta2_xx_xy_xxxz_0[i] = ta2_xx_x_xxxz_0[i] * pa_y[i] - ta2_xx_x_xxxz_1[i] * pc_y[i];

        ta2_xx_xy_xxyy_0[i] =
            2.0 * ta2_xx_x_xxy_0[i] * fe_0 - 2.0 * ta2_xx_x_xxy_1[i] * fe_0 + ta2_xx_x_xxyy_0[i] * pa_y[i] - ta2_xx_x_xxyy_1[i] * pc_y[i];

        ta2_xx_xy_xxyz_0[i] = ta2_xx_x_xxz_0[i] * fe_0 - ta2_xx_x_xxz_1[i] * fe_0 + ta2_xx_x_xxyz_0[i] * pa_y[i] - ta2_xx_x_xxyz_1[i] * pc_y[i];

        ta2_xx_xy_xxzz_0[i] = ta2_xx_x_xxzz_0[i] * pa_y[i] - ta2_xx_x_xxzz_1[i] * pc_y[i];

        ta2_xx_xy_xyyy_0[i] =
            3.0 * ta2_xx_x_xyy_0[i] * fe_0 - 3.0 * ta2_xx_x_xyy_1[i] * fe_0 + ta2_xx_x_xyyy_0[i] * pa_y[i] - ta2_xx_x_xyyy_1[i] * pc_y[i];

        ta2_xx_xy_xyyz_0[i] =
            2.0 * ta2_xx_x_xyz_0[i] * fe_0 - 2.0 * ta2_xx_x_xyz_1[i] * fe_0 + ta2_xx_x_xyyz_0[i] * pa_y[i] - ta2_xx_x_xyyz_1[i] * pc_y[i];

        ta2_xx_xy_xyzz_0[i] = ta2_xx_x_xzz_0[i] * fe_0 - ta2_xx_x_xzz_1[i] * fe_0 + ta2_xx_x_xyzz_0[i] * pa_y[i] - ta2_xx_x_xyzz_1[i] * pc_y[i];

        ta2_xx_xy_xzzz_0[i] = ta2_xx_x_xzzz_0[i] * pa_y[i] - ta2_xx_x_xzzz_1[i] * pc_y[i];

        ta2_xx_xy_yyyy_0[i] = 2.0 * ta1_x_y_yyyy_1[i] + ta2_xx_y_yyyy_0[i] * pa_x[i] - ta2_xx_y_yyyy_1[i] * pc_x[i];

        ta2_xx_xy_yyyz_0[i] = 2.0 * ta1_x_y_yyyz_1[i] + ta2_xx_y_yyyz_0[i] * pa_x[i] - ta2_xx_y_yyyz_1[i] * pc_x[i];

        ta2_xx_xy_yyzz_0[i] = 2.0 * ta1_x_y_yyzz_1[i] + ta2_xx_y_yyzz_0[i] * pa_x[i] - ta2_xx_y_yyzz_1[i] * pc_x[i];

        ta2_xx_xy_yzzz_0[i] = 2.0 * ta1_x_y_yzzz_1[i] + ta2_xx_y_yzzz_0[i] * pa_x[i] - ta2_xx_y_yzzz_1[i] * pc_x[i];

        ta2_xx_xy_zzzz_0[i] = ta2_xx_x_zzzz_0[i] * pa_y[i] - ta2_xx_x_zzzz_1[i] * pc_y[i];
    }

    // Set up 30-45 components of targeted buffer : DG

    auto ta2_xx_xz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 30);

    auto ta2_xx_xz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 31);

    auto ta2_xx_xz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 32);

    auto ta2_xx_xz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 33);

    auto ta2_xx_xz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 34);

    auto ta2_xx_xz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 35);

    auto ta2_xx_xz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 36);

    auto ta2_xx_xz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 37);

    auto ta2_xx_xz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 38);

    auto ta2_xx_xz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 39);

    auto ta2_xx_xz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 40);

    auto ta2_xx_xz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 41);

    auto ta2_xx_xz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 42);

    auto ta2_xx_xz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 43);

    auto ta2_xx_xz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 44);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_x_z_yyyz_1,   \
                             ta1_x_z_yyzz_1,   \
                             ta1_x_z_yzzz_1,   \
                             ta1_x_z_zzzz_1,   \
                             ta2_xx_x_xxx_0,   \
                             ta2_xx_x_xxx_1,   \
                             ta2_xx_x_xxxx_0,  \
                             ta2_xx_x_xxxx_1,  \
                             ta2_xx_x_xxxy_0,  \
                             ta2_xx_x_xxxy_1,  \
                             ta2_xx_x_xxxz_0,  \
                             ta2_xx_x_xxxz_1,  \
                             ta2_xx_x_xxy_0,   \
                             ta2_xx_x_xxy_1,   \
                             ta2_xx_x_xxyy_0,  \
                             ta2_xx_x_xxyy_1,  \
                             ta2_xx_x_xxyz_0,  \
                             ta2_xx_x_xxyz_1,  \
                             ta2_xx_x_xxz_0,   \
                             ta2_xx_x_xxz_1,   \
                             ta2_xx_x_xxzz_0,  \
                             ta2_xx_x_xxzz_1,  \
                             ta2_xx_x_xyy_0,   \
                             ta2_xx_x_xyy_1,   \
                             ta2_xx_x_xyyy_0,  \
                             ta2_xx_x_xyyy_1,  \
                             ta2_xx_x_xyyz_0,  \
                             ta2_xx_x_xyyz_1,  \
                             ta2_xx_x_xyz_0,   \
                             ta2_xx_x_xyz_1,   \
                             ta2_xx_x_xyzz_0,  \
                             ta2_xx_x_xyzz_1,  \
                             ta2_xx_x_xzz_0,   \
                             ta2_xx_x_xzz_1,   \
                             ta2_xx_x_xzzz_0,  \
                             ta2_xx_x_xzzz_1,  \
                             ta2_xx_x_yyyy_0,  \
                             ta2_xx_x_yyyy_1,  \
                             ta2_xx_xz_xxxx_0, \
                             ta2_xx_xz_xxxy_0, \
                             ta2_xx_xz_xxxz_0, \
                             ta2_xx_xz_xxyy_0, \
                             ta2_xx_xz_xxyz_0, \
                             ta2_xx_xz_xxzz_0, \
                             ta2_xx_xz_xyyy_0, \
                             ta2_xx_xz_xyyz_0, \
                             ta2_xx_xz_xyzz_0, \
                             ta2_xx_xz_xzzz_0, \
                             ta2_xx_xz_yyyy_0, \
                             ta2_xx_xz_yyyz_0, \
                             ta2_xx_xz_yyzz_0, \
                             ta2_xx_xz_yzzz_0, \
                             ta2_xx_xz_zzzz_0, \
                             ta2_xx_z_yyyz_0,  \
                             ta2_xx_z_yyyz_1,  \
                             ta2_xx_z_yyzz_0,  \
                             ta2_xx_z_yyzz_1,  \
                             ta2_xx_z_yzzz_0,  \
                             ta2_xx_z_yzzz_1,  \
                             ta2_xx_z_zzzz_0,  \
                             ta2_xx_z_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xz_xxxx_0[i] = ta2_xx_x_xxxx_0[i] * pa_z[i] - ta2_xx_x_xxxx_1[i] * pc_z[i];

        ta2_xx_xz_xxxy_0[i] = ta2_xx_x_xxxy_0[i] * pa_z[i] - ta2_xx_x_xxxy_1[i] * pc_z[i];

        ta2_xx_xz_xxxz_0[i] = ta2_xx_x_xxx_0[i] * fe_0 - ta2_xx_x_xxx_1[i] * fe_0 + ta2_xx_x_xxxz_0[i] * pa_z[i] - ta2_xx_x_xxxz_1[i] * pc_z[i];

        ta2_xx_xz_xxyy_0[i] = ta2_xx_x_xxyy_0[i] * pa_z[i] - ta2_xx_x_xxyy_1[i] * pc_z[i];

        ta2_xx_xz_xxyz_0[i] = ta2_xx_x_xxy_0[i] * fe_0 - ta2_xx_x_xxy_1[i] * fe_0 + ta2_xx_x_xxyz_0[i] * pa_z[i] - ta2_xx_x_xxyz_1[i] * pc_z[i];

        ta2_xx_xz_xxzz_0[i] =
            2.0 * ta2_xx_x_xxz_0[i] * fe_0 - 2.0 * ta2_xx_x_xxz_1[i] * fe_0 + ta2_xx_x_xxzz_0[i] * pa_z[i] - ta2_xx_x_xxzz_1[i] * pc_z[i];

        ta2_xx_xz_xyyy_0[i] = ta2_xx_x_xyyy_0[i] * pa_z[i] - ta2_xx_x_xyyy_1[i] * pc_z[i];

        ta2_xx_xz_xyyz_0[i] = ta2_xx_x_xyy_0[i] * fe_0 - ta2_xx_x_xyy_1[i] * fe_0 + ta2_xx_x_xyyz_0[i] * pa_z[i] - ta2_xx_x_xyyz_1[i] * pc_z[i];

        ta2_xx_xz_xyzz_0[i] =
            2.0 * ta2_xx_x_xyz_0[i] * fe_0 - 2.0 * ta2_xx_x_xyz_1[i] * fe_0 + ta2_xx_x_xyzz_0[i] * pa_z[i] - ta2_xx_x_xyzz_1[i] * pc_z[i];

        ta2_xx_xz_xzzz_0[i] =
            3.0 * ta2_xx_x_xzz_0[i] * fe_0 - 3.0 * ta2_xx_x_xzz_1[i] * fe_0 + ta2_xx_x_xzzz_0[i] * pa_z[i] - ta2_xx_x_xzzz_1[i] * pc_z[i];

        ta2_xx_xz_yyyy_0[i] = ta2_xx_x_yyyy_0[i] * pa_z[i] - ta2_xx_x_yyyy_1[i] * pc_z[i];

        ta2_xx_xz_yyyz_0[i] = 2.0 * ta1_x_z_yyyz_1[i] + ta2_xx_z_yyyz_0[i] * pa_x[i] - ta2_xx_z_yyyz_1[i] * pc_x[i];

        ta2_xx_xz_yyzz_0[i] = 2.0 * ta1_x_z_yyzz_1[i] + ta2_xx_z_yyzz_0[i] * pa_x[i] - ta2_xx_z_yyzz_1[i] * pc_x[i];

        ta2_xx_xz_yzzz_0[i] = 2.0 * ta1_x_z_yzzz_1[i] + ta2_xx_z_yzzz_0[i] * pa_x[i] - ta2_xx_z_yzzz_1[i] * pc_x[i];

        ta2_xx_xz_zzzz_0[i] = 2.0 * ta1_x_z_zzzz_1[i] + ta2_xx_z_zzzz_0[i] * pa_x[i] - ta2_xx_z_zzzz_1[i] * pc_x[i];
    }

    // Set up 45-60 components of targeted buffer : DG

    auto ta2_xx_yy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 45);

    auto ta2_xx_yy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 46);

    auto ta2_xx_yy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 47);

    auto ta2_xx_yy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 48);

    auto ta2_xx_yy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 49);

    auto ta2_xx_yy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 50);

    auto ta2_xx_yy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 51);

    auto ta2_xx_yy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 52);

    auto ta2_xx_yy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 53);

    auto ta2_xx_yy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 54);

    auto ta2_xx_yy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 55);

    auto ta2_xx_yy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 56);

    auto ta2_xx_yy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 57);

    auto ta2_xx_yy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 58);

    auto ta2_xx_yy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 59);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta2_xx_0_xxxx_0,  \
                             ta2_xx_0_xxxx_1,  \
                             ta2_xx_0_xxxy_0,  \
                             ta2_xx_0_xxxy_1,  \
                             ta2_xx_0_xxxz_0,  \
                             ta2_xx_0_xxxz_1,  \
                             ta2_xx_0_xxyy_0,  \
                             ta2_xx_0_xxyy_1,  \
                             ta2_xx_0_xxyz_0,  \
                             ta2_xx_0_xxyz_1,  \
                             ta2_xx_0_xxzz_0,  \
                             ta2_xx_0_xxzz_1,  \
                             ta2_xx_0_xyyy_0,  \
                             ta2_xx_0_xyyy_1,  \
                             ta2_xx_0_xyyz_0,  \
                             ta2_xx_0_xyyz_1,  \
                             ta2_xx_0_xyzz_0,  \
                             ta2_xx_0_xyzz_1,  \
                             ta2_xx_0_xzzz_0,  \
                             ta2_xx_0_xzzz_1,  \
                             ta2_xx_0_yyyy_0,  \
                             ta2_xx_0_yyyy_1,  \
                             ta2_xx_0_yyyz_0,  \
                             ta2_xx_0_yyyz_1,  \
                             ta2_xx_0_yyzz_0,  \
                             ta2_xx_0_yyzz_1,  \
                             ta2_xx_0_yzzz_0,  \
                             ta2_xx_0_yzzz_1,  \
                             ta2_xx_0_zzzz_0,  \
                             ta2_xx_0_zzzz_1,  \
                             ta2_xx_y_xxx_0,   \
                             ta2_xx_y_xxx_1,   \
                             ta2_xx_y_xxxx_0,  \
                             ta2_xx_y_xxxx_1,  \
                             ta2_xx_y_xxxy_0,  \
                             ta2_xx_y_xxxy_1,  \
                             ta2_xx_y_xxxz_0,  \
                             ta2_xx_y_xxxz_1,  \
                             ta2_xx_y_xxy_0,   \
                             ta2_xx_y_xxy_1,   \
                             ta2_xx_y_xxyy_0,  \
                             ta2_xx_y_xxyy_1,  \
                             ta2_xx_y_xxyz_0,  \
                             ta2_xx_y_xxyz_1,  \
                             ta2_xx_y_xxz_0,   \
                             ta2_xx_y_xxz_1,   \
                             ta2_xx_y_xxzz_0,  \
                             ta2_xx_y_xxzz_1,  \
                             ta2_xx_y_xyy_0,   \
                             ta2_xx_y_xyy_1,   \
                             ta2_xx_y_xyyy_0,  \
                             ta2_xx_y_xyyy_1,  \
                             ta2_xx_y_xyyz_0,  \
                             ta2_xx_y_xyyz_1,  \
                             ta2_xx_y_xyz_0,   \
                             ta2_xx_y_xyz_1,   \
                             ta2_xx_y_xyzz_0,  \
                             ta2_xx_y_xyzz_1,  \
                             ta2_xx_y_xzz_0,   \
                             ta2_xx_y_xzz_1,   \
                             ta2_xx_y_xzzz_0,  \
                             ta2_xx_y_xzzz_1,  \
                             ta2_xx_y_yyy_0,   \
                             ta2_xx_y_yyy_1,   \
                             ta2_xx_y_yyyy_0,  \
                             ta2_xx_y_yyyy_1,  \
                             ta2_xx_y_yyyz_0,  \
                             ta2_xx_y_yyyz_1,  \
                             ta2_xx_y_yyz_0,   \
                             ta2_xx_y_yyz_1,   \
                             ta2_xx_y_yyzz_0,  \
                             ta2_xx_y_yyzz_1,  \
                             ta2_xx_y_yzz_0,   \
                             ta2_xx_y_yzz_1,   \
                             ta2_xx_y_yzzz_0,  \
                             ta2_xx_y_yzzz_1,  \
                             ta2_xx_y_zzz_0,   \
                             ta2_xx_y_zzz_1,   \
                             ta2_xx_y_zzzz_0,  \
                             ta2_xx_y_zzzz_1,  \
                             ta2_xx_yy_xxxx_0, \
                             ta2_xx_yy_xxxy_0, \
                             ta2_xx_yy_xxxz_0, \
                             ta2_xx_yy_xxyy_0, \
                             ta2_xx_yy_xxyz_0, \
                             ta2_xx_yy_xxzz_0, \
                             ta2_xx_yy_xyyy_0, \
                             ta2_xx_yy_xyyz_0, \
                             ta2_xx_yy_xyzz_0, \
                             ta2_xx_yy_xzzz_0, \
                             ta2_xx_yy_yyyy_0, \
                             ta2_xx_yy_yyyz_0, \
                             ta2_xx_yy_yyzz_0, \
                             ta2_xx_yy_yzzz_0, \
                             ta2_xx_yy_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yy_xxxx_0[i] = ta2_xx_0_xxxx_0[i] * fe_0 - ta2_xx_0_xxxx_1[i] * fe_0 + ta2_xx_y_xxxx_0[i] * pa_y[i] - ta2_xx_y_xxxx_1[i] * pc_y[i];

        ta2_xx_yy_xxxy_0[i] = ta2_xx_0_xxxy_0[i] * fe_0 - ta2_xx_0_xxxy_1[i] * fe_0 + ta2_xx_y_xxx_0[i] * fe_0 - ta2_xx_y_xxx_1[i] * fe_0 +
                              ta2_xx_y_xxxy_0[i] * pa_y[i] - ta2_xx_y_xxxy_1[i] * pc_y[i];

        ta2_xx_yy_xxxz_0[i] = ta2_xx_0_xxxz_0[i] * fe_0 - ta2_xx_0_xxxz_1[i] * fe_0 + ta2_xx_y_xxxz_0[i] * pa_y[i] - ta2_xx_y_xxxz_1[i] * pc_y[i];

        ta2_xx_yy_xxyy_0[i] = ta2_xx_0_xxyy_0[i] * fe_0 - ta2_xx_0_xxyy_1[i] * fe_0 + 2.0 * ta2_xx_y_xxy_0[i] * fe_0 -
                              2.0 * ta2_xx_y_xxy_1[i] * fe_0 + ta2_xx_y_xxyy_0[i] * pa_y[i] - ta2_xx_y_xxyy_1[i] * pc_y[i];

        ta2_xx_yy_xxyz_0[i] = ta2_xx_0_xxyz_0[i] * fe_0 - ta2_xx_0_xxyz_1[i] * fe_0 + ta2_xx_y_xxz_0[i] * fe_0 - ta2_xx_y_xxz_1[i] * fe_0 +
                              ta2_xx_y_xxyz_0[i] * pa_y[i] - ta2_xx_y_xxyz_1[i] * pc_y[i];

        ta2_xx_yy_xxzz_0[i] = ta2_xx_0_xxzz_0[i] * fe_0 - ta2_xx_0_xxzz_1[i] * fe_0 + ta2_xx_y_xxzz_0[i] * pa_y[i] - ta2_xx_y_xxzz_1[i] * pc_y[i];

        ta2_xx_yy_xyyy_0[i] = ta2_xx_0_xyyy_0[i] * fe_0 - ta2_xx_0_xyyy_1[i] * fe_0 + 3.0 * ta2_xx_y_xyy_0[i] * fe_0 -
                              3.0 * ta2_xx_y_xyy_1[i] * fe_0 + ta2_xx_y_xyyy_0[i] * pa_y[i] - ta2_xx_y_xyyy_1[i] * pc_y[i];

        ta2_xx_yy_xyyz_0[i] = ta2_xx_0_xyyz_0[i] * fe_0 - ta2_xx_0_xyyz_1[i] * fe_0 + 2.0 * ta2_xx_y_xyz_0[i] * fe_0 -
                              2.0 * ta2_xx_y_xyz_1[i] * fe_0 + ta2_xx_y_xyyz_0[i] * pa_y[i] - ta2_xx_y_xyyz_1[i] * pc_y[i];

        ta2_xx_yy_xyzz_0[i] = ta2_xx_0_xyzz_0[i] * fe_0 - ta2_xx_0_xyzz_1[i] * fe_0 + ta2_xx_y_xzz_0[i] * fe_0 - ta2_xx_y_xzz_1[i] * fe_0 +
                              ta2_xx_y_xyzz_0[i] * pa_y[i] - ta2_xx_y_xyzz_1[i] * pc_y[i];

        ta2_xx_yy_xzzz_0[i] = ta2_xx_0_xzzz_0[i] * fe_0 - ta2_xx_0_xzzz_1[i] * fe_0 + ta2_xx_y_xzzz_0[i] * pa_y[i] - ta2_xx_y_xzzz_1[i] * pc_y[i];

        ta2_xx_yy_yyyy_0[i] = ta2_xx_0_yyyy_0[i] * fe_0 - ta2_xx_0_yyyy_1[i] * fe_0 + 4.0 * ta2_xx_y_yyy_0[i] * fe_0 -
                              4.0 * ta2_xx_y_yyy_1[i] * fe_0 + ta2_xx_y_yyyy_0[i] * pa_y[i] - ta2_xx_y_yyyy_1[i] * pc_y[i];

        ta2_xx_yy_yyyz_0[i] = ta2_xx_0_yyyz_0[i] * fe_0 - ta2_xx_0_yyyz_1[i] * fe_0 + 3.0 * ta2_xx_y_yyz_0[i] * fe_0 -
                              3.0 * ta2_xx_y_yyz_1[i] * fe_0 + ta2_xx_y_yyyz_0[i] * pa_y[i] - ta2_xx_y_yyyz_1[i] * pc_y[i];

        ta2_xx_yy_yyzz_0[i] = ta2_xx_0_yyzz_0[i] * fe_0 - ta2_xx_0_yyzz_1[i] * fe_0 + 2.0 * ta2_xx_y_yzz_0[i] * fe_0 -
                              2.0 * ta2_xx_y_yzz_1[i] * fe_0 + ta2_xx_y_yyzz_0[i] * pa_y[i] - ta2_xx_y_yyzz_1[i] * pc_y[i];

        ta2_xx_yy_yzzz_0[i] = ta2_xx_0_yzzz_0[i] * fe_0 - ta2_xx_0_yzzz_1[i] * fe_0 + ta2_xx_y_zzz_0[i] * fe_0 - ta2_xx_y_zzz_1[i] * fe_0 +
                              ta2_xx_y_yzzz_0[i] * pa_y[i] - ta2_xx_y_yzzz_1[i] * pc_y[i];

        ta2_xx_yy_zzzz_0[i] = ta2_xx_0_zzzz_0[i] * fe_0 - ta2_xx_0_zzzz_1[i] * fe_0 + ta2_xx_y_zzzz_0[i] * pa_y[i] - ta2_xx_y_zzzz_1[i] * pc_y[i];
    }

    // Set up 60-75 components of targeted buffer : DG

    auto ta2_xx_yz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 60);

    auto ta2_xx_yz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 61);

    auto ta2_xx_yz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 62);

    auto ta2_xx_yz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 63);

    auto ta2_xx_yz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 64);

    auto ta2_xx_yz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 65);

    auto ta2_xx_yz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 66);

    auto ta2_xx_yz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 67);

    auto ta2_xx_yz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 68);

    auto ta2_xx_yz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 69);

    auto ta2_xx_yz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 70);

    auto ta2_xx_yz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 71);

    auto ta2_xx_yz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 72);

    auto ta2_xx_yz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 73);

    auto ta2_xx_yz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 74);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta2_xx_y_xxxy_0,  \
                             ta2_xx_y_xxxy_1,  \
                             ta2_xx_y_xxyy_0,  \
                             ta2_xx_y_xxyy_1,  \
                             ta2_xx_y_xyyy_0,  \
                             ta2_xx_y_xyyy_1,  \
                             ta2_xx_y_yyyy_0,  \
                             ta2_xx_y_yyyy_1,  \
                             ta2_xx_yz_xxxx_0, \
                             ta2_xx_yz_xxxy_0, \
                             ta2_xx_yz_xxxz_0, \
                             ta2_xx_yz_xxyy_0, \
                             ta2_xx_yz_xxyz_0, \
                             ta2_xx_yz_xxzz_0, \
                             ta2_xx_yz_xyyy_0, \
                             ta2_xx_yz_xyyz_0, \
                             ta2_xx_yz_xyzz_0, \
                             ta2_xx_yz_xzzz_0, \
                             ta2_xx_yz_yyyy_0, \
                             ta2_xx_yz_yyyz_0, \
                             ta2_xx_yz_yyzz_0, \
                             ta2_xx_yz_yzzz_0, \
                             ta2_xx_yz_zzzz_0, \
                             ta2_xx_z_xxxx_0,  \
                             ta2_xx_z_xxxx_1,  \
                             ta2_xx_z_xxxz_0,  \
                             ta2_xx_z_xxxz_1,  \
                             ta2_xx_z_xxyz_0,  \
                             ta2_xx_z_xxyz_1,  \
                             ta2_xx_z_xxz_0,   \
                             ta2_xx_z_xxz_1,   \
                             ta2_xx_z_xxzz_0,  \
                             ta2_xx_z_xxzz_1,  \
                             ta2_xx_z_xyyz_0,  \
                             ta2_xx_z_xyyz_1,  \
                             ta2_xx_z_xyz_0,   \
                             ta2_xx_z_xyz_1,   \
                             ta2_xx_z_xyzz_0,  \
                             ta2_xx_z_xyzz_1,  \
                             ta2_xx_z_xzz_0,   \
                             ta2_xx_z_xzz_1,   \
                             ta2_xx_z_xzzz_0,  \
                             ta2_xx_z_xzzz_1,  \
                             ta2_xx_z_yyyz_0,  \
                             ta2_xx_z_yyyz_1,  \
                             ta2_xx_z_yyz_0,   \
                             ta2_xx_z_yyz_1,   \
                             ta2_xx_z_yyzz_0,  \
                             ta2_xx_z_yyzz_1,  \
                             ta2_xx_z_yzz_0,   \
                             ta2_xx_z_yzz_1,   \
                             ta2_xx_z_yzzz_0,  \
                             ta2_xx_z_yzzz_1,  \
                             ta2_xx_z_zzz_0,   \
                             ta2_xx_z_zzz_1,   \
                             ta2_xx_z_zzzz_0,  \
                             ta2_xx_z_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yz_xxxx_0[i] = ta2_xx_z_xxxx_0[i] * pa_y[i] - ta2_xx_z_xxxx_1[i] * pc_y[i];

        ta2_xx_yz_xxxy_0[i] = ta2_xx_y_xxxy_0[i] * pa_z[i] - ta2_xx_y_xxxy_1[i] * pc_z[i];

        ta2_xx_yz_xxxz_0[i] = ta2_xx_z_xxxz_0[i] * pa_y[i] - ta2_xx_z_xxxz_1[i] * pc_y[i];

        ta2_xx_yz_xxyy_0[i] = ta2_xx_y_xxyy_0[i] * pa_z[i] - ta2_xx_y_xxyy_1[i] * pc_z[i];

        ta2_xx_yz_xxyz_0[i] = ta2_xx_z_xxz_0[i] * fe_0 - ta2_xx_z_xxz_1[i] * fe_0 + ta2_xx_z_xxyz_0[i] * pa_y[i] - ta2_xx_z_xxyz_1[i] * pc_y[i];

        ta2_xx_yz_xxzz_0[i] = ta2_xx_z_xxzz_0[i] * pa_y[i] - ta2_xx_z_xxzz_1[i] * pc_y[i];

        ta2_xx_yz_xyyy_0[i] = ta2_xx_y_xyyy_0[i] * pa_z[i] - ta2_xx_y_xyyy_1[i] * pc_z[i];

        ta2_xx_yz_xyyz_0[i] =
            2.0 * ta2_xx_z_xyz_0[i] * fe_0 - 2.0 * ta2_xx_z_xyz_1[i] * fe_0 + ta2_xx_z_xyyz_0[i] * pa_y[i] - ta2_xx_z_xyyz_1[i] * pc_y[i];

        ta2_xx_yz_xyzz_0[i] = ta2_xx_z_xzz_0[i] * fe_0 - ta2_xx_z_xzz_1[i] * fe_0 + ta2_xx_z_xyzz_0[i] * pa_y[i] - ta2_xx_z_xyzz_1[i] * pc_y[i];

        ta2_xx_yz_xzzz_0[i] = ta2_xx_z_xzzz_0[i] * pa_y[i] - ta2_xx_z_xzzz_1[i] * pc_y[i];

        ta2_xx_yz_yyyy_0[i] = ta2_xx_y_yyyy_0[i] * pa_z[i] - ta2_xx_y_yyyy_1[i] * pc_z[i];

        ta2_xx_yz_yyyz_0[i] =
            3.0 * ta2_xx_z_yyz_0[i] * fe_0 - 3.0 * ta2_xx_z_yyz_1[i] * fe_0 + ta2_xx_z_yyyz_0[i] * pa_y[i] - ta2_xx_z_yyyz_1[i] * pc_y[i];

        ta2_xx_yz_yyzz_0[i] =
            2.0 * ta2_xx_z_yzz_0[i] * fe_0 - 2.0 * ta2_xx_z_yzz_1[i] * fe_0 + ta2_xx_z_yyzz_0[i] * pa_y[i] - ta2_xx_z_yyzz_1[i] * pc_y[i];

        ta2_xx_yz_yzzz_0[i] = ta2_xx_z_zzz_0[i] * fe_0 - ta2_xx_z_zzz_1[i] * fe_0 + ta2_xx_z_yzzz_0[i] * pa_y[i] - ta2_xx_z_yzzz_1[i] * pc_y[i];

        ta2_xx_yz_zzzz_0[i] = ta2_xx_z_zzzz_0[i] * pa_y[i] - ta2_xx_z_zzzz_1[i] * pc_y[i];
    }

    // Set up 75-90 components of targeted buffer : DG

    auto ta2_xx_zz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 75);

    auto ta2_xx_zz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 76);

    auto ta2_xx_zz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 77);

    auto ta2_xx_zz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 78);

    auto ta2_xx_zz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 79);

    auto ta2_xx_zz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 80);

    auto ta2_xx_zz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 81);

    auto ta2_xx_zz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 82);

    auto ta2_xx_zz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 83);

    auto ta2_xx_zz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 84);

    auto ta2_xx_zz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 85);

    auto ta2_xx_zz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 86);

    auto ta2_xx_zz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 87);

    auto ta2_xx_zz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 88);

    auto ta2_xx_zz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 89);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta2_xx_0_xxxx_0,  \
                             ta2_xx_0_xxxx_1,  \
                             ta2_xx_0_xxxy_0,  \
                             ta2_xx_0_xxxy_1,  \
                             ta2_xx_0_xxxz_0,  \
                             ta2_xx_0_xxxz_1,  \
                             ta2_xx_0_xxyy_0,  \
                             ta2_xx_0_xxyy_1,  \
                             ta2_xx_0_xxyz_0,  \
                             ta2_xx_0_xxyz_1,  \
                             ta2_xx_0_xxzz_0,  \
                             ta2_xx_0_xxzz_1,  \
                             ta2_xx_0_xyyy_0,  \
                             ta2_xx_0_xyyy_1,  \
                             ta2_xx_0_xyyz_0,  \
                             ta2_xx_0_xyyz_1,  \
                             ta2_xx_0_xyzz_0,  \
                             ta2_xx_0_xyzz_1,  \
                             ta2_xx_0_xzzz_0,  \
                             ta2_xx_0_xzzz_1,  \
                             ta2_xx_0_yyyy_0,  \
                             ta2_xx_0_yyyy_1,  \
                             ta2_xx_0_yyyz_0,  \
                             ta2_xx_0_yyyz_1,  \
                             ta2_xx_0_yyzz_0,  \
                             ta2_xx_0_yyzz_1,  \
                             ta2_xx_0_yzzz_0,  \
                             ta2_xx_0_yzzz_1,  \
                             ta2_xx_0_zzzz_0,  \
                             ta2_xx_0_zzzz_1,  \
                             ta2_xx_z_xxx_0,   \
                             ta2_xx_z_xxx_1,   \
                             ta2_xx_z_xxxx_0,  \
                             ta2_xx_z_xxxx_1,  \
                             ta2_xx_z_xxxy_0,  \
                             ta2_xx_z_xxxy_1,  \
                             ta2_xx_z_xxxz_0,  \
                             ta2_xx_z_xxxz_1,  \
                             ta2_xx_z_xxy_0,   \
                             ta2_xx_z_xxy_1,   \
                             ta2_xx_z_xxyy_0,  \
                             ta2_xx_z_xxyy_1,  \
                             ta2_xx_z_xxyz_0,  \
                             ta2_xx_z_xxyz_1,  \
                             ta2_xx_z_xxz_0,   \
                             ta2_xx_z_xxz_1,   \
                             ta2_xx_z_xxzz_0,  \
                             ta2_xx_z_xxzz_1,  \
                             ta2_xx_z_xyy_0,   \
                             ta2_xx_z_xyy_1,   \
                             ta2_xx_z_xyyy_0,  \
                             ta2_xx_z_xyyy_1,  \
                             ta2_xx_z_xyyz_0,  \
                             ta2_xx_z_xyyz_1,  \
                             ta2_xx_z_xyz_0,   \
                             ta2_xx_z_xyz_1,   \
                             ta2_xx_z_xyzz_0,  \
                             ta2_xx_z_xyzz_1,  \
                             ta2_xx_z_xzz_0,   \
                             ta2_xx_z_xzz_1,   \
                             ta2_xx_z_xzzz_0,  \
                             ta2_xx_z_xzzz_1,  \
                             ta2_xx_z_yyy_0,   \
                             ta2_xx_z_yyy_1,   \
                             ta2_xx_z_yyyy_0,  \
                             ta2_xx_z_yyyy_1,  \
                             ta2_xx_z_yyyz_0,  \
                             ta2_xx_z_yyyz_1,  \
                             ta2_xx_z_yyz_0,   \
                             ta2_xx_z_yyz_1,   \
                             ta2_xx_z_yyzz_0,  \
                             ta2_xx_z_yyzz_1,  \
                             ta2_xx_z_yzz_0,   \
                             ta2_xx_z_yzz_1,   \
                             ta2_xx_z_yzzz_0,  \
                             ta2_xx_z_yzzz_1,  \
                             ta2_xx_z_zzz_0,   \
                             ta2_xx_z_zzz_1,   \
                             ta2_xx_z_zzzz_0,  \
                             ta2_xx_z_zzzz_1,  \
                             ta2_xx_zz_xxxx_0, \
                             ta2_xx_zz_xxxy_0, \
                             ta2_xx_zz_xxxz_0, \
                             ta2_xx_zz_xxyy_0, \
                             ta2_xx_zz_xxyz_0, \
                             ta2_xx_zz_xxzz_0, \
                             ta2_xx_zz_xyyy_0, \
                             ta2_xx_zz_xyyz_0, \
                             ta2_xx_zz_xyzz_0, \
                             ta2_xx_zz_xzzz_0, \
                             ta2_xx_zz_yyyy_0, \
                             ta2_xx_zz_yyyz_0, \
                             ta2_xx_zz_yyzz_0, \
                             ta2_xx_zz_yzzz_0, \
                             ta2_xx_zz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_zz_xxxx_0[i] = ta2_xx_0_xxxx_0[i] * fe_0 - ta2_xx_0_xxxx_1[i] * fe_0 + ta2_xx_z_xxxx_0[i] * pa_z[i] - ta2_xx_z_xxxx_1[i] * pc_z[i];

        ta2_xx_zz_xxxy_0[i] = ta2_xx_0_xxxy_0[i] * fe_0 - ta2_xx_0_xxxy_1[i] * fe_0 + ta2_xx_z_xxxy_0[i] * pa_z[i] - ta2_xx_z_xxxy_1[i] * pc_z[i];

        ta2_xx_zz_xxxz_0[i] = ta2_xx_0_xxxz_0[i] * fe_0 - ta2_xx_0_xxxz_1[i] * fe_0 + ta2_xx_z_xxx_0[i] * fe_0 - ta2_xx_z_xxx_1[i] * fe_0 +
                              ta2_xx_z_xxxz_0[i] * pa_z[i] - ta2_xx_z_xxxz_1[i] * pc_z[i];

        ta2_xx_zz_xxyy_0[i] = ta2_xx_0_xxyy_0[i] * fe_0 - ta2_xx_0_xxyy_1[i] * fe_0 + ta2_xx_z_xxyy_0[i] * pa_z[i] - ta2_xx_z_xxyy_1[i] * pc_z[i];

        ta2_xx_zz_xxyz_0[i] = ta2_xx_0_xxyz_0[i] * fe_0 - ta2_xx_0_xxyz_1[i] * fe_0 + ta2_xx_z_xxy_0[i] * fe_0 - ta2_xx_z_xxy_1[i] * fe_0 +
                              ta2_xx_z_xxyz_0[i] * pa_z[i] - ta2_xx_z_xxyz_1[i] * pc_z[i];

        ta2_xx_zz_xxzz_0[i] = ta2_xx_0_xxzz_0[i] * fe_0 - ta2_xx_0_xxzz_1[i] * fe_0 + 2.0 * ta2_xx_z_xxz_0[i] * fe_0 -
                              2.0 * ta2_xx_z_xxz_1[i] * fe_0 + ta2_xx_z_xxzz_0[i] * pa_z[i] - ta2_xx_z_xxzz_1[i] * pc_z[i];

        ta2_xx_zz_xyyy_0[i] = ta2_xx_0_xyyy_0[i] * fe_0 - ta2_xx_0_xyyy_1[i] * fe_0 + ta2_xx_z_xyyy_0[i] * pa_z[i] - ta2_xx_z_xyyy_1[i] * pc_z[i];

        ta2_xx_zz_xyyz_0[i] = ta2_xx_0_xyyz_0[i] * fe_0 - ta2_xx_0_xyyz_1[i] * fe_0 + ta2_xx_z_xyy_0[i] * fe_0 - ta2_xx_z_xyy_1[i] * fe_0 +
                              ta2_xx_z_xyyz_0[i] * pa_z[i] - ta2_xx_z_xyyz_1[i] * pc_z[i];

        ta2_xx_zz_xyzz_0[i] = ta2_xx_0_xyzz_0[i] * fe_0 - ta2_xx_0_xyzz_1[i] * fe_0 + 2.0 * ta2_xx_z_xyz_0[i] * fe_0 -
                              2.0 * ta2_xx_z_xyz_1[i] * fe_0 + ta2_xx_z_xyzz_0[i] * pa_z[i] - ta2_xx_z_xyzz_1[i] * pc_z[i];

        ta2_xx_zz_xzzz_0[i] = ta2_xx_0_xzzz_0[i] * fe_0 - ta2_xx_0_xzzz_1[i] * fe_0 + 3.0 * ta2_xx_z_xzz_0[i] * fe_0 -
                              3.0 * ta2_xx_z_xzz_1[i] * fe_0 + ta2_xx_z_xzzz_0[i] * pa_z[i] - ta2_xx_z_xzzz_1[i] * pc_z[i];

        ta2_xx_zz_yyyy_0[i] = ta2_xx_0_yyyy_0[i] * fe_0 - ta2_xx_0_yyyy_1[i] * fe_0 + ta2_xx_z_yyyy_0[i] * pa_z[i] - ta2_xx_z_yyyy_1[i] * pc_z[i];

        ta2_xx_zz_yyyz_0[i] = ta2_xx_0_yyyz_0[i] * fe_0 - ta2_xx_0_yyyz_1[i] * fe_0 + ta2_xx_z_yyy_0[i] * fe_0 - ta2_xx_z_yyy_1[i] * fe_0 +
                              ta2_xx_z_yyyz_0[i] * pa_z[i] - ta2_xx_z_yyyz_1[i] * pc_z[i];

        ta2_xx_zz_yyzz_0[i] = ta2_xx_0_yyzz_0[i] * fe_0 - ta2_xx_0_yyzz_1[i] * fe_0 + 2.0 * ta2_xx_z_yyz_0[i] * fe_0 -
                              2.0 * ta2_xx_z_yyz_1[i] * fe_0 + ta2_xx_z_yyzz_0[i] * pa_z[i] - ta2_xx_z_yyzz_1[i] * pc_z[i];

        ta2_xx_zz_yzzz_0[i] = ta2_xx_0_yzzz_0[i] * fe_0 - ta2_xx_0_yzzz_1[i] * fe_0 + 3.0 * ta2_xx_z_yzz_0[i] * fe_0 -
                              3.0 * ta2_xx_z_yzz_1[i] * fe_0 + ta2_xx_z_yzzz_0[i] * pa_z[i] - ta2_xx_z_yzzz_1[i] * pc_z[i];

        ta2_xx_zz_zzzz_0[i] = ta2_xx_0_zzzz_0[i] * fe_0 - ta2_xx_0_zzzz_1[i] * fe_0 + 4.0 * ta2_xx_z_zzz_0[i] * fe_0 -
                              4.0 * ta2_xx_z_zzz_1[i] * fe_0 + ta2_xx_z_zzzz_0[i] * pa_z[i] - ta2_xx_z_zzzz_1[i] * pc_z[i];
    }

    // Set up 90-105 components of targeted buffer : DG

    auto ta2_xy_xx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 90);

    auto ta2_xy_xx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 91);

    auto ta2_xy_xx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 92);

    auto ta2_xy_xx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 93);

    auto ta2_xy_xx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 94);

    auto ta2_xy_xx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 95);

    auto ta2_xy_xx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 96);

    auto ta2_xy_xx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 97);

    auto ta2_xy_xx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 98);

    auto ta2_xy_xx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 99);

    auto ta2_xy_xx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 100);

    auto ta2_xy_xx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 101);

    auto ta2_xy_xx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 102);

    auto ta2_xy_xx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 103);

    auto ta2_xy_xx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 104);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_y_x_xxxx_1,   \
                             ta1_y_x_xxxy_1,   \
                             ta1_y_x_xxxz_1,   \
                             ta1_y_x_xxyy_1,   \
                             ta1_y_x_xxyz_1,   \
                             ta1_y_x_xxzz_1,   \
                             ta1_y_x_xyyy_1,   \
                             ta1_y_x_xyyz_1,   \
                             ta1_y_x_xyzz_1,   \
                             ta1_y_x_xzzz_1,   \
                             ta1_y_x_yyyy_1,   \
                             ta1_y_x_yyyz_1,   \
                             ta1_y_x_yyzz_1,   \
                             ta1_y_x_yzzz_1,   \
                             ta1_y_x_zzzz_1,   \
                             ta2_xy_0_xxxx_0,  \
                             ta2_xy_0_xxxx_1,  \
                             ta2_xy_0_xxxy_0,  \
                             ta2_xy_0_xxxy_1,  \
                             ta2_xy_0_xxxz_0,  \
                             ta2_xy_0_xxxz_1,  \
                             ta2_xy_0_xxyy_0,  \
                             ta2_xy_0_xxyy_1,  \
                             ta2_xy_0_xxyz_0,  \
                             ta2_xy_0_xxyz_1,  \
                             ta2_xy_0_xxzz_0,  \
                             ta2_xy_0_xxzz_1,  \
                             ta2_xy_0_xyyy_0,  \
                             ta2_xy_0_xyyy_1,  \
                             ta2_xy_0_xyyz_0,  \
                             ta2_xy_0_xyyz_1,  \
                             ta2_xy_0_xyzz_0,  \
                             ta2_xy_0_xyzz_1,  \
                             ta2_xy_0_xzzz_0,  \
                             ta2_xy_0_xzzz_1,  \
                             ta2_xy_0_yyyy_0,  \
                             ta2_xy_0_yyyy_1,  \
                             ta2_xy_0_yyyz_0,  \
                             ta2_xy_0_yyyz_1,  \
                             ta2_xy_0_yyzz_0,  \
                             ta2_xy_0_yyzz_1,  \
                             ta2_xy_0_yzzz_0,  \
                             ta2_xy_0_yzzz_1,  \
                             ta2_xy_0_zzzz_0,  \
                             ta2_xy_0_zzzz_1,  \
                             ta2_xy_x_xxx_0,   \
                             ta2_xy_x_xxx_1,   \
                             ta2_xy_x_xxxx_0,  \
                             ta2_xy_x_xxxx_1,  \
                             ta2_xy_x_xxxy_0,  \
                             ta2_xy_x_xxxy_1,  \
                             ta2_xy_x_xxxz_0,  \
                             ta2_xy_x_xxxz_1,  \
                             ta2_xy_x_xxy_0,   \
                             ta2_xy_x_xxy_1,   \
                             ta2_xy_x_xxyy_0,  \
                             ta2_xy_x_xxyy_1,  \
                             ta2_xy_x_xxyz_0,  \
                             ta2_xy_x_xxyz_1,  \
                             ta2_xy_x_xxz_0,   \
                             ta2_xy_x_xxz_1,   \
                             ta2_xy_x_xxzz_0,  \
                             ta2_xy_x_xxzz_1,  \
                             ta2_xy_x_xyy_0,   \
                             ta2_xy_x_xyy_1,   \
                             ta2_xy_x_xyyy_0,  \
                             ta2_xy_x_xyyy_1,  \
                             ta2_xy_x_xyyz_0,  \
                             ta2_xy_x_xyyz_1,  \
                             ta2_xy_x_xyz_0,   \
                             ta2_xy_x_xyz_1,   \
                             ta2_xy_x_xyzz_0,  \
                             ta2_xy_x_xyzz_1,  \
                             ta2_xy_x_xzz_0,   \
                             ta2_xy_x_xzz_1,   \
                             ta2_xy_x_xzzz_0,  \
                             ta2_xy_x_xzzz_1,  \
                             ta2_xy_x_yyy_0,   \
                             ta2_xy_x_yyy_1,   \
                             ta2_xy_x_yyyy_0,  \
                             ta2_xy_x_yyyy_1,  \
                             ta2_xy_x_yyyz_0,  \
                             ta2_xy_x_yyyz_1,  \
                             ta2_xy_x_yyz_0,   \
                             ta2_xy_x_yyz_1,   \
                             ta2_xy_x_yyzz_0,  \
                             ta2_xy_x_yyzz_1,  \
                             ta2_xy_x_yzz_0,   \
                             ta2_xy_x_yzz_1,   \
                             ta2_xy_x_yzzz_0,  \
                             ta2_xy_x_yzzz_1,  \
                             ta2_xy_x_zzz_0,   \
                             ta2_xy_x_zzz_1,   \
                             ta2_xy_x_zzzz_0,  \
                             ta2_xy_x_zzzz_1,  \
                             ta2_xy_xx_xxxx_0, \
                             ta2_xy_xx_xxxy_0, \
                             ta2_xy_xx_xxxz_0, \
                             ta2_xy_xx_xxyy_0, \
                             ta2_xy_xx_xxyz_0, \
                             ta2_xy_xx_xxzz_0, \
                             ta2_xy_xx_xyyy_0, \
                             ta2_xy_xx_xyyz_0, \
                             ta2_xy_xx_xyzz_0, \
                             ta2_xy_xx_xzzz_0, \
                             ta2_xy_xx_yyyy_0, \
                             ta2_xy_xx_yyyz_0, \
                             ta2_xy_xx_yyzz_0, \
                             ta2_xy_xx_yzzz_0, \
                             ta2_xy_xx_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xx_xxxx_0[i] = ta2_xy_0_xxxx_0[i] * fe_0 - ta2_xy_0_xxxx_1[i] * fe_0 + 4.0 * ta2_xy_x_xxx_0[i] * fe_0 -
                              4.0 * ta2_xy_x_xxx_1[i] * fe_0 + ta1_y_x_xxxx_1[i] + ta2_xy_x_xxxx_0[i] * pa_x[i] - ta2_xy_x_xxxx_1[i] * pc_x[i];

        ta2_xy_xx_xxxy_0[i] = ta2_xy_0_xxxy_0[i] * fe_0 - ta2_xy_0_xxxy_1[i] * fe_0 + 3.0 * ta2_xy_x_xxy_0[i] * fe_0 -
                              3.0 * ta2_xy_x_xxy_1[i] * fe_0 + ta1_y_x_xxxy_1[i] + ta2_xy_x_xxxy_0[i] * pa_x[i] - ta2_xy_x_xxxy_1[i] * pc_x[i];

        ta2_xy_xx_xxxz_0[i] = ta2_xy_0_xxxz_0[i] * fe_0 - ta2_xy_0_xxxz_1[i] * fe_0 + 3.0 * ta2_xy_x_xxz_0[i] * fe_0 -
                              3.0 * ta2_xy_x_xxz_1[i] * fe_0 + ta1_y_x_xxxz_1[i] + ta2_xy_x_xxxz_0[i] * pa_x[i] - ta2_xy_x_xxxz_1[i] * pc_x[i];

        ta2_xy_xx_xxyy_0[i] = ta2_xy_0_xxyy_0[i] * fe_0 - ta2_xy_0_xxyy_1[i] * fe_0 + 2.0 * ta2_xy_x_xyy_0[i] * fe_0 -
                              2.0 * ta2_xy_x_xyy_1[i] * fe_0 + ta1_y_x_xxyy_1[i] + ta2_xy_x_xxyy_0[i] * pa_x[i] - ta2_xy_x_xxyy_1[i] * pc_x[i];

        ta2_xy_xx_xxyz_0[i] = ta2_xy_0_xxyz_0[i] * fe_0 - ta2_xy_0_xxyz_1[i] * fe_0 + 2.0 * ta2_xy_x_xyz_0[i] * fe_0 -
                              2.0 * ta2_xy_x_xyz_1[i] * fe_0 + ta1_y_x_xxyz_1[i] + ta2_xy_x_xxyz_0[i] * pa_x[i] - ta2_xy_x_xxyz_1[i] * pc_x[i];

        ta2_xy_xx_xxzz_0[i] = ta2_xy_0_xxzz_0[i] * fe_0 - ta2_xy_0_xxzz_1[i] * fe_0 + 2.0 * ta2_xy_x_xzz_0[i] * fe_0 -
                              2.0 * ta2_xy_x_xzz_1[i] * fe_0 + ta1_y_x_xxzz_1[i] + ta2_xy_x_xxzz_0[i] * pa_x[i] - ta2_xy_x_xxzz_1[i] * pc_x[i];

        ta2_xy_xx_xyyy_0[i] = ta2_xy_0_xyyy_0[i] * fe_0 - ta2_xy_0_xyyy_1[i] * fe_0 + ta2_xy_x_yyy_0[i] * fe_0 - ta2_xy_x_yyy_1[i] * fe_0 +
                              ta1_y_x_xyyy_1[i] + ta2_xy_x_xyyy_0[i] * pa_x[i] - ta2_xy_x_xyyy_1[i] * pc_x[i];

        ta2_xy_xx_xyyz_0[i] = ta2_xy_0_xyyz_0[i] * fe_0 - ta2_xy_0_xyyz_1[i] * fe_0 + ta2_xy_x_yyz_0[i] * fe_0 - ta2_xy_x_yyz_1[i] * fe_0 +
                              ta1_y_x_xyyz_1[i] + ta2_xy_x_xyyz_0[i] * pa_x[i] - ta2_xy_x_xyyz_1[i] * pc_x[i];

        ta2_xy_xx_xyzz_0[i] = ta2_xy_0_xyzz_0[i] * fe_0 - ta2_xy_0_xyzz_1[i] * fe_0 + ta2_xy_x_yzz_0[i] * fe_0 - ta2_xy_x_yzz_1[i] * fe_0 +
                              ta1_y_x_xyzz_1[i] + ta2_xy_x_xyzz_0[i] * pa_x[i] - ta2_xy_x_xyzz_1[i] * pc_x[i];

        ta2_xy_xx_xzzz_0[i] = ta2_xy_0_xzzz_0[i] * fe_0 - ta2_xy_0_xzzz_1[i] * fe_0 + ta2_xy_x_zzz_0[i] * fe_0 - ta2_xy_x_zzz_1[i] * fe_0 +
                              ta1_y_x_xzzz_1[i] + ta2_xy_x_xzzz_0[i] * pa_x[i] - ta2_xy_x_xzzz_1[i] * pc_x[i];

        ta2_xy_xx_yyyy_0[i] =
            ta2_xy_0_yyyy_0[i] * fe_0 - ta2_xy_0_yyyy_1[i] * fe_0 + ta1_y_x_yyyy_1[i] + ta2_xy_x_yyyy_0[i] * pa_x[i] - ta2_xy_x_yyyy_1[i] * pc_x[i];

        ta2_xy_xx_yyyz_0[i] =
            ta2_xy_0_yyyz_0[i] * fe_0 - ta2_xy_0_yyyz_1[i] * fe_0 + ta1_y_x_yyyz_1[i] + ta2_xy_x_yyyz_0[i] * pa_x[i] - ta2_xy_x_yyyz_1[i] * pc_x[i];

        ta2_xy_xx_yyzz_0[i] =
            ta2_xy_0_yyzz_0[i] * fe_0 - ta2_xy_0_yyzz_1[i] * fe_0 + ta1_y_x_yyzz_1[i] + ta2_xy_x_yyzz_0[i] * pa_x[i] - ta2_xy_x_yyzz_1[i] * pc_x[i];

        ta2_xy_xx_yzzz_0[i] =
            ta2_xy_0_yzzz_0[i] * fe_0 - ta2_xy_0_yzzz_1[i] * fe_0 + ta1_y_x_yzzz_1[i] + ta2_xy_x_yzzz_0[i] * pa_x[i] - ta2_xy_x_yzzz_1[i] * pc_x[i];

        ta2_xy_xx_zzzz_0[i] =
            ta2_xy_0_zzzz_0[i] * fe_0 - ta2_xy_0_zzzz_1[i] * fe_0 + ta1_y_x_zzzz_1[i] + ta2_xy_x_zzzz_0[i] * pa_x[i] - ta2_xy_x_zzzz_1[i] * pc_x[i];
    }

    // Set up 105-120 components of targeted buffer : DG

    auto ta2_xy_xy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 105);

    auto ta2_xy_xy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 106);

    auto ta2_xy_xy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 107);

    auto ta2_xy_xy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 108);

    auto ta2_xy_xy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 109);

    auto ta2_xy_xy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 110);

    auto ta2_xy_xy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 111);

    auto ta2_xy_xy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 112);

    auto ta2_xy_xy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 113);

    auto ta2_xy_xy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 114);

    auto ta2_xy_xy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 115);

    auto ta2_xy_xy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 116);

    auto ta2_xy_xy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 117);

    auto ta2_xy_xy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 118);

    auto ta2_xy_xy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 119);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_x_x_xxxx_1,   \
                             ta1_x_x_xxxz_1,   \
                             ta1_x_x_xxzz_1,   \
                             ta1_x_x_xzzz_1,   \
                             ta1_y_y_xxxy_1,   \
                             ta1_y_y_xxyy_1,   \
                             ta1_y_y_xxyz_1,   \
                             ta1_y_y_xyyy_1,   \
                             ta1_y_y_xyyz_1,   \
                             ta1_y_y_xyzz_1,   \
                             ta1_y_y_yyyy_1,   \
                             ta1_y_y_yyyz_1,   \
                             ta1_y_y_yyzz_1,   \
                             ta1_y_y_yzzz_1,   \
                             ta1_y_y_zzzz_1,   \
                             ta2_xy_x_xxxx_0,  \
                             ta2_xy_x_xxxx_1,  \
                             ta2_xy_x_xxxz_0,  \
                             ta2_xy_x_xxxz_1,  \
                             ta2_xy_x_xxzz_0,  \
                             ta2_xy_x_xxzz_1,  \
                             ta2_xy_x_xzzz_0,  \
                             ta2_xy_x_xzzz_1,  \
                             ta2_xy_xy_xxxx_0, \
                             ta2_xy_xy_xxxy_0, \
                             ta2_xy_xy_xxxz_0, \
                             ta2_xy_xy_xxyy_0, \
                             ta2_xy_xy_xxyz_0, \
                             ta2_xy_xy_xxzz_0, \
                             ta2_xy_xy_xyyy_0, \
                             ta2_xy_xy_xyyz_0, \
                             ta2_xy_xy_xyzz_0, \
                             ta2_xy_xy_xzzz_0, \
                             ta2_xy_xy_yyyy_0, \
                             ta2_xy_xy_yyyz_0, \
                             ta2_xy_xy_yyzz_0, \
                             ta2_xy_xy_yzzz_0, \
                             ta2_xy_xy_zzzz_0, \
                             ta2_xy_y_xxxy_0,  \
                             ta2_xy_y_xxxy_1,  \
                             ta2_xy_y_xxy_0,   \
                             ta2_xy_y_xxy_1,   \
                             ta2_xy_y_xxyy_0,  \
                             ta2_xy_y_xxyy_1,  \
                             ta2_xy_y_xxyz_0,  \
                             ta2_xy_y_xxyz_1,  \
                             ta2_xy_y_xyy_0,   \
                             ta2_xy_y_xyy_1,   \
                             ta2_xy_y_xyyy_0,  \
                             ta2_xy_y_xyyy_1,  \
                             ta2_xy_y_xyyz_0,  \
                             ta2_xy_y_xyyz_1,  \
                             ta2_xy_y_xyz_0,   \
                             ta2_xy_y_xyz_1,   \
                             ta2_xy_y_xyzz_0,  \
                             ta2_xy_y_xyzz_1,  \
                             ta2_xy_y_yyy_0,   \
                             ta2_xy_y_yyy_1,   \
                             ta2_xy_y_yyyy_0,  \
                             ta2_xy_y_yyyy_1,  \
                             ta2_xy_y_yyyz_0,  \
                             ta2_xy_y_yyyz_1,  \
                             ta2_xy_y_yyz_0,   \
                             ta2_xy_y_yyz_1,   \
                             ta2_xy_y_yyzz_0,  \
                             ta2_xy_y_yyzz_1,  \
                             ta2_xy_y_yzz_0,   \
                             ta2_xy_y_yzz_1,   \
                             ta2_xy_y_yzzz_0,  \
                             ta2_xy_y_yzzz_1,  \
                             ta2_xy_y_zzzz_0,  \
                             ta2_xy_y_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xy_xxxx_0[i] = ta1_x_x_xxxx_1[i] + ta2_xy_x_xxxx_0[i] * pa_y[i] - ta2_xy_x_xxxx_1[i] * pc_y[i];

        ta2_xy_xy_xxxy_0[i] = 3.0 * ta2_xy_y_xxy_0[i] * fe_0 - 3.0 * ta2_xy_y_xxy_1[i] * fe_0 + ta1_y_y_xxxy_1[i] + ta2_xy_y_xxxy_0[i] * pa_x[i] -
                              ta2_xy_y_xxxy_1[i] * pc_x[i];

        ta2_xy_xy_xxxz_0[i] = ta1_x_x_xxxz_1[i] + ta2_xy_x_xxxz_0[i] * pa_y[i] - ta2_xy_x_xxxz_1[i] * pc_y[i];

        ta2_xy_xy_xxyy_0[i] = 2.0 * ta2_xy_y_xyy_0[i] * fe_0 - 2.0 * ta2_xy_y_xyy_1[i] * fe_0 + ta1_y_y_xxyy_1[i] + ta2_xy_y_xxyy_0[i] * pa_x[i] -
                              ta2_xy_y_xxyy_1[i] * pc_x[i];

        ta2_xy_xy_xxyz_0[i] = 2.0 * ta2_xy_y_xyz_0[i] * fe_0 - 2.0 * ta2_xy_y_xyz_1[i] * fe_0 + ta1_y_y_xxyz_1[i] + ta2_xy_y_xxyz_0[i] * pa_x[i] -
                              ta2_xy_y_xxyz_1[i] * pc_x[i];

        ta2_xy_xy_xxzz_0[i] = ta1_x_x_xxzz_1[i] + ta2_xy_x_xxzz_0[i] * pa_y[i] - ta2_xy_x_xxzz_1[i] * pc_y[i];

        ta2_xy_xy_xyyy_0[i] =
            ta2_xy_y_yyy_0[i] * fe_0 - ta2_xy_y_yyy_1[i] * fe_0 + ta1_y_y_xyyy_1[i] + ta2_xy_y_xyyy_0[i] * pa_x[i] - ta2_xy_y_xyyy_1[i] * pc_x[i];

        ta2_xy_xy_xyyz_0[i] =
            ta2_xy_y_yyz_0[i] * fe_0 - ta2_xy_y_yyz_1[i] * fe_0 + ta1_y_y_xyyz_1[i] + ta2_xy_y_xyyz_0[i] * pa_x[i] - ta2_xy_y_xyyz_1[i] * pc_x[i];

        ta2_xy_xy_xyzz_0[i] =
            ta2_xy_y_yzz_0[i] * fe_0 - ta2_xy_y_yzz_1[i] * fe_0 + ta1_y_y_xyzz_1[i] + ta2_xy_y_xyzz_0[i] * pa_x[i] - ta2_xy_y_xyzz_1[i] * pc_x[i];

        ta2_xy_xy_xzzz_0[i] = ta1_x_x_xzzz_1[i] + ta2_xy_x_xzzz_0[i] * pa_y[i] - ta2_xy_x_xzzz_1[i] * pc_y[i];

        ta2_xy_xy_yyyy_0[i] = ta1_y_y_yyyy_1[i] + ta2_xy_y_yyyy_0[i] * pa_x[i] - ta2_xy_y_yyyy_1[i] * pc_x[i];

        ta2_xy_xy_yyyz_0[i] = ta1_y_y_yyyz_1[i] + ta2_xy_y_yyyz_0[i] * pa_x[i] - ta2_xy_y_yyyz_1[i] * pc_x[i];

        ta2_xy_xy_yyzz_0[i] = ta1_y_y_yyzz_1[i] + ta2_xy_y_yyzz_0[i] * pa_x[i] - ta2_xy_y_yyzz_1[i] * pc_x[i];

        ta2_xy_xy_yzzz_0[i] = ta1_y_y_yzzz_1[i] + ta2_xy_y_yzzz_0[i] * pa_x[i] - ta2_xy_y_yzzz_1[i] * pc_x[i];

        ta2_xy_xy_zzzz_0[i] = ta1_y_y_zzzz_1[i] + ta2_xy_y_zzzz_0[i] * pa_x[i] - ta2_xy_y_zzzz_1[i] * pc_x[i];
    }

    // Set up 120-135 components of targeted buffer : DG

    auto ta2_xy_xz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 120);

    auto ta2_xy_xz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 121);

    auto ta2_xy_xz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 122);

    auto ta2_xy_xz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 123);

    auto ta2_xy_xz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 124);

    auto ta2_xy_xz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 125);

    auto ta2_xy_xz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 126);

    auto ta2_xy_xz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 127);

    auto ta2_xy_xz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 128);

    auto ta2_xy_xz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 129);

    auto ta2_xy_xz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 130);

    auto ta2_xy_xz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 131);

    auto ta2_xy_xz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 132);

    auto ta2_xy_xz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 133);

    auto ta2_xy_xz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 134);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_y_z_yyyz_1,   \
                             ta1_y_z_yyzz_1,   \
                             ta1_y_z_yzzz_1,   \
                             ta1_y_z_zzzz_1,   \
                             ta2_xy_x_xxx_0,   \
                             ta2_xy_x_xxx_1,   \
                             ta2_xy_x_xxxx_0,  \
                             ta2_xy_x_xxxx_1,  \
                             ta2_xy_x_xxxy_0,  \
                             ta2_xy_x_xxxy_1,  \
                             ta2_xy_x_xxxz_0,  \
                             ta2_xy_x_xxxz_1,  \
                             ta2_xy_x_xxy_0,   \
                             ta2_xy_x_xxy_1,   \
                             ta2_xy_x_xxyy_0,  \
                             ta2_xy_x_xxyy_1,  \
                             ta2_xy_x_xxyz_0,  \
                             ta2_xy_x_xxyz_1,  \
                             ta2_xy_x_xxz_0,   \
                             ta2_xy_x_xxz_1,   \
                             ta2_xy_x_xxzz_0,  \
                             ta2_xy_x_xxzz_1,  \
                             ta2_xy_x_xyy_0,   \
                             ta2_xy_x_xyy_1,   \
                             ta2_xy_x_xyyy_0,  \
                             ta2_xy_x_xyyy_1,  \
                             ta2_xy_x_xyyz_0,  \
                             ta2_xy_x_xyyz_1,  \
                             ta2_xy_x_xyz_0,   \
                             ta2_xy_x_xyz_1,   \
                             ta2_xy_x_xyzz_0,  \
                             ta2_xy_x_xyzz_1,  \
                             ta2_xy_x_xzz_0,   \
                             ta2_xy_x_xzz_1,   \
                             ta2_xy_x_xzzz_0,  \
                             ta2_xy_x_xzzz_1,  \
                             ta2_xy_x_yyyy_0,  \
                             ta2_xy_x_yyyy_1,  \
                             ta2_xy_xz_xxxx_0, \
                             ta2_xy_xz_xxxy_0, \
                             ta2_xy_xz_xxxz_0, \
                             ta2_xy_xz_xxyy_0, \
                             ta2_xy_xz_xxyz_0, \
                             ta2_xy_xz_xxzz_0, \
                             ta2_xy_xz_xyyy_0, \
                             ta2_xy_xz_xyyz_0, \
                             ta2_xy_xz_xyzz_0, \
                             ta2_xy_xz_xzzz_0, \
                             ta2_xy_xz_yyyy_0, \
                             ta2_xy_xz_yyyz_0, \
                             ta2_xy_xz_yyzz_0, \
                             ta2_xy_xz_yzzz_0, \
                             ta2_xy_xz_zzzz_0, \
                             ta2_xy_z_yyyz_0,  \
                             ta2_xy_z_yyyz_1,  \
                             ta2_xy_z_yyzz_0,  \
                             ta2_xy_z_yyzz_1,  \
                             ta2_xy_z_yzzz_0,  \
                             ta2_xy_z_yzzz_1,  \
                             ta2_xy_z_zzzz_0,  \
                             ta2_xy_z_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xz_xxxx_0[i] = ta2_xy_x_xxxx_0[i] * pa_z[i] - ta2_xy_x_xxxx_1[i] * pc_z[i];

        ta2_xy_xz_xxxy_0[i] = ta2_xy_x_xxxy_0[i] * pa_z[i] - ta2_xy_x_xxxy_1[i] * pc_z[i];

        ta2_xy_xz_xxxz_0[i] = ta2_xy_x_xxx_0[i] * fe_0 - ta2_xy_x_xxx_1[i] * fe_0 + ta2_xy_x_xxxz_0[i] * pa_z[i] - ta2_xy_x_xxxz_1[i] * pc_z[i];

        ta2_xy_xz_xxyy_0[i] = ta2_xy_x_xxyy_0[i] * pa_z[i] - ta2_xy_x_xxyy_1[i] * pc_z[i];

        ta2_xy_xz_xxyz_0[i] = ta2_xy_x_xxy_0[i] * fe_0 - ta2_xy_x_xxy_1[i] * fe_0 + ta2_xy_x_xxyz_0[i] * pa_z[i] - ta2_xy_x_xxyz_1[i] * pc_z[i];

        ta2_xy_xz_xxzz_0[i] =
            2.0 * ta2_xy_x_xxz_0[i] * fe_0 - 2.0 * ta2_xy_x_xxz_1[i] * fe_0 + ta2_xy_x_xxzz_0[i] * pa_z[i] - ta2_xy_x_xxzz_1[i] * pc_z[i];

        ta2_xy_xz_xyyy_0[i] = ta2_xy_x_xyyy_0[i] * pa_z[i] - ta2_xy_x_xyyy_1[i] * pc_z[i];

        ta2_xy_xz_xyyz_0[i] = ta2_xy_x_xyy_0[i] * fe_0 - ta2_xy_x_xyy_1[i] * fe_0 + ta2_xy_x_xyyz_0[i] * pa_z[i] - ta2_xy_x_xyyz_1[i] * pc_z[i];

        ta2_xy_xz_xyzz_0[i] =
            2.0 * ta2_xy_x_xyz_0[i] * fe_0 - 2.0 * ta2_xy_x_xyz_1[i] * fe_0 + ta2_xy_x_xyzz_0[i] * pa_z[i] - ta2_xy_x_xyzz_1[i] * pc_z[i];

        ta2_xy_xz_xzzz_0[i] =
            3.0 * ta2_xy_x_xzz_0[i] * fe_0 - 3.0 * ta2_xy_x_xzz_1[i] * fe_0 + ta2_xy_x_xzzz_0[i] * pa_z[i] - ta2_xy_x_xzzz_1[i] * pc_z[i];

        ta2_xy_xz_yyyy_0[i] = ta2_xy_x_yyyy_0[i] * pa_z[i] - ta2_xy_x_yyyy_1[i] * pc_z[i];

        ta2_xy_xz_yyyz_0[i] = ta1_y_z_yyyz_1[i] + ta2_xy_z_yyyz_0[i] * pa_x[i] - ta2_xy_z_yyyz_1[i] * pc_x[i];

        ta2_xy_xz_yyzz_0[i] = ta1_y_z_yyzz_1[i] + ta2_xy_z_yyzz_0[i] * pa_x[i] - ta2_xy_z_yyzz_1[i] * pc_x[i];

        ta2_xy_xz_yzzz_0[i] = ta1_y_z_yzzz_1[i] + ta2_xy_z_yzzz_0[i] * pa_x[i] - ta2_xy_z_yzzz_1[i] * pc_x[i];

        ta2_xy_xz_zzzz_0[i] = ta1_y_z_zzzz_1[i] + ta2_xy_z_zzzz_0[i] * pa_x[i] - ta2_xy_z_zzzz_1[i] * pc_x[i];
    }

    // Set up 135-150 components of targeted buffer : DG

    auto ta2_xy_yy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 135);

    auto ta2_xy_yy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 136);

    auto ta2_xy_yy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 137);

    auto ta2_xy_yy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 138);

    auto ta2_xy_yy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 139);

    auto ta2_xy_yy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 140);

    auto ta2_xy_yy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 141);

    auto ta2_xy_yy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 142);

    auto ta2_xy_yy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 143);

    auto ta2_xy_yy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 144);

    auto ta2_xy_yy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 145);

    auto ta2_xy_yy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 146);

    auto ta2_xy_yy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 147);

    auto ta2_xy_yy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 148);

    auto ta2_xy_yy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 149);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_x_y_xxxx_1,   \
                             ta1_x_y_xxxy_1,   \
                             ta1_x_y_xxxz_1,   \
                             ta1_x_y_xxyy_1,   \
                             ta1_x_y_xxyz_1,   \
                             ta1_x_y_xxzz_1,   \
                             ta1_x_y_xyyy_1,   \
                             ta1_x_y_xyyz_1,   \
                             ta1_x_y_xyzz_1,   \
                             ta1_x_y_xzzz_1,   \
                             ta1_x_y_yyyy_1,   \
                             ta1_x_y_yyyz_1,   \
                             ta1_x_y_yyzz_1,   \
                             ta1_x_y_yzzz_1,   \
                             ta1_x_y_zzzz_1,   \
                             ta2_xy_0_xxxx_0,  \
                             ta2_xy_0_xxxx_1,  \
                             ta2_xy_0_xxxy_0,  \
                             ta2_xy_0_xxxy_1,  \
                             ta2_xy_0_xxxz_0,  \
                             ta2_xy_0_xxxz_1,  \
                             ta2_xy_0_xxyy_0,  \
                             ta2_xy_0_xxyy_1,  \
                             ta2_xy_0_xxyz_0,  \
                             ta2_xy_0_xxyz_1,  \
                             ta2_xy_0_xxzz_0,  \
                             ta2_xy_0_xxzz_1,  \
                             ta2_xy_0_xyyy_0,  \
                             ta2_xy_0_xyyy_1,  \
                             ta2_xy_0_xyyz_0,  \
                             ta2_xy_0_xyyz_1,  \
                             ta2_xy_0_xyzz_0,  \
                             ta2_xy_0_xyzz_1,  \
                             ta2_xy_0_xzzz_0,  \
                             ta2_xy_0_xzzz_1,  \
                             ta2_xy_0_yyyy_0,  \
                             ta2_xy_0_yyyy_1,  \
                             ta2_xy_0_yyyz_0,  \
                             ta2_xy_0_yyyz_1,  \
                             ta2_xy_0_yyzz_0,  \
                             ta2_xy_0_yyzz_1,  \
                             ta2_xy_0_yzzz_0,  \
                             ta2_xy_0_yzzz_1,  \
                             ta2_xy_0_zzzz_0,  \
                             ta2_xy_0_zzzz_1,  \
                             ta2_xy_y_xxx_0,   \
                             ta2_xy_y_xxx_1,   \
                             ta2_xy_y_xxxx_0,  \
                             ta2_xy_y_xxxx_1,  \
                             ta2_xy_y_xxxy_0,  \
                             ta2_xy_y_xxxy_1,  \
                             ta2_xy_y_xxxz_0,  \
                             ta2_xy_y_xxxz_1,  \
                             ta2_xy_y_xxy_0,   \
                             ta2_xy_y_xxy_1,   \
                             ta2_xy_y_xxyy_0,  \
                             ta2_xy_y_xxyy_1,  \
                             ta2_xy_y_xxyz_0,  \
                             ta2_xy_y_xxyz_1,  \
                             ta2_xy_y_xxz_0,   \
                             ta2_xy_y_xxz_1,   \
                             ta2_xy_y_xxzz_0,  \
                             ta2_xy_y_xxzz_1,  \
                             ta2_xy_y_xyy_0,   \
                             ta2_xy_y_xyy_1,   \
                             ta2_xy_y_xyyy_0,  \
                             ta2_xy_y_xyyy_1,  \
                             ta2_xy_y_xyyz_0,  \
                             ta2_xy_y_xyyz_1,  \
                             ta2_xy_y_xyz_0,   \
                             ta2_xy_y_xyz_1,   \
                             ta2_xy_y_xyzz_0,  \
                             ta2_xy_y_xyzz_1,  \
                             ta2_xy_y_xzz_0,   \
                             ta2_xy_y_xzz_1,   \
                             ta2_xy_y_xzzz_0,  \
                             ta2_xy_y_xzzz_1,  \
                             ta2_xy_y_yyy_0,   \
                             ta2_xy_y_yyy_1,   \
                             ta2_xy_y_yyyy_0,  \
                             ta2_xy_y_yyyy_1,  \
                             ta2_xy_y_yyyz_0,  \
                             ta2_xy_y_yyyz_1,  \
                             ta2_xy_y_yyz_0,   \
                             ta2_xy_y_yyz_1,   \
                             ta2_xy_y_yyzz_0,  \
                             ta2_xy_y_yyzz_1,  \
                             ta2_xy_y_yzz_0,   \
                             ta2_xy_y_yzz_1,   \
                             ta2_xy_y_yzzz_0,  \
                             ta2_xy_y_yzzz_1,  \
                             ta2_xy_y_zzz_0,   \
                             ta2_xy_y_zzz_1,   \
                             ta2_xy_y_zzzz_0,  \
                             ta2_xy_y_zzzz_1,  \
                             ta2_xy_yy_xxxx_0, \
                             ta2_xy_yy_xxxy_0, \
                             ta2_xy_yy_xxxz_0, \
                             ta2_xy_yy_xxyy_0, \
                             ta2_xy_yy_xxyz_0, \
                             ta2_xy_yy_xxzz_0, \
                             ta2_xy_yy_xyyy_0, \
                             ta2_xy_yy_xyyz_0, \
                             ta2_xy_yy_xyzz_0, \
                             ta2_xy_yy_xzzz_0, \
                             ta2_xy_yy_yyyy_0, \
                             ta2_xy_yy_yyyz_0, \
                             ta2_xy_yy_yyzz_0, \
                             ta2_xy_yy_yzzz_0, \
                             ta2_xy_yy_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yy_xxxx_0[i] =
            ta2_xy_0_xxxx_0[i] * fe_0 - ta2_xy_0_xxxx_1[i] * fe_0 + ta1_x_y_xxxx_1[i] + ta2_xy_y_xxxx_0[i] * pa_y[i] - ta2_xy_y_xxxx_1[i] * pc_y[i];

        ta2_xy_yy_xxxy_0[i] = ta2_xy_0_xxxy_0[i] * fe_0 - ta2_xy_0_xxxy_1[i] * fe_0 + ta2_xy_y_xxx_0[i] * fe_0 - ta2_xy_y_xxx_1[i] * fe_0 +
                              ta1_x_y_xxxy_1[i] + ta2_xy_y_xxxy_0[i] * pa_y[i] - ta2_xy_y_xxxy_1[i] * pc_y[i];

        ta2_xy_yy_xxxz_0[i] =
            ta2_xy_0_xxxz_0[i] * fe_0 - ta2_xy_0_xxxz_1[i] * fe_0 + ta1_x_y_xxxz_1[i] + ta2_xy_y_xxxz_0[i] * pa_y[i] - ta2_xy_y_xxxz_1[i] * pc_y[i];

        ta2_xy_yy_xxyy_0[i] = ta2_xy_0_xxyy_0[i] * fe_0 - ta2_xy_0_xxyy_1[i] * fe_0 + 2.0 * ta2_xy_y_xxy_0[i] * fe_0 -
                              2.0 * ta2_xy_y_xxy_1[i] * fe_0 + ta1_x_y_xxyy_1[i] + ta2_xy_y_xxyy_0[i] * pa_y[i] - ta2_xy_y_xxyy_1[i] * pc_y[i];

        ta2_xy_yy_xxyz_0[i] = ta2_xy_0_xxyz_0[i] * fe_0 - ta2_xy_0_xxyz_1[i] * fe_0 + ta2_xy_y_xxz_0[i] * fe_0 - ta2_xy_y_xxz_1[i] * fe_0 +
                              ta1_x_y_xxyz_1[i] + ta2_xy_y_xxyz_0[i] * pa_y[i] - ta2_xy_y_xxyz_1[i] * pc_y[i];

        ta2_xy_yy_xxzz_0[i] =
            ta2_xy_0_xxzz_0[i] * fe_0 - ta2_xy_0_xxzz_1[i] * fe_0 + ta1_x_y_xxzz_1[i] + ta2_xy_y_xxzz_0[i] * pa_y[i] - ta2_xy_y_xxzz_1[i] * pc_y[i];

        ta2_xy_yy_xyyy_0[i] = ta2_xy_0_xyyy_0[i] * fe_0 - ta2_xy_0_xyyy_1[i] * fe_0 + 3.0 * ta2_xy_y_xyy_0[i] * fe_0 -
                              3.0 * ta2_xy_y_xyy_1[i] * fe_0 + ta1_x_y_xyyy_1[i] + ta2_xy_y_xyyy_0[i] * pa_y[i] - ta2_xy_y_xyyy_1[i] * pc_y[i];

        ta2_xy_yy_xyyz_0[i] = ta2_xy_0_xyyz_0[i] * fe_0 - ta2_xy_0_xyyz_1[i] * fe_0 + 2.0 * ta2_xy_y_xyz_0[i] * fe_0 -
                              2.0 * ta2_xy_y_xyz_1[i] * fe_0 + ta1_x_y_xyyz_1[i] + ta2_xy_y_xyyz_0[i] * pa_y[i] - ta2_xy_y_xyyz_1[i] * pc_y[i];

        ta2_xy_yy_xyzz_0[i] = ta2_xy_0_xyzz_0[i] * fe_0 - ta2_xy_0_xyzz_1[i] * fe_0 + ta2_xy_y_xzz_0[i] * fe_0 - ta2_xy_y_xzz_1[i] * fe_0 +
                              ta1_x_y_xyzz_1[i] + ta2_xy_y_xyzz_0[i] * pa_y[i] - ta2_xy_y_xyzz_1[i] * pc_y[i];

        ta2_xy_yy_xzzz_0[i] =
            ta2_xy_0_xzzz_0[i] * fe_0 - ta2_xy_0_xzzz_1[i] * fe_0 + ta1_x_y_xzzz_1[i] + ta2_xy_y_xzzz_0[i] * pa_y[i] - ta2_xy_y_xzzz_1[i] * pc_y[i];

        ta2_xy_yy_yyyy_0[i] = ta2_xy_0_yyyy_0[i] * fe_0 - ta2_xy_0_yyyy_1[i] * fe_0 + 4.0 * ta2_xy_y_yyy_0[i] * fe_0 -
                              4.0 * ta2_xy_y_yyy_1[i] * fe_0 + ta1_x_y_yyyy_1[i] + ta2_xy_y_yyyy_0[i] * pa_y[i] - ta2_xy_y_yyyy_1[i] * pc_y[i];

        ta2_xy_yy_yyyz_0[i] = ta2_xy_0_yyyz_0[i] * fe_0 - ta2_xy_0_yyyz_1[i] * fe_0 + 3.0 * ta2_xy_y_yyz_0[i] * fe_0 -
                              3.0 * ta2_xy_y_yyz_1[i] * fe_0 + ta1_x_y_yyyz_1[i] + ta2_xy_y_yyyz_0[i] * pa_y[i] - ta2_xy_y_yyyz_1[i] * pc_y[i];

        ta2_xy_yy_yyzz_0[i] = ta2_xy_0_yyzz_0[i] * fe_0 - ta2_xy_0_yyzz_1[i] * fe_0 + 2.0 * ta2_xy_y_yzz_0[i] * fe_0 -
                              2.0 * ta2_xy_y_yzz_1[i] * fe_0 + ta1_x_y_yyzz_1[i] + ta2_xy_y_yyzz_0[i] * pa_y[i] - ta2_xy_y_yyzz_1[i] * pc_y[i];

        ta2_xy_yy_yzzz_0[i] = ta2_xy_0_yzzz_0[i] * fe_0 - ta2_xy_0_yzzz_1[i] * fe_0 + ta2_xy_y_zzz_0[i] * fe_0 - ta2_xy_y_zzz_1[i] * fe_0 +
                              ta1_x_y_yzzz_1[i] + ta2_xy_y_yzzz_0[i] * pa_y[i] - ta2_xy_y_yzzz_1[i] * pc_y[i];

        ta2_xy_yy_zzzz_0[i] =
            ta2_xy_0_zzzz_0[i] * fe_0 - ta2_xy_0_zzzz_1[i] * fe_0 + ta1_x_y_zzzz_1[i] + ta2_xy_y_zzzz_0[i] * pa_y[i] - ta2_xy_y_zzzz_1[i] * pc_y[i];
    }

    // Set up 150-165 components of targeted buffer : DG

    auto ta2_xy_yz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 150);

    auto ta2_xy_yz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 151);

    auto ta2_xy_yz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 152);

    auto ta2_xy_yz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 153);

    auto ta2_xy_yz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 154);

    auto ta2_xy_yz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 155);

    auto ta2_xy_yz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 156);

    auto ta2_xy_yz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 157);

    auto ta2_xy_yz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 158);

    auto ta2_xy_yz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 159);

    auto ta2_xy_yz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 160);

    auto ta2_xy_yz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 161);

    auto ta2_xy_yz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 162);

    auto ta2_xy_yz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 163);

    auto ta2_xy_yz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 164);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_x_z_xxxz_1,   \
                             ta1_x_z_xxzz_1,   \
                             ta1_x_z_xzzz_1,   \
                             ta1_x_z_zzzz_1,   \
                             ta2_xy_y_xxxx_0,  \
                             ta2_xy_y_xxxx_1,  \
                             ta2_xy_y_xxxy_0,  \
                             ta2_xy_y_xxxy_1,  \
                             ta2_xy_y_xxy_0,   \
                             ta2_xy_y_xxy_1,   \
                             ta2_xy_y_xxyy_0,  \
                             ta2_xy_y_xxyy_1,  \
                             ta2_xy_y_xxyz_0,  \
                             ta2_xy_y_xxyz_1,  \
                             ta2_xy_y_xyy_0,   \
                             ta2_xy_y_xyy_1,   \
                             ta2_xy_y_xyyy_0,  \
                             ta2_xy_y_xyyy_1,  \
                             ta2_xy_y_xyyz_0,  \
                             ta2_xy_y_xyyz_1,  \
                             ta2_xy_y_xyz_0,   \
                             ta2_xy_y_xyz_1,   \
                             ta2_xy_y_xyzz_0,  \
                             ta2_xy_y_xyzz_1,  \
                             ta2_xy_y_yyy_0,   \
                             ta2_xy_y_yyy_1,   \
                             ta2_xy_y_yyyy_0,  \
                             ta2_xy_y_yyyy_1,  \
                             ta2_xy_y_yyyz_0,  \
                             ta2_xy_y_yyyz_1,  \
                             ta2_xy_y_yyz_0,   \
                             ta2_xy_y_yyz_1,   \
                             ta2_xy_y_yyzz_0,  \
                             ta2_xy_y_yyzz_1,  \
                             ta2_xy_y_yzz_0,   \
                             ta2_xy_y_yzz_1,   \
                             ta2_xy_y_yzzz_0,  \
                             ta2_xy_y_yzzz_1,  \
                             ta2_xy_yz_xxxx_0, \
                             ta2_xy_yz_xxxy_0, \
                             ta2_xy_yz_xxxz_0, \
                             ta2_xy_yz_xxyy_0, \
                             ta2_xy_yz_xxyz_0, \
                             ta2_xy_yz_xxzz_0, \
                             ta2_xy_yz_xyyy_0, \
                             ta2_xy_yz_xyyz_0, \
                             ta2_xy_yz_xyzz_0, \
                             ta2_xy_yz_xzzz_0, \
                             ta2_xy_yz_yyyy_0, \
                             ta2_xy_yz_yyyz_0, \
                             ta2_xy_yz_yyzz_0, \
                             ta2_xy_yz_yzzz_0, \
                             ta2_xy_yz_zzzz_0, \
                             ta2_xy_z_xxxz_0,  \
                             ta2_xy_z_xxxz_1,  \
                             ta2_xy_z_xxzz_0,  \
                             ta2_xy_z_xxzz_1,  \
                             ta2_xy_z_xzzz_0,  \
                             ta2_xy_z_xzzz_1,  \
                             ta2_xy_z_zzzz_0,  \
                             ta2_xy_z_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yz_xxxx_0[i] = ta2_xy_y_xxxx_0[i] * pa_z[i] - ta2_xy_y_xxxx_1[i] * pc_z[i];

        ta2_xy_yz_xxxy_0[i] = ta2_xy_y_xxxy_0[i] * pa_z[i] - ta2_xy_y_xxxy_1[i] * pc_z[i];

        ta2_xy_yz_xxxz_0[i] = ta1_x_z_xxxz_1[i] + ta2_xy_z_xxxz_0[i] * pa_y[i] - ta2_xy_z_xxxz_1[i] * pc_y[i];

        ta2_xy_yz_xxyy_0[i] = ta2_xy_y_xxyy_0[i] * pa_z[i] - ta2_xy_y_xxyy_1[i] * pc_z[i];

        ta2_xy_yz_xxyz_0[i] = ta2_xy_y_xxy_0[i] * fe_0 - ta2_xy_y_xxy_1[i] * fe_0 + ta2_xy_y_xxyz_0[i] * pa_z[i] - ta2_xy_y_xxyz_1[i] * pc_z[i];

        ta2_xy_yz_xxzz_0[i] = ta1_x_z_xxzz_1[i] + ta2_xy_z_xxzz_0[i] * pa_y[i] - ta2_xy_z_xxzz_1[i] * pc_y[i];

        ta2_xy_yz_xyyy_0[i] = ta2_xy_y_xyyy_0[i] * pa_z[i] - ta2_xy_y_xyyy_1[i] * pc_z[i];

        ta2_xy_yz_xyyz_0[i] = ta2_xy_y_xyy_0[i] * fe_0 - ta2_xy_y_xyy_1[i] * fe_0 + ta2_xy_y_xyyz_0[i] * pa_z[i] - ta2_xy_y_xyyz_1[i] * pc_z[i];

        ta2_xy_yz_xyzz_0[i] =
            2.0 * ta2_xy_y_xyz_0[i] * fe_0 - 2.0 * ta2_xy_y_xyz_1[i] * fe_0 + ta2_xy_y_xyzz_0[i] * pa_z[i] - ta2_xy_y_xyzz_1[i] * pc_z[i];

        ta2_xy_yz_xzzz_0[i] = ta1_x_z_xzzz_1[i] + ta2_xy_z_xzzz_0[i] * pa_y[i] - ta2_xy_z_xzzz_1[i] * pc_y[i];

        ta2_xy_yz_yyyy_0[i] = ta2_xy_y_yyyy_0[i] * pa_z[i] - ta2_xy_y_yyyy_1[i] * pc_z[i];

        ta2_xy_yz_yyyz_0[i] = ta2_xy_y_yyy_0[i] * fe_0 - ta2_xy_y_yyy_1[i] * fe_0 + ta2_xy_y_yyyz_0[i] * pa_z[i] - ta2_xy_y_yyyz_1[i] * pc_z[i];

        ta2_xy_yz_yyzz_0[i] =
            2.0 * ta2_xy_y_yyz_0[i] * fe_0 - 2.0 * ta2_xy_y_yyz_1[i] * fe_0 + ta2_xy_y_yyzz_0[i] * pa_z[i] - ta2_xy_y_yyzz_1[i] * pc_z[i];

        ta2_xy_yz_yzzz_0[i] =
            3.0 * ta2_xy_y_yzz_0[i] * fe_0 - 3.0 * ta2_xy_y_yzz_1[i] * fe_0 + ta2_xy_y_yzzz_0[i] * pa_z[i] - ta2_xy_y_yzzz_1[i] * pc_z[i];

        ta2_xy_yz_zzzz_0[i] = ta1_x_z_zzzz_1[i] + ta2_xy_z_zzzz_0[i] * pa_y[i] - ta2_xy_z_zzzz_1[i] * pc_y[i];
    }

    // Set up 165-180 components of targeted buffer : DG

    auto ta2_xy_zz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 165);

    auto ta2_xy_zz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 166);

    auto ta2_xy_zz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 167);

    auto ta2_xy_zz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 168);

    auto ta2_xy_zz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 169);

    auto ta2_xy_zz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 170);

    auto ta2_xy_zz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 171);

    auto ta2_xy_zz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 172);

    auto ta2_xy_zz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 173);

    auto ta2_xy_zz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 174);

    auto ta2_xy_zz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 175);

    auto ta2_xy_zz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 176);

    auto ta2_xy_zz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 177);

    auto ta2_xy_zz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 178);

    auto ta2_xy_zz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 179);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta2_xy_0_xxxx_0,  \
                             ta2_xy_0_xxxx_1,  \
                             ta2_xy_0_xxxy_0,  \
                             ta2_xy_0_xxxy_1,  \
                             ta2_xy_0_xxxz_0,  \
                             ta2_xy_0_xxxz_1,  \
                             ta2_xy_0_xxyy_0,  \
                             ta2_xy_0_xxyy_1,  \
                             ta2_xy_0_xxyz_0,  \
                             ta2_xy_0_xxyz_1,  \
                             ta2_xy_0_xxzz_0,  \
                             ta2_xy_0_xxzz_1,  \
                             ta2_xy_0_xyyy_0,  \
                             ta2_xy_0_xyyy_1,  \
                             ta2_xy_0_xyyz_0,  \
                             ta2_xy_0_xyyz_1,  \
                             ta2_xy_0_xyzz_0,  \
                             ta2_xy_0_xyzz_1,  \
                             ta2_xy_0_xzzz_0,  \
                             ta2_xy_0_xzzz_1,  \
                             ta2_xy_0_yyyy_0,  \
                             ta2_xy_0_yyyy_1,  \
                             ta2_xy_0_yyyz_0,  \
                             ta2_xy_0_yyyz_1,  \
                             ta2_xy_0_yyzz_0,  \
                             ta2_xy_0_yyzz_1,  \
                             ta2_xy_0_yzzz_0,  \
                             ta2_xy_0_yzzz_1,  \
                             ta2_xy_0_zzzz_0,  \
                             ta2_xy_0_zzzz_1,  \
                             ta2_xy_z_xxx_0,   \
                             ta2_xy_z_xxx_1,   \
                             ta2_xy_z_xxxx_0,  \
                             ta2_xy_z_xxxx_1,  \
                             ta2_xy_z_xxxy_0,  \
                             ta2_xy_z_xxxy_1,  \
                             ta2_xy_z_xxxz_0,  \
                             ta2_xy_z_xxxz_1,  \
                             ta2_xy_z_xxy_0,   \
                             ta2_xy_z_xxy_1,   \
                             ta2_xy_z_xxyy_0,  \
                             ta2_xy_z_xxyy_1,  \
                             ta2_xy_z_xxyz_0,  \
                             ta2_xy_z_xxyz_1,  \
                             ta2_xy_z_xxz_0,   \
                             ta2_xy_z_xxz_1,   \
                             ta2_xy_z_xxzz_0,  \
                             ta2_xy_z_xxzz_1,  \
                             ta2_xy_z_xyy_0,   \
                             ta2_xy_z_xyy_1,   \
                             ta2_xy_z_xyyy_0,  \
                             ta2_xy_z_xyyy_1,  \
                             ta2_xy_z_xyyz_0,  \
                             ta2_xy_z_xyyz_1,  \
                             ta2_xy_z_xyz_0,   \
                             ta2_xy_z_xyz_1,   \
                             ta2_xy_z_xyzz_0,  \
                             ta2_xy_z_xyzz_1,  \
                             ta2_xy_z_xzz_0,   \
                             ta2_xy_z_xzz_1,   \
                             ta2_xy_z_xzzz_0,  \
                             ta2_xy_z_xzzz_1,  \
                             ta2_xy_z_yyy_0,   \
                             ta2_xy_z_yyy_1,   \
                             ta2_xy_z_yyyy_0,  \
                             ta2_xy_z_yyyy_1,  \
                             ta2_xy_z_yyyz_0,  \
                             ta2_xy_z_yyyz_1,  \
                             ta2_xy_z_yyz_0,   \
                             ta2_xy_z_yyz_1,   \
                             ta2_xy_z_yyzz_0,  \
                             ta2_xy_z_yyzz_1,  \
                             ta2_xy_z_yzz_0,   \
                             ta2_xy_z_yzz_1,   \
                             ta2_xy_z_yzzz_0,  \
                             ta2_xy_z_yzzz_1,  \
                             ta2_xy_z_zzz_0,   \
                             ta2_xy_z_zzz_1,   \
                             ta2_xy_z_zzzz_0,  \
                             ta2_xy_z_zzzz_1,  \
                             ta2_xy_zz_xxxx_0, \
                             ta2_xy_zz_xxxy_0, \
                             ta2_xy_zz_xxxz_0, \
                             ta2_xy_zz_xxyy_0, \
                             ta2_xy_zz_xxyz_0, \
                             ta2_xy_zz_xxzz_0, \
                             ta2_xy_zz_xyyy_0, \
                             ta2_xy_zz_xyyz_0, \
                             ta2_xy_zz_xyzz_0, \
                             ta2_xy_zz_xzzz_0, \
                             ta2_xy_zz_yyyy_0, \
                             ta2_xy_zz_yyyz_0, \
                             ta2_xy_zz_yyzz_0, \
                             ta2_xy_zz_yzzz_0, \
                             ta2_xy_zz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_zz_xxxx_0[i] = ta2_xy_0_xxxx_0[i] * fe_0 - ta2_xy_0_xxxx_1[i] * fe_0 + ta2_xy_z_xxxx_0[i] * pa_z[i] - ta2_xy_z_xxxx_1[i] * pc_z[i];

        ta2_xy_zz_xxxy_0[i] = ta2_xy_0_xxxy_0[i] * fe_0 - ta2_xy_0_xxxy_1[i] * fe_0 + ta2_xy_z_xxxy_0[i] * pa_z[i] - ta2_xy_z_xxxy_1[i] * pc_z[i];

        ta2_xy_zz_xxxz_0[i] = ta2_xy_0_xxxz_0[i] * fe_0 - ta2_xy_0_xxxz_1[i] * fe_0 + ta2_xy_z_xxx_0[i] * fe_0 - ta2_xy_z_xxx_1[i] * fe_0 +
                              ta2_xy_z_xxxz_0[i] * pa_z[i] - ta2_xy_z_xxxz_1[i] * pc_z[i];

        ta2_xy_zz_xxyy_0[i] = ta2_xy_0_xxyy_0[i] * fe_0 - ta2_xy_0_xxyy_1[i] * fe_0 + ta2_xy_z_xxyy_0[i] * pa_z[i] - ta2_xy_z_xxyy_1[i] * pc_z[i];

        ta2_xy_zz_xxyz_0[i] = ta2_xy_0_xxyz_0[i] * fe_0 - ta2_xy_0_xxyz_1[i] * fe_0 + ta2_xy_z_xxy_0[i] * fe_0 - ta2_xy_z_xxy_1[i] * fe_0 +
                              ta2_xy_z_xxyz_0[i] * pa_z[i] - ta2_xy_z_xxyz_1[i] * pc_z[i];

        ta2_xy_zz_xxzz_0[i] = ta2_xy_0_xxzz_0[i] * fe_0 - ta2_xy_0_xxzz_1[i] * fe_0 + 2.0 * ta2_xy_z_xxz_0[i] * fe_0 -
                              2.0 * ta2_xy_z_xxz_1[i] * fe_0 + ta2_xy_z_xxzz_0[i] * pa_z[i] - ta2_xy_z_xxzz_1[i] * pc_z[i];

        ta2_xy_zz_xyyy_0[i] = ta2_xy_0_xyyy_0[i] * fe_0 - ta2_xy_0_xyyy_1[i] * fe_0 + ta2_xy_z_xyyy_0[i] * pa_z[i] - ta2_xy_z_xyyy_1[i] * pc_z[i];

        ta2_xy_zz_xyyz_0[i] = ta2_xy_0_xyyz_0[i] * fe_0 - ta2_xy_0_xyyz_1[i] * fe_0 + ta2_xy_z_xyy_0[i] * fe_0 - ta2_xy_z_xyy_1[i] * fe_0 +
                              ta2_xy_z_xyyz_0[i] * pa_z[i] - ta2_xy_z_xyyz_1[i] * pc_z[i];

        ta2_xy_zz_xyzz_0[i] = ta2_xy_0_xyzz_0[i] * fe_0 - ta2_xy_0_xyzz_1[i] * fe_0 + 2.0 * ta2_xy_z_xyz_0[i] * fe_0 -
                              2.0 * ta2_xy_z_xyz_1[i] * fe_0 + ta2_xy_z_xyzz_0[i] * pa_z[i] - ta2_xy_z_xyzz_1[i] * pc_z[i];

        ta2_xy_zz_xzzz_0[i] = ta2_xy_0_xzzz_0[i] * fe_0 - ta2_xy_0_xzzz_1[i] * fe_0 + 3.0 * ta2_xy_z_xzz_0[i] * fe_0 -
                              3.0 * ta2_xy_z_xzz_1[i] * fe_0 + ta2_xy_z_xzzz_0[i] * pa_z[i] - ta2_xy_z_xzzz_1[i] * pc_z[i];

        ta2_xy_zz_yyyy_0[i] = ta2_xy_0_yyyy_0[i] * fe_0 - ta2_xy_0_yyyy_1[i] * fe_0 + ta2_xy_z_yyyy_0[i] * pa_z[i] - ta2_xy_z_yyyy_1[i] * pc_z[i];

        ta2_xy_zz_yyyz_0[i] = ta2_xy_0_yyyz_0[i] * fe_0 - ta2_xy_0_yyyz_1[i] * fe_0 + ta2_xy_z_yyy_0[i] * fe_0 - ta2_xy_z_yyy_1[i] * fe_0 +
                              ta2_xy_z_yyyz_0[i] * pa_z[i] - ta2_xy_z_yyyz_1[i] * pc_z[i];

        ta2_xy_zz_yyzz_0[i] = ta2_xy_0_yyzz_0[i] * fe_0 - ta2_xy_0_yyzz_1[i] * fe_0 + 2.0 * ta2_xy_z_yyz_0[i] * fe_0 -
                              2.0 * ta2_xy_z_yyz_1[i] * fe_0 + ta2_xy_z_yyzz_0[i] * pa_z[i] - ta2_xy_z_yyzz_1[i] * pc_z[i];

        ta2_xy_zz_yzzz_0[i] = ta2_xy_0_yzzz_0[i] * fe_0 - ta2_xy_0_yzzz_1[i] * fe_0 + 3.0 * ta2_xy_z_yzz_0[i] * fe_0 -
                              3.0 * ta2_xy_z_yzz_1[i] * fe_0 + ta2_xy_z_yzzz_0[i] * pa_z[i] - ta2_xy_z_yzzz_1[i] * pc_z[i];

        ta2_xy_zz_zzzz_0[i] = ta2_xy_0_zzzz_0[i] * fe_0 - ta2_xy_0_zzzz_1[i] * fe_0 + 4.0 * ta2_xy_z_zzz_0[i] * fe_0 -
                              4.0 * ta2_xy_z_zzz_1[i] * fe_0 + ta2_xy_z_zzzz_0[i] * pa_z[i] - ta2_xy_z_zzzz_1[i] * pc_z[i];
    }

    // Set up 180-195 components of targeted buffer : DG

    auto ta2_xz_xx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 180);

    auto ta2_xz_xx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 181);

    auto ta2_xz_xx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 182);

    auto ta2_xz_xx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 183);

    auto ta2_xz_xx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 184);

    auto ta2_xz_xx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 185);

    auto ta2_xz_xx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 186);

    auto ta2_xz_xx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 187);

    auto ta2_xz_xx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 188);

    auto ta2_xz_xx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 189);

    auto ta2_xz_xx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 190);

    auto ta2_xz_xx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 191);

    auto ta2_xz_xx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 192);

    auto ta2_xz_xx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 193);

    auto ta2_xz_xx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 194);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_z_x_xxxx_1,   \
                             ta1_z_x_xxxy_1,   \
                             ta1_z_x_xxxz_1,   \
                             ta1_z_x_xxyy_1,   \
                             ta1_z_x_xxyz_1,   \
                             ta1_z_x_xxzz_1,   \
                             ta1_z_x_xyyy_1,   \
                             ta1_z_x_xyyz_1,   \
                             ta1_z_x_xyzz_1,   \
                             ta1_z_x_xzzz_1,   \
                             ta1_z_x_yyyy_1,   \
                             ta1_z_x_yyyz_1,   \
                             ta1_z_x_yyzz_1,   \
                             ta1_z_x_yzzz_1,   \
                             ta1_z_x_zzzz_1,   \
                             ta2_xz_0_xxxx_0,  \
                             ta2_xz_0_xxxx_1,  \
                             ta2_xz_0_xxxy_0,  \
                             ta2_xz_0_xxxy_1,  \
                             ta2_xz_0_xxxz_0,  \
                             ta2_xz_0_xxxz_1,  \
                             ta2_xz_0_xxyy_0,  \
                             ta2_xz_0_xxyy_1,  \
                             ta2_xz_0_xxyz_0,  \
                             ta2_xz_0_xxyz_1,  \
                             ta2_xz_0_xxzz_0,  \
                             ta2_xz_0_xxzz_1,  \
                             ta2_xz_0_xyyy_0,  \
                             ta2_xz_0_xyyy_1,  \
                             ta2_xz_0_xyyz_0,  \
                             ta2_xz_0_xyyz_1,  \
                             ta2_xz_0_xyzz_0,  \
                             ta2_xz_0_xyzz_1,  \
                             ta2_xz_0_xzzz_0,  \
                             ta2_xz_0_xzzz_1,  \
                             ta2_xz_0_yyyy_0,  \
                             ta2_xz_0_yyyy_1,  \
                             ta2_xz_0_yyyz_0,  \
                             ta2_xz_0_yyyz_1,  \
                             ta2_xz_0_yyzz_0,  \
                             ta2_xz_0_yyzz_1,  \
                             ta2_xz_0_yzzz_0,  \
                             ta2_xz_0_yzzz_1,  \
                             ta2_xz_0_zzzz_0,  \
                             ta2_xz_0_zzzz_1,  \
                             ta2_xz_x_xxx_0,   \
                             ta2_xz_x_xxx_1,   \
                             ta2_xz_x_xxxx_0,  \
                             ta2_xz_x_xxxx_1,  \
                             ta2_xz_x_xxxy_0,  \
                             ta2_xz_x_xxxy_1,  \
                             ta2_xz_x_xxxz_0,  \
                             ta2_xz_x_xxxz_1,  \
                             ta2_xz_x_xxy_0,   \
                             ta2_xz_x_xxy_1,   \
                             ta2_xz_x_xxyy_0,  \
                             ta2_xz_x_xxyy_1,  \
                             ta2_xz_x_xxyz_0,  \
                             ta2_xz_x_xxyz_1,  \
                             ta2_xz_x_xxz_0,   \
                             ta2_xz_x_xxz_1,   \
                             ta2_xz_x_xxzz_0,  \
                             ta2_xz_x_xxzz_1,  \
                             ta2_xz_x_xyy_0,   \
                             ta2_xz_x_xyy_1,   \
                             ta2_xz_x_xyyy_0,  \
                             ta2_xz_x_xyyy_1,  \
                             ta2_xz_x_xyyz_0,  \
                             ta2_xz_x_xyyz_1,  \
                             ta2_xz_x_xyz_0,   \
                             ta2_xz_x_xyz_1,   \
                             ta2_xz_x_xyzz_0,  \
                             ta2_xz_x_xyzz_1,  \
                             ta2_xz_x_xzz_0,   \
                             ta2_xz_x_xzz_1,   \
                             ta2_xz_x_xzzz_0,  \
                             ta2_xz_x_xzzz_1,  \
                             ta2_xz_x_yyy_0,   \
                             ta2_xz_x_yyy_1,   \
                             ta2_xz_x_yyyy_0,  \
                             ta2_xz_x_yyyy_1,  \
                             ta2_xz_x_yyyz_0,  \
                             ta2_xz_x_yyyz_1,  \
                             ta2_xz_x_yyz_0,   \
                             ta2_xz_x_yyz_1,   \
                             ta2_xz_x_yyzz_0,  \
                             ta2_xz_x_yyzz_1,  \
                             ta2_xz_x_yzz_0,   \
                             ta2_xz_x_yzz_1,   \
                             ta2_xz_x_yzzz_0,  \
                             ta2_xz_x_yzzz_1,  \
                             ta2_xz_x_zzz_0,   \
                             ta2_xz_x_zzz_1,   \
                             ta2_xz_x_zzzz_0,  \
                             ta2_xz_x_zzzz_1,  \
                             ta2_xz_xx_xxxx_0, \
                             ta2_xz_xx_xxxy_0, \
                             ta2_xz_xx_xxxz_0, \
                             ta2_xz_xx_xxyy_0, \
                             ta2_xz_xx_xxyz_0, \
                             ta2_xz_xx_xxzz_0, \
                             ta2_xz_xx_xyyy_0, \
                             ta2_xz_xx_xyyz_0, \
                             ta2_xz_xx_xyzz_0, \
                             ta2_xz_xx_xzzz_0, \
                             ta2_xz_xx_yyyy_0, \
                             ta2_xz_xx_yyyz_0, \
                             ta2_xz_xx_yyzz_0, \
                             ta2_xz_xx_yzzz_0, \
                             ta2_xz_xx_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xx_xxxx_0[i] = ta2_xz_0_xxxx_0[i] * fe_0 - ta2_xz_0_xxxx_1[i] * fe_0 + 4.0 * ta2_xz_x_xxx_0[i] * fe_0 -
                              4.0 * ta2_xz_x_xxx_1[i] * fe_0 + ta1_z_x_xxxx_1[i] + ta2_xz_x_xxxx_0[i] * pa_x[i] - ta2_xz_x_xxxx_1[i] * pc_x[i];

        ta2_xz_xx_xxxy_0[i] = ta2_xz_0_xxxy_0[i] * fe_0 - ta2_xz_0_xxxy_1[i] * fe_0 + 3.0 * ta2_xz_x_xxy_0[i] * fe_0 -
                              3.0 * ta2_xz_x_xxy_1[i] * fe_0 + ta1_z_x_xxxy_1[i] + ta2_xz_x_xxxy_0[i] * pa_x[i] - ta2_xz_x_xxxy_1[i] * pc_x[i];

        ta2_xz_xx_xxxz_0[i] = ta2_xz_0_xxxz_0[i] * fe_0 - ta2_xz_0_xxxz_1[i] * fe_0 + 3.0 * ta2_xz_x_xxz_0[i] * fe_0 -
                              3.0 * ta2_xz_x_xxz_1[i] * fe_0 + ta1_z_x_xxxz_1[i] + ta2_xz_x_xxxz_0[i] * pa_x[i] - ta2_xz_x_xxxz_1[i] * pc_x[i];

        ta2_xz_xx_xxyy_0[i] = ta2_xz_0_xxyy_0[i] * fe_0 - ta2_xz_0_xxyy_1[i] * fe_0 + 2.0 * ta2_xz_x_xyy_0[i] * fe_0 -
                              2.0 * ta2_xz_x_xyy_1[i] * fe_0 + ta1_z_x_xxyy_1[i] + ta2_xz_x_xxyy_0[i] * pa_x[i] - ta2_xz_x_xxyy_1[i] * pc_x[i];

        ta2_xz_xx_xxyz_0[i] = ta2_xz_0_xxyz_0[i] * fe_0 - ta2_xz_0_xxyz_1[i] * fe_0 + 2.0 * ta2_xz_x_xyz_0[i] * fe_0 -
                              2.0 * ta2_xz_x_xyz_1[i] * fe_0 + ta1_z_x_xxyz_1[i] + ta2_xz_x_xxyz_0[i] * pa_x[i] - ta2_xz_x_xxyz_1[i] * pc_x[i];

        ta2_xz_xx_xxzz_0[i] = ta2_xz_0_xxzz_0[i] * fe_0 - ta2_xz_0_xxzz_1[i] * fe_0 + 2.0 * ta2_xz_x_xzz_0[i] * fe_0 -
                              2.0 * ta2_xz_x_xzz_1[i] * fe_0 + ta1_z_x_xxzz_1[i] + ta2_xz_x_xxzz_0[i] * pa_x[i] - ta2_xz_x_xxzz_1[i] * pc_x[i];

        ta2_xz_xx_xyyy_0[i] = ta2_xz_0_xyyy_0[i] * fe_0 - ta2_xz_0_xyyy_1[i] * fe_0 + ta2_xz_x_yyy_0[i] * fe_0 - ta2_xz_x_yyy_1[i] * fe_0 +
                              ta1_z_x_xyyy_1[i] + ta2_xz_x_xyyy_0[i] * pa_x[i] - ta2_xz_x_xyyy_1[i] * pc_x[i];

        ta2_xz_xx_xyyz_0[i] = ta2_xz_0_xyyz_0[i] * fe_0 - ta2_xz_0_xyyz_1[i] * fe_0 + ta2_xz_x_yyz_0[i] * fe_0 - ta2_xz_x_yyz_1[i] * fe_0 +
                              ta1_z_x_xyyz_1[i] + ta2_xz_x_xyyz_0[i] * pa_x[i] - ta2_xz_x_xyyz_1[i] * pc_x[i];

        ta2_xz_xx_xyzz_0[i] = ta2_xz_0_xyzz_0[i] * fe_0 - ta2_xz_0_xyzz_1[i] * fe_0 + ta2_xz_x_yzz_0[i] * fe_0 - ta2_xz_x_yzz_1[i] * fe_0 +
                              ta1_z_x_xyzz_1[i] + ta2_xz_x_xyzz_0[i] * pa_x[i] - ta2_xz_x_xyzz_1[i] * pc_x[i];

        ta2_xz_xx_xzzz_0[i] = ta2_xz_0_xzzz_0[i] * fe_0 - ta2_xz_0_xzzz_1[i] * fe_0 + ta2_xz_x_zzz_0[i] * fe_0 - ta2_xz_x_zzz_1[i] * fe_0 +
                              ta1_z_x_xzzz_1[i] + ta2_xz_x_xzzz_0[i] * pa_x[i] - ta2_xz_x_xzzz_1[i] * pc_x[i];

        ta2_xz_xx_yyyy_0[i] =
            ta2_xz_0_yyyy_0[i] * fe_0 - ta2_xz_0_yyyy_1[i] * fe_0 + ta1_z_x_yyyy_1[i] + ta2_xz_x_yyyy_0[i] * pa_x[i] - ta2_xz_x_yyyy_1[i] * pc_x[i];

        ta2_xz_xx_yyyz_0[i] =
            ta2_xz_0_yyyz_0[i] * fe_0 - ta2_xz_0_yyyz_1[i] * fe_0 + ta1_z_x_yyyz_1[i] + ta2_xz_x_yyyz_0[i] * pa_x[i] - ta2_xz_x_yyyz_1[i] * pc_x[i];

        ta2_xz_xx_yyzz_0[i] =
            ta2_xz_0_yyzz_0[i] * fe_0 - ta2_xz_0_yyzz_1[i] * fe_0 + ta1_z_x_yyzz_1[i] + ta2_xz_x_yyzz_0[i] * pa_x[i] - ta2_xz_x_yyzz_1[i] * pc_x[i];

        ta2_xz_xx_yzzz_0[i] =
            ta2_xz_0_yzzz_0[i] * fe_0 - ta2_xz_0_yzzz_1[i] * fe_0 + ta1_z_x_yzzz_1[i] + ta2_xz_x_yzzz_0[i] * pa_x[i] - ta2_xz_x_yzzz_1[i] * pc_x[i];

        ta2_xz_xx_zzzz_0[i] =
            ta2_xz_0_zzzz_0[i] * fe_0 - ta2_xz_0_zzzz_1[i] * fe_0 + ta1_z_x_zzzz_1[i] + ta2_xz_x_zzzz_0[i] * pa_x[i] - ta2_xz_x_zzzz_1[i] * pc_x[i];
    }

    // Set up 195-210 components of targeted buffer : DG

    auto ta2_xz_xy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 195);

    auto ta2_xz_xy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 196);

    auto ta2_xz_xy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 197);

    auto ta2_xz_xy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 198);

    auto ta2_xz_xy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 199);

    auto ta2_xz_xy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 200);

    auto ta2_xz_xy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 201);

    auto ta2_xz_xy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 202);

    auto ta2_xz_xy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 203);

    auto ta2_xz_xy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 204);

    auto ta2_xz_xy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 205);

    auto ta2_xz_xy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 206);

    auto ta2_xz_xy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 207);

    auto ta2_xz_xy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 208);

    auto ta2_xz_xy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 209);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_z_y_yyyy_1,   \
                             ta1_z_y_yyyz_1,   \
                             ta1_z_y_yyzz_1,   \
                             ta1_z_y_yzzz_1,   \
                             ta2_xz_x_xxx_0,   \
                             ta2_xz_x_xxx_1,   \
                             ta2_xz_x_xxxx_0,  \
                             ta2_xz_x_xxxx_1,  \
                             ta2_xz_x_xxxy_0,  \
                             ta2_xz_x_xxxy_1,  \
                             ta2_xz_x_xxxz_0,  \
                             ta2_xz_x_xxxz_1,  \
                             ta2_xz_x_xxy_0,   \
                             ta2_xz_x_xxy_1,   \
                             ta2_xz_x_xxyy_0,  \
                             ta2_xz_x_xxyy_1,  \
                             ta2_xz_x_xxyz_0,  \
                             ta2_xz_x_xxyz_1,  \
                             ta2_xz_x_xxz_0,   \
                             ta2_xz_x_xxz_1,   \
                             ta2_xz_x_xxzz_0,  \
                             ta2_xz_x_xxzz_1,  \
                             ta2_xz_x_xyy_0,   \
                             ta2_xz_x_xyy_1,   \
                             ta2_xz_x_xyyy_0,  \
                             ta2_xz_x_xyyy_1,  \
                             ta2_xz_x_xyyz_0,  \
                             ta2_xz_x_xyyz_1,  \
                             ta2_xz_x_xyz_0,   \
                             ta2_xz_x_xyz_1,   \
                             ta2_xz_x_xyzz_0,  \
                             ta2_xz_x_xyzz_1,  \
                             ta2_xz_x_xzz_0,   \
                             ta2_xz_x_xzz_1,   \
                             ta2_xz_x_xzzz_0,  \
                             ta2_xz_x_xzzz_1,  \
                             ta2_xz_x_zzzz_0,  \
                             ta2_xz_x_zzzz_1,  \
                             ta2_xz_xy_xxxx_0, \
                             ta2_xz_xy_xxxy_0, \
                             ta2_xz_xy_xxxz_0, \
                             ta2_xz_xy_xxyy_0, \
                             ta2_xz_xy_xxyz_0, \
                             ta2_xz_xy_xxzz_0, \
                             ta2_xz_xy_xyyy_0, \
                             ta2_xz_xy_xyyz_0, \
                             ta2_xz_xy_xyzz_0, \
                             ta2_xz_xy_xzzz_0, \
                             ta2_xz_xy_yyyy_0, \
                             ta2_xz_xy_yyyz_0, \
                             ta2_xz_xy_yyzz_0, \
                             ta2_xz_xy_yzzz_0, \
                             ta2_xz_xy_zzzz_0, \
                             ta2_xz_y_yyyy_0,  \
                             ta2_xz_y_yyyy_1,  \
                             ta2_xz_y_yyyz_0,  \
                             ta2_xz_y_yyyz_1,  \
                             ta2_xz_y_yyzz_0,  \
                             ta2_xz_y_yyzz_1,  \
                             ta2_xz_y_yzzz_0,  \
                             ta2_xz_y_yzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xy_xxxx_0[i] = ta2_xz_x_xxxx_0[i] * pa_y[i] - ta2_xz_x_xxxx_1[i] * pc_y[i];

        ta2_xz_xy_xxxy_0[i] = ta2_xz_x_xxx_0[i] * fe_0 - ta2_xz_x_xxx_1[i] * fe_0 + ta2_xz_x_xxxy_0[i] * pa_y[i] - ta2_xz_x_xxxy_1[i] * pc_y[i];

        ta2_xz_xy_xxxz_0[i] = ta2_xz_x_xxxz_0[i] * pa_y[i] - ta2_xz_x_xxxz_1[i] * pc_y[i];

        ta2_xz_xy_xxyy_0[i] =
            2.0 * ta2_xz_x_xxy_0[i] * fe_0 - 2.0 * ta2_xz_x_xxy_1[i] * fe_0 + ta2_xz_x_xxyy_0[i] * pa_y[i] - ta2_xz_x_xxyy_1[i] * pc_y[i];

        ta2_xz_xy_xxyz_0[i] = ta2_xz_x_xxz_0[i] * fe_0 - ta2_xz_x_xxz_1[i] * fe_0 + ta2_xz_x_xxyz_0[i] * pa_y[i] - ta2_xz_x_xxyz_1[i] * pc_y[i];

        ta2_xz_xy_xxzz_0[i] = ta2_xz_x_xxzz_0[i] * pa_y[i] - ta2_xz_x_xxzz_1[i] * pc_y[i];

        ta2_xz_xy_xyyy_0[i] =
            3.0 * ta2_xz_x_xyy_0[i] * fe_0 - 3.0 * ta2_xz_x_xyy_1[i] * fe_0 + ta2_xz_x_xyyy_0[i] * pa_y[i] - ta2_xz_x_xyyy_1[i] * pc_y[i];

        ta2_xz_xy_xyyz_0[i] =
            2.0 * ta2_xz_x_xyz_0[i] * fe_0 - 2.0 * ta2_xz_x_xyz_1[i] * fe_0 + ta2_xz_x_xyyz_0[i] * pa_y[i] - ta2_xz_x_xyyz_1[i] * pc_y[i];

        ta2_xz_xy_xyzz_0[i] = ta2_xz_x_xzz_0[i] * fe_0 - ta2_xz_x_xzz_1[i] * fe_0 + ta2_xz_x_xyzz_0[i] * pa_y[i] - ta2_xz_x_xyzz_1[i] * pc_y[i];

        ta2_xz_xy_xzzz_0[i] = ta2_xz_x_xzzz_0[i] * pa_y[i] - ta2_xz_x_xzzz_1[i] * pc_y[i];

        ta2_xz_xy_yyyy_0[i] = ta1_z_y_yyyy_1[i] + ta2_xz_y_yyyy_0[i] * pa_x[i] - ta2_xz_y_yyyy_1[i] * pc_x[i];

        ta2_xz_xy_yyyz_0[i] = ta1_z_y_yyyz_1[i] + ta2_xz_y_yyyz_0[i] * pa_x[i] - ta2_xz_y_yyyz_1[i] * pc_x[i];

        ta2_xz_xy_yyzz_0[i] = ta1_z_y_yyzz_1[i] + ta2_xz_y_yyzz_0[i] * pa_x[i] - ta2_xz_y_yyzz_1[i] * pc_x[i];

        ta2_xz_xy_yzzz_0[i] = ta1_z_y_yzzz_1[i] + ta2_xz_y_yzzz_0[i] * pa_x[i] - ta2_xz_y_yzzz_1[i] * pc_x[i];

        ta2_xz_xy_zzzz_0[i] = ta2_xz_x_zzzz_0[i] * pa_y[i] - ta2_xz_x_zzzz_1[i] * pc_y[i];
    }

    // Set up 210-225 components of targeted buffer : DG

    auto ta2_xz_xz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 210);

    auto ta2_xz_xz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 211);

    auto ta2_xz_xz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 212);

    auto ta2_xz_xz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 213);

    auto ta2_xz_xz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 214);

    auto ta2_xz_xz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 215);

    auto ta2_xz_xz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 216);

    auto ta2_xz_xz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 217);

    auto ta2_xz_xz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 218);

    auto ta2_xz_xz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 219);

    auto ta2_xz_xz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 220);

    auto ta2_xz_xz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 221);

    auto ta2_xz_xz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 222);

    auto ta2_xz_xz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 223);

    auto ta2_xz_xz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 224);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_x_x_xxxx_1,   \
                             ta1_x_x_xxxy_1,   \
                             ta1_x_x_xxyy_1,   \
                             ta1_x_x_xyyy_1,   \
                             ta1_z_z_xxxz_1,   \
                             ta1_z_z_xxyz_1,   \
                             ta1_z_z_xxzz_1,   \
                             ta1_z_z_xyyz_1,   \
                             ta1_z_z_xyzz_1,   \
                             ta1_z_z_xzzz_1,   \
                             ta1_z_z_yyyy_1,   \
                             ta1_z_z_yyyz_1,   \
                             ta1_z_z_yyzz_1,   \
                             ta1_z_z_yzzz_1,   \
                             ta1_z_z_zzzz_1,   \
                             ta2_xz_x_xxxx_0,  \
                             ta2_xz_x_xxxx_1,  \
                             ta2_xz_x_xxxy_0,  \
                             ta2_xz_x_xxxy_1,  \
                             ta2_xz_x_xxyy_0,  \
                             ta2_xz_x_xxyy_1,  \
                             ta2_xz_x_xyyy_0,  \
                             ta2_xz_x_xyyy_1,  \
                             ta2_xz_xz_xxxx_0, \
                             ta2_xz_xz_xxxy_0, \
                             ta2_xz_xz_xxxz_0, \
                             ta2_xz_xz_xxyy_0, \
                             ta2_xz_xz_xxyz_0, \
                             ta2_xz_xz_xxzz_0, \
                             ta2_xz_xz_xyyy_0, \
                             ta2_xz_xz_xyyz_0, \
                             ta2_xz_xz_xyzz_0, \
                             ta2_xz_xz_xzzz_0, \
                             ta2_xz_xz_yyyy_0, \
                             ta2_xz_xz_yyyz_0, \
                             ta2_xz_xz_yyzz_0, \
                             ta2_xz_xz_yzzz_0, \
                             ta2_xz_xz_zzzz_0, \
                             ta2_xz_z_xxxz_0,  \
                             ta2_xz_z_xxxz_1,  \
                             ta2_xz_z_xxyz_0,  \
                             ta2_xz_z_xxyz_1,  \
                             ta2_xz_z_xxz_0,   \
                             ta2_xz_z_xxz_1,   \
                             ta2_xz_z_xxzz_0,  \
                             ta2_xz_z_xxzz_1,  \
                             ta2_xz_z_xyyz_0,  \
                             ta2_xz_z_xyyz_1,  \
                             ta2_xz_z_xyz_0,   \
                             ta2_xz_z_xyz_1,   \
                             ta2_xz_z_xyzz_0,  \
                             ta2_xz_z_xyzz_1,  \
                             ta2_xz_z_xzz_0,   \
                             ta2_xz_z_xzz_1,   \
                             ta2_xz_z_xzzz_0,  \
                             ta2_xz_z_xzzz_1,  \
                             ta2_xz_z_yyyy_0,  \
                             ta2_xz_z_yyyy_1,  \
                             ta2_xz_z_yyyz_0,  \
                             ta2_xz_z_yyyz_1,  \
                             ta2_xz_z_yyz_0,   \
                             ta2_xz_z_yyz_1,   \
                             ta2_xz_z_yyzz_0,  \
                             ta2_xz_z_yyzz_1,  \
                             ta2_xz_z_yzz_0,   \
                             ta2_xz_z_yzz_1,   \
                             ta2_xz_z_yzzz_0,  \
                             ta2_xz_z_yzzz_1,  \
                             ta2_xz_z_zzz_0,   \
                             ta2_xz_z_zzz_1,   \
                             ta2_xz_z_zzzz_0,  \
                             ta2_xz_z_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xz_xxxx_0[i] = ta1_x_x_xxxx_1[i] + ta2_xz_x_xxxx_0[i] * pa_z[i] - ta2_xz_x_xxxx_1[i] * pc_z[i];

        ta2_xz_xz_xxxy_0[i] = ta1_x_x_xxxy_1[i] + ta2_xz_x_xxxy_0[i] * pa_z[i] - ta2_xz_x_xxxy_1[i] * pc_z[i];

        ta2_xz_xz_xxxz_0[i] = 3.0 * ta2_xz_z_xxz_0[i] * fe_0 - 3.0 * ta2_xz_z_xxz_1[i] * fe_0 + ta1_z_z_xxxz_1[i] + ta2_xz_z_xxxz_0[i] * pa_x[i] -
                              ta2_xz_z_xxxz_1[i] * pc_x[i];

        ta2_xz_xz_xxyy_0[i] = ta1_x_x_xxyy_1[i] + ta2_xz_x_xxyy_0[i] * pa_z[i] - ta2_xz_x_xxyy_1[i] * pc_z[i];

        ta2_xz_xz_xxyz_0[i] = 2.0 * ta2_xz_z_xyz_0[i] * fe_0 - 2.0 * ta2_xz_z_xyz_1[i] * fe_0 + ta1_z_z_xxyz_1[i] + ta2_xz_z_xxyz_0[i] * pa_x[i] -
                              ta2_xz_z_xxyz_1[i] * pc_x[i];

        ta2_xz_xz_xxzz_0[i] = 2.0 * ta2_xz_z_xzz_0[i] * fe_0 - 2.0 * ta2_xz_z_xzz_1[i] * fe_0 + ta1_z_z_xxzz_1[i] + ta2_xz_z_xxzz_0[i] * pa_x[i] -
                              ta2_xz_z_xxzz_1[i] * pc_x[i];

        ta2_xz_xz_xyyy_0[i] = ta1_x_x_xyyy_1[i] + ta2_xz_x_xyyy_0[i] * pa_z[i] - ta2_xz_x_xyyy_1[i] * pc_z[i];

        ta2_xz_xz_xyyz_0[i] =
            ta2_xz_z_yyz_0[i] * fe_0 - ta2_xz_z_yyz_1[i] * fe_0 + ta1_z_z_xyyz_1[i] + ta2_xz_z_xyyz_0[i] * pa_x[i] - ta2_xz_z_xyyz_1[i] * pc_x[i];

        ta2_xz_xz_xyzz_0[i] =
            ta2_xz_z_yzz_0[i] * fe_0 - ta2_xz_z_yzz_1[i] * fe_0 + ta1_z_z_xyzz_1[i] + ta2_xz_z_xyzz_0[i] * pa_x[i] - ta2_xz_z_xyzz_1[i] * pc_x[i];

        ta2_xz_xz_xzzz_0[i] =
            ta2_xz_z_zzz_0[i] * fe_0 - ta2_xz_z_zzz_1[i] * fe_0 + ta1_z_z_xzzz_1[i] + ta2_xz_z_xzzz_0[i] * pa_x[i] - ta2_xz_z_xzzz_1[i] * pc_x[i];

        ta2_xz_xz_yyyy_0[i] = ta1_z_z_yyyy_1[i] + ta2_xz_z_yyyy_0[i] * pa_x[i] - ta2_xz_z_yyyy_1[i] * pc_x[i];

        ta2_xz_xz_yyyz_0[i] = ta1_z_z_yyyz_1[i] + ta2_xz_z_yyyz_0[i] * pa_x[i] - ta2_xz_z_yyyz_1[i] * pc_x[i];

        ta2_xz_xz_yyzz_0[i] = ta1_z_z_yyzz_1[i] + ta2_xz_z_yyzz_0[i] * pa_x[i] - ta2_xz_z_yyzz_1[i] * pc_x[i];

        ta2_xz_xz_yzzz_0[i] = ta1_z_z_yzzz_1[i] + ta2_xz_z_yzzz_0[i] * pa_x[i] - ta2_xz_z_yzzz_1[i] * pc_x[i];

        ta2_xz_xz_zzzz_0[i] = ta1_z_z_zzzz_1[i] + ta2_xz_z_zzzz_0[i] * pa_x[i] - ta2_xz_z_zzzz_1[i] * pc_x[i];
    }

    // Set up 225-240 components of targeted buffer : DG

    auto ta2_xz_yy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 225);

    auto ta2_xz_yy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 226);

    auto ta2_xz_yy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 227);

    auto ta2_xz_yy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 228);

    auto ta2_xz_yy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 229);

    auto ta2_xz_yy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 230);

    auto ta2_xz_yy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 231);

    auto ta2_xz_yy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 232);

    auto ta2_xz_yy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 233);

    auto ta2_xz_yy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 234);

    auto ta2_xz_yy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 235);

    auto ta2_xz_yy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 236);

    auto ta2_xz_yy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 237);

    auto ta2_xz_yy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 238);

    auto ta2_xz_yy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 239);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta2_xz_0_xxxx_0,  \
                             ta2_xz_0_xxxx_1,  \
                             ta2_xz_0_xxxy_0,  \
                             ta2_xz_0_xxxy_1,  \
                             ta2_xz_0_xxxz_0,  \
                             ta2_xz_0_xxxz_1,  \
                             ta2_xz_0_xxyy_0,  \
                             ta2_xz_0_xxyy_1,  \
                             ta2_xz_0_xxyz_0,  \
                             ta2_xz_0_xxyz_1,  \
                             ta2_xz_0_xxzz_0,  \
                             ta2_xz_0_xxzz_1,  \
                             ta2_xz_0_xyyy_0,  \
                             ta2_xz_0_xyyy_1,  \
                             ta2_xz_0_xyyz_0,  \
                             ta2_xz_0_xyyz_1,  \
                             ta2_xz_0_xyzz_0,  \
                             ta2_xz_0_xyzz_1,  \
                             ta2_xz_0_xzzz_0,  \
                             ta2_xz_0_xzzz_1,  \
                             ta2_xz_0_yyyy_0,  \
                             ta2_xz_0_yyyy_1,  \
                             ta2_xz_0_yyyz_0,  \
                             ta2_xz_0_yyyz_1,  \
                             ta2_xz_0_yyzz_0,  \
                             ta2_xz_0_yyzz_1,  \
                             ta2_xz_0_yzzz_0,  \
                             ta2_xz_0_yzzz_1,  \
                             ta2_xz_0_zzzz_0,  \
                             ta2_xz_0_zzzz_1,  \
                             ta2_xz_y_xxx_0,   \
                             ta2_xz_y_xxx_1,   \
                             ta2_xz_y_xxxx_0,  \
                             ta2_xz_y_xxxx_1,  \
                             ta2_xz_y_xxxy_0,  \
                             ta2_xz_y_xxxy_1,  \
                             ta2_xz_y_xxxz_0,  \
                             ta2_xz_y_xxxz_1,  \
                             ta2_xz_y_xxy_0,   \
                             ta2_xz_y_xxy_1,   \
                             ta2_xz_y_xxyy_0,  \
                             ta2_xz_y_xxyy_1,  \
                             ta2_xz_y_xxyz_0,  \
                             ta2_xz_y_xxyz_1,  \
                             ta2_xz_y_xxz_0,   \
                             ta2_xz_y_xxz_1,   \
                             ta2_xz_y_xxzz_0,  \
                             ta2_xz_y_xxzz_1,  \
                             ta2_xz_y_xyy_0,   \
                             ta2_xz_y_xyy_1,   \
                             ta2_xz_y_xyyy_0,  \
                             ta2_xz_y_xyyy_1,  \
                             ta2_xz_y_xyyz_0,  \
                             ta2_xz_y_xyyz_1,  \
                             ta2_xz_y_xyz_0,   \
                             ta2_xz_y_xyz_1,   \
                             ta2_xz_y_xyzz_0,  \
                             ta2_xz_y_xyzz_1,  \
                             ta2_xz_y_xzz_0,   \
                             ta2_xz_y_xzz_1,   \
                             ta2_xz_y_xzzz_0,  \
                             ta2_xz_y_xzzz_1,  \
                             ta2_xz_y_yyy_0,   \
                             ta2_xz_y_yyy_1,   \
                             ta2_xz_y_yyyy_0,  \
                             ta2_xz_y_yyyy_1,  \
                             ta2_xz_y_yyyz_0,  \
                             ta2_xz_y_yyyz_1,  \
                             ta2_xz_y_yyz_0,   \
                             ta2_xz_y_yyz_1,   \
                             ta2_xz_y_yyzz_0,  \
                             ta2_xz_y_yyzz_1,  \
                             ta2_xz_y_yzz_0,   \
                             ta2_xz_y_yzz_1,   \
                             ta2_xz_y_yzzz_0,  \
                             ta2_xz_y_yzzz_1,  \
                             ta2_xz_y_zzz_0,   \
                             ta2_xz_y_zzz_1,   \
                             ta2_xz_y_zzzz_0,  \
                             ta2_xz_y_zzzz_1,  \
                             ta2_xz_yy_xxxx_0, \
                             ta2_xz_yy_xxxy_0, \
                             ta2_xz_yy_xxxz_0, \
                             ta2_xz_yy_xxyy_0, \
                             ta2_xz_yy_xxyz_0, \
                             ta2_xz_yy_xxzz_0, \
                             ta2_xz_yy_xyyy_0, \
                             ta2_xz_yy_xyyz_0, \
                             ta2_xz_yy_xyzz_0, \
                             ta2_xz_yy_xzzz_0, \
                             ta2_xz_yy_yyyy_0, \
                             ta2_xz_yy_yyyz_0, \
                             ta2_xz_yy_yyzz_0, \
                             ta2_xz_yy_yzzz_0, \
                             ta2_xz_yy_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yy_xxxx_0[i] = ta2_xz_0_xxxx_0[i] * fe_0 - ta2_xz_0_xxxx_1[i] * fe_0 + ta2_xz_y_xxxx_0[i] * pa_y[i] - ta2_xz_y_xxxx_1[i] * pc_y[i];

        ta2_xz_yy_xxxy_0[i] = ta2_xz_0_xxxy_0[i] * fe_0 - ta2_xz_0_xxxy_1[i] * fe_0 + ta2_xz_y_xxx_0[i] * fe_0 - ta2_xz_y_xxx_1[i] * fe_0 +
                              ta2_xz_y_xxxy_0[i] * pa_y[i] - ta2_xz_y_xxxy_1[i] * pc_y[i];

        ta2_xz_yy_xxxz_0[i] = ta2_xz_0_xxxz_0[i] * fe_0 - ta2_xz_0_xxxz_1[i] * fe_0 + ta2_xz_y_xxxz_0[i] * pa_y[i] - ta2_xz_y_xxxz_1[i] * pc_y[i];

        ta2_xz_yy_xxyy_0[i] = ta2_xz_0_xxyy_0[i] * fe_0 - ta2_xz_0_xxyy_1[i] * fe_0 + 2.0 * ta2_xz_y_xxy_0[i] * fe_0 -
                              2.0 * ta2_xz_y_xxy_1[i] * fe_0 + ta2_xz_y_xxyy_0[i] * pa_y[i] - ta2_xz_y_xxyy_1[i] * pc_y[i];

        ta2_xz_yy_xxyz_0[i] = ta2_xz_0_xxyz_0[i] * fe_0 - ta2_xz_0_xxyz_1[i] * fe_0 + ta2_xz_y_xxz_0[i] * fe_0 - ta2_xz_y_xxz_1[i] * fe_0 +
                              ta2_xz_y_xxyz_0[i] * pa_y[i] - ta2_xz_y_xxyz_1[i] * pc_y[i];

        ta2_xz_yy_xxzz_0[i] = ta2_xz_0_xxzz_0[i] * fe_0 - ta2_xz_0_xxzz_1[i] * fe_0 + ta2_xz_y_xxzz_0[i] * pa_y[i] - ta2_xz_y_xxzz_1[i] * pc_y[i];

        ta2_xz_yy_xyyy_0[i] = ta2_xz_0_xyyy_0[i] * fe_0 - ta2_xz_0_xyyy_1[i] * fe_0 + 3.0 * ta2_xz_y_xyy_0[i] * fe_0 -
                              3.0 * ta2_xz_y_xyy_1[i] * fe_0 + ta2_xz_y_xyyy_0[i] * pa_y[i] - ta2_xz_y_xyyy_1[i] * pc_y[i];

        ta2_xz_yy_xyyz_0[i] = ta2_xz_0_xyyz_0[i] * fe_0 - ta2_xz_0_xyyz_1[i] * fe_0 + 2.0 * ta2_xz_y_xyz_0[i] * fe_0 -
                              2.0 * ta2_xz_y_xyz_1[i] * fe_0 + ta2_xz_y_xyyz_0[i] * pa_y[i] - ta2_xz_y_xyyz_1[i] * pc_y[i];

        ta2_xz_yy_xyzz_0[i] = ta2_xz_0_xyzz_0[i] * fe_0 - ta2_xz_0_xyzz_1[i] * fe_0 + ta2_xz_y_xzz_0[i] * fe_0 - ta2_xz_y_xzz_1[i] * fe_0 +
                              ta2_xz_y_xyzz_0[i] * pa_y[i] - ta2_xz_y_xyzz_1[i] * pc_y[i];

        ta2_xz_yy_xzzz_0[i] = ta2_xz_0_xzzz_0[i] * fe_0 - ta2_xz_0_xzzz_1[i] * fe_0 + ta2_xz_y_xzzz_0[i] * pa_y[i] - ta2_xz_y_xzzz_1[i] * pc_y[i];

        ta2_xz_yy_yyyy_0[i] = ta2_xz_0_yyyy_0[i] * fe_0 - ta2_xz_0_yyyy_1[i] * fe_0 + 4.0 * ta2_xz_y_yyy_0[i] * fe_0 -
                              4.0 * ta2_xz_y_yyy_1[i] * fe_0 + ta2_xz_y_yyyy_0[i] * pa_y[i] - ta2_xz_y_yyyy_1[i] * pc_y[i];

        ta2_xz_yy_yyyz_0[i] = ta2_xz_0_yyyz_0[i] * fe_0 - ta2_xz_0_yyyz_1[i] * fe_0 + 3.0 * ta2_xz_y_yyz_0[i] * fe_0 -
                              3.0 * ta2_xz_y_yyz_1[i] * fe_0 + ta2_xz_y_yyyz_0[i] * pa_y[i] - ta2_xz_y_yyyz_1[i] * pc_y[i];

        ta2_xz_yy_yyzz_0[i] = ta2_xz_0_yyzz_0[i] * fe_0 - ta2_xz_0_yyzz_1[i] * fe_0 + 2.0 * ta2_xz_y_yzz_0[i] * fe_0 -
                              2.0 * ta2_xz_y_yzz_1[i] * fe_0 + ta2_xz_y_yyzz_0[i] * pa_y[i] - ta2_xz_y_yyzz_1[i] * pc_y[i];

        ta2_xz_yy_yzzz_0[i] = ta2_xz_0_yzzz_0[i] * fe_0 - ta2_xz_0_yzzz_1[i] * fe_0 + ta2_xz_y_zzz_0[i] * fe_0 - ta2_xz_y_zzz_1[i] * fe_0 +
                              ta2_xz_y_yzzz_0[i] * pa_y[i] - ta2_xz_y_yzzz_1[i] * pc_y[i];

        ta2_xz_yy_zzzz_0[i] = ta2_xz_0_zzzz_0[i] * fe_0 - ta2_xz_0_zzzz_1[i] * fe_0 + ta2_xz_y_zzzz_0[i] * pa_y[i] - ta2_xz_y_zzzz_1[i] * pc_y[i];
    }

    // Set up 240-255 components of targeted buffer : DG

    auto ta2_xz_yz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 240);

    auto ta2_xz_yz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 241);

    auto ta2_xz_yz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 242);

    auto ta2_xz_yz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 243);

    auto ta2_xz_yz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 244);

    auto ta2_xz_yz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 245);

    auto ta2_xz_yz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 246);

    auto ta2_xz_yz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 247);

    auto ta2_xz_yz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 248);

    auto ta2_xz_yz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 249);

    auto ta2_xz_yz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 250);

    auto ta2_xz_yz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 251);

    auto ta2_xz_yz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 252);

    auto ta2_xz_yz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 253);

    auto ta2_xz_yz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 254);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_x_y_xxxy_1,   \
                             ta1_x_y_xxyy_1,   \
                             ta1_x_y_xyyy_1,   \
                             ta1_x_y_yyyy_1,   \
                             ta2_xz_y_xxxy_0,  \
                             ta2_xz_y_xxxy_1,  \
                             ta2_xz_y_xxyy_0,  \
                             ta2_xz_y_xxyy_1,  \
                             ta2_xz_y_xyyy_0,  \
                             ta2_xz_y_xyyy_1,  \
                             ta2_xz_y_yyyy_0,  \
                             ta2_xz_y_yyyy_1,  \
                             ta2_xz_yz_xxxx_0, \
                             ta2_xz_yz_xxxy_0, \
                             ta2_xz_yz_xxxz_0, \
                             ta2_xz_yz_xxyy_0, \
                             ta2_xz_yz_xxyz_0, \
                             ta2_xz_yz_xxzz_0, \
                             ta2_xz_yz_xyyy_0, \
                             ta2_xz_yz_xyyz_0, \
                             ta2_xz_yz_xyzz_0, \
                             ta2_xz_yz_xzzz_0, \
                             ta2_xz_yz_yyyy_0, \
                             ta2_xz_yz_yyyz_0, \
                             ta2_xz_yz_yyzz_0, \
                             ta2_xz_yz_yzzz_0, \
                             ta2_xz_yz_zzzz_0, \
                             ta2_xz_z_xxxx_0,  \
                             ta2_xz_z_xxxx_1,  \
                             ta2_xz_z_xxxz_0,  \
                             ta2_xz_z_xxxz_1,  \
                             ta2_xz_z_xxyz_0,  \
                             ta2_xz_z_xxyz_1,  \
                             ta2_xz_z_xxz_0,   \
                             ta2_xz_z_xxz_1,   \
                             ta2_xz_z_xxzz_0,  \
                             ta2_xz_z_xxzz_1,  \
                             ta2_xz_z_xyyz_0,  \
                             ta2_xz_z_xyyz_1,  \
                             ta2_xz_z_xyz_0,   \
                             ta2_xz_z_xyz_1,   \
                             ta2_xz_z_xyzz_0,  \
                             ta2_xz_z_xyzz_1,  \
                             ta2_xz_z_xzz_0,   \
                             ta2_xz_z_xzz_1,   \
                             ta2_xz_z_xzzz_0,  \
                             ta2_xz_z_xzzz_1,  \
                             ta2_xz_z_yyyz_0,  \
                             ta2_xz_z_yyyz_1,  \
                             ta2_xz_z_yyz_0,   \
                             ta2_xz_z_yyz_1,   \
                             ta2_xz_z_yyzz_0,  \
                             ta2_xz_z_yyzz_1,  \
                             ta2_xz_z_yzz_0,   \
                             ta2_xz_z_yzz_1,   \
                             ta2_xz_z_yzzz_0,  \
                             ta2_xz_z_yzzz_1,  \
                             ta2_xz_z_zzz_0,   \
                             ta2_xz_z_zzz_1,   \
                             ta2_xz_z_zzzz_0,  \
                             ta2_xz_z_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yz_xxxx_0[i] = ta2_xz_z_xxxx_0[i] * pa_y[i] - ta2_xz_z_xxxx_1[i] * pc_y[i];

        ta2_xz_yz_xxxy_0[i] = ta1_x_y_xxxy_1[i] + ta2_xz_y_xxxy_0[i] * pa_z[i] - ta2_xz_y_xxxy_1[i] * pc_z[i];

        ta2_xz_yz_xxxz_0[i] = ta2_xz_z_xxxz_0[i] * pa_y[i] - ta2_xz_z_xxxz_1[i] * pc_y[i];

        ta2_xz_yz_xxyy_0[i] = ta1_x_y_xxyy_1[i] + ta2_xz_y_xxyy_0[i] * pa_z[i] - ta2_xz_y_xxyy_1[i] * pc_z[i];

        ta2_xz_yz_xxyz_0[i] = ta2_xz_z_xxz_0[i] * fe_0 - ta2_xz_z_xxz_1[i] * fe_0 + ta2_xz_z_xxyz_0[i] * pa_y[i] - ta2_xz_z_xxyz_1[i] * pc_y[i];

        ta2_xz_yz_xxzz_0[i] = ta2_xz_z_xxzz_0[i] * pa_y[i] - ta2_xz_z_xxzz_1[i] * pc_y[i];

        ta2_xz_yz_xyyy_0[i] = ta1_x_y_xyyy_1[i] + ta2_xz_y_xyyy_0[i] * pa_z[i] - ta2_xz_y_xyyy_1[i] * pc_z[i];

        ta2_xz_yz_xyyz_0[i] =
            2.0 * ta2_xz_z_xyz_0[i] * fe_0 - 2.0 * ta2_xz_z_xyz_1[i] * fe_0 + ta2_xz_z_xyyz_0[i] * pa_y[i] - ta2_xz_z_xyyz_1[i] * pc_y[i];

        ta2_xz_yz_xyzz_0[i] = ta2_xz_z_xzz_0[i] * fe_0 - ta2_xz_z_xzz_1[i] * fe_0 + ta2_xz_z_xyzz_0[i] * pa_y[i] - ta2_xz_z_xyzz_1[i] * pc_y[i];

        ta2_xz_yz_xzzz_0[i] = ta2_xz_z_xzzz_0[i] * pa_y[i] - ta2_xz_z_xzzz_1[i] * pc_y[i];

        ta2_xz_yz_yyyy_0[i] = ta1_x_y_yyyy_1[i] + ta2_xz_y_yyyy_0[i] * pa_z[i] - ta2_xz_y_yyyy_1[i] * pc_z[i];

        ta2_xz_yz_yyyz_0[i] =
            3.0 * ta2_xz_z_yyz_0[i] * fe_0 - 3.0 * ta2_xz_z_yyz_1[i] * fe_0 + ta2_xz_z_yyyz_0[i] * pa_y[i] - ta2_xz_z_yyyz_1[i] * pc_y[i];

        ta2_xz_yz_yyzz_0[i] =
            2.0 * ta2_xz_z_yzz_0[i] * fe_0 - 2.0 * ta2_xz_z_yzz_1[i] * fe_0 + ta2_xz_z_yyzz_0[i] * pa_y[i] - ta2_xz_z_yyzz_1[i] * pc_y[i];

        ta2_xz_yz_yzzz_0[i] = ta2_xz_z_zzz_0[i] * fe_0 - ta2_xz_z_zzz_1[i] * fe_0 + ta2_xz_z_yzzz_0[i] * pa_y[i] - ta2_xz_z_yzzz_1[i] * pc_y[i];

        ta2_xz_yz_zzzz_0[i] = ta2_xz_z_zzzz_0[i] * pa_y[i] - ta2_xz_z_zzzz_1[i] * pc_y[i];
    }

    // Set up 255-270 components of targeted buffer : DG

    auto ta2_xz_zz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 255);

    auto ta2_xz_zz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 256);

    auto ta2_xz_zz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 257);

    auto ta2_xz_zz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 258);

    auto ta2_xz_zz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 259);

    auto ta2_xz_zz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 260);

    auto ta2_xz_zz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 261);

    auto ta2_xz_zz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 262);

    auto ta2_xz_zz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 263);

    auto ta2_xz_zz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 264);

    auto ta2_xz_zz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 265);

    auto ta2_xz_zz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 266);

    auto ta2_xz_zz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 267);

    auto ta2_xz_zz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 268);

    auto ta2_xz_zz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 269);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta1_x_z_xxxx_1,   \
                             ta1_x_z_xxxy_1,   \
                             ta1_x_z_xxxz_1,   \
                             ta1_x_z_xxyy_1,   \
                             ta1_x_z_xxyz_1,   \
                             ta1_x_z_xxzz_1,   \
                             ta1_x_z_xyyy_1,   \
                             ta1_x_z_xyyz_1,   \
                             ta1_x_z_xyzz_1,   \
                             ta1_x_z_xzzz_1,   \
                             ta1_x_z_yyyy_1,   \
                             ta1_x_z_yyyz_1,   \
                             ta1_x_z_yyzz_1,   \
                             ta1_x_z_yzzz_1,   \
                             ta1_x_z_zzzz_1,   \
                             ta2_xz_0_xxxx_0,  \
                             ta2_xz_0_xxxx_1,  \
                             ta2_xz_0_xxxy_0,  \
                             ta2_xz_0_xxxy_1,  \
                             ta2_xz_0_xxxz_0,  \
                             ta2_xz_0_xxxz_1,  \
                             ta2_xz_0_xxyy_0,  \
                             ta2_xz_0_xxyy_1,  \
                             ta2_xz_0_xxyz_0,  \
                             ta2_xz_0_xxyz_1,  \
                             ta2_xz_0_xxzz_0,  \
                             ta2_xz_0_xxzz_1,  \
                             ta2_xz_0_xyyy_0,  \
                             ta2_xz_0_xyyy_1,  \
                             ta2_xz_0_xyyz_0,  \
                             ta2_xz_0_xyyz_1,  \
                             ta2_xz_0_xyzz_0,  \
                             ta2_xz_0_xyzz_1,  \
                             ta2_xz_0_xzzz_0,  \
                             ta2_xz_0_xzzz_1,  \
                             ta2_xz_0_yyyy_0,  \
                             ta2_xz_0_yyyy_1,  \
                             ta2_xz_0_yyyz_0,  \
                             ta2_xz_0_yyyz_1,  \
                             ta2_xz_0_yyzz_0,  \
                             ta2_xz_0_yyzz_1,  \
                             ta2_xz_0_yzzz_0,  \
                             ta2_xz_0_yzzz_1,  \
                             ta2_xz_0_zzzz_0,  \
                             ta2_xz_0_zzzz_1,  \
                             ta2_xz_z_xxx_0,   \
                             ta2_xz_z_xxx_1,   \
                             ta2_xz_z_xxxx_0,  \
                             ta2_xz_z_xxxx_1,  \
                             ta2_xz_z_xxxy_0,  \
                             ta2_xz_z_xxxy_1,  \
                             ta2_xz_z_xxxz_0,  \
                             ta2_xz_z_xxxz_1,  \
                             ta2_xz_z_xxy_0,   \
                             ta2_xz_z_xxy_1,   \
                             ta2_xz_z_xxyy_0,  \
                             ta2_xz_z_xxyy_1,  \
                             ta2_xz_z_xxyz_0,  \
                             ta2_xz_z_xxyz_1,  \
                             ta2_xz_z_xxz_0,   \
                             ta2_xz_z_xxz_1,   \
                             ta2_xz_z_xxzz_0,  \
                             ta2_xz_z_xxzz_1,  \
                             ta2_xz_z_xyy_0,   \
                             ta2_xz_z_xyy_1,   \
                             ta2_xz_z_xyyy_0,  \
                             ta2_xz_z_xyyy_1,  \
                             ta2_xz_z_xyyz_0,  \
                             ta2_xz_z_xyyz_1,  \
                             ta2_xz_z_xyz_0,   \
                             ta2_xz_z_xyz_1,   \
                             ta2_xz_z_xyzz_0,  \
                             ta2_xz_z_xyzz_1,  \
                             ta2_xz_z_xzz_0,   \
                             ta2_xz_z_xzz_1,   \
                             ta2_xz_z_xzzz_0,  \
                             ta2_xz_z_xzzz_1,  \
                             ta2_xz_z_yyy_0,   \
                             ta2_xz_z_yyy_1,   \
                             ta2_xz_z_yyyy_0,  \
                             ta2_xz_z_yyyy_1,  \
                             ta2_xz_z_yyyz_0,  \
                             ta2_xz_z_yyyz_1,  \
                             ta2_xz_z_yyz_0,   \
                             ta2_xz_z_yyz_1,   \
                             ta2_xz_z_yyzz_0,  \
                             ta2_xz_z_yyzz_1,  \
                             ta2_xz_z_yzz_0,   \
                             ta2_xz_z_yzz_1,   \
                             ta2_xz_z_yzzz_0,  \
                             ta2_xz_z_yzzz_1,  \
                             ta2_xz_z_zzz_0,   \
                             ta2_xz_z_zzz_1,   \
                             ta2_xz_z_zzzz_0,  \
                             ta2_xz_z_zzzz_1,  \
                             ta2_xz_zz_xxxx_0, \
                             ta2_xz_zz_xxxy_0, \
                             ta2_xz_zz_xxxz_0, \
                             ta2_xz_zz_xxyy_0, \
                             ta2_xz_zz_xxyz_0, \
                             ta2_xz_zz_xxzz_0, \
                             ta2_xz_zz_xyyy_0, \
                             ta2_xz_zz_xyyz_0, \
                             ta2_xz_zz_xyzz_0, \
                             ta2_xz_zz_xzzz_0, \
                             ta2_xz_zz_yyyy_0, \
                             ta2_xz_zz_yyyz_0, \
                             ta2_xz_zz_yyzz_0, \
                             ta2_xz_zz_yzzz_0, \
                             ta2_xz_zz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_zz_xxxx_0[i] =
            ta2_xz_0_xxxx_0[i] * fe_0 - ta2_xz_0_xxxx_1[i] * fe_0 + ta1_x_z_xxxx_1[i] + ta2_xz_z_xxxx_0[i] * pa_z[i] - ta2_xz_z_xxxx_1[i] * pc_z[i];

        ta2_xz_zz_xxxy_0[i] =
            ta2_xz_0_xxxy_0[i] * fe_0 - ta2_xz_0_xxxy_1[i] * fe_0 + ta1_x_z_xxxy_1[i] + ta2_xz_z_xxxy_0[i] * pa_z[i] - ta2_xz_z_xxxy_1[i] * pc_z[i];

        ta2_xz_zz_xxxz_0[i] = ta2_xz_0_xxxz_0[i] * fe_0 - ta2_xz_0_xxxz_1[i] * fe_0 + ta2_xz_z_xxx_0[i] * fe_0 - ta2_xz_z_xxx_1[i] * fe_0 +
                              ta1_x_z_xxxz_1[i] + ta2_xz_z_xxxz_0[i] * pa_z[i] - ta2_xz_z_xxxz_1[i] * pc_z[i];

        ta2_xz_zz_xxyy_0[i] =
            ta2_xz_0_xxyy_0[i] * fe_0 - ta2_xz_0_xxyy_1[i] * fe_0 + ta1_x_z_xxyy_1[i] + ta2_xz_z_xxyy_0[i] * pa_z[i] - ta2_xz_z_xxyy_1[i] * pc_z[i];

        ta2_xz_zz_xxyz_0[i] = ta2_xz_0_xxyz_0[i] * fe_0 - ta2_xz_0_xxyz_1[i] * fe_0 + ta2_xz_z_xxy_0[i] * fe_0 - ta2_xz_z_xxy_1[i] * fe_0 +
                              ta1_x_z_xxyz_1[i] + ta2_xz_z_xxyz_0[i] * pa_z[i] - ta2_xz_z_xxyz_1[i] * pc_z[i];

        ta2_xz_zz_xxzz_0[i] = ta2_xz_0_xxzz_0[i] * fe_0 - ta2_xz_0_xxzz_1[i] * fe_0 + 2.0 * ta2_xz_z_xxz_0[i] * fe_0 -
                              2.0 * ta2_xz_z_xxz_1[i] * fe_0 + ta1_x_z_xxzz_1[i] + ta2_xz_z_xxzz_0[i] * pa_z[i] - ta2_xz_z_xxzz_1[i] * pc_z[i];

        ta2_xz_zz_xyyy_0[i] =
            ta2_xz_0_xyyy_0[i] * fe_0 - ta2_xz_0_xyyy_1[i] * fe_0 + ta1_x_z_xyyy_1[i] + ta2_xz_z_xyyy_0[i] * pa_z[i] - ta2_xz_z_xyyy_1[i] * pc_z[i];

        ta2_xz_zz_xyyz_0[i] = ta2_xz_0_xyyz_0[i] * fe_0 - ta2_xz_0_xyyz_1[i] * fe_0 + ta2_xz_z_xyy_0[i] * fe_0 - ta2_xz_z_xyy_1[i] * fe_0 +
                              ta1_x_z_xyyz_1[i] + ta2_xz_z_xyyz_0[i] * pa_z[i] - ta2_xz_z_xyyz_1[i] * pc_z[i];

        ta2_xz_zz_xyzz_0[i] = ta2_xz_0_xyzz_0[i] * fe_0 - ta2_xz_0_xyzz_1[i] * fe_0 + 2.0 * ta2_xz_z_xyz_0[i] * fe_0 -
                              2.0 * ta2_xz_z_xyz_1[i] * fe_0 + ta1_x_z_xyzz_1[i] + ta2_xz_z_xyzz_0[i] * pa_z[i] - ta2_xz_z_xyzz_1[i] * pc_z[i];

        ta2_xz_zz_xzzz_0[i] = ta2_xz_0_xzzz_0[i] * fe_0 - ta2_xz_0_xzzz_1[i] * fe_0 + 3.0 * ta2_xz_z_xzz_0[i] * fe_0 -
                              3.0 * ta2_xz_z_xzz_1[i] * fe_0 + ta1_x_z_xzzz_1[i] + ta2_xz_z_xzzz_0[i] * pa_z[i] - ta2_xz_z_xzzz_1[i] * pc_z[i];

        ta2_xz_zz_yyyy_0[i] =
            ta2_xz_0_yyyy_0[i] * fe_0 - ta2_xz_0_yyyy_1[i] * fe_0 + ta1_x_z_yyyy_1[i] + ta2_xz_z_yyyy_0[i] * pa_z[i] - ta2_xz_z_yyyy_1[i] * pc_z[i];

        ta2_xz_zz_yyyz_0[i] = ta2_xz_0_yyyz_0[i] * fe_0 - ta2_xz_0_yyyz_1[i] * fe_0 + ta2_xz_z_yyy_0[i] * fe_0 - ta2_xz_z_yyy_1[i] * fe_0 +
                              ta1_x_z_yyyz_1[i] + ta2_xz_z_yyyz_0[i] * pa_z[i] - ta2_xz_z_yyyz_1[i] * pc_z[i];

        ta2_xz_zz_yyzz_0[i] = ta2_xz_0_yyzz_0[i] * fe_0 - ta2_xz_0_yyzz_1[i] * fe_0 + 2.0 * ta2_xz_z_yyz_0[i] * fe_0 -
                              2.0 * ta2_xz_z_yyz_1[i] * fe_0 + ta1_x_z_yyzz_1[i] + ta2_xz_z_yyzz_0[i] * pa_z[i] - ta2_xz_z_yyzz_1[i] * pc_z[i];

        ta2_xz_zz_yzzz_0[i] = ta2_xz_0_yzzz_0[i] * fe_0 - ta2_xz_0_yzzz_1[i] * fe_0 + 3.0 * ta2_xz_z_yzz_0[i] * fe_0 -
                              3.0 * ta2_xz_z_yzz_1[i] * fe_0 + ta1_x_z_yzzz_1[i] + ta2_xz_z_yzzz_0[i] * pa_z[i] - ta2_xz_z_yzzz_1[i] * pc_z[i];

        ta2_xz_zz_zzzz_0[i] = ta2_xz_0_zzzz_0[i] * fe_0 - ta2_xz_0_zzzz_1[i] * fe_0 + 4.0 * ta2_xz_z_zzz_0[i] * fe_0 -
                              4.0 * ta2_xz_z_zzz_1[i] * fe_0 + ta1_x_z_zzzz_1[i] + ta2_xz_z_zzzz_0[i] * pa_z[i] - ta2_xz_z_zzzz_1[i] * pc_z[i];
    }

    // Set up 270-285 components of targeted buffer : DG

    auto ta2_yy_xx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 270);

    auto ta2_yy_xx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 271);

    auto ta2_yy_xx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 272);

    auto ta2_yy_xx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 273);

    auto ta2_yy_xx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 274);

    auto ta2_yy_xx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 275);

    auto ta2_yy_xx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 276);

    auto ta2_yy_xx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 277);

    auto ta2_yy_xx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 278);

    auto ta2_yy_xx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 279);

    auto ta2_yy_xx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 280);

    auto ta2_yy_xx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 281);

    auto ta2_yy_xx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 282);

    auto ta2_yy_xx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 283);

    auto ta2_yy_xx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 284);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta2_yy_0_xxxx_0,  \
                             ta2_yy_0_xxxx_1,  \
                             ta2_yy_0_xxxy_0,  \
                             ta2_yy_0_xxxy_1,  \
                             ta2_yy_0_xxxz_0,  \
                             ta2_yy_0_xxxz_1,  \
                             ta2_yy_0_xxyy_0,  \
                             ta2_yy_0_xxyy_1,  \
                             ta2_yy_0_xxyz_0,  \
                             ta2_yy_0_xxyz_1,  \
                             ta2_yy_0_xxzz_0,  \
                             ta2_yy_0_xxzz_1,  \
                             ta2_yy_0_xyyy_0,  \
                             ta2_yy_0_xyyy_1,  \
                             ta2_yy_0_xyyz_0,  \
                             ta2_yy_0_xyyz_1,  \
                             ta2_yy_0_xyzz_0,  \
                             ta2_yy_0_xyzz_1,  \
                             ta2_yy_0_xzzz_0,  \
                             ta2_yy_0_xzzz_1,  \
                             ta2_yy_0_yyyy_0,  \
                             ta2_yy_0_yyyy_1,  \
                             ta2_yy_0_yyyz_0,  \
                             ta2_yy_0_yyyz_1,  \
                             ta2_yy_0_yyzz_0,  \
                             ta2_yy_0_yyzz_1,  \
                             ta2_yy_0_yzzz_0,  \
                             ta2_yy_0_yzzz_1,  \
                             ta2_yy_0_zzzz_0,  \
                             ta2_yy_0_zzzz_1,  \
                             ta2_yy_x_xxx_0,   \
                             ta2_yy_x_xxx_1,   \
                             ta2_yy_x_xxxx_0,  \
                             ta2_yy_x_xxxx_1,  \
                             ta2_yy_x_xxxy_0,  \
                             ta2_yy_x_xxxy_1,  \
                             ta2_yy_x_xxxz_0,  \
                             ta2_yy_x_xxxz_1,  \
                             ta2_yy_x_xxy_0,   \
                             ta2_yy_x_xxy_1,   \
                             ta2_yy_x_xxyy_0,  \
                             ta2_yy_x_xxyy_1,  \
                             ta2_yy_x_xxyz_0,  \
                             ta2_yy_x_xxyz_1,  \
                             ta2_yy_x_xxz_0,   \
                             ta2_yy_x_xxz_1,   \
                             ta2_yy_x_xxzz_0,  \
                             ta2_yy_x_xxzz_1,  \
                             ta2_yy_x_xyy_0,   \
                             ta2_yy_x_xyy_1,   \
                             ta2_yy_x_xyyy_0,  \
                             ta2_yy_x_xyyy_1,  \
                             ta2_yy_x_xyyz_0,  \
                             ta2_yy_x_xyyz_1,  \
                             ta2_yy_x_xyz_0,   \
                             ta2_yy_x_xyz_1,   \
                             ta2_yy_x_xyzz_0,  \
                             ta2_yy_x_xyzz_1,  \
                             ta2_yy_x_xzz_0,   \
                             ta2_yy_x_xzz_1,   \
                             ta2_yy_x_xzzz_0,  \
                             ta2_yy_x_xzzz_1,  \
                             ta2_yy_x_yyy_0,   \
                             ta2_yy_x_yyy_1,   \
                             ta2_yy_x_yyyy_0,  \
                             ta2_yy_x_yyyy_1,  \
                             ta2_yy_x_yyyz_0,  \
                             ta2_yy_x_yyyz_1,  \
                             ta2_yy_x_yyz_0,   \
                             ta2_yy_x_yyz_1,   \
                             ta2_yy_x_yyzz_0,  \
                             ta2_yy_x_yyzz_1,  \
                             ta2_yy_x_yzz_0,   \
                             ta2_yy_x_yzz_1,   \
                             ta2_yy_x_yzzz_0,  \
                             ta2_yy_x_yzzz_1,  \
                             ta2_yy_x_zzz_0,   \
                             ta2_yy_x_zzz_1,   \
                             ta2_yy_x_zzzz_0,  \
                             ta2_yy_x_zzzz_1,  \
                             ta2_yy_xx_xxxx_0, \
                             ta2_yy_xx_xxxy_0, \
                             ta2_yy_xx_xxxz_0, \
                             ta2_yy_xx_xxyy_0, \
                             ta2_yy_xx_xxyz_0, \
                             ta2_yy_xx_xxzz_0, \
                             ta2_yy_xx_xyyy_0, \
                             ta2_yy_xx_xyyz_0, \
                             ta2_yy_xx_xyzz_0, \
                             ta2_yy_xx_xzzz_0, \
                             ta2_yy_xx_yyyy_0, \
                             ta2_yy_xx_yyyz_0, \
                             ta2_yy_xx_yyzz_0, \
                             ta2_yy_xx_yzzz_0, \
                             ta2_yy_xx_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xx_xxxx_0[i] = ta2_yy_0_xxxx_0[i] * fe_0 - ta2_yy_0_xxxx_1[i] * fe_0 + 4.0 * ta2_yy_x_xxx_0[i] * fe_0 -
                              4.0 * ta2_yy_x_xxx_1[i] * fe_0 + ta2_yy_x_xxxx_0[i] * pa_x[i] - ta2_yy_x_xxxx_1[i] * pc_x[i];

        ta2_yy_xx_xxxy_0[i] = ta2_yy_0_xxxy_0[i] * fe_0 - ta2_yy_0_xxxy_1[i] * fe_0 + 3.0 * ta2_yy_x_xxy_0[i] * fe_0 -
                              3.0 * ta2_yy_x_xxy_1[i] * fe_0 + ta2_yy_x_xxxy_0[i] * pa_x[i] - ta2_yy_x_xxxy_1[i] * pc_x[i];

        ta2_yy_xx_xxxz_0[i] = ta2_yy_0_xxxz_0[i] * fe_0 - ta2_yy_0_xxxz_1[i] * fe_0 + 3.0 * ta2_yy_x_xxz_0[i] * fe_0 -
                              3.0 * ta2_yy_x_xxz_1[i] * fe_0 + ta2_yy_x_xxxz_0[i] * pa_x[i] - ta2_yy_x_xxxz_1[i] * pc_x[i];

        ta2_yy_xx_xxyy_0[i] = ta2_yy_0_xxyy_0[i] * fe_0 - ta2_yy_0_xxyy_1[i] * fe_0 + 2.0 * ta2_yy_x_xyy_0[i] * fe_0 -
                              2.0 * ta2_yy_x_xyy_1[i] * fe_0 + ta2_yy_x_xxyy_0[i] * pa_x[i] - ta2_yy_x_xxyy_1[i] * pc_x[i];

        ta2_yy_xx_xxyz_0[i] = ta2_yy_0_xxyz_0[i] * fe_0 - ta2_yy_0_xxyz_1[i] * fe_0 + 2.0 * ta2_yy_x_xyz_0[i] * fe_0 -
                              2.0 * ta2_yy_x_xyz_1[i] * fe_0 + ta2_yy_x_xxyz_0[i] * pa_x[i] - ta2_yy_x_xxyz_1[i] * pc_x[i];

        ta2_yy_xx_xxzz_0[i] = ta2_yy_0_xxzz_0[i] * fe_0 - ta2_yy_0_xxzz_1[i] * fe_0 + 2.0 * ta2_yy_x_xzz_0[i] * fe_0 -
                              2.0 * ta2_yy_x_xzz_1[i] * fe_0 + ta2_yy_x_xxzz_0[i] * pa_x[i] - ta2_yy_x_xxzz_1[i] * pc_x[i];

        ta2_yy_xx_xyyy_0[i] = ta2_yy_0_xyyy_0[i] * fe_0 - ta2_yy_0_xyyy_1[i] * fe_0 + ta2_yy_x_yyy_0[i] * fe_0 - ta2_yy_x_yyy_1[i] * fe_0 +
                              ta2_yy_x_xyyy_0[i] * pa_x[i] - ta2_yy_x_xyyy_1[i] * pc_x[i];

        ta2_yy_xx_xyyz_0[i] = ta2_yy_0_xyyz_0[i] * fe_0 - ta2_yy_0_xyyz_1[i] * fe_0 + ta2_yy_x_yyz_0[i] * fe_0 - ta2_yy_x_yyz_1[i] * fe_0 +
                              ta2_yy_x_xyyz_0[i] * pa_x[i] - ta2_yy_x_xyyz_1[i] * pc_x[i];

        ta2_yy_xx_xyzz_0[i] = ta2_yy_0_xyzz_0[i] * fe_0 - ta2_yy_0_xyzz_1[i] * fe_0 + ta2_yy_x_yzz_0[i] * fe_0 - ta2_yy_x_yzz_1[i] * fe_0 +
                              ta2_yy_x_xyzz_0[i] * pa_x[i] - ta2_yy_x_xyzz_1[i] * pc_x[i];

        ta2_yy_xx_xzzz_0[i] = ta2_yy_0_xzzz_0[i] * fe_0 - ta2_yy_0_xzzz_1[i] * fe_0 + ta2_yy_x_zzz_0[i] * fe_0 - ta2_yy_x_zzz_1[i] * fe_0 +
                              ta2_yy_x_xzzz_0[i] * pa_x[i] - ta2_yy_x_xzzz_1[i] * pc_x[i];

        ta2_yy_xx_yyyy_0[i] = ta2_yy_0_yyyy_0[i] * fe_0 - ta2_yy_0_yyyy_1[i] * fe_0 + ta2_yy_x_yyyy_0[i] * pa_x[i] - ta2_yy_x_yyyy_1[i] * pc_x[i];

        ta2_yy_xx_yyyz_0[i] = ta2_yy_0_yyyz_0[i] * fe_0 - ta2_yy_0_yyyz_1[i] * fe_0 + ta2_yy_x_yyyz_0[i] * pa_x[i] - ta2_yy_x_yyyz_1[i] * pc_x[i];

        ta2_yy_xx_yyzz_0[i] = ta2_yy_0_yyzz_0[i] * fe_0 - ta2_yy_0_yyzz_1[i] * fe_0 + ta2_yy_x_yyzz_0[i] * pa_x[i] - ta2_yy_x_yyzz_1[i] * pc_x[i];

        ta2_yy_xx_yzzz_0[i] = ta2_yy_0_yzzz_0[i] * fe_0 - ta2_yy_0_yzzz_1[i] * fe_0 + ta2_yy_x_yzzz_0[i] * pa_x[i] - ta2_yy_x_yzzz_1[i] * pc_x[i];

        ta2_yy_xx_zzzz_0[i] = ta2_yy_0_zzzz_0[i] * fe_0 - ta2_yy_0_zzzz_1[i] * fe_0 + ta2_yy_x_zzzz_0[i] * pa_x[i] - ta2_yy_x_zzzz_1[i] * pc_x[i];
    }

    // Set up 285-300 components of targeted buffer : DG

    auto ta2_yy_xy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 285);

    auto ta2_yy_xy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 286);

    auto ta2_yy_xy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 287);

    auto ta2_yy_xy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 288);

    auto ta2_yy_xy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 289);

    auto ta2_yy_xy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 290);

    auto ta2_yy_xy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 291);

    auto ta2_yy_xy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 292);

    auto ta2_yy_xy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 293);

    auto ta2_yy_xy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 294);

    auto ta2_yy_xy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 295);

    auto ta2_yy_xy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 296);

    auto ta2_yy_xy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 297);

    auto ta2_yy_xy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 298);

    auto ta2_yy_xy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 299);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_y_x_xxxx_1,   \
                             ta1_y_x_xxxz_1,   \
                             ta1_y_x_xxzz_1,   \
                             ta1_y_x_xzzz_1,   \
                             ta2_yy_x_xxxx_0,  \
                             ta2_yy_x_xxxx_1,  \
                             ta2_yy_x_xxxz_0,  \
                             ta2_yy_x_xxxz_1,  \
                             ta2_yy_x_xxzz_0,  \
                             ta2_yy_x_xxzz_1,  \
                             ta2_yy_x_xzzz_0,  \
                             ta2_yy_x_xzzz_1,  \
                             ta2_yy_xy_xxxx_0, \
                             ta2_yy_xy_xxxy_0, \
                             ta2_yy_xy_xxxz_0, \
                             ta2_yy_xy_xxyy_0, \
                             ta2_yy_xy_xxyz_0, \
                             ta2_yy_xy_xxzz_0, \
                             ta2_yy_xy_xyyy_0, \
                             ta2_yy_xy_xyyz_0, \
                             ta2_yy_xy_xyzz_0, \
                             ta2_yy_xy_xzzz_0, \
                             ta2_yy_xy_yyyy_0, \
                             ta2_yy_xy_yyyz_0, \
                             ta2_yy_xy_yyzz_0, \
                             ta2_yy_xy_yzzz_0, \
                             ta2_yy_xy_zzzz_0, \
                             ta2_yy_y_xxxy_0,  \
                             ta2_yy_y_xxxy_1,  \
                             ta2_yy_y_xxy_0,   \
                             ta2_yy_y_xxy_1,   \
                             ta2_yy_y_xxyy_0,  \
                             ta2_yy_y_xxyy_1,  \
                             ta2_yy_y_xxyz_0,  \
                             ta2_yy_y_xxyz_1,  \
                             ta2_yy_y_xyy_0,   \
                             ta2_yy_y_xyy_1,   \
                             ta2_yy_y_xyyy_0,  \
                             ta2_yy_y_xyyy_1,  \
                             ta2_yy_y_xyyz_0,  \
                             ta2_yy_y_xyyz_1,  \
                             ta2_yy_y_xyz_0,   \
                             ta2_yy_y_xyz_1,   \
                             ta2_yy_y_xyzz_0,  \
                             ta2_yy_y_xyzz_1,  \
                             ta2_yy_y_yyy_0,   \
                             ta2_yy_y_yyy_1,   \
                             ta2_yy_y_yyyy_0,  \
                             ta2_yy_y_yyyy_1,  \
                             ta2_yy_y_yyyz_0,  \
                             ta2_yy_y_yyyz_1,  \
                             ta2_yy_y_yyz_0,   \
                             ta2_yy_y_yyz_1,   \
                             ta2_yy_y_yyzz_0,  \
                             ta2_yy_y_yyzz_1,  \
                             ta2_yy_y_yzz_0,   \
                             ta2_yy_y_yzz_1,   \
                             ta2_yy_y_yzzz_0,  \
                             ta2_yy_y_yzzz_1,  \
                             ta2_yy_y_zzzz_0,  \
                             ta2_yy_y_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xy_xxxx_0[i] = 2.0 * ta1_y_x_xxxx_1[i] + ta2_yy_x_xxxx_0[i] * pa_y[i] - ta2_yy_x_xxxx_1[i] * pc_y[i];

        ta2_yy_xy_xxxy_0[i] =
            3.0 * ta2_yy_y_xxy_0[i] * fe_0 - 3.0 * ta2_yy_y_xxy_1[i] * fe_0 + ta2_yy_y_xxxy_0[i] * pa_x[i] - ta2_yy_y_xxxy_1[i] * pc_x[i];

        ta2_yy_xy_xxxz_0[i] = 2.0 * ta1_y_x_xxxz_1[i] + ta2_yy_x_xxxz_0[i] * pa_y[i] - ta2_yy_x_xxxz_1[i] * pc_y[i];

        ta2_yy_xy_xxyy_0[i] =
            2.0 * ta2_yy_y_xyy_0[i] * fe_0 - 2.0 * ta2_yy_y_xyy_1[i] * fe_0 + ta2_yy_y_xxyy_0[i] * pa_x[i] - ta2_yy_y_xxyy_1[i] * pc_x[i];

        ta2_yy_xy_xxyz_0[i] =
            2.0 * ta2_yy_y_xyz_0[i] * fe_0 - 2.0 * ta2_yy_y_xyz_1[i] * fe_0 + ta2_yy_y_xxyz_0[i] * pa_x[i] - ta2_yy_y_xxyz_1[i] * pc_x[i];

        ta2_yy_xy_xxzz_0[i] = 2.0 * ta1_y_x_xxzz_1[i] + ta2_yy_x_xxzz_0[i] * pa_y[i] - ta2_yy_x_xxzz_1[i] * pc_y[i];

        ta2_yy_xy_xyyy_0[i] = ta2_yy_y_yyy_0[i] * fe_0 - ta2_yy_y_yyy_1[i] * fe_0 + ta2_yy_y_xyyy_0[i] * pa_x[i] - ta2_yy_y_xyyy_1[i] * pc_x[i];

        ta2_yy_xy_xyyz_0[i] = ta2_yy_y_yyz_0[i] * fe_0 - ta2_yy_y_yyz_1[i] * fe_0 + ta2_yy_y_xyyz_0[i] * pa_x[i] - ta2_yy_y_xyyz_1[i] * pc_x[i];

        ta2_yy_xy_xyzz_0[i] = ta2_yy_y_yzz_0[i] * fe_0 - ta2_yy_y_yzz_1[i] * fe_0 + ta2_yy_y_xyzz_0[i] * pa_x[i] - ta2_yy_y_xyzz_1[i] * pc_x[i];

        ta2_yy_xy_xzzz_0[i] = 2.0 * ta1_y_x_xzzz_1[i] + ta2_yy_x_xzzz_0[i] * pa_y[i] - ta2_yy_x_xzzz_1[i] * pc_y[i];

        ta2_yy_xy_yyyy_0[i] = ta2_yy_y_yyyy_0[i] * pa_x[i] - ta2_yy_y_yyyy_1[i] * pc_x[i];

        ta2_yy_xy_yyyz_0[i] = ta2_yy_y_yyyz_0[i] * pa_x[i] - ta2_yy_y_yyyz_1[i] * pc_x[i];

        ta2_yy_xy_yyzz_0[i] = ta2_yy_y_yyzz_0[i] * pa_x[i] - ta2_yy_y_yyzz_1[i] * pc_x[i];

        ta2_yy_xy_yzzz_0[i] = ta2_yy_y_yzzz_0[i] * pa_x[i] - ta2_yy_y_yzzz_1[i] * pc_x[i];

        ta2_yy_xy_zzzz_0[i] = ta2_yy_y_zzzz_0[i] * pa_x[i] - ta2_yy_y_zzzz_1[i] * pc_x[i];
    }

    // Set up 300-315 components of targeted buffer : DG

    auto ta2_yy_xz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 300);

    auto ta2_yy_xz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 301);

    auto ta2_yy_xz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 302);

    auto ta2_yy_xz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 303);

    auto ta2_yy_xz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 304);

    auto ta2_yy_xz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 305);

    auto ta2_yy_xz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 306);

    auto ta2_yy_xz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 307);

    auto ta2_yy_xz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 308);

    auto ta2_yy_xz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 309);

    auto ta2_yy_xz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 310);

    auto ta2_yy_xz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 311);

    auto ta2_yy_xz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 312);

    auto ta2_yy_xz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 313);

    auto ta2_yy_xz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 314);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta2_yy_x_xxxx_0,  \
                             ta2_yy_x_xxxx_1,  \
                             ta2_yy_x_xxxy_0,  \
                             ta2_yy_x_xxxy_1,  \
                             ta2_yy_x_xxyy_0,  \
                             ta2_yy_x_xxyy_1,  \
                             ta2_yy_x_xyyy_0,  \
                             ta2_yy_x_xyyy_1,  \
                             ta2_yy_xz_xxxx_0, \
                             ta2_yy_xz_xxxy_0, \
                             ta2_yy_xz_xxxz_0, \
                             ta2_yy_xz_xxyy_0, \
                             ta2_yy_xz_xxyz_0, \
                             ta2_yy_xz_xxzz_0, \
                             ta2_yy_xz_xyyy_0, \
                             ta2_yy_xz_xyyz_0, \
                             ta2_yy_xz_xyzz_0, \
                             ta2_yy_xz_xzzz_0, \
                             ta2_yy_xz_yyyy_0, \
                             ta2_yy_xz_yyyz_0, \
                             ta2_yy_xz_yyzz_0, \
                             ta2_yy_xz_yzzz_0, \
                             ta2_yy_xz_zzzz_0, \
                             ta2_yy_z_xxxz_0,  \
                             ta2_yy_z_xxxz_1,  \
                             ta2_yy_z_xxyz_0,  \
                             ta2_yy_z_xxyz_1,  \
                             ta2_yy_z_xxz_0,   \
                             ta2_yy_z_xxz_1,   \
                             ta2_yy_z_xxzz_0,  \
                             ta2_yy_z_xxzz_1,  \
                             ta2_yy_z_xyyz_0,  \
                             ta2_yy_z_xyyz_1,  \
                             ta2_yy_z_xyz_0,   \
                             ta2_yy_z_xyz_1,   \
                             ta2_yy_z_xyzz_0,  \
                             ta2_yy_z_xyzz_1,  \
                             ta2_yy_z_xzz_0,   \
                             ta2_yy_z_xzz_1,   \
                             ta2_yy_z_xzzz_0,  \
                             ta2_yy_z_xzzz_1,  \
                             ta2_yy_z_yyyy_0,  \
                             ta2_yy_z_yyyy_1,  \
                             ta2_yy_z_yyyz_0,  \
                             ta2_yy_z_yyyz_1,  \
                             ta2_yy_z_yyz_0,   \
                             ta2_yy_z_yyz_1,   \
                             ta2_yy_z_yyzz_0,  \
                             ta2_yy_z_yyzz_1,  \
                             ta2_yy_z_yzz_0,   \
                             ta2_yy_z_yzz_1,   \
                             ta2_yy_z_yzzz_0,  \
                             ta2_yy_z_yzzz_1,  \
                             ta2_yy_z_zzz_0,   \
                             ta2_yy_z_zzz_1,   \
                             ta2_yy_z_zzzz_0,  \
                             ta2_yy_z_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xz_xxxx_0[i] = ta2_yy_x_xxxx_0[i] * pa_z[i] - ta2_yy_x_xxxx_1[i] * pc_z[i];

        ta2_yy_xz_xxxy_0[i] = ta2_yy_x_xxxy_0[i] * pa_z[i] - ta2_yy_x_xxxy_1[i] * pc_z[i];

        ta2_yy_xz_xxxz_0[i] =
            3.0 * ta2_yy_z_xxz_0[i] * fe_0 - 3.0 * ta2_yy_z_xxz_1[i] * fe_0 + ta2_yy_z_xxxz_0[i] * pa_x[i] - ta2_yy_z_xxxz_1[i] * pc_x[i];

        ta2_yy_xz_xxyy_0[i] = ta2_yy_x_xxyy_0[i] * pa_z[i] - ta2_yy_x_xxyy_1[i] * pc_z[i];

        ta2_yy_xz_xxyz_0[i] =
            2.0 * ta2_yy_z_xyz_0[i] * fe_0 - 2.0 * ta2_yy_z_xyz_1[i] * fe_0 + ta2_yy_z_xxyz_0[i] * pa_x[i] - ta2_yy_z_xxyz_1[i] * pc_x[i];

        ta2_yy_xz_xxzz_0[i] =
            2.0 * ta2_yy_z_xzz_0[i] * fe_0 - 2.0 * ta2_yy_z_xzz_1[i] * fe_0 + ta2_yy_z_xxzz_0[i] * pa_x[i] - ta2_yy_z_xxzz_1[i] * pc_x[i];

        ta2_yy_xz_xyyy_0[i] = ta2_yy_x_xyyy_0[i] * pa_z[i] - ta2_yy_x_xyyy_1[i] * pc_z[i];

        ta2_yy_xz_xyyz_0[i] = ta2_yy_z_yyz_0[i] * fe_0 - ta2_yy_z_yyz_1[i] * fe_0 + ta2_yy_z_xyyz_0[i] * pa_x[i] - ta2_yy_z_xyyz_1[i] * pc_x[i];

        ta2_yy_xz_xyzz_0[i] = ta2_yy_z_yzz_0[i] * fe_0 - ta2_yy_z_yzz_1[i] * fe_0 + ta2_yy_z_xyzz_0[i] * pa_x[i] - ta2_yy_z_xyzz_1[i] * pc_x[i];

        ta2_yy_xz_xzzz_0[i] = ta2_yy_z_zzz_0[i] * fe_0 - ta2_yy_z_zzz_1[i] * fe_0 + ta2_yy_z_xzzz_0[i] * pa_x[i] - ta2_yy_z_xzzz_1[i] * pc_x[i];

        ta2_yy_xz_yyyy_0[i] = ta2_yy_z_yyyy_0[i] * pa_x[i] - ta2_yy_z_yyyy_1[i] * pc_x[i];

        ta2_yy_xz_yyyz_0[i] = ta2_yy_z_yyyz_0[i] * pa_x[i] - ta2_yy_z_yyyz_1[i] * pc_x[i];

        ta2_yy_xz_yyzz_0[i] = ta2_yy_z_yyzz_0[i] * pa_x[i] - ta2_yy_z_yyzz_1[i] * pc_x[i];

        ta2_yy_xz_yzzz_0[i] = ta2_yy_z_yzzz_0[i] * pa_x[i] - ta2_yy_z_yzzz_1[i] * pc_x[i];

        ta2_yy_xz_zzzz_0[i] = ta2_yy_z_zzzz_0[i] * pa_x[i] - ta2_yy_z_zzzz_1[i] * pc_x[i];
    }

    // Set up 315-330 components of targeted buffer : DG

    auto ta2_yy_yy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 315);

    auto ta2_yy_yy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 316);

    auto ta2_yy_yy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 317);

    auto ta2_yy_yy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 318);

    auto ta2_yy_yy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 319);

    auto ta2_yy_yy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 320);

    auto ta2_yy_yy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 321);

    auto ta2_yy_yy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 322);

    auto ta2_yy_yy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 323);

    auto ta2_yy_yy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 324);

    auto ta2_yy_yy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 325);

    auto ta2_yy_yy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 326);

    auto ta2_yy_yy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 327);

    auto ta2_yy_yy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 328);

    auto ta2_yy_yy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 329);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_y_y_xxxx_1,   \
                             ta1_y_y_xxxy_1,   \
                             ta1_y_y_xxxz_1,   \
                             ta1_y_y_xxyy_1,   \
                             ta1_y_y_xxyz_1,   \
                             ta1_y_y_xxzz_1,   \
                             ta1_y_y_xyyy_1,   \
                             ta1_y_y_xyyz_1,   \
                             ta1_y_y_xyzz_1,   \
                             ta1_y_y_xzzz_1,   \
                             ta1_y_y_yyyy_1,   \
                             ta1_y_y_yyyz_1,   \
                             ta1_y_y_yyzz_1,   \
                             ta1_y_y_yzzz_1,   \
                             ta1_y_y_zzzz_1,   \
                             ta2_yy_0_xxxx_0,  \
                             ta2_yy_0_xxxx_1,  \
                             ta2_yy_0_xxxy_0,  \
                             ta2_yy_0_xxxy_1,  \
                             ta2_yy_0_xxxz_0,  \
                             ta2_yy_0_xxxz_1,  \
                             ta2_yy_0_xxyy_0,  \
                             ta2_yy_0_xxyy_1,  \
                             ta2_yy_0_xxyz_0,  \
                             ta2_yy_0_xxyz_1,  \
                             ta2_yy_0_xxzz_0,  \
                             ta2_yy_0_xxzz_1,  \
                             ta2_yy_0_xyyy_0,  \
                             ta2_yy_0_xyyy_1,  \
                             ta2_yy_0_xyyz_0,  \
                             ta2_yy_0_xyyz_1,  \
                             ta2_yy_0_xyzz_0,  \
                             ta2_yy_0_xyzz_1,  \
                             ta2_yy_0_xzzz_0,  \
                             ta2_yy_0_xzzz_1,  \
                             ta2_yy_0_yyyy_0,  \
                             ta2_yy_0_yyyy_1,  \
                             ta2_yy_0_yyyz_0,  \
                             ta2_yy_0_yyyz_1,  \
                             ta2_yy_0_yyzz_0,  \
                             ta2_yy_0_yyzz_1,  \
                             ta2_yy_0_yzzz_0,  \
                             ta2_yy_0_yzzz_1,  \
                             ta2_yy_0_zzzz_0,  \
                             ta2_yy_0_zzzz_1,  \
                             ta2_yy_y_xxx_0,   \
                             ta2_yy_y_xxx_1,   \
                             ta2_yy_y_xxxx_0,  \
                             ta2_yy_y_xxxx_1,  \
                             ta2_yy_y_xxxy_0,  \
                             ta2_yy_y_xxxy_1,  \
                             ta2_yy_y_xxxz_0,  \
                             ta2_yy_y_xxxz_1,  \
                             ta2_yy_y_xxy_0,   \
                             ta2_yy_y_xxy_1,   \
                             ta2_yy_y_xxyy_0,  \
                             ta2_yy_y_xxyy_1,  \
                             ta2_yy_y_xxyz_0,  \
                             ta2_yy_y_xxyz_1,  \
                             ta2_yy_y_xxz_0,   \
                             ta2_yy_y_xxz_1,   \
                             ta2_yy_y_xxzz_0,  \
                             ta2_yy_y_xxzz_1,  \
                             ta2_yy_y_xyy_0,   \
                             ta2_yy_y_xyy_1,   \
                             ta2_yy_y_xyyy_0,  \
                             ta2_yy_y_xyyy_1,  \
                             ta2_yy_y_xyyz_0,  \
                             ta2_yy_y_xyyz_1,  \
                             ta2_yy_y_xyz_0,   \
                             ta2_yy_y_xyz_1,   \
                             ta2_yy_y_xyzz_0,  \
                             ta2_yy_y_xyzz_1,  \
                             ta2_yy_y_xzz_0,   \
                             ta2_yy_y_xzz_1,   \
                             ta2_yy_y_xzzz_0,  \
                             ta2_yy_y_xzzz_1,  \
                             ta2_yy_y_yyy_0,   \
                             ta2_yy_y_yyy_1,   \
                             ta2_yy_y_yyyy_0,  \
                             ta2_yy_y_yyyy_1,  \
                             ta2_yy_y_yyyz_0,  \
                             ta2_yy_y_yyyz_1,  \
                             ta2_yy_y_yyz_0,   \
                             ta2_yy_y_yyz_1,   \
                             ta2_yy_y_yyzz_0,  \
                             ta2_yy_y_yyzz_1,  \
                             ta2_yy_y_yzz_0,   \
                             ta2_yy_y_yzz_1,   \
                             ta2_yy_y_yzzz_0,  \
                             ta2_yy_y_yzzz_1,  \
                             ta2_yy_y_zzz_0,   \
                             ta2_yy_y_zzz_1,   \
                             ta2_yy_y_zzzz_0,  \
                             ta2_yy_y_zzzz_1,  \
                             ta2_yy_yy_xxxx_0, \
                             ta2_yy_yy_xxxy_0, \
                             ta2_yy_yy_xxxz_0, \
                             ta2_yy_yy_xxyy_0, \
                             ta2_yy_yy_xxyz_0, \
                             ta2_yy_yy_xxzz_0, \
                             ta2_yy_yy_xyyy_0, \
                             ta2_yy_yy_xyyz_0, \
                             ta2_yy_yy_xyzz_0, \
                             ta2_yy_yy_xzzz_0, \
                             ta2_yy_yy_yyyy_0, \
                             ta2_yy_yy_yyyz_0, \
                             ta2_yy_yy_yyzz_0, \
                             ta2_yy_yy_yzzz_0, \
                             ta2_yy_yy_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yy_xxxx_0[i] = ta2_yy_0_xxxx_0[i] * fe_0 - ta2_yy_0_xxxx_1[i] * fe_0 + 2.0 * ta1_y_y_xxxx_1[i] + ta2_yy_y_xxxx_0[i] * pa_y[i] -
                              ta2_yy_y_xxxx_1[i] * pc_y[i];

        ta2_yy_yy_xxxy_0[i] = ta2_yy_0_xxxy_0[i] * fe_0 - ta2_yy_0_xxxy_1[i] * fe_0 + ta2_yy_y_xxx_0[i] * fe_0 - ta2_yy_y_xxx_1[i] * fe_0 +
                              2.0 * ta1_y_y_xxxy_1[i] + ta2_yy_y_xxxy_0[i] * pa_y[i] - ta2_yy_y_xxxy_1[i] * pc_y[i];

        ta2_yy_yy_xxxz_0[i] = ta2_yy_0_xxxz_0[i] * fe_0 - ta2_yy_0_xxxz_1[i] * fe_0 + 2.0 * ta1_y_y_xxxz_1[i] + ta2_yy_y_xxxz_0[i] * pa_y[i] -
                              ta2_yy_y_xxxz_1[i] * pc_y[i];

        ta2_yy_yy_xxyy_0[i] = ta2_yy_0_xxyy_0[i] * fe_0 - ta2_yy_0_xxyy_1[i] * fe_0 + 2.0 * ta2_yy_y_xxy_0[i] * fe_0 -
                              2.0 * ta2_yy_y_xxy_1[i] * fe_0 + 2.0 * ta1_y_y_xxyy_1[i] + ta2_yy_y_xxyy_0[i] * pa_y[i] - ta2_yy_y_xxyy_1[i] * pc_y[i];

        ta2_yy_yy_xxyz_0[i] = ta2_yy_0_xxyz_0[i] * fe_0 - ta2_yy_0_xxyz_1[i] * fe_0 + ta2_yy_y_xxz_0[i] * fe_0 - ta2_yy_y_xxz_1[i] * fe_0 +
                              2.0 * ta1_y_y_xxyz_1[i] + ta2_yy_y_xxyz_0[i] * pa_y[i] - ta2_yy_y_xxyz_1[i] * pc_y[i];

        ta2_yy_yy_xxzz_0[i] = ta2_yy_0_xxzz_0[i] * fe_0 - ta2_yy_0_xxzz_1[i] * fe_0 + 2.0 * ta1_y_y_xxzz_1[i] + ta2_yy_y_xxzz_0[i] * pa_y[i] -
                              ta2_yy_y_xxzz_1[i] * pc_y[i];

        ta2_yy_yy_xyyy_0[i] = ta2_yy_0_xyyy_0[i] * fe_0 - ta2_yy_0_xyyy_1[i] * fe_0 + 3.0 * ta2_yy_y_xyy_0[i] * fe_0 -
                              3.0 * ta2_yy_y_xyy_1[i] * fe_0 + 2.0 * ta1_y_y_xyyy_1[i] + ta2_yy_y_xyyy_0[i] * pa_y[i] - ta2_yy_y_xyyy_1[i] * pc_y[i];

        ta2_yy_yy_xyyz_0[i] = ta2_yy_0_xyyz_0[i] * fe_0 - ta2_yy_0_xyyz_1[i] * fe_0 + 2.0 * ta2_yy_y_xyz_0[i] * fe_0 -
                              2.0 * ta2_yy_y_xyz_1[i] * fe_0 + 2.0 * ta1_y_y_xyyz_1[i] + ta2_yy_y_xyyz_0[i] * pa_y[i] - ta2_yy_y_xyyz_1[i] * pc_y[i];

        ta2_yy_yy_xyzz_0[i] = ta2_yy_0_xyzz_0[i] * fe_0 - ta2_yy_0_xyzz_1[i] * fe_0 + ta2_yy_y_xzz_0[i] * fe_0 - ta2_yy_y_xzz_1[i] * fe_0 +
                              2.0 * ta1_y_y_xyzz_1[i] + ta2_yy_y_xyzz_0[i] * pa_y[i] - ta2_yy_y_xyzz_1[i] * pc_y[i];

        ta2_yy_yy_xzzz_0[i] = ta2_yy_0_xzzz_0[i] * fe_0 - ta2_yy_0_xzzz_1[i] * fe_0 + 2.0 * ta1_y_y_xzzz_1[i] + ta2_yy_y_xzzz_0[i] * pa_y[i] -
                              ta2_yy_y_xzzz_1[i] * pc_y[i];

        ta2_yy_yy_yyyy_0[i] = ta2_yy_0_yyyy_0[i] * fe_0 - ta2_yy_0_yyyy_1[i] * fe_0 + 4.0 * ta2_yy_y_yyy_0[i] * fe_0 -
                              4.0 * ta2_yy_y_yyy_1[i] * fe_0 + 2.0 * ta1_y_y_yyyy_1[i] + ta2_yy_y_yyyy_0[i] * pa_y[i] - ta2_yy_y_yyyy_1[i] * pc_y[i];

        ta2_yy_yy_yyyz_0[i] = ta2_yy_0_yyyz_0[i] * fe_0 - ta2_yy_0_yyyz_1[i] * fe_0 + 3.0 * ta2_yy_y_yyz_0[i] * fe_0 -
                              3.0 * ta2_yy_y_yyz_1[i] * fe_0 + 2.0 * ta1_y_y_yyyz_1[i] + ta2_yy_y_yyyz_0[i] * pa_y[i] - ta2_yy_y_yyyz_1[i] * pc_y[i];

        ta2_yy_yy_yyzz_0[i] = ta2_yy_0_yyzz_0[i] * fe_0 - ta2_yy_0_yyzz_1[i] * fe_0 + 2.0 * ta2_yy_y_yzz_0[i] * fe_0 -
                              2.0 * ta2_yy_y_yzz_1[i] * fe_0 + 2.0 * ta1_y_y_yyzz_1[i] + ta2_yy_y_yyzz_0[i] * pa_y[i] - ta2_yy_y_yyzz_1[i] * pc_y[i];

        ta2_yy_yy_yzzz_0[i] = ta2_yy_0_yzzz_0[i] * fe_0 - ta2_yy_0_yzzz_1[i] * fe_0 + ta2_yy_y_zzz_0[i] * fe_0 - ta2_yy_y_zzz_1[i] * fe_0 +
                              2.0 * ta1_y_y_yzzz_1[i] + ta2_yy_y_yzzz_0[i] * pa_y[i] - ta2_yy_y_yzzz_1[i] * pc_y[i];

        ta2_yy_yy_zzzz_0[i] = ta2_yy_0_zzzz_0[i] * fe_0 - ta2_yy_0_zzzz_1[i] * fe_0 + 2.0 * ta1_y_y_zzzz_1[i] + ta2_yy_y_zzzz_0[i] * pa_y[i] -
                              ta2_yy_y_zzzz_1[i] * pc_y[i];
    }

    // Set up 330-345 components of targeted buffer : DG

    auto ta2_yy_yz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 330);

    auto ta2_yy_yz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 331);

    auto ta2_yy_yz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 332);

    auto ta2_yy_yz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 333);

    auto ta2_yy_yz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 334);

    auto ta2_yy_yz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 335);

    auto ta2_yy_yz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 336);

    auto ta2_yy_yz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 337);

    auto ta2_yy_yz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 338);

    auto ta2_yy_yz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 339);

    auto ta2_yy_yz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 340);

    auto ta2_yy_yz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 341);

    auto ta2_yy_yz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 342);

    auto ta2_yy_yz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 343);

    auto ta2_yy_yz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 344);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_y_z_xxxz_1,   \
                             ta1_y_z_xxzz_1,   \
                             ta1_y_z_xzzz_1,   \
                             ta1_y_z_zzzz_1,   \
                             ta2_yy_y_xxxx_0,  \
                             ta2_yy_y_xxxx_1,  \
                             ta2_yy_y_xxxy_0,  \
                             ta2_yy_y_xxxy_1,  \
                             ta2_yy_y_xxy_0,   \
                             ta2_yy_y_xxy_1,   \
                             ta2_yy_y_xxyy_0,  \
                             ta2_yy_y_xxyy_1,  \
                             ta2_yy_y_xxyz_0,  \
                             ta2_yy_y_xxyz_1,  \
                             ta2_yy_y_xyy_0,   \
                             ta2_yy_y_xyy_1,   \
                             ta2_yy_y_xyyy_0,  \
                             ta2_yy_y_xyyy_1,  \
                             ta2_yy_y_xyyz_0,  \
                             ta2_yy_y_xyyz_1,  \
                             ta2_yy_y_xyz_0,   \
                             ta2_yy_y_xyz_1,   \
                             ta2_yy_y_xyzz_0,  \
                             ta2_yy_y_xyzz_1,  \
                             ta2_yy_y_yyy_0,   \
                             ta2_yy_y_yyy_1,   \
                             ta2_yy_y_yyyy_0,  \
                             ta2_yy_y_yyyy_1,  \
                             ta2_yy_y_yyyz_0,  \
                             ta2_yy_y_yyyz_1,  \
                             ta2_yy_y_yyz_0,   \
                             ta2_yy_y_yyz_1,   \
                             ta2_yy_y_yyzz_0,  \
                             ta2_yy_y_yyzz_1,  \
                             ta2_yy_y_yzz_0,   \
                             ta2_yy_y_yzz_1,   \
                             ta2_yy_y_yzzz_0,  \
                             ta2_yy_y_yzzz_1,  \
                             ta2_yy_yz_xxxx_0, \
                             ta2_yy_yz_xxxy_0, \
                             ta2_yy_yz_xxxz_0, \
                             ta2_yy_yz_xxyy_0, \
                             ta2_yy_yz_xxyz_0, \
                             ta2_yy_yz_xxzz_0, \
                             ta2_yy_yz_xyyy_0, \
                             ta2_yy_yz_xyyz_0, \
                             ta2_yy_yz_xyzz_0, \
                             ta2_yy_yz_xzzz_0, \
                             ta2_yy_yz_yyyy_0, \
                             ta2_yy_yz_yyyz_0, \
                             ta2_yy_yz_yyzz_0, \
                             ta2_yy_yz_yzzz_0, \
                             ta2_yy_yz_zzzz_0, \
                             ta2_yy_z_xxxz_0,  \
                             ta2_yy_z_xxxz_1,  \
                             ta2_yy_z_xxzz_0,  \
                             ta2_yy_z_xxzz_1,  \
                             ta2_yy_z_xzzz_0,  \
                             ta2_yy_z_xzzz_1,  \
                             ta2_yy_z_zzzz_0,  \
                             ta2_yy_z_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yz_xxxx_0[i] = ta2_yy_y_xxxx_0[i] * pa_z[i] - ta2_yy_y_xxxx_1[i] * pc_z[i];

        ta2_yy_yz_xxxy_0[i] = ta2_yy_y_xxxy_0[i] * pa_z[i] - ta2_yy_y_xxxy_1[i] * pc_z[i];

        ta2_yy_yz_xxxz_0[i] = 2.0 * ta1_y_z_xxxz_1[i] + ta2_yy_z_xxxz_0[i] * pa_y[i] - ta2_yy_z_xxxz_1[i] * pc_y[i];

        ta2_yy_yz_xxyy_0[i] = ta2_yy_y_xxyy_0[i] * pa_z[i] - ta2_yy_y_xxyy_1[i] * pc_z[i];

        ta2_yy_yz_xxyz_0[i] = ta2_yy_y_xxy_0[i] * fe_0 - ta2_yy_y_xxy_1[i] * fe_0 + ta2_yy_y_xxyz_0[i] * pa_z[i] - ta2_yy_y_xxyz_1[i] * pc_z[i];

        ta2_yy_yz_xxzz_0[i] = 2.0 * ta1_y_z_xxzz_1[i] + ta2_yy_z_xxzz_0[i] * pa_y[i] - ta2_yy_z_xxzz_1[i] * pc_y[i];

        ta2_yy_yz_xyyy_0[i] = ta2_yy_y_xyyy_0[i] * pa_z[i] - ta2_yy_y_xyyy_1[i] * pc_z[i];

        ta2_yy_yz_xyyz_0[i] = ta2_yy_y_xyy_0[i] * fe_0 - ta2_yy_y_xyy_1[i] * fe_0 + ta2_yy_y_xyyz_0[i] * pa_z[i] - ta2_yy_y_xyyz_1[i] * pc_z[i];

        ta2_yy_yz_xyzz_0[i] =
            2.0 * ta2_yy_y_xyz_0[i] * fe_0 - 2.0 * ta2_yy_y_xyz_1[i] * fe_0 + ta2_yy_y_xyzz_0[i] * pa_z[i] - ta2_yy_y_xyzz_1[i] * pc_z[i];

        ta2_yy_yz_xzzz_0[i] = 2.0 * ta1_y_z_xzzz_1[i] + ta2_yy_z_xzzz_0[i] * pa_y[i] - ta2_yy_z_xzzz_1[i] * pc_y[i];

        ta2_yy_yz_yyyy_0[i] = ta2_yy_y_yyyy_0[i] * pa_z[i] - ta2_yy_y_yyyy_1[i] * pc_z[i];

        ta2_yy_yz_yyyz_0[i] = ta2_yy_y_yyy_0[i] * fe_0 - ta2_yy_y_yyy_1[i] * fe_0 + ta2_yy_y_yyyz_0[i] * pa_z[i] - ta2_yy_y_yyyz_1[i] * pc_z[i];

        ta2_yy_yz_yyzz_0[i] =
            2.0 * ta2_yy_y_yyz_0[i] * fe_0 - 2.0 * ta2_yy_y_yyz_1[i] * fe_0 + ta2_yy_y_yyzz_0[i] * pa_z[i] - ta2_yy_y_yyzz_1[i] * pc_z[i];

        ta2_yy_yz_yzzz_0[i] =
            3.0 * ta2_yy_y_yzz_0[i] * fe_0 - 3.0 * ta2_yy_y_yzz_1[i] * fe_0 + ta2_yy_y_yzzz_0[i] * pa_z[i] - ta2_yy_y_yzzz_1[i] * pc_z[i];

        ta2_yy_yz_zzzz_0[i] = 2.0 * ta1_y_z_zzzz_1[i] + ta2_yy_z_zzzz_0[i] * pa_y[i] - ta2_yy_z_zzzz_1[i] * pc_y[i];
    }

    // Set up 345-360 components of targeted buffer : DG

    auto ta2_yy_zz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 345);

    auto ta2_yy_zz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 346);

    auto ta2_yy_zz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 347);

    auto ta2_yy_zz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 348);

    auto ta2_yy_zz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 349);

    auto ta2_yy_zz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 350);

    auto ta2_yy_zz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 351);

    auto ta2_yy_zz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 352);

    auto ta2_yy_zz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 353);

    auto ta2_yy_zz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 354);

    auto ta2_yy_zz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 355);

    auto ta2_yy_zz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 356);

    auto ta2_yy_zz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 357);

    auto ta2_yy_zz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 358);

    auto ta2_yy_zz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 359);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta2_yy_0_xxxx_0,  \
                             ta2_yy_0_xxxx_1,  \
                             ta2_yy_0_xxxy_0,  \
                             ta2_yy_0_xxxy_1,  \
                             ta2_yy_0_xxxz_0,  \
                             ta2_yy_0_xxxz_1,  \
                             ta2_yy_0_xxyy_0,  \
                             ta2_yy_0_xxyy_1,  \
                             ta2_yy_0_xxyz_0,  \
                             ta2_yy_0_xxyz_1,  \
                             ta2_yy_0_xxzz_0,  \
                             ta2_yy_0_xxzz_1,  \
                             ta2_yy_0_xyyy_0,  \
                             ta2_yy_0_xyyy_1,  \
                             ta2_yy_0_xyyz_0,  \
                             ta2_yy_0_xyyz_1,  \
                             ta2_yy_0_xyzz_0,  \
                             ta2_yy_0_xyzz_1,  \
                             ta2_yy_0_xzzz_0,  \
                             ta2_yy_0_xzzz_1,  \
                             ta2_yy_0_yyyy_0,  \
                             ta2_yy_0_yyyy_1,  \
                             ta2_yy_0_yyyz_0,  \
                             ta2_yy_0_yyyz_1,  \
                             ta2_yy_0_yyzz_0,  \
                             ta2_yy_0_yyzz_1,  \
                             ta2_yy_0_yzzz_0,  \
                             ta2_yy_0_yzzz_1,  \
                             ta2_yy_0_zzzz_0,  \
                             ta2_yy_0_zzzz_1,  \
                             ta2_yy_z_xxx_0,   \
                             ta2_yy_z_xxx_1,   \
                             ta2_yy_z_xxxx_0,  \
                             ta2_yy_z_xxxx_1,  \
                             ta2_yy_z_xxxy_0,  \
                             ta2_yy_z_xxxy_1,  \
                             ta2_yy_z_xxxz_0,  \
                             ta2_yy_z_xxxz_1,  \
                             ta2_yy_z_xxy_0,   \
                             ta2_yy_z_xxy_1,   \
                             ta2_yy_z_xxyy_0,  \
                             ta2_yy_z_xxyy_1,  \
                             ta2_yy_z_xxyz_0,  \
                             ta2_yy_z_xxyz_1,  \
                             ta2_yy_z_xxz_0,   \
                             ta2_yy_z_xxz_1,   \
                             ta2_yy_z_xxzz_0,  \
                             ta2_yy_z_xxzz_1,  \
                             ta2_yy_z_xyy_0,   \
                             ta2_yy_z_xyy_1,   \
                             ta2_yy_z_xyyy_0,  \
                             ta2_yy_z_xyyy_1,  \
                             ta2_yy_z_xyyz_0,  \
                             ta2_yy_z_xyyz_1,  \
                             ta2_yy_z_xyz_0,   \
                             ta2_yy_z_xyz_1,   \
                             ta2_yy_z_xyzz_0,  \
                             ta2_yy_z_xyzz_1,  \
                             ta2_yy_z_xzz_0,   \
                             ta2_yy_z_xzz_1,   \
                             ta2_yy_z_xzzz_0,  \
                             ta2_yy_z_xzzz_1,  \
                             ta2_yy_z_yyy_0,   \
                             ta2_yy_z_yyy_1,   \
                             ta2_yy_z_yyyy_0,  \
                             ta2_yy_z_yyyy_1,  \
                             ta2_yy_z_yyyz_0,  \
                             ta2_yy_z_yyyz_1,  \
                             ta2_yy_z_yyz_0,   \
                             ta2_yy_z_yyz_1,   \
                             ta2_yy_z_yyzz_0,  \
                             ta2_yy_z_yyzz_1,  \
                             ta2_yy_z_yzz_0,   \
                             ta2_yy_z_yzz_1,   \
                             ta2_yy_z_yzzz_0,  \
                             ta2_yy_z_yzzz_1,  \
                             ta2_yy_z_zzz_0,   \
                             ta2_yy_z_zzz_1,   \
                             ta2_yy_z_zzzz_0,  \
                             ta2_yy_z_zzzz_1,  \
                             ta2_yy_zz_xxxx_0, \
                             ta2_yy_zz_xxxy_0, \
                             ta2_yy_zz_xxxz_0, \
                             ta2_yy_zz_xxyy_0, \
                             ta2_yy_zz_xxyz_0, \
                             ta2_yy_zz_xxzz_0, \
                             ta2_yy_zz_xyyy_0, \
                             ta2_yy_zz_xyyz_0, \
                             ta2_yy_zz_xyzz_0, \
                             ta2_yy_zz_xzzz_0, \
                             ta2_yy_zz_yyyy_0, \
                             ta2_yy_zz_yyyz_0, \
                             ta2_yy_zz_yyzz_0, \
                             ta2_yy_zz_yzzz_0, \
                             ta2_yy_zz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_zz_xxxx_0[i] = ta2_yy_0_xxxx_0[i] * fe_0 - ta2_yy_0_xxxx_1[i] * fe_0 + ta2_yy_z_xxxx_0[i] * pa_z[i] - ta2_yy_z_xxxx_1[i] * pc_z[i];

        ta2_yy_zz_xxxy_0[i] = ta2_yy_0_xxxy_0[i] * fe_0 - ta2_yy_0_xxxy_1[i] * fe_0 + ta2_yy_z_xxxy_0[i] * pa_z[i] - ta2_yy_z_xxxy_1[i] * pc_z[i];

        ta2_yy_zz_xxxz_0[i] = ta2_yy_0_xxxz_0[i] * fe_0 - ta2_yy_0_xxxz_1[i] * fe_0 + ta2_yy_z_xxx_0[i] * fe_0 - ta2_yy_z_xxx_1[i] * fe_0 +
                              ta2_yy_z_xxxz_0[i] * pa_z[i] - ta2_yy_z_xxxz_1[i] * pc_z[i];

        ta2_yy_zz_xxyy_0[i] = ta2_yy_0_xxyy_0[i] * fe_0 - ta2_yy_0_xxyy_1[i] * fe_0 + ta2_yy_z_xxyy_0[i] * pa_z[i] - ta2_yy_z_xxyy_1[i] * pc_z[i];

        ta2_yy_zz_xxyz_0[i] = ta2_yy_0_xxyz_0[i] * fe_0 - ta2_yy_0_xxyz_1[i] * fe_0 + ta2_yy_z_xxy_0[i] * fe_0 - ta2_yy_z_xxy_1[i] * fe_0 +
                              ta2_yy_z_xxyz_0[i] * pa_z[i] - ta2_yy_z_xxyz_1[i] * pc_z[i];

        ta2_yy_zz_xxzz_0[i] = ta2_yy_0_xxzz_0[i] * fe_0 - ta2_yy_0_xxzz_1[i] * fe_0 + 2.0 * ta2_yy_z_xxz_0[i] * fe_0 -
                              2.0 * ta2_yy_z_xxz_1[i] * fe_0 + ta2_yy_z_xxzz_0[i] * pa_z[i] - ta2_yy_z_xxzz_1[i] * pc_z[i];

        ta2_yy_zz_xyyy_0[i] = ta2_yy_0_xyyy_0[i] * fe_0 - ta2_yy_0_xyyy_1[i] * fe_0 + ta2_yy_z_xyyy_0[i] * pa_z[i] - ta2_yy_z_xyyy_1[i] * pc_z[i];

        ta2_yy_zz_xyyz_0[i] = ta2_yy_0_xyyz_0[i] * fe_0 - ta2_yy_0_xyyz_1[i] * fe_0 + ta2_yy_z_xyy_0[i] * fe_0 - ta2_yy_z_xyy_1[i] * fe_0 +
                              ta2_yy_z_xyyz_0[i] * pa_z[i] - ta2_yy_z_xyyz_1[i] * pc_z[i];

        ta2_yy_zz_xyzz_0[i] = ta2_yy_0_xyzz_0[i] * fe_0 - ta2_yy_0_xyzz_1[i] * fe_0 + 2.0 * ta2_yy_z_xyz_0[i] * fe_0 -
                              2.0 * ta2_yy_z_xyz_1[i] * fe_0 + ta2_yy_z_xyzz_0[i] * pa_z[i] - ta2_yy_z_xyzz_1[i] * pc_z[i];

        ta2_yy_zz_xzzz_0[i] = ta2_yy_0_xzzz_0[i] * fe_0 - ta2_yy_0_xzzz_1[i] * fe_0 + 3.0 * ta2_yy_z_xzz_0[i] * fe_0 -
                              3.0 * ta2_yy_z_xzz_1[i] * fe_0 + ta2_yy_z_xzzz_0[i] * pa_z[i] - ta2_yy_z_xzzz_1[i] * pc_z[i];

        ta2_yy_zz_yyyy_0[i] = ta2_yy_0_yyyy_0[i] * fe_0 - ta2_yy_0_yyyy_1[i] * fe_0 + ta2_yy_z_yyyy_0[i] * pa_z[i] - ta2_yy_z_yyyy_1[i] * pc_z[i];

        ta2_yy_zz_yyyz_0[i] = ta2_yy_0_yyyz_0[i] * fe_0 - ta2_yy_0_yyyz_1[i] * fe_0 + ta2_yy_z_yyy_0[i] * fe_0 - ta2_yy_z_yyy_1[i] * fe_0 +
                              ta2_yy_z_yyyz_0[i] * pa_z[i] - ta2_yy_z_yyyz_1[i] * pc_z[i];

        ta2_yy_zz_yyzz_0[i] = ta2_yy_0_yyzz_0[i] * fe_0 - ta2_yy_0_yyzz_1[i] * fe_0 + 2.0 * ta2_yy_z_yyz_0[i] * fe_0 -
                              2.0 * ta2_yy_z_yyz_1[i] * fe_0 + ta2_yy_z_yyzz_0[i] * pa_z[i] - ta2_yy_z_yyzz_1[i] * pc_z[i];

        ta2_yy_zz_yzzz_0[i] = ta2_yy_0_yzzz_0[i] * fe_0 - ta2_yy_0_yzzz_1[i] * fe_0 + 3.0 * ta2_yy_z_yzz_0[i] * fe_0 -
                              3.0 * ta2_yy_z_yzz_1[i] * fe_0 + ta2_yy_z_yzzz_0[i] * pa_z[i] - ta2_yy_z_yzzz_1[i] * pc_z[i];

        ta2_yy_zz_zzzz_0[i] = ta2_yy_0_zzzz_0[i] * fe_0 - ta2_yy_0_zzzz_1[i] * fe_0 + 4.0 * ta2_yy_z_zzz_0[i] * fe_0 -
                              4.0 * ta2_yy_z_zzz_1[i] * fe_0 + ta2_yy_z_zzzz_0[i] * pa_z[i] - ta2_yy_z_zzzz_1[i] * pc_z[i];
    }

    // Set up 360-375 components of targeted buffer : DG

    auto ta2_yz_xx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 360);

    auto ta2_yz_xx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 361);

    auto ta2_yz_xx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 362);

    auto ta2_yz_xx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 363);

    auto ta2_yz_xx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 364);

    auto ta2_yz_xx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 365);

    auto ta2_yz_xx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 366);

    auto ta2_yz_xx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 367);

    auto ta2_yz_xx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 368);

    auto ta2_yz_xx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 369);

    auto ta2_yz_xx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 370);

    auto ta2_yz_xx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 371);

    auto ta2_yz_xx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 372);

    auto ta2_yz_xx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 373);

    auto ta2_yz_xx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 374);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta2_yz_0_xxxx_0,  \
                             ta2_yz_0_xxxx_1,  \
                             ta2_yz_0_xxxy_0,  \
                             ta2_yz_0_xxxy_1,  \
                             ta2_yz_0_xxxz_0,  \
                             ta2_yz_0_xxxz_1,  \
                             ta2_yz_0_xxyy_0,  \
                             ta2_yz_0_xxyy_1,  \
                             ta2_yz_0_xxyz_0,  \
                             ta2_yz_0_xxyz_1,  \
                             ta2_yz_0_xxzz_0,  \
                             ta2_yz_0_xxzz_1,  \
                             ta2_yz_0_xyyy_0,  \
                             ta2_yz_0_xyyy_1,  \
                             ta2_yz_0_xyyz_0,  \
                             ta2_yz_0_xyyz_1,  \
                             ta2_yz_0_xyzz_0,  \
                             ta2_yz_0_xyzz_1,  \
                             ta2_yz_0_xzzz_0,  \
                             ta2_yz_0_xzzz_1,  \
                             ta2_yz_0_yyyy_0,  \
                             ta2_yz_0_yyyy_1,  \
                             ta2_yz_0_yyyz_0,  \
                             ta2_yz_0_yyyz_1,  \
                             ta2_yz_0_yyzz_0,  \
                             ta2_yz_0_yyzz_1,  \
                             ta2_yz_0_yzzz_0,  \
                             ta2_yz_0_yzzz_1,  \
                             ta2_yz_0_zzzz_0,  \
                             ta2_yz_0_zzzz_1,  \
                             ta2_yz_x_xxx_0,   \
                             ta2_yz_x_xxx_1,   \
                             ta2_yz_x_xxxx_0,  \
                             ta2_yz_x_xxxx_1,  \
                             ta2_yz_x_xxxy_0,  \
                             ta2_yz_x_xxxy_1,  \
                             ta2_yz_x_xxxz_0,  \
                             ta2_yz_x_xxxz_1,  \
                             ta2_yz_x_xxy_0,   \
                             ta2_yz_x_xxy_1,   \
                             ta2_yz_x_xxyy_0,  \
                             ta2_yz_x_xxyy_1,  \
                             ta2_yz_x_xxyz_0,  \
                             ta2_yz_x_xxyz_1,  \
                             ta2_yz_x_xxz_0,   \
                             ta2_yz_x_xxz_1,   \
                             ta2_yz_x_xxzz_0,  \
                             ta2_yz_x_xxzz_1,  \
                             ta2_yz_x_xyy_0,   \
                             ta2_yz_x_xyy_1,   \
                             ta2_yz_x_xyyy_0,  \
                             ta2_yz_x_xyyy_1,  \
                             ta2_yz_x_xyyz_0,  \
                             ta2_yz_x_xyyz_1,  \
                             ta2_yz_x_xyz_0,   \
                             ta2_yz_x_xyz_1,   \
                             ta2_yz_x_xyzz_0,  \
                             ta2_yz_x_xyzz_1,  \
                             ta2_yz_x_xzz_0,   \
                             ta2_yz_x_xzz_1,   \
                             ta2_yz_x_xzzz_0,  \
                             ta2_yz_x_xzzz_1,  \
                             ta2_yz_x_yyy_0,   \
                             ta2_yz_x_yyy_1,   \
                             ta2_yz_x_yyyy_0,  \
                             ta2_yz_x_yyyy_1,  \
                             ta2_yz_x_yyyz_0,  \
                             ta2_yz_x_yyyz_1,  \
                             ta2_yz_x_yyz_0,   \
                             ta2_yz_x_yyz_1,   \
                             ta2_yz_x_yyzz_0,  \
                             ta2_yz_x_yyzz_1,  \
                             ta2_yz_x_yzz_0,   \
                             ta2_yz_x_yzz_1,   \
                             ta2_yz_x_yzzz_0,  \
                             ta2_yz_x_yzzz_1,  \
                             ta2_yz_x_zzz_0,   \
                             ta2_yz_x_zzz_1,   \
                             ta2_yz_x_zzzz_0,  \
                             ta2_yz_x_zzzz_1,  \
                             ta2_yz_xx_xxxx_0, \
                             ta2_yz_xx_xxxy_0, \
                             ta2_yz_xx_xxxz_0, \
                             ta2_yz_xx_xxyy_0, \
                             ta2_yz_xx_xxyz_0, \
                             ta2_yz_xx_xxzz_0, \
                             ta2_yz_xx_xyyy_0, \
                             ta2_yz_xx_xyyz_0, \
                             ta2_yz_xx_xyzz_0, \
                             ta2_yz_xx_xzzz_0, \
                             ta2_yz_xx_yyyy_0, \
                             ta2_yz_xx_yyyz_0, \
                             ta2_yz_xx_yyzz_0, \
                             ta2_yz_xx_yzzz_0, \
                             ta2_yz_xx_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xx_xxxx_0[i] = ta2_yz_0_xxxx_0[i] * fe_0 - ta2_yz_0_xxxx_1[i] * fe_0 + 4.0 * ta2_yz_x_xxx_0[i] * fe_0 -
                              4.0 * ta2_yz_x_xxx_1[i] * fe_0 + ta2_yz_x_xxxx_0[i] * pa_x[i] - ta2_yz_x_xxxx_1[i] * pc_x[i];

        ta2_yz_xx_xxxy_0[i] = ta2_yz_0_xxxy_0[i] * fe_0 - ta2_yz_0_xxxy_1[i] * fe_0 + 3.0 * ta2_yz_x_xxy_0[i] * fe_0 -
                              3.0 * ta2_yz_x_xxy_1[i] * fe_0 + ta2_yz_x_xxxy_0[i] * pa_x[i] - ta2_yz_x_xxxy_1[i] * pc_x[i];

        ta2_yz_xx_xxxz_0[i] = ta2_yz_0_xxxz_0[i] * fe_0 - ta2_yz_0_xxxz_1[i] * fe_0 + 3.0 * ta2_yz_x_xxz_0[i] * fe_0 -
                              3.0 * ta2_yz_x_xxz_1[i] * fe_0 + ta2_yz_x_xxxz_0[i] * pa_x[i] - ta2_yz_x_xxxz_1[i] * pc_x[i];

        ta2_yz_xx_xxyy_0[i] = ta2_yz_0_xxyy_0[i] * fe_0 - ta2_yz_0_xxyy_1[i] * fe_0 + 2.0 * ta2_yz_x_xyy_0[i] * fe_0 -
                              2.0 * ta2_yz_x_xyy_1[i] * fe_0 + ta2_yz_x_xxyy_0[i] * pa_x[i] - ta2_yz_x_xxyy_1[i] * pc_x[i];

        ta2_yz_xx_xxyz_0[i] = ta2_yz_0_xxyz_0[i] * fe_0 - ta2_yz_0_xxyz_1[i] * fe_0 + 2.0 * ta2_yz_x_xyz_0[i] * fe_0 -
                              2.0 * ta2_yz_x_xyz_1[i] * fe_0 + ta2_yz_x_xxyz_0[i] * pa_x[i] - ta2_yz_x_xxyz_1[i] * pc_x[i];

        ta2_yz_xx_xxzz_0[i] = ta2_yz_0_xxzz_0[i] * fe_0 - ta2_yz_0_xxzz_1[i] * fe_0 + 2.0 * ta2_yz_x_xzz_0[i] * fe_0 -
                              2.0 * ta2_yz_x_xzz_1[i] * fe_0 + ta2_yz_x_xxzz_0[i] * pa_x[i] - ta2_yz_x_xxzz_1[i] * pc_x[i];

        ta2_yz_xx_xyyy_0[i] = ta2_yz_0_xyyy_0[i] * fe_0 - ta2_yz_0_xyyy_1[i] * fe_0 + ta2_yz_x_yyy_0[i] * fe_0 - ta2_yz_x_yyy_1[i] * fe_0 +
                              ta2_yz_x_xyyy_0[i] * pa_x[i] - ta2_yz_x_xyyy_1[i] * pc_x[i];

        ta2_yz_xx_xyyz_0[i] = ta2_yz_0_xyyz_0[i] * fe_0 - ta2_yz_0_xyyz_1[i] * fe_0 + ta2_yz_x_yyz_0[i] * fe_0 - ta2_yz_x_yyz_1[i] * fe_0 +
                              ta2_yz_x_xyyz_0[i] * pa_x[i] - ta2_yz_x_xyyz_1[i] * pc_x[i];

        ta2_yz_xx_xyzz_0[i] = ta2_yz_0_xyzz_0[i] * fe_0 - ta2_yz_0_xyzz_1[i] * fe_0 + ta2_yz_x_yzz_0[i] * fe_0 - ta2_yz_x_yzz_1[i] * fe_0 +
                              ta2_yz_x_xyzz_0[i] * pa_x[i] - ta2_yz_x_xyzz_1[i] * pc_x[i];

        ta2_yz_xx_xzzz_0[i] = ta2_yz_0_xzzz_0[i] * fe_0 - ta2_yz_0_xzzz_1[i] * fe_0 + ta2_yz_x_zzz_0[i] * fe_0 - ta2_yz_x_zzz_1[i] * fe_0 +
                              ta2_yz_x_xzzz_0[i] * pa_x[i] - ta2_yz_x_xzzz_1[i] * pc_x[i];

        ta2_yz_xx_yyyy_0[i] = ta2_yz_0_yyyy_0[i] * fe_0 - ta2_yz_0_yyyy_1[i] * fe_0 + ta2_yz_x_yyyy_0[i] * pa_x[i] - ta2_yz_x_yyyy_1[i] * pc_x[i];

        ta2_yz_xx_yyyz_0[i] = ta2_yz_0_yyyz_0[i] * fe_0 - ta2_yz_0_yyyz_1[i] * fe_0 + ta2_yz_x_yyyz_0[i] * pa_x[i] - ta2_yz_x_yyyz_1[i] * pc_x[i];

        ta2_yz_xx_yyzz_0[i] = ta2_yz_0_yyzz_0[i] * fe_0 - ta2_yz_0_yyzz_1[i] * fe_0 + ta2_yz_x_yyzz_0[i] * pa_x[i] - ta2_yz_x_yyzz_1[i] * pc_x[i];

        ta2_yz_xx_yzzz_0[i] = ta2_yz_0_yzzz_0[i] * fe_0 - ta2_yz_0_yzzz_1[i] * fe_0 + ta2_yz_x_yzzz_0[i] * pa_x[i] - ta2_yz_x_yzzz_1[i] * pc_x[i];

        ta2_yz_xx_zzzz_0[i] = ta2_yz_0_zzzz_0[i] * fe_0 - ta2_yz_0_zzzz_1[i] * fe_0 + ta2_yz_x_zzzz_0[i] * pa_x[i] - ta2_yz_x_zzzz_1[i] * pc_x[i];
    }

    // Set up 375-390 components of targeted buffer : DG

    auto ta2_yz_xy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 375);

    auto ta2_yz_xy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 376);

    auto ta2_yz_xy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 377);

    auto ta2_yz_xy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 378);

    auto ta2_yz_xy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 379);

    auto ta2_yz_xy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 380);

    auto ta2_yz_xy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 381);

    auto ta2_yz_xy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 382);

    auto ta2_yz_xy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 383);

    auto ta2_yz_xy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 384);

    auto ta2_yz_xy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 385);

    auto ta2_yz_xy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 386);

    auto ta2_yz_xy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 387);

    auto ta2_yz_xy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 388);

    auto ta2_yz_xy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 389);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_z_x_xxxx_1,   \
                             ta1_z_x_xxxz_1,   \
                             ta1_z_x_xxzz_1,   \
                             ta1_z_x_xzzz_1,   \
                             ta2_yz_x_xxxx_0,  \
                             ta2_yz_x_xxxx_1,  \
                             ta2_yz_x_xxxz_0,  \
                             ta2_yz_x_xxxz_1,  \
                             ta2_yz_x_xxzz_0,  \
                             ta2_yz_x_xxzz_1,  \
                             ta2_yz_x_xzzz_0,  \
                             ta2_yz_x_xzzz_1,  \
                             ta2_yz_xy_xxxx_0, \
                             ta2_yz_xy_xxxy_0, \
                             ta2_yz_xy_xxxz_0, \
                             ta2_yz_xy_xxyy_0, \
                             ta2_yz_xy_xxyz_0, \
                             ta2_yz_xy_xxzz_0, \
                             ta2_yz_xy_xyyy_0, \
                             ta2_yz_xy_xyyz_0, \
                             ta2_yz_xy_xyzz_0, \
                             ta2_yz_xy_xzzz_0, \
                             ta2_yz_xy_yyyy_0, \
                             ta2_yz_xy_yyyz_0, \
                             ta2_yz_xy_yyzz_0, \
                             ta2_yz_xy_yzzz_0, \
                             ta2_yz_xy_zzzz_0, \
                             ta2_yz_y_xxxy_0,  \
                             ta2_yz_y_xxxy_1,  \
                             ta2_yz_y_xxy_0,   \
                             ta2_yz_y_xxy_1,   \
                             ta2_yz_y_xxyy_0,  \
                             ta2_yz_y_xxyy_1,  \
                             ta2_yz_y_xxyz_0,  \
                             ta2_yz_y_xxyz_1,  \
                             ta2_yz_y_xyy_0,   \
                             ta2_yz_y_xyy_1,   \
                             ta2_yz_y_xyyy_0,  \
                             ta2_yz_y_xyyy_1,  \
                             ta2_yz_y_xyyz_0,  \
                             ta2_yz_y_xyyz_1,  \
                             ta2_yz_y_xyz_0,   \
                             ta2_yz_y_xyz_1,   \
                             ta2_yz_y_xyzz_0,  \
                             ta2_yz_y_xyzz_1,  \
                             ta2_yz_y_yyy_0,   \
                             ta2_yz_y_yyy_1,   \
                             ta2_yz_y_yyyy_0,  \
                             ta2_yz_y_yyyy_1,  \
                             ta2_yz_y_yyyz_0,  \
                             ta2_yz_y_yyyz_1,  \
                             ta2_yz_y_yyz_0,   \
                             ta2_yz_y_yyz_1,   \
                             ta2_yz_y_yyzz_0,  \
                             ta2_yz_y_yyzz_1,  \
                             ta2_yz_y_yzz_0,   \
                             ta2_yz_y_yzz_1,   \
                             ta2_yz_y_yzzz_0,  \
                             ta2_yz_y_yzzz_1,  \
                             ta2_yz_y_zzzz_0,  \
                             ta2_yz_y_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xy_xxxx_0[i] = ta1_z_x_xxxx_1[i] + ta2_yz_x_xxxx_0[i] * pa_y[i] - ta2_yz_x_xxxx_1[i] * pc_y[i];

        ta2_yz_xy_xxxy_0[i] =
            3.0 * ta2_yz_y_xxy_0[i] * fe_0 - 3.0 * ta2_yz_y_xxy_1[i] * fe_0 + ta2_yz_y_xxxy_0[i] * pa_x[i] - ta2_yz_y_xxxy_1[i] * pc_x[i];

        ta2_yz_xy_xxxz_0[i] = ta1_z_x_xxxz_1[i] + ta2_yz_x_xxxz_0[i] * pa_y[i] - ta2_yz_x_xxxz_1[i] * pc_y[i];

        ta2_yz_xy_xxyy_0[i] =
            2.0 * ta2_yz_y_xyy_0[i] * fe_0 - 2.0 * ta2_yz_y_xyy_1[i] * fe_0 + ta2_yz_y_xxyy_0[i] * pa_x[i] - ta2_yz_y_xxyy_1[i] * pc_x[i];

        ta2_yz_xy_xxyz_0[i] =
            2.0 * ta2_yz_y_xyz_0[i] * fe_0 - 2.0 * ta2_yz_y_xyz_1[i] * fe_0 + ta2_yz_y_xxyz_0[i] * pa_x[i] - ta2_yz_y_xxyz_1[i] * pc_x[i];

        ta2_yz_xy_xxzz_0[i] = ta1_z_x_xxzz_1[i] + ta2_yz_x_xxzz_0[i] * pa_y[i] - ta2_yz_x_xxzz_1[i] * pc_y[i];

        ta2_yz_xy_xyyy_0[i] = ta2_yz_y_yyy_0[i] * fe_0 - ta2_yz_y_yyy_1[i] * fe_0 + ta2_yz_y_xyyy_0[i] * pa_x[i] - ta2_yz_y_xyyy_1[i] * pc_x[i];

        ta2_yz_xy_xyyz_0[i] = ta2_yz_y_yyz_0[i] * fe_0 - ta2_yz_y_yyz_1[i] * fe_0 + ta2_yz_y_xyyz_0[i] * pa_x[i] - ta2_yz_y_xyyz_1[i] * pc_x[i];

        ta2_yz_xy_xyzz_0[i] = ta2_yz_y_yzz_0[i] * fe_0 - ta2_yz_y_yzz_1[i] * fe_0 + ta2_yz_y_xyzz_0[i] * pa_x[i] - ta2_yz_y_xyzz_1[i] * pc_x[i];

        ta2_yz_xy_xzzz_0[i] = ta1_z_x_xzzz_1[i] + ta2_yz_x_xzzz_0[i] * pa_y[i] - ta2_yz_x_xzzz_1[i] * pc_y[i];

        ta2_yz_xy_yyyy_0[i] = ta2_yz_y_yyyy_0[i] * pa_x[i] - ta2_yz_y_yyyy_1[i] * pc_x[i];

        ta2_yz_xy_yyyz_0[i] = ta2_yz_y_yyyz_0[i] * pa_x[i] - ta2_yz_y_yyyz_1[i] * pc_x[i];

        ta2_yz_xy_yyzz_0[i] = ta2_yz_y_yyzz_0[i] * pa_x[i] - ta2_yz_y_yyzz_1[i] * pc_x[i];

        ta2_yz_xy_yzzz_0[i] = ta2_yz_y_yzzz_0[i] * pa_x[i] - ta2_yz_y_yzzz_1[i] * pc_x[i];

        ta2_yz_xy_zzzz_0[i] = ta2_yz_y_zzzz_0[i] * pa_x[i] - ta2_yz_y_zzzz_1[i] * pc_x[i];
    }

    // Set up 390-405 components of targeted buffer : DG

    auto ta2_yz_xz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 390);

    auto ta2_yz_xz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 391);

    auto ta2_yz_xz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 392);

    auto ta2_yz_xz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 393);

    auto ta2_yz_xz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 394);

    auto ta2_yz_xz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 395);

    auto ta2_yz_xz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 396);

    auto ta2_yz_xz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 397);

    auto ta2_yz_xz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 398);

    auto ta2_yz_xz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 399);

    auto ta2_yz_xz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 400);

    auto ta2_yz_xz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 401);

    auto ta2_yz_xz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 402);

    auto ta2_yz_xz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 403);

    auto ta2_yz_xz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 404);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_y_x_xxxx_1,   \
                             ta1_y_x_xxxy_1,   \
                             ta1_y_x_xxyy_1,   \
                             ta1_y_x_xyyy_1,   \
                             ta2_yz_x_xxxx_0,  \
                             ta2_yz_x_xxxx_1,  \
                             ta2_yz_x_xxxy_0,  \
                             ta2_yz_x_xxxy_1,  \
                             ta2_yz_x_xxyy_0,  \
                             ta2_yz_x_xxyy_1,  \
                             ta2_yz_x_xyyy_0,  \
                             ta2_yz_x_xyyy_1,  \
                             ta2_yz_xz_xxxx_0, \
                             ta2_yz_xz_xxxy_0, \
                             ta2_yz_xz_xxxz_0, \
                             ta2_yz_xz_xxyy_0, \
                             ta2_yz_xz_xxyz_0, \
                             ta2_yz_xz_xxzz_0, \
                             ta2_yz_xz_xyyy_0, \
                             ta2_yz_xz_xyyz_0, \
                             ta2_yz_xz_xyzz_0, \
                             ta2_yz_xz_xzzz_0, \
                             ta2_yz_xz_yyyy_0, \
                             ta2_yz_xz_yyyz_0, \
                             ta2_yz_xz_yyzz_0, \
                             ta2_yz_xz_yzzz_0, \
                             ta2_yz_xz_zzzz_0, \
                             ta2_yz_z_xxxz_0,  \
                             ta2_yz_z_xxxz_1,  \
                             ta2_yz_z_xxyz_0,  \
                             ta2_yz_z_xxyz_1,  \
                             ta2_yz_z_xxz_0,   \
                             ta2_yz_z_xxz_1,   \
                             ta2_yz_z_xxzz_0,  \
                             ta2_yz_z_xxzz_1,  \
                             ta2_yz_z_xyyz_0,  \
                             ta2_yz_z_xyyz_1,  \
                             ta2_yz_z_xyz_0,   \
                             ta2_yz_z_xyz_1,   \
                             ta2_yz_z_xyzz_0,  \
                             ta2_yz_z_xyzz_1,  \
                             ta2_yz_z_xzz_0,   \
                             ta2_yz_z_xzz_1,   \
                             ta2_yz_z_xzzz_0,  \
                             ta2_yz_z_xzzz_1,  \
                             ta2_yz_z_yyyy_0,  \
                             ta2_yz_z_yyyy_1,  \
                             ta2_yz_z_yyyz_0,  \
                             ta2_yz_z_yyyz_1,  \
                             ta2_yz_z_yyz_0,   \
                             ta2_yz_z_yyz_1,   \
                             ta2_yz_z_yyzz_0,  \
                             ta2_yz_z_yyzz_1,  \
                             ta2_yz_z_yzz_0,   \
                             ta2_yz_z_yzz_1,   \
                             ta2_yz_z_yzzz_0,  \
                             ta2_yz_z_yzzz_1,  \
                             ta2_yz_z_zzz_0,   \
                             ta2_yz_z_zzz_1,   \
                             ta2_yz_z_zzzz_0,  \
                             ta2_yz_z_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xz_xxxx_0[i] = ta1_y_x_xxxx_1[i] + ta2_yz_x_xxxx_0[i] * pa_z[i] - ta2_yz_x_xxxx_1[i] * pc_z[i];

        ta2_yz_xz_xxxy_0[i] = ta1_y_x_xxxy_1[i] + ta2_yz_x_xxxy_0[i] * pa_z[i] - ta2_yz_x_xxxy_1[i] * pc_z[i];

        ta2_yz_xz_xxxz_0[i] =
            3.0 * ta2_yz_z_xxz_0[i] * fe_0 - 3.0 * ta2_yz_z_xxz_1[i] * fe_0 + ta2_yz_z_xxxz_0[i] * pa_x[i] - ta2_yz_z_xxxz_1[i] * pc_x[i];

        ta2_yz_xz_xxyy_0[i] = ta1_y_x_xxyy_1[i] + ta2_yz_x_xxyy_0[i] * pa_z[i] - ta2_yz_x_xxyy_1[i] * pc_z[i];

        ta2_yz_xz_xxyz_0[i] =
            2.0 * ta2_yz_z_xyz_0[i] * fe_0 - 2.0 * ta2_yz_z_xyz_1[i] * fe_0 + ta2_yz_z_xxyz_0[i] * pa_x[i] - ta2_yz_z_xxyz_1[i] * pc_x[i];

        ta2_yz_xz_xxzz_0[i] =
            2.0 * ta2_yz_z_xzz_0[i] * fe_0 - 2.0 * ta2_yz_z_xzz_1[i] * fe_0 + ta2_yz_z_xxzz_0[i] * pa_x[i] - ta2_yz_z_xxzz_1[i] * pc_x[i];

        ta2_yz_xz_xyyy_0[i] = ta1_y_x_xyyy_1[i] + ta2_yz_x_xyyy_0[i] * pa_z[i] - ta2_yz_x_xyyy_1[i] * pc_z[i];

        ta2_yz_xz_xyyz_0[i] = ta2_yz_z_yyz_0[i] * fe_0 - ta2_yz_z_yyz_1[i] * fe_0 + ta2_yz_z_xyyz_0[i] * pa_x[i] - ta2_yz_z_xyyz_1[i] * pc_x[i];

        ta2_yz_xz_xyzz_0[i] = ta2_yz_z_yzz_0[i] * fe_0 - ta2_yz_z_yzz_1[i] * fe_0 + ta2_yz_z_xyzz_0[i] * pa_x[i] - ta2_yz_z_xyzz_1[i] * pc_x[i];

        ta2_yz_xz_xzzz_0[i] = ta2_yz_z_zzz_0[i] * fe_0 - ta2_yz_z_zzz_1[i] * fe_0 + ta2_yz_z_xzzz_0[i] * pa_x[i] - ta2_yz_z_xzzz_1[i] * pc_x[i];

        ta2_yz_xz_yyyy_0[i] = ta2_yz_z_yyyy_0[i] * pa_x[i] - ta2_yz_z_yyyy_1[i] * pc_x[i];

        ta2_yz_xz_yyyz_0[i] = ta2_yz_z_yyyz_0[i] * pa_x[i] - ta2_yz_z_yyyz_1[i] * pc_x[i];

        ta2_yz_xz_yyzz_0[i] = ta2_yz_z_yyzz_0[i] * pa_x[i] - ta2_yz_z_yyzz_1[i] * pc_x[i];

        ta2_yz_xz_yzzz_0[i] = ta2_yz_z_yzzz_0[i] * pa_x[i] - ta2_yz_z_yzzz_1[i] * pc_x[i];

        ta2_yz_xz_zzzz_0[i] = ta2_yz_z_zzzz_0[i] * pa_x[i] - ta2_yz_z_zzzz_1[i] * pc_x[i];
    }

    // Set up 405-420 components of targeted buffer : DG

    auto ta2_yz_yy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 405);

    auto ta2_yz_yy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 406);

    auto ta2_yz_yy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 407);

    auto ta2_yz_yy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 408);

    auto ta2_yz_yy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 409);

    auto ta2_yz_yy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 410);

    auto ta2_yz_yy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 411);

    auto ta2_yz_yy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 412);

    auto ta2_yz_yy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 413);

    auto ta2_yz_yy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 414);

    auto ta2_yz_yy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 415);

    auto ta2_yz_yy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 416);

    auto ta2_yz_yy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 417);

    auto ta2_yz_yy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 418);

    auto ta2_yz_yy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 419);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_z_y_xxxx_1,   \
                             ta1_z_y_xxxy_1,   \
                             ta1_z_y_xxxz_1,   \
                             ta1_z_y_xxyy_1,   \
                             ta1_z_y_xxyz_1,   \
                             ta1_z_y_xxzz_1,   \
                             ta1_z_y_xyyy_1,   \
                             ta1_z_y_xyyz_1,   \
                             ta1_z_y_xyzz_1,   \
                             ta1_z_y_xzzz_1,   \
                             ta1_z_y_yyyy_1,   \
                             ta1_z_y_yyyz_1,   \
                             ta1_z_y_yyzz_1,   \
                             ta1_z_y_yzzz_1,   \
                             ta1_z_y_zzzz_1,   \
                             ta2_yz_0_xxxx_0,  \
                             ta2_yz_0_xxxx_1,  \
                             ta2_yz_0_xxxy_0,  \
                             ta2_yz_0_xxxy_1,  \
                             ta2_yz_0_xxxz_0,  \
                             ta2_yz_0_xxxz_1,  \
                             ta2_yz_0_xxyy_0,  \
                             ta2_yz_0_xxyy_1,  \
                             ta2_yz_0_xxyz_0,  \
                             ta2_yz_0_xxyz_1,  \
                             ta2_yz_0_xxzz_0,  \
                             ta2_yz_0_xxzz_1,  \
                             ta2_yz_0_xyyy_0,  \
                             ta2_yz_0_xyyy_1,  \
                             ta2_yz_0_xyyz_0,  \
                             ta2_yz_0_xyyz_1,  \
                             ta2_yz_0_xyzz_0,  \
                             ta2_yz_0_xyzz_1,  \
                             ta2_yz_0_xzzz_0,  \
                             ta2_yz_0_xzzz_1,  \
                             ta2_yz_0_yyyy_0,  \
                             ta2_yz_0_yyyy_1,  \
                             ta2_yz_0_yyyz_0,  \
                             ta2_yz_0_yyyz_1,  \
                             ta2_yz_0_yyzz_0,  \
                             ta2_yz_0_yyzz_1,  \
                             ta2_yz_0_yzzz_0,  \
                             ta2_yz_0_yzzz_1,  \
                             ta2_yz_0_zzzz_0,  \
                             ta2_yz_0_zzzz_1,  \
                             ta2_yz_y_xxx_0,   \
                             ta2_yz_y_xxx_1,   \
                             ta2_yz_y_xxxx_0,  \
                             ta2_yz_y_xxxx_1,  \
                             ta2_yz_y_xxxy_0,  \
                             ta2_yz_y_xxxy_1,  \
                             ta2_yz_y_xxxz_0,  \
                             ta2_yz_y_xxxz_1,  \
                             ta2_yz_y_xxy_0,   \
                             ta2_yz_y_xxy_1,   \
                             ta2_yz_y_xxyy_0,  \
                             ta2_yz_y_xxyy_1,  \
                             ta2_yz_y_xxyz_0,  \
                             ta2_yz_y_xxyz_1,  \
                             ta2_yz_y_xxz_0,   \
                             ta2_yz_y_xxz_1,   \
                             ta2_yz_y_xxzz_0,  \
                             ta2_yz_y_xxzz_1,  \
                             ta2_yz_y_xyy_0,   \
                             ta2_yz_y_xyy_1,   \
                             ta2_yz_y_xyyy_0,  \
                             ta2_yz_y_xyyy_1,  \
                             ta2_yz_y_xyyz_0,  \
                             ta2_yz_y_xyyz_1,  \
                             ta2_yz_y_xyz_0,   \
                             ta2_yz_y_xyz_1,   \
                             ta2_yz_y_xyzz_0,  \
                             ta2_yz_y_xyzz_1,  \
                             ta2_yz_y_xzz_0,   \
                             ta2_yz_y_xzz_1,   \
                             ta2_yz_y_xzzz_0,  \
                             ta2_yz_y_xzzz_1,  \
                             ta2_yz_y_yyy_0,   \
                             ta2_yz_y_yyy_1,   \
                             ta2_yz_y_yyyy_0,  \
                             ta2_yz_y_yyyy_1,  \
                             ta2_yz_y_yyyz_0,  \
                             ta2_yz_y_yyyz_1,  \
                             ta2_yz_y_yyz_0,   \
                             ta2_yz_y_yyz_1,   \
                             ta2_yz_y_yyzz_0,  \
                             ta2_yz_y_yyzz_1,  \
                             ta2_yz_y_yzz_0,   \
                             ta2_yz_y_yzz_1,   \
                             ta2_yz_y_yzzz_0,  \
                             ta2_yz_y_yzzz_1,  \
                             ta2_yz_y_zzz_0,   \
                             ta2_yz_y_zzz_1,   \
                             ta2_yz_y_zzzz_0,  \
                             ta2_yz_y_zzzz_1,  \
                             ta2_yz_yy_xxxx_0, \
                             ta2_yz_yy_xxxy_0, \
                             ta2_yz_yy_xxxz_0, \
                             ta2_yz_yy_xxyy_0, \
                             ta2_yz_yy_xxyz_0, \
                             ta2_yz_yy_xxzz_0, \
                             ta2_yz_yy_xyyy_0, \
                             ta2_yz_yy_xyyz_0, \
                             ta2_yz_yy_xyzz_0, \
                             ta2_yz_yy_xzzz_0, \
                             ta2_yz_yy_yyyy_0, \
                             ta2_yz_yy_yyyz_0, \
                             ta2_yz_yy_yyzz_0, \
                             ta2_yz_yy_yzzz_0, \
                             ta2_yz_yy_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yy_xxxx_0[i] =
            ta2_yz_0_xxxx_0[i] * fe_0 - ta2_yz_0_xxxx_1[i] * fe_0 + ta1_z_y_xxxx_1[i] + ta2_yz_y_xxxx_0[i] * pa_y[i] - ta2_yz_y_xxxx_1[i] * pc_y[i];

        ta2_yz_yy_xxxy_0[i] = ta2_yz_0_xxxy_0[i] * fe_0 - ta2_yz_0_xxxy_1[i] * fe_0 + ta2_yz_y_xxx_0[i] * fe_0 - ta2_yz_y_xxx_1[i] * fe_0 +
                              ta1_z_y_xxxy_1[i] + ta2_yz_y_xxxy_0[i] * pa_y[i] - ta2_yz_y_xxxy_1[i] * pc_y[i];

        ta2_yz_yy_xxxz_0[i] =
            ta2_yz_0_xxxz_0[i] * fe_0 - ta2_yz_0_xxxz_1[i] * fe_0 + ta1_z_y_xxxz_1[i] + ta2_yz_y_xxxz_0[i] * pa_y[i] - ta2_yz_y_xxxz_1[i] * pc_y[i];

        ta2_yz_yy_xxyy_0[i] = ta2_yz_0_xxyy_0[i] * fe_0 - ta2_yz_0_xxyy_1[i] * fe_0 + 2.0 * ta2_yz_y_xxy_0[i] * fe_0 -
                              2.0 * ta2_yz_y_xxy_1[i] * fe_0 + ta1_z_y_xxyy_1[i] + ta2_yz_y_xxyy_0[i] * pa_y[i] - ta2_yz_y_xxyy_1[i] * pc_y[i];

        ta2_yz_yy_xxyz_0[i] = ta2_yz_0_xxyz_0[i] * fe_0 - ta2_yz_0_xxyz_1[i] * fe_0 + ta2_yz_y_xxz_0[i] * fe_0 - ta2_yz_y_xxz_1[i] * fe_0 +
                              ta1_z_y_xxyz_1[i] + ta2_yz_y_xxyz_0[i] * pa_y[i] - ta2_yz_y_xxyz_1[i] * pc_y[i];

        ta2_yz_yy_xxzz_0[i] =
            ta2_yz_0_xxzz_0[i] * fe_0 - ta2_yz_0_xxzz_1[i] * fe_0 + ta1_z_y_xxzz_1[i] + ta2_yz_y_xxzz_0[i] * pa_y[i] - ta2_yz_y_xxzz_1[i] * pc_y[i];

        ta2_yz_yy_xyyy_0[i] = ta2_yz_0_xyyy_0[i] * fe_0 - ta2_yz_0_xyyy_1[i] * fe_0 + 3.0 * ta2_yz_y_xyy_0[i] * fe_0 -
                              3.0 * ta2_yz_y_xyy_1[i] * fe_0 + ta1_z_y_xyyy_1[i] + ta2_yz_y_xyyy_0[i] * pa_y[i] - ta2_yz_y_xyyy_1[i] * pc_y[i];

        ta2_yz_yy_xyyz_0[i] = ta2_yz_0_xyyz_0[i] * fe_0 - ta2_yz_0_xyyz_1[i] * fe_0 + 2.0 * ta2_yz_y_xyz_0[i] * fe_0 -
                              2.0 * ta2_yz_y_xyz_1[i] * fe_0 + ta1_z_y_xyyz_1[i] + ta2_yz_y_xyyz_0[i] * pa_y[i] - ta2_yz_y_xyyz_1[i] * pc_y[i];

        ta2_yz_yy_xyzz_0[i] = ta2_yz_0_xyzz_0[i] * fe_0 - ta2_yz_0_xyzz_1[i] * fe_0 + ta2_yz_y_xzz_0[i] * fe_0 - ta2_yz_y_xzz_1[i] * fe_0 +
                              ta1_z_y_xyzz_1[i] + ta2_yz_y_xyzz_0[i] * pa_y[i] - ta2_yz_y_xyzz_1[i] * pc_y[i];

        ta2_yz_yy_xzzz_0[i] =
            ta2_yz_0_xzzz_0[i] * fe_0 - ta2_yz_0_xzzz_1[i] * fe_0 + ta1_z_y_xzzz_1[i] + ta2_yz_y_xzzz_0[i] * pa_y[i] - ta2_yz_y_xzzz_1[i] * pc_y[i];

        ta2_yz_yy_yyyy_0[i] = ta2_yz_0_yyyy_0[i] * fe_0 - ta2_yz_0_yyyy_1[i] * fe_0 + 4.0 * ta2_yz_y_yyy_0[i] * fe_0 -
                              4.0 * ta2_yz_y_yyy_1[i] * fe_0 + ta1_z_y_yyyy_1[i] + ta2_yz_y_yyyy_0[i] * pa_y[i] - ta2_yz_y_yyyy_1[i] * pc_y[i];

        ta2_yz_yy_yyyz_0[i] = ta2_yz_0_yyyz_0[i] * fe_0 - ta2_yz_0_yyyz_1[i] * fe_0 + 3.0 * ta2_yz_y_yyz_0[i] * fe_0 -
                              3.0 * ta2_yz_y_yyz_1[i] * fe_0 + ta1_z_y_yyyz_1[i] + ta2_yz_y_yyyz_0[i] * pa_y[i] - ta2_yz_y_yyyz_1[i] * pc_y[i];

        ta2_yz_yy_yyzz_0[i] = ta2_yz_0_yyzz_0[i] * fe_0 - ta2_yz_0_yyzz_1[i] * fe_0 + 2.0 * ta2_yz_y_yzz_0[i] * fe_0 -
                              2.0 * ta2_yz_y_yzz_1[i] * fe_0 + ta1_z_y_yyzz_1[i] + ta2_yz_y_yyzz_0[i] * pa_y[i] - ta2_yz_y_yyzz_1[i] * pc_y[i];

        ta2_yz_yy_yzzz_0[i] = ta2_yz_0_yzzz_0[i] * fe_0 - ta2_yz_0_yzzz_1[i] * fe_0 + ta2_yz_y_zzz_0[i] * fe_0 - ta2_yz_y_zzz_1[i] * fe_0 +
                              ta1_z_y_yzzz_1[i] + ta2_yz_y_yzzz_0[i] * pa_y[i] - ta2_yz_y_yzzz_1[i] * pc_y[i];

        ta2_yz_yy_zzzz_0[i] =
            ta2_yz_0_zzzz_0[i] * fe_0 - ta2_yz_0_zzzz_1[i] * fe_0 + ta1_z_y_zzzz_1[i] + ta2_yz_y_zzzz_0[i] * pa_y[i] - ta2_yz_y_zzzz_1[i] * pc_y[i];
    }

    // Set up 420-435 components of targeted buffer : DG

    auto ta2_yz_yz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 420);

    auto ta2_yz_yz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 421);

    auto ta2_yz_yz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 422);

    auto ta2_yz_yz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 423);

    auto ta2_yz_yz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 424);

    auto ta2_yz_yz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 425);

    auto ta2_yz_yz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 426);

    auto ta2_yz_yz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 427);

    auto ta2_yz_yz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 428);

    auto ta2_yz_yz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 429);

    auto ta2_yz_yz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 430);

    auto ta2_yz_yz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 431);

    auto ta2_yz_yz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 432);

    auto ta2_yz_yz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 433);

    auto ta2_yz_yz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 434);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_y_y_xxxy_1,   \
                             ta1_y_y_xxyy_1,   \
                             ta1_y_y_xyyy_1,   \
                             ta1_y_y_yyyy_1,   \
                             ta1_z_z_xxxx_1,   \
                             ta1_z_z_xxxz_1,   \
                             ta1_z_z_xxyz_1,   \
                             ta1_z_z_xxzz_1,   \
                             ta1_z_z_xyyz_1,   \
                             ta1_z_z_xyzz_1,   \
                             ta1_z_z_xzzz_1,   \
                             ta1_z_z_yyyz_1,   \
                             ta1_z_z_yyzz_1,   \
                             ta1_z_z_yzzz_1,   \
                             ta1_z_z_zzzz_1,   \
                             ta2_yz_y_xxxy_0,  \
                             ta2_yz_y_xxxy_1,  \
                             ta2_yz_y_xxyy_0,  \
                             ta2_yz_y_xxyy_1,  \
                             ta2_yz_y_xyyy_0,  \
                             ta2_yz_y_xyyy_1,  \
                             ta2_yz_y_yyyy_0,  \
                             ta2_yz_y_yyyy_1,  \
                             ta2_yz_yz_xxxx_0, \
                             ta2_yz_yz_xxxy_0, \
                             ta2_yz_yz_xxxz_0, \
                             ta2_yz_yz_xxyy_0, \
                             ta2_yz_yz_xxyz_0, \
                             ta2_yz_yz_xxzz_0, \
                             ta2_yz_yz_xyyy_0, \
                             ta2_yz_yz_xyyz_0, \
                             ta2_yz_yz_xyzz_0, \
                             ta2_yz_yz_xzzz_0, \
                             ta2_yz_yz_yyyy_0, \
                             ta2_yz_yz_yyyz_0, \
                             ta2_yz_yz_yyzz_0, \
                             ta2_yz_yz_yzzz_0, \
                             ta2_yz_yz_zzzz_0, \
                             ta2_yz_z_xxxx_0,  \
                             ta2_yz_z_xxxx_1,  \
                             ta2_yz_z_xxxz_0,  \
                             ta2_yz_z_xxxz_1,  \
                             ta2_yz_z_xxyz_0,  \
                             ta2_yz_z_xxyz_1,  \
                             ta2_yz_z_xxz_0,   \
                             ta2_yz_z_xxz_1,   \
                             ta2_yz_z_xxzz_0,  \
                             ta2_yz_z_xxzz_1,  \
                             ta2_yz_z_xyyz_0,  \
                             ta2_yz_z_xyyz_1,  \
                             ta2_yz_z_xyz_0,   \
                             ta2_yz_z_xyz_1,   \
                             ta2_yz_z_xyzz_0,  \
                             ta2_yz_z_xyzz_1,  \
                             ta2_yz_z_xzz_0,   \
                             ta2_yz_z_xzz_1,   \
                             ta2_yz_z_xzzz_0,  \
                             ta2_yz_z_xzzz_1,  \
                             ta2_yz_z_yyyz_0,  \
                             ta2_yz_z_yyyz_1,  \
                             ta2_yz_z_yyz_0,   \
                             ta2_yz_z_yyz_1,   \
                             ta2_yz_z_yyzz_0,  \
                             ta2_yz_z_yyzz_1,  \
                             ta2_yz_z_yzz_0,   \
                             ta2_yz_z_yzz_1,   \
                             ta2_yz_z_yzzz_0,  \
                             ta2_yz_z_yzzz_1,  \
                             ta2_yz_z_zzz_0,   \
                             ta2_yz_z_zzz_1,   \
                             ta2_yz_z_zzzz_0,  \
                             ta2_yz_z_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yz_xxxx_0[i] = ta1_z_z_xxxx_1[i] + ta2_yz_z_xxxx_0[i] * pa_y[i] - ta2_yz_z_xxxx_1[i] * pc_y[i];

        ta2_yz_yz_xxxy_0[i] = ta1_y_y_xxxy_1[i] + ta2_yz_y_xxxy_0[i] * pa_z[i] - ta2_yz_y_xxxy_1[i] * pc_z[i];

        ta2_yz_yz_xxxz_0[i] = ta1_z_z_xxxz_1[i] + ta2_yz_z_xxxz_0[i] * pa_y[i] - ta2_yz_z_xxxz_1[i] * pc_y[i];

        ta2_yz_yz_xxyy_0[i] = ta1_y_y_xxyy_1[i] + ta2_yz_y_xxyy_0[i] * pa_z[i] - ta2_yz_y_xxyy_1[i] * pc_z[i];

        ta2_yz_yz_xxyz_0[i] =
            ta2_yz_z_xxz_0[i] * fe_0 - ta2_yz_z_xxz_1[i] * fe_0 + ta1_z_z_xxyz_1[i] + ta2_yz_z_xxyz_0[i] * pa_y[i] - ta2_yz_z_xxyz_1[i] * pc_y[i];

        ta2_yz_yz_xxzz_0[i] = ta1_z_z_xxzz_1[i] + ta2_yz_z_xxzz_0[i] * pa_y[i] - ta2_yz_z_xxzz_1[i] * pc_y[i];

        ta2_yz_yz_xyyy_0[i] = ta1_y_y_xyyy_1[i] + ta2_yz_y_xyyy_0[i] * pa_z[i] - ta2_yz_y_xyyy_1[i] * pc_z[i];

        ta2_yz_yz_xyyz_0[i] = 2.0 * ta2_yz_z_xyz_0[i] * fe_0 - 2.0 * ta2_yz_z_xyz_1[i] * fe_0 + ta1_z_z_xyyz_1[i] + ta2_yz_z_xyyz_0[i] * pa_y[i] -
                              ta2_yz_z_xyyz_1[i] * pc_y[i];

        ta2_yz_yz_xyzz_0[i] =
            ta2_yz_z_xzz_0[i] * fe_0 - ta2_yz_z_xzz_1[i] * fe_0 + ta1_z_z_xyzz_1[i] + ta2_yz_z_xyzz_0[i] * pa_y[i] - ta2_yz_z_xyzz_1[i] * pc_y[i];

        ta2_yz_yz_xzzz_0[i] = ta1_z_z_xzzz_1[i] + ta2_yz_z_xzzz_0[i] * pa_y[i] - ta2_yz_z_xzzz_1[i] * pc_y[i];

        ta2_yz_yz_yyyy_0[i] = ta1_y_y_yyyy_1[i] + ta2_yz_y_yyyy_0[i] * pa_z[i] - ta2_yz_y_yyyy_1[i] * pc_z[i];

        ta2_yz_yz_yyyz_0[i] = 3.0 * ta2_yz_z_yyz_0[i] * fe_0 - 3.0 * ta2_yz_z_yyz_1[i] * fe_0 + ta1_z_z_yyyz_1[i] + ta2_yz_z_yyyz_0[i] * pa_y[i] -
                              ta2_yz_z_yyyz_1[i] * pc_y[i];

        ta2_yz_yz_yyzz_0[i] = 2.0 * ta2_yz_z_yzz_0[i] * fe_0 - 2.0 * ta2_yz_z_yzz_1[i] * fe_0 + ta1_z_z_yyzz_1[i] + ta2_yz_z_yyzz_0[i] * pa_y[i] -
                              ta2_yz_z_yyzz_1[i] * pc_y[i];

        ta2_yz_yz_yzzz_0[i] =
            ta2_yz_z_zzz_0[i] * fe_0 - ta2_yz_z_zzz_1[i] * fe_0 + ta1_z_z_yzzz_1[i] + ta2_yz_z_yzzz_0[i] * pa_y[i] - ta2_yz_z_yzzz_1[i] * pc_y[i];

        ta2_yz_yz_zzzz_0[i] = ta1_z_z_zzzz_1[i] + ta2_yz_z_zzzz_0[i] * pa_y[i] - ta2_yz_z_zzzz_1[i] * pc_y[i];
    }

    // Set up 435-450 components of targeted buffer : DG

    auto ta2_yz_zz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 435);

    auto ta2_yz_zz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 436);

    auto ta2_yz_zz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 437);

    auto ta2_yz_zz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 438);

    auto ta2_yz_zz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 439);

    auto ta2_yz_zz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 440);

    auto ta2_yz_zz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 441);

    auto ta2_yz_zz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 442);

    auto ta2_yz_zz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 443);

    auto ta2_yz_zz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 444);

    auto ta2_yz_zz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 445);

    auto ta2_yz_zz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 446);

    auto ta2_yz_zz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 447);

    auto ta2_yz_zz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 448);

    auto ta2_yz_zz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 449);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta1_y_z_xxxx_1,   \
                             ta1_y_z_xxxy_1,   \
                             ta1_y_z_xxxz_1,   \
                             ta1_y_z_xxyy_1,   \
                             ta1_y_z_xxyz_1,   \
                             ta1_y_z_xxzz_1,   \
                             ta1_y_z_xyyy_1,   \
                             ta1_y_z_xyyz_1,   \
                             ta1_y_z_xyzz_1,   \
                             ta1_y_z_xzzz_1,   \
                             ta1_y_z_yyyy_1,   \
                             ta1_y_z_yyyz_1,   \
                             ta1_y_z_yyzz_1,   \
                             ta1_y_z_yzzz_1,   \
                             ta1_y_z_zzzz_1,   \
                             ta2_yz_0_xxxx_0,  \
                             ta2_yz_0_xxxx_1,  \
                             ta2_yz_0_xxxy_0,  \
                             ta2_yz_0_xxxy_1,  \
                             ta2_yz_0_xxxz_0,  \
                             ta2_yz_0_xxxz_1,  \
                             ta2_yz_0_xxyy_0,  \
                             ta2_yz_0_xxyy_1,  \
                             ta2_yz_0_xxyz_0,  \
                             ta2_yz_0_xxyz_1,  \
                             ta2_yz_0_xxzz_0,  \
                             ta2_yz_0_xxzz_1,  \
                             ta2_yz_0_xyyy_0,  \
                             ta2_yz_0_xyyy_1,  \
                             ta2_yz_0_xyyz_0,  \
                             ta2_yz_0_xyyz_1,  \
                             ta2_yz_0_xyzz_0,  \
                             ta2_yz_0_xyzz_1,  \
                             ta2_yz_0_xzzz_0,  \
                             ta2_yz_0_xzzz_1,  \
                             ta2_yz_0_yyyy_0,  \
                             ta2_yz_0_yyyy_1,  \
                             ta2_yz_0_yyyz_0,  \
                             ta2_yz_0_yyyz_1,  \
                             ta2_yz_0_yyzz_0,  \
                             ta2_yz_0_yyzz_1,  \
                             ta2_yz_0_yzzz_0,  \
                             ta2_yz_0_yzzz_1,  \
                             ta2_yz_0_zzzz_0,  \
                             ta2_yz_0_zzzz_1,  \
                             ta2_yz_z_xxx_0,   \
                             ta2_yz_z_xxx_1,   \
                             ta2_yz_z_xxxx_0,  \
                             ta2_yz_z_xxxx_1,  \
                             ta2_yz_z_xxxy_0,  \
                             ta2_yz_z_xxxy_1,  \
                             ta2_yz_z_xxxz_0,  \
                             ta2_yz_z_xxxz_1,  \
                             ta2_yz_z_xxy_0,   \
                             ta2_yz_z_xxy_1,   \
                             ta2_yz_z_xxyy_0,  \
                             ta2_yz_z_xxyy_1,  \
                             ta2_yz_z_xxyz_0,  \
                             ta2_yz_z_xxyz_1,  \
                             ta2_yz_z_xxz_0,   \
                             ta2_yz_z_xxz_1,   \
                             ta2_yz_z_xxzz_0,  \
                             ta2_yz_z_xxzz_1,  \
                             ta2_yz_z_xyy_0,   \
                             ta2_yz_z_xyy_1,   \
                             ta2_yz_z_xyyy_0,  \
                             ta2_yz_z_xyyy_1,  \
                             ta2_yz_z_xyyz_0,  \
                             ta2_yz_z_xyyz_1,  \
                             ta2_yz_z_xyz_0,   \
                             ta2_yz_z_xyz_1,   \
                             ta2_yz_z_xyzz_0,  \
                             ta2_yz_z_xyzz_1,  \
                             ta2_yz_z_xzz_0,   \
                             ta2_yz_z_xzz_1,   \
                             ta2_yz_z_xzzz_0,  \
                             ta2_yz_z_xzzz_1,  \
                             ta2_yz_z_yyy_0,   \
                             ta2_yz_z_yyy_1,   \
                             ta2_yz_z_yyyy_0,  \
                             ta2_yz_z_yyyy_1,  \
                             ta2_yz_z_yyyz_0,  \
                             ta2_yz_z_yyyz_1,  \
                             ta2_yz_z_yyz_0,   \
                             ta2_yz_z_yyz_1,   \
                             ta2_yz_z_yyzz_0,  \
                             ta2_yz_z_yyzz_1,  \
                             ta2_yz_z_yzz_0,   \
                             ta2_yz_z_yzz_1,   \
                             ta2_yz_z_yzzz_0,  \
                             ta2_yz_z_yzzz_1,  \
                             ta2_yz_z_zzz_0,   \
                             ta2_yz_z_zzz_1,   \
                             ta2_yz_z_zzzz_0,  \
                             ta2_yz_z_zzzz_1,  \
                             ta2_yz_zz_xxxx_0, \
                             ta2_yz_zz_xxxy_0, \
                             ta2_yz_zz_xxxz_0, \
                             ta2_yz_zz_xxyy_0, \
                             ta2_yz_zz_xxyz_0, \
                             ta2_yz_zz_xxzz_0, \
                             ta2_yz_zz_xyyy_0, \
                             ta2_yz_zz_xyyz_0, \
                             ta2_yz_zz_xyzz_0, \
                             ta2_yz_zz_xzzz_0, \
                             ta2_yz_zz_yyyy_0, \
                             ta2_yz_zz_yyyz_0, \
                             ta2_yz_zz_yyzz_0, \
                             ta2_yz_zz_yzzz_0, \
                             ta2_yz_zz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_zz_xxxx_0[i] =
            ta2_yz_0_xxxx_0[i] * fe_0 - ta2_yz_0_xxxx_1[i] * fe_0 + ta1_y_z_xxxx_1[i] + ta2_yz_z_xxxx_0[i] * pa_z[i] - ta2_yz_z_xxxx_1[i] * pc_z[i];

        ta2_yz_zz_xxxy_0[i] =
            ta2_yz_0_xxxy_0[i] * fe_0 - ta2_yz_0_xxxy_1[i] * fe_0 + ta1_y_z_xxxy_1[i] + ta2_yz_z_xxxy_0[i] * pa_z[i] - ta2_yz_z_xxxy_1[i] * pc_z[i];

        ta2_yz_zz_xxxz_0[i] = ta2_yz_0_xxxz_0[i] * fe_0 - ta2_yz_0_xxxz_1[i] * fe_0 + ta2_yz_z_xxx_0[i] * fe_0 - ta2_yz_z_xxx_1[i] * fe_0 +
                              ta1_y_z_xxxz_1[i] + ta2_yz_z_xxxz_0[i] * pa_z[i] - ta2_yz_z_xxxz_1[i] * pc_z[i];

        ta2_yz_zz_xxyy_0[i] =
            ta2_yz_0_xxyy_0[i] * fe_0 - ta2_yz_0_xxyy_1[i] * fe_0 + ta1_y_z_xxyy_1[i] + ta2_yz_z_xxyy_0[i] * pa_z[i] - ta2_yz_z_xxyy_1[i] * pc_z[i];

        ta2_yz_zz_xxyz_0[i] = ta2_yz_0_xxyz_0[i] * fe_0 - ta2_yz_0_xxyz_1[i] * fe_0 + ta2_yz_z_xxy_0[i] * fe_0 - ta2_yz_z_xxy_1[i] * fe_0 +
                              ta1_y_z_xxyz_1[i] + ta2_yz_z_xxyz_0[i] * pa_z[i] - ta2_yz_z_xxyz_1[i] * pc_z[i];

        ta2_yz_zz_xxzz_0[i] = ta2_yz_0_xxzz_0[i] * fe_0 - ta2_yz_0_xxzz_1[i] * fe_0 + 2.0 * ta2_yz_z_xxz_0[i] * fe_0 -
                              2.0 * ta2_yz_z_xxz_1[i] * fe_0 + ta1_y_z_xxzz_1[i] + ta2_yz_z_xxzz_0[i] * pa_z[i] - ta2_yz_z_xxzz_1[i] * pc_z[i];

        ta2_yz_zz_xyyy_0[i] =
            ta2_yz_0_xyyy_0[i] * fe_0 - ta2_yz_0_xyyy_1[i] * fe_0 + ta1_y_z_xyyy_1[i] + ta2_yz_z_xyyy_0[i] * pa_z[i] - ta2_yz_z_xyyy_1[i] * pc_z[i];

        ta2_yz_zz_xyyz_0[i] = ta2_yz_0_xyyz_0[i] * fe_0 - ta2_yz_0_xyyz_1[i] * fe_0 + ta2_yz_z_xyy_0[i] * fe_0 - ta2_yz_z_xyy_1[i] * fe_0 +
                              ta1_y_z_xyyz_1[i] + ta2_yz_z_xyyz_0[i] * pa_z[i] - ta2_yz_z_xyyz_1[i] * pc_z[i];

        ta2_yz_zz_xyzz_0[i] = ta2_yz_0_xyzz_0[i] * fe_0 - ta2_yz_0_xyzz_1[i] * fe_0 + 2.0 * ta2_yz_z_xyz_0[i] * fe_0 -
                              2.0 * ta2_yz_z_xyz_1[i] * fe_0 + ta1_y_z_xyzz_1[i] + ta2_yz_z_xyzz_0[i] * pa_z[i] - ta2_yz_z_xyzz_1[i] * pc_z[i];

        ta2_yz_zz_xzzz_0[i] = ta2_yz_0_xzzz_0[i] * fe_0 - ta2_yz_0_xzzz_1[i] * fe_0 + 3.0 * ta2_yz_z_xzz_0[i] * fe_0 -
                              3.0 * ta2_yz_z_xzz_1[i] * fe_0 + ta1_y_z_xzzz_1[i] + ta2_yz_z_xzzz_0[i] * pa_z[i] - ta2_yz_z_xzzz_1[i] * pc_z[i];

        ta2_yz_zz_yyyy_0[i] =
            ta2_yz_0_yyyy_0[i] * fe_0 - ta2_yz_0_yyyy_1[i] * fe_0 + ta1_y_z_yyyy_1[i] + ta2_yz_z_yyyy_0[i] * pa_z[i] - ta2_yz_z_yyyy_1[i] * pc_z[i];

        ta2_yz_zz_yyyz_0[i] = ta2_yz_0_yyyz_0[i] * fe_0 - ta2_yz_0_yyyz_1[i] * fe_0 + ta2_yz_z_yyy_0[i] * fe_0 - ta2_yz_z_yyy_1[i] * fe_0 +
                              ta1_y_z_yyyz_1[i] + ta2_yz_z_yyyz_0[i] * pa_z[i] - ta2_yz_z_yyyz_1[i] * pc_z[i];

        ta2_yz_zz_yyzz_0[i] = ta2_yz_0_yyzz_0[i] * fe_0 - ta2_yz_0_yyzz_1[i] * fe_0 + 2.0 * ta2_yz_z_yyz_0[i] * fe_0 -
                              2.0 * ta2_yz_z_yyz_1[i] * fe_0 + ta1_y_z_yyzz_1[i] + ta2_yz_z_yyzz_0[i] * pa_z[i] - ta2_yz_z_yyzz_1[i] * pc_z[i];

        ta2_yz_zz_yzzz_0[i] = ta2_yz_0_yzzz_0[i] * fe_0 - ta2_yz_0_yzzz_1[i] * fe_0 + 3.0 * ta2_yz_z_yzz_0[i] * fe_0 -
                              3.0 * ta2_yz_z_yzz_1[i] * fe_0 + ta1_y_z_yzzz_1[i] + ta2_yz_z_yzzz_0[i] * pa_z[i] - ta2_yz_z_yzzz_1[i] * pc_z[i];

        ta2_yz_zz_zzzz_0[i] = ta2_yz_0_zzzz_0[i] * fe_0 - ta2_yz_0_zzzz_1[i] * fe_0 + 4.0 * ta2_yz_z_zzz_0[i] * fe_0 -
                              4.0 * ta2_yz_z_zzz_1[i] * fe_0 + ta1_y_z_zzzz_1[i] + ta2_yz_z_zzzz_0[i] * pa_z[i] - ta2_yz_z_zzzz_1[i] * pc_z[i];
    }

    // Set up 450-465 components of targeted buffer : DG

    auto ta2_zz_xx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 450);

    auto ta2_zz_xx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 451);

    auto ta2_zz_xx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 452);

    auto ta2_zz_xx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 453);

    auto ta2_zz_xx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 454);

    auto ta2_zz_xx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 455);

    auto ta2_zz_xx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 456);

    auto ta2_zz_xx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 457);

    auto ta2_zz_xx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 458);

    auto ta2_zz_xx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 459);

    auto ta2_zz_xx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 460);

    auto ta2_zz_xx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 461);

    auto ta2_zz_xx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 462);

    auto ta2_zz_xx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 463);

    auto ta2_zz_xx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 464);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta2_zz_0_xxxx_0,  \
                             ta2_zz_0_xxxx_1,  \
                             ta2_zz_0_xxxy_0,  \
                             ta2_zz_0_xxxy_1,  \
                             ta2_zz_0_xxxz_0,  \
                             ta2_zz_0_xxxz_1,  \
                             ta2_zz_0_xxyy_0,  \
                             ta2_zz_0_xxyy_1,  \
                             ta2_zz_0_xxyz_0,  \
                             ta2_zz_0_xxyz_1,  \
                             ta2_zz_0_xxzz_0,  \
                             ta2_zz_0_xxzz_1,  \
                             ta2_zz_0_xyyy_0,  \
                             ta2_zz_0_xyyy_1,  \
                             ta2_zz_0_xyyz_0,  \
                             ta2_zz_0_xyyz_1,  \
                             ta2_zz_0_xyzz_0,  \
                             ta2_zz_0_xyzz_1,  \
                             ta2_zz_0_xzzz_0,  \
                             ta2_zz_0_xzzz_1,  \
                             ta2_zz_0_yyyy_0,  \
                             ta2_zz_0_yyyy_1,  \
                             ta2_zz_0_yyyz_0,  \
                             ta2_zz_0_yyyz_1,  \
                             ta2_zz_0_yyzz_0,  \
                             ta2_zz_0_yyzz_1,  \
                             ta2_zz_0_yzzz_0,  \
                             ta2_zz_0_yzzz_1,  \
                             ta2_zz_0_zzzz_0,  \
                             ta2_zz_0_zzzz_1,  \
                             ta2_zz_x_xxx_0,   \
                             ta2_zz_x_xxx_1,   \
                             ta2_zz_x_xxxx_0,  \
                             ta2_zz_x_xxxx_1,  \
                             ta2_zz_x_xxxy_0,  \
                             ta2_zz_x_xxxy_1,  \
                             ta2_zz_x_xxxz_0,  \
                             ta2_zz_x_xxxz_1,  \
                             ta2_zz_x_xxy_0,   \
                             ta2_zz_x_xxy_1,   \
                             ta2_zz_x_xxyy_0,  \
                             ta2_zz_x_xxyy_1,  \
                             ta2_zz_x_xxyz_0,  \
                             ta2_zz_x_xxyz_1,  \
                             ta2_zz_x_xxz_0,   \
                             ta2_zz_x_xxz_1,   \
                             ta2_zz_x_xxzz_0,  \
                             ta2_zz_x_xxzz_1,  \
                             ta2_zz_x_xyy_0,   \
                             ta2_zz_x_xyy_1,   \
                             ta2_zz_x_xyyy_0,  \
                             ta2_zz_x_xyyy_1,  \
                             ta2_zz_x_xyyz_0,  \
                             ta2_zz_x_xyyz_1,  \
                             ta2_zz_x_xyz_0,   \
                             ta2_zz_x_xyz_1,   \
                             ta2_zz_x_xyzz_0,  \
                             ta2_zz_x_xyzz_1,  \
                             ta2_zz_x_xzz_0,   \
                             ta2_zz_x_xzz_1,   \
                             ta2_zz_x_xzzz_0,  \
                             ta2_zz_x_xzzz_1,  \
                             ta2_zz_x_yyy_0,   \
                             ta2_zz_x_yyy_1,   \
                             ta2_zz_x_yyyy_0,  \
                             ta2_zz_x_yyyy_1,  \
                             ta2_zz_x_yyyz_0,  \
                             ta2_zz_x_yyyz_1,  \
                             ta2_zz_x_yyz_0,   \
                             ta2_zz_x_yyz_1,   \
                             ta2_zz_x_yyzz_0,  \
                             ta2_zz_x_yyzz_1,  \
                             ta2_zz_x_yzz_0,   \
                             ta2_zz_x_yzz_1,   \
                             ta2_zz_x_yzzz_0,  \
                             ta2_zz_x_yzzz_1,  \
                             ta2_zz_x_zzz_0,   \
                             ta2_zz_x_zzz_1,   \
                             ta2_zz_x_zzzz_0,  \
                             ta2_zz_x_zzzz_1,  \
                             ta2_zz_xx_xxxx_0, \
                             ta2_zz_xx_xxxy_0, \
                             ta2_zz_xx_xxxz_0, \
                             ta2_zz_xx_xxyy_0, \
                             ta2_zz_xx_xxyz_0, \
                             ta2_zz_xx_xxzz_0, \
                             ta2_zz_xx_xyyy_0, \
                             ta2_zz_xx_xyyz_0, \
                             ta2_zz_xx_xyzz_0, \
                             ta2_zz_xx_xzzz_0, \
                             ta2_zz_xx_yyyy_0, \
                             ta2_zz_xx_yyyz_0, \
                             ta2_zz_xx_yyzz_0, \
                             ta2_zz_xx_yzzz_0, \
                             ta2_zz_xx_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xx_xxxx_0[i] = ta2_zz_0_xxxx_0[i] * fe_0 - ta2_zz_0_xxxx_1[i] * fe_0 + 4.0 * ta2_zz_x_xxx_0[i] * fe_0 -
                              4.0 * ta2_zz_x_xxx_1[i] * fe_0 + ta2_zz_x_xxxx_0[i] * pa_x[i] - ta2_zz_x_xxxx_1[i] * pc_x[i];

        ta2_zz_xx_xxxy_0[i] = ta2_zz_0_xxxy_0[i] * fe_0 - ta2_zz_0_xxxy_1[i] * fe_0 + 3.0 * ta2_zz_x_xxy_0[i] * fe_0 -
                              3.0 * ta2_zz_x_xxy_1[i] * fe_0 + ta2_zz_x_xxxy_0[i] * pa_x[i] - ta2_zz_x_xxxy_1[i] * pc_x[i];

        ta2_zz_xx_xxxz_0[i] = ta2_zz_0_xxxz_0[i] * fe_0 - ta2_zz_0_xxxz_1[i] * fe_0 + 3.0 * ta2_zz_x_xxz_0[i] * fe_0 -
                              3.0 * ta2_zz_x_xxz_1[i] * fe_0 + ta2_zz_x_xxxz_0[i] * pa_x[i] - ta2_zz_x_xxxz_1[i] * pc_x[i];

        ta2_zz_xx_xxyy_0[i] = ta2_zz_0_xxyy_0[i] * fe_0 - ta2_zz_0_xxyy_1[i] * fe_0 + 2.0 * ta2_zz_x_xyy_0[i] * fe_0 -
                              2.0 * ta2_zz_x_xyy_1[i] * fe_0 + ta2_zz_x_xxyy_0[i] * pa_x[i] - ta2_zz_x_xxyy_1[i] * pc_x[i];

        ta2_zz_xx_xxyz_0[i] = ta2_zz_0_xxyz_0[i] * fe_0 - ta2_zz_0_xxyz_1[i] * fe_0 + 2.0 * ta2_zz_x_xyz_0[i] * fe_0 -
                              2.0 * ta2_zz_x_xyz_1[i] * fe_0 + ta2_zz_x_xxyz_0[i] * pa_x[i] - ta2_zz_x_xxyz_1[i] * pc_x[i];

        ta2_zz_xx_xxzz_0[i] = ta2_zz_0_xxzz_0[i] * fe_0 - ta2_zz_0_xxzz_1[i] * fe_0 + 2.0 * ta2_zz_x_xzz_0[i] * fe_0 -
                              2.0 * ta2_zz_x_xzz_1[i] * fe_0 + ta2_zz_x_xxzz_0[i] * pa_x[i] - ta2_zz_x_xxzz_1[i] * pc_x[i];

        ta2_zz_xx_xyyy_0[i] = ta2_zz_0_xyyy_0[i] * fe_0 - ta2_zz_0_xyyy_1[i] * fe_0 + ta2_zz_x_yyy_0[i] * fe_0 - ta2_zz_x_yyy_1[i] * fe_0 +
                              ta2_zz_x_xyyy_0[i] * pa_x[i] - ta2_zz_x_xyyy_1[i] * pc_x[i];

        ta2_zz_xx_xyyz_0[i] = ta2_zz_0_xyyz_0[i] * fe_0 - ta2_zz_0_xyyz_1[i] * fe_0 + ta2_zz_x_yyz_0[i] * fe_0 - ta2_zz_x_yyz_1[i] * fe_0 +
                              ta2_zz_x_xyyz_0[i] * pa_x[i] - ta2_zz_x_xyyz_1[i] * pc_x[i];

        ta2_zz_xx_xyzz_0[i] = ta2_zz_0_xyzz_0[i] * fe_0 - ta2_zz_0_xyzz_1[i] * fe_0 + ta2_zz_x_yzz_0[i] * fe_0 - ta2_zz_x_yzz_1[i] * fe_0 +
                              ta2_zz_x_xyzz_0[i] * pa_x[i] - ta2_zz_x_xyzz_1[i] * pc_x[i];

        ta2_zz_xx_xzzz_0[i] = ta2_zz_0_xzzz_0[i] * fe_0 - ta2_zz_0_xzzz_1[i] * fe_0 + ta2_zz_x_zzz_0[i] * fe_0 - ta2_zz_x_zzz_1[i] * fe_0 +
                              ta2_zz_x_xzzz_0[i] * pa_x[i] - ta2_zz_x_xzzz_1[i] * pc_x[i];

        ta2_zz_xx_yyyy_0[i] = ta2_zz_0_yyyy_0[i] * fe_0 - ta2_zz_0_yyyy_1[i] * fe_0 + ta2_zz_x_yyyy_0[i] * pa_x[i] - ta2_zz_x_yyyy_1[i] * pc_x[i];

        ta2_zz_xx_yyyz_0[i] = ta2_zz_0_yyyz_0[i] * fe_0 - ta2_zz_0_yyyz_1[i] * fe_0 + ta2_zz_x_yyyz_0[i] * pa_x[i] - ta2_zz_x_yyyz_1[i] * pc_x[i];

        ta2_zz_xx_yyzz_0[i] = ta2_zz_0_yyzz_0[i] * fe_0 - ta2_zz_0_yyzz_1[i] * fe_0 + ta2_zz_x_yyzz_0[i] * pa_x[i] - ta2_zz_x_yyzz_1[i] * pc_x[i];

        ta2_zz_xx_yzzz_0[i] = ta2_zz_0_yzzz_0[i] * fe_0 - ta2_zz_0_yzzz_1[i] * fe_0 + ta2_zz_x_yzzz_0[i] * pa_x[i] - ta2_zz_x_yzzz_1[i] * pc_x[i];

        ta2_zz_xx_zzzz_0[i] = ta2_zz_0_zzzz_0[i] * fe_0 - ta2_zz_0_zzzz_1[i] * fe_0 + ta2_zz_x_zzzz_0[i] * pa_x[i] - ta2_zz_x_zzzz_1[i] * pc_x[i];
    }

    // Set up 465-480 components of targeted buffer : DG

    auto ta2_zz_xy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 465);

    auto ta2_zz_xy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 466);

    auto ta2_zz_xy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 467);

    auto ta2_zz_xy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 468);

    auto ta2_zz_xy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 469);

    auto ta2_zz_xy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 470);

    auto ta2_zz_xy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 471);

    auto ta2_zz_xy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 472);

    auto ta2_zz_xy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 473);

    auto ta2_zz_xy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 474);

    auto ta2_zz_xy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 475);

    auto ta2_zz_xy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 476);

    auto ta2_zz_xy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 477);

    auto ta2_zz_xy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 478);

    auto ta2_zz_xy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 479);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta2_zz_x_xxxx_0,  \
                             ta2_zz_x_xxxx_1,  \
                             ta2_zz_x_xxxz_0,  \
                             ta2_zz_x_xxxz_1,  \
                             ta2_zz_x_xxzz_0,  \
                             ta2_zz_x_xxzz_1,  \
                             ta2_zz_x_xzzz_0,  \
                             ta2_zz_x_xzzz_1,  \
                             ta2_zz_xy_xxxx_0, \
                             ta2_zz_xy_xxxy_0, \
                             ta2_zz_xy_xxxz_0, \
                             ta2_zz_xy_xxyy_0, \
                             ta2_zz_xy_xxyz_0, \
                             ta2_zz_xy_xxzz_0, \
                             ta2_zz_xy_xyyy_0, \
                             ta2_zz_xy_xyyz_0, \
                             ta2_zz_xy_xyzz_0, \
                             ta2_zz_xy_xzzz_0, \
                             ta2_zz_xy_yyyy_0, \
                             ta2_zz_xy_yyyz_0, \
                             ta2_zz_xy_yyzz_0, \
                             ta2_zz_xy_yzzz_0, \
                             ta2_zz_xy_zzzz_0, \
                             ta2_zz_y_xxxy_0,  \
                             ta2_zz_y_xxxy_1,  \
                             ta2_zz_y_xxy_0,   \
                             ta2_zz_y_xxy_1,   \
                             ta2_zz_y_xxyy_0,  \
                             ta2_zz_y_xxyy_1,  \
                             ta2_zz_y_xxyz_0,  \
                             ta2_zz_y_xxyz_1,  \
                             ta2_zz_y_xyy_0,   \
                             ta2_zz_y_xyy_1,   \
                             ta2_zz_y_xyyy_0,  \
                             ta2_zz_y_xyyy_1,  \
                             ta2_zz_y_xyyz_0,  \
                             ta2_zz_y_xyyz_1,  \
                             ta2_zz_y_xyz_0,   \
                             ta2_zz_y_xyz_1,   \
                             ta2_zz_y_xyzz_0,  \
                             ta2_zz_y_xyzz_1,  \
                             ta2_zz_y_yyy_0,   \
                             ta2_zz_y_yyy_1,   \
                             ta2_zz_y_yyyy_0,  \
                             ta2_zz_y_yyyy_1,  \
                             ta2_zz_y_yyyz_0,  \
                             ta2_zz_y_yyyz_1,  \
                             ta2_zz_y_yyz_0,   \
                             ta2_zz_y_yyz_1,   \
                             ta2_zz_y_yyzz_0,  \
                             ta2_zz_y_yyzz_1,  \
                             ta2_zz_y_yzz_0,   \
                             ta2_zz_y_yzz_1,   \
                             ta2_zz_y_yzzz_0,  \
                             ta2_zz_y_yzzz_1,  \
                             ta2_zz_y_zzzz_0,  \
                             ta2_zz_y_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xy_xxxx_0[i] = ta2_zz_x_xxxx_0[i] * pa_y[i] - ta2_zz_x_xxxx_1[i] * pc_y[i];

        ta2_zz_xy_xxxy_0[i] =
            3.0 * ta2_zz_y_xxy_0[i] * fe_0 - 3.0 * ta2_zz_y_xxy_1[i] * fe_0 + ta2_zz_y_xxxy_0[i] * pa_x[i] - ta2_zz_y_xxxy_1[i] * pc_x[i];

        ta2_zz_xy_xxxz_0[i] = ta2_zz_x_xxxz_0[i] * pa_y[i] - ta2_zz_x_xxxz_1[i] * pc_y[i];

        ta2_zz_xy_xxyy_0[i] =
            2.0 * ta2_zz_y_xyy_0[i] * fe_0 - 2.0 * ta2_zz_y_xyy_1[i] * fe_0 + ta2_zz_y_xxyy_0[i] * pa_x[i] - ta2_zz_y_xxyy_1[i] * pc_x[i];

        ta2_zz_xy_xxyz_0[i] =
            2.0 * ta2_zz_y_xyz_0[i] * fe_0 - 2.0 * ta2_zz_y_xyz_1[i] * fe_0 + ta2_zz_y_xxyz_0[i] * pa_x[i] - ta2_zz_y_xxyz_1[i] * pc_x[i];

        ta2_zz_xy_xxzz_0[i] = ta2_zz_x_xxzz_0[i] * pa_y[i] - ta2_zz_x_xxzz_1[i] * pc_y[i];

        ta2_zz_xy_xyyy_0[i] = ta2_zz_y_yyy_0[i] * fe_0 - ta2_zz_y_yyy_1[i] * fe_0 + ta2_zz_y_xyyy_0[i] * pa_x[i] - ta2_zz_y_xyyy_1[i] * pc_x[i];

        ta2_zz_xy_xyyz_0[i] = ta2_zz_y_yyz_0[i] * fe_0 - ta2_zz_y_yyz_1[i] * fe_0 + ta2_zz_y_xyyz_0[i] * pa_x[i] - ta2_zz_y_xyyz_1[i] * pc_x[i];

        ta2_zz_xy_xyzz_0[i] = ta2_zz_y_yzz_0[i] * fe_0 - ta2_zz_y_yzz_1[i] * fe_0 + ta2_zz_y_xyzz_0[i] * pa_x[i] - ta2_zz_y_xyzz_1[i] * pc_x[i];

        ta2_zz_xy_xzzz_0[i] = ta2_zz_x_xzzz_0[i] * pa_y[i] - ta2_zz_x_xzzz_1[i] * pc_y[i];

        ta2_zz_xy_yyyy_0[i] = ta2_zz_y_yyyy_0[i] * pa_x[i] - ta2_zz_y_yyyy_1[i] * pc_x[i];

        ta2_zz_xy_yyyz_0[i] = ta2_zz_y_yyyz_0[i] * pa_x[i] - ta2_zz_y_yyyz_1[i] * pc_x[i];

        ta2_zz_xy_yyzz_0[i] = ta2_zz_y_yyzz_0[i] * pa_x[i] - ta2_zz_y_yyzz_1[i] * pc_x[i];

        ta2_zz_xy_yzzz_0[i] = ta2_zz_y_yzzz_0[i] * pa_x[i] - ta2_zz_y_yzzz_1[i] * pc_x[i];

        ta2_zz_xy_zzzz_0[i] = ta2_zz_y_zzzz_0[i] * pa_x[i] - ta2_zz_y_zzzz_1[i] * pc_x[i];
    }

    // Set up 480-495 components of targeted buffer : DG

    auto ta2_zz_xz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 480);

    auto ta2_zz_xz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 481);

    auto ta2_zz_xz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 482);

    auto ta2_zz_xz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 483);

    auto ta2_zz_xz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 484);

    auto ta2_zz_xz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 485);

    auto ta2_zz_xz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 486);

    auto ta2_zz_xz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 487);

    auto ta2_zz_xz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 488);

    auto ta2_zz_xz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 489);

    auto ta2_zz_xz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 490);

    auto ta2_zz_xz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 491);

    auto ta2_zz_xz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 492);

    auto ta2_zz_xz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 493);

    auto ta2_zz_xz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 494);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_z_x_xxxx_1,   \
                             ta1_z_x_xxxy_1,   \
                             ta1_z_x_xxyy_1,   \
                             ta1_z_x_xyyy_1,   \
                             ta2_zz_x_xxxx_0,  \
                             ta2_zz_x_xxxx_1,  \
                             ta2_zz_x_xxxy_0,  \
                             ta2_zz_x_xxxy_1,  \
                             ta2_zz_x_xxyy_0,  \
                             ta2_zz_x_xxyy_1,  \
                             ta2_zz_x_xyyy_0,  \
                             ta2_zz_x_xyyy_1,  \
                             ta2_zz_xz_xxxx_0, \
                             ta2_zz_xz_xxxy_0, \
                             ta2_zz_xz_xxxz_0, \
                             ta2_zz_xz_xxyy_0, \
                             ta2_zz_xz_xxyz_0, \
                             ta2_zz_xz_xxzz_0, \
                             ta2_zz_xz_xyyy_0, \
                             ta2_zz_xz_xyyz_0, \
                             ta2_zz_xz_xyzz_0, \
                             ta2_zz_xz_xzzz_0, \
                             ta2_zz_xz_yyyy_0, \
                             ta2_zz_xz_yyyz_0, \
                             ta2_zz_xz_yyzz_0, \
                             ta2_zz_xz_yzzz_0, \
                             ta2_zz_xz_zzzz_0, \
                             ta2_zz_z_xxxz_0,  \
                             ta2_zz_z_xxxz_1,  \
                             ta2_zz_z_xxyz_0,  \
                             ta2_zz_z_xxyz_1,  \
                             ta2_zz_z_xxz_0,   \
                             ta2_zz_z_xxz_1,   \
                             ta2_zz_z_xxzz_0,  \
                             ta2_zz_z_xxzz_1,  \
                             ta2_zz_z_xyyz_0,  \
                             ta2_zz_z_xyyz_1,  \
                             ta2_zz_z_xyz_0,   \
                             ta2_zz_z_xyz_1,   \
                             ta2_zz_z_xyzz_0,  \
                             ta2_zz_z_xyzz_1,  \
                             ta2_zz_z_xzz_0,   \
                             ta2_zz_z_xzz_1,   \
                             ta2_zz_z_xzzz_0,  \
                             ta2_zz_z_xzzz_1,  \
                             ta2_zz_z_yyyy_0,  \
                             ta2_zz_z_yyyy_1,  \
                             ta2_zz_z_yyyz_0,  \
                             ta2_zz_z_yyyz_1,  \
                             ta2_zz_z_yyz_0,   \
                             ta2_zz_z_yyz_1,   \
                             ta2_zz_z_yyzz_0,  \
                             ta2_zz_z_yyzz_1,  \
                             ta2_zz_z_yzz_0,   \
                             ta2_zz_z_yzz_1,   \
                             ta2_zz_z_yzzz_0,  \
                             ta2_zz_z_yzzz_1,  \
                             ta2_zz_z_zzz_0,   \
                             ta2_zz_z_zzz_1,   \
                             ta2_zz_z_zzzz_0,  \
                             ta2_zz_z_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xz_xxxx_0[i] = 2.0 * ta1_z_x_xxxx_1[i] + ta2_zz_x_xxxx_0[i] * pa_z[i] - ta2_zz_x_xxxx_1[i] * pc_z[i];

        ta2_zz_xz_xxxy_0[i] = 2.0 * ta1_z_x_xxxy_1[i] + ta2_zz_x_xxxy_0[i] * pa_z[i] - ta2_zz_x_xxxy_1[i] * pc_z[i];

        ta2_zz_xz_xxxz_0[i] =
            3.0 * ta2_zz_z_xxz_0[i] * fe_0 - 3.0 * ta2_zz_z_xxz_1[i] * fe_0 + ta2_zz_z_xxxz_0[i] * pa_x[i] - ta2_zz_z_xxxz_1[i] * pc_x[i];

        ta2_zz_xz_xxyy_0[i] = 2.0 * ta1_z_x_xxyy_1[i] + ta2_zz_x_xxyy_0[i] * pa_z[i] - ta2_zz_x_xxyy_1[i] * pc_z[i];

        ta2_zz_xz_xxyz_0[i] =
            2.0 * ta2_zz_z_xyz_0[i] * fe_0 - 2.0 * ta2_zz_z_xyz_1[i] * fe_0 + ta2_zz_z_xxyz_0[i] * pa_x[i] - ta2_zz_z_xxyz_1[i] * pc_x[i];

        ta2_zz_xz_xxzz_0[i] =
            2.0 * ta2_zz_z_xzz_0[i] * fe_0 - 2.0 * ta2_zz_z_xzz_1[i] * fe_0 + ta2_zz_z_xxzz_0[i] * pa_x[i] - ta2_zz_z_xxzz_1[i] * pc_x[i];

        ta2_zz_xz_xyyy_0[i] = 2.0 * ta1_z_x_xyyy_1[i] + ta2_zz_x_xyyy_0[i] * pa_z[i] - ta2_zz_x_xyyy_1[i] * pc_z[i];

        ta2_zz_xz_xyyz_0[i] = ta2_zz_z_yyz_0[i] * fe_0 - ta2_zz_z_yyz_1[i] * fe_0 + ta2_zz_z_xyyz_0[i] * pa_x[i] - ta2_zz_z_xyyz_1[i] * pc_x[i];

        ta2_zz_xz_xyzz_0[i] = ta2_zz_z_yzz_0[i] * fe_0 - ta2_zz_z_yzz_1[i] * fe_0 + ta2_zz_z_xyzz_0[i] * pa_x[i] - ta2_zz_z_xyzz_1[i] * pc_x[i];

        ta2_zz_xz_xzzz_0[i] = ta2_zz_z_zzz_0[i] * fe_0 - ta2_zz_z_zzz_1[i] * fe_0 + ta2_zz_z_xzzz_0[i] * pa_x[i] - ta2_zz_z_xzzz_1[i] * pc_x[i];

        ta2_zz_xz_yyyy_0[i] = ta2_zz_z_yyyy_0[i] * pa_x[i] - ta2_zz_z_yyyy_1[i] * pc_x[i];

        ta2_zz_xz_yyyz_0[i] = ta2_zz_z_yyyz_0[i] * pa_x[i] - ta2_zz_z_yyyz_1[i] * pc_x[i];

        ta2_zz_xz_yyzz_0[i] = ta2_zz_z_yyzz_0[i] * pa_x[i] - ta2_zz_z_yyzz_1[i] * pc_x[i];

        ta2_zz_xz_yzzz_0[i] = ta2_zz_z_yzzz_0[i] * pa_x[i] - ta2_zz_z_yzzz_1[i] * pc_x[i];

        ta2_zz_xz_zzzz_0[i] = ta2_zz_z_zzzz_0[i] * pa_x[i] - ta2_zz_z_zzzz_1[i] * pc_x[i];
    }

    // Set up 495-510 components of targeted buffer : DG

    auto ta2_zz_yy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 495);

    auto ta2_zz_yy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 496);

    auto ta2_zz_yy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 497);

    auto ta2_zz_yy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 498);

    auto ta2_zz_yy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 499);

    auto ta2_zz_yy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 500);

    auto ta2_zz_yy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 501);

    auto ta2_zz_yy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 502);

    auto ta2_zz_yy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 503);

    auto ta2_zz_yy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 504);

    auto ta2_zz_yy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 505);

    auto ta2_zz_yy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 506);

    auto ta2_zz_yy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 507);

    auto ta2_zz_yy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 508);

    auto ta2_zz_yy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 509);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta2_zz_0_xxxx_0,  \
                             ta2_zz_0_xxxx_1,  \
                             ta2_zz_0_xxxy_0,  \
                             ta2_zz_0_xxxy_1,  \
                             ta2_zz_0_xxxz_0,  \
                             ta2_zz_0_xxxz_1,  \
                             ta2_zz_0_xxyy_0,  \
                             ta2_zz_0_xxyy_1,  \
                             ta2_zz_0_xxyz_0,  \
                             ta2_zz_0_xxyz_1,  \
                             ta2_zz_0_xxzz_0,  \
                             ta2_zz_0_xxzz_1,  \
                             ta2_zz_0_xyyy_0,  \
                             ta2_zz_0_xyyy_1,  \
                             ta2_zz_0_xyyz_0,  \
                             ta2_zz_0_xyyz_1,  \
                             ta2_zz_0_xyzz_0,  \
                             ta2_zz_0_xyzz_1,  \
                             ta2_zz_0_xzzz_0,  \
                             ta2_zz_0_xzzz_1,  \
                             ta2_zz_0_yyyy_0,  \
                             ta2_zz_0_yyyy_1,  \
                             ta2_zz_0_yyyz_0,  \
                             ta2_zz_0_yyyz_1,  \
                             ta2_zz_0_yyzz_0,  \
                             ta2_zz_0_yyzz_1,  \
                             ta2_zz_0_yzzz_0,  \
                             ta2_zz_0_yzzz_1,  \
                             ta2_zz_0_zzzz_0,  \
                             ta2_zz_0_zzzz_1,  \
                             ta2_zz_y_xxx_0,   \
                             ta2_zz_y_xxx_1,   \
                             ta2_zz_y_xxxx_0,  \
                             ta2_zz_y_xxxx_1,  \
                             ta2_zz_y_xxxy_0,  \
                             ta2_zz_y_xxxy_1,  \
                             ta2_zz_y_xxxz_0,  \
                             ta2_zz_y_xxxz_1,  \
                             ta2_zz_y_xxy_0,   \
                             ta2_zz_y_xxy_1,   \
                             ta2_zz_y_xxyy_0,  \
                             ta2_zz_y_xxyy_1,  \
                             ta2_zz_y_xxyz_0,  \
                             ta2_zz_y_xxyz_1,  \
                             ta2_zz_y_xxz_0,   \
                             ta2_zz_y_xxz_1,   \
                             ta2_zz_y_xxzz_0,  \
                             ta2_zz_y_xxzz_1,  \
                             ta2_zz_y_xyy_0,   \
                             ta2_zz_y_xyy_1,   \
                             ta2_zz_y_xyyy_0,  \
                             ta2_zz_y_xyyy_1,  \
                             ta2_zz_y_xyyz_0,  \
                             ta2_zz_y_xyyz_1,  \
                             ta2_zz_y_xyz_0,   \
                             ta2_zz_y_xyz_1,   \
                             ta2_zz_y_xyzz_0,  \
                             ta2_zz_y_xyzz_1,  \
                             ta2_zz_y_xzz_0,   \
                             ta2_zz_y_xzz_1,   \
                             ta2_zz_y_xzzz_0,  \
                             ta2_zz_y_xzzz_1,  \
                             ta2_zz_y_yyy_0,   \
                             ta2_zz_y_yyy_1,   \
                             ta2_zz_y_yyyy_0,  \
                             ta2_zz_y_yyyy_1,  \
                             ta2_zz_y_yyyz_0,  \
                             ta2_zz_y_yyyz_1,  \
                             ta2_zz_y_yyz_0,   \
                             ta2_zz_y_yyz_1,   \
                             ta2_zz_y_yyzz_0,  \
                             ta2_zz_y_yyzz_1,  \
                             ta2_zz_y_yzz_0,   \
                             ta2_zz_y_yzz_1,   \
                             ta2_zz_y_yzzz_0,  \
                             ta2_zz_y_yzzz_1,  \
                             ta2_zz_y_zzz_0,   \
                             ta2_zz_y_zzz_1,   \
                             ta2_zz_y_zzzz_0,  \
                             ta2_zz_y_zzzz_1,  \
                             ta2_zz_yy_xxxx_0, \
                             ta2_zz_yy_xxxy_0, \
                             ta2_zz_yy_xxxz_0, \
                             ta2_zz_yy_xxyy_0, \
                             ta2_zz_yy_xxyz_0, \
                             ta2_zz_yy_xxzz_0, \
                             ta2_zz_yy_xyyy_0, \
                             ta2_zz_yy_xyyz_0, \
                             ta2_zz_yy_xyzz_0, \
                             ta2_zz_yy_xzzz_0, \
                             ta2_zz_yy_yyyy_0, \
                             ta2_zz_yy_yyyz_0, \
                             ta2_zz_yy_yyzz_0, \
                             ta2_zz_yy_yzzz_0, \
                             ta2_zz_yy_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yy_xxxx_0[i] = ta2_zz_0_xxxx_0[i] * fe_0 - ta2_zz_0_xxxx_1[i] * fe_0 + ta2_zz_y_xxxx_0[i] * pa_y[i] - ta2_zz_y_xxxx_1[i] * pc_y[i];

        ta2_zz_yy_xxxy_0[i] = ta2_zz_0_xxxy_0[i] * fe_0 - ta2_zz_0_xxxy_1[i] * fe_0 + ta2_zz_y_xxx_0[i] * fe_0 - ta2_zz_y_xxx_1[i] * fe_0 +
                              ta2_zz_y_xxxy_0[i] * pa_y[i] - ta2_zz_y_xxxy_1[i] * pc_y[i];

        ta2_zz_yy_xxxz_0[i] = ta2_zz_0_xxxz_0[i] * fe_0 - ta2_zz_0_xxxz_1[i] * fe_0 + ta2_zz_y_xxxz_0[i] * pa_y[i] - ta2_zz_y_xxxz_1[i] * pc_y[i];

        ta2_zz_yy_xxyy_0[i] = ta2_zz_0_xxyy_0[i] * fe_0 - ta2_zz_0_xxyy_1[i] * fe_0 + 2.0 * ta2_zz_y_xxy_0[i] * fe_0 -
                              2.0 * ta2_zz_y_xxy_1[i] * fe_0 + ta2_zz_y_xxyy_0[i] * pa_y[i] - ta2_zz_y_xxyy_1[i] * pc_y[i];

        ta2_zz_yy_xxyz_0[i] = ta2_zz_0_xxyz_0[i] * fe_0 - ta2_zz_0_xxyz_1[i] * fe_0 + ta2_zz_y_xxz_0[i] * fe_0 - ta2_zz_y_xxz_1[i] * fe_0 +
                              ta2_zz_y_xxyz_0[i] * pa_y[i] - ta2_zz_y_xxyz_1[i] * pc_y[i];

        ta2_zz_yy_xxzz_0[i] = ta2_zz_0_xxzz_0[i] * fe_0 - ta2_zz_0_xxzz_1[i] * fe_0 + ta2_zz_y_xxzz_0[i] * pa_y[i] - ta2_zz_y_xxzz_1[i] * pc_y[i];

        ta2_zz_yy_xyyy_0[i] = ta2_zz_0_xyyy_0[i] * fe_0 - ta2_zz_0_xyyy_1[i] * fe_0 + 3.0 * ta2_zz_y_xyy_0[i] * fe_0 -
                              3.0 * ta2_zz_y_xyy_1[i] * fe_0 + ta2_zz_y_xyyy_0[i] * pa_y[i] - ta2_zz_y_xyyy_1[i] * pc_y[i];

        ta2_zz_yy_xyyz_0[i] = ta2_zz_0_xyyz_0[i] * fe_0 - ta2_zz_0_xyyz_1[i] * fe_0 + 2.0 * ta2_zz_y_xyz_0[i] * fe_0 -
                              2.0 * ta2_zz_y_xyz_1[i] * fe_0 + ta2_zz_y_xyyz_0[i] * pa_y[i] - ta2_zz_y_xyyz_1[i] * pc_y[i];

        ta2_zz_yy_xyzz_0[i] = ta2_zz_0_xyzz_0[i] * fe_0 - ta2_zz_0_xyzz_1[i] * fe_0 + ta2_zz_y_xzz_0[i] * fe_0 - ta2_zz_y_xzz_1[i] * fe_0 +
                              ta2_zz_y_xyzz_0[i] * pa_y[i] - ta2_zz_y_xyzz_1[i] * pc_y[i];

        ta2_zz_yy_xzzz_0[i] = ta2_zz_0_xzzz_0[i] * fe_0 - ta2_zz_0_xzzz_1[i] * fe_0 + ta2_zz_y_xzzz_0[i] * pa_y[i] - ta2_zz_y_xzzz_1[i] * pc_y[i];

        ta2_zz_yy_yyyy_0[i] = ta2_zz_0_yyyy_0[i] * fe_0 - ta2_zz_0_yyyy_1[i] * fe_0 + 4.0 * ta2_zz_y_yyy_0[i] * fe_0 -
                              4.0 * ta2_zz_y_yyy_1[i] * fe_0 + ta2_zz_y_yyyy_0[i] * pa_y[i] - ta2_zz_y_yyyy_1[i] * pc_y[i];

        ta2_zz_yy_yyyz_0[i] = ta2_zz_0_yyyz_0[i] * fe_0 - ta2_zz_0_yyyz_1[i] * fe_0 + 3.0 * ta2_zz_y_yyz_0[i] * fe_0 -
                              3.0 * ta2_zz_y_yyz_1[i] * fe_0 + ta2_zz_y_yyyz_0[i] * pa_y[i] - ta2_zz_y_yyyz_1[i] * pc_y[i];

        ta2_zz_yy_yyzz_0[i] = ta2_zz_0_yyzz_0[i] * fe_0 - ta2_zz_0_yyzz_1[i] * fe_0 + 2.0 * ta2_zz_y_yzz_0[i] * fe_0 -
                              2.0 * ta2_zz_y_yzz_1[i] * fe_0 + ta2_zz_y_yyzz_0[i] * pa_y[i] - ta2_zz_y_yyzz_1[i] * pc_y[i];

        ta2_zz_yy_yzzz_0[i] = ta2_zz_0_yzzz_0[i] * fe_0 - ta2_zz_0_yzzz_1[i] * fe_0 + ta2_zz_y_zzz_0[i] * fe_0 - ta2_zz_y_zzz_1[i] * fe_0 +
                              ta2_zz_y_yzzz_0[i] * pa_y[i] - ta2_zz_y_yzzz_1[i] * pc_y[i];

        ta2_zz_yy_zzzz_0[i] = ta2_zz_0_zzzz_0[i] * fe_0 - ta2_zz_0_zzzz_1[i] * fe_0 + ta2_zz_y_zzzz_0[i] * pa_y[i] - ta2_zz_y_zzzz_1[i] * pc_y[i];
    }

    // Set up 510-525 components of targeted buffer : DG

    auto ta2_zz_yz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 510);

    auto ta2_zz_yz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 511);

    auto ta2_zz_yz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 512);

    auto ta2_zz_yz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 513);

    auto ta2_zz_yz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 514);

    auto ta2_zz_yz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 515);

    auto ta2_zz_yz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 516);

    auto ta2_zz_yz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 517);

    auto ta2_zz_yz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 518);

    auto ta2_zz_yz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 519);

    auto ta2_zz_yz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 520);

    auto ta2_zz_yz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 521);

    auto ta2_zz_yz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 522);

    auto ta2_zz_yz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 523);

    auto ta2_zz_yz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 524);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_z_y_xxxy_1,   \
                             ta1_z_y_xxyy_1,   \
                             ta1_z_y_xyyy_1,   \
                             ta1_z_y_yyyy_1,   \
                             ta2_zz_y_xxxy_0,  \
                             ta2_zz_y_xxxy_1,  \
                             ta2_zz_y_xxyy_0,  \
                             ta2_zz_y_xxyy_1,  \
                             ta2_zz_y_xyyy_0,  \
                             ta2_zz_y_xyyy_1,  \
                             ta2_zz_y_yyyy_0,  \
                             ta2_zz_y_yyyy_1,  \
                             ta2_zz_yz_xxxx_0, \
                             ta2_zz_yz_xxxy_0, \
                             ta2_zz_yz_xxxz_0, \
                             ta2_zz_yz_xxyy_0, \
                             ta2_zz_yz_xxyz_0, \
                             ta2_zz_yz_xxzz_0, \
                             ta2_zz_yz_xyyy_0, \
                             ta2_zz_yz_xyyz_0, \
                             ta2_zz_yz_xyzz_0, \
                             ta2_zz_yz_xzzz_0, \
                             ta2_zz_yz_yyyy_0, \
                             ta2_zz_yz_yyyz_0, \
                             ta2_zz_yz_yyzz_0, \
                             ta2_zz_yz_yzzz_0, \
                             ta2_zz_yz_zzzz_0, \
                             ta2_zz_z_xxxx_0,  \
                             ta2_zz_z_xxxx_1,  \
                             ta2_zz_z_xxxz_0,  \
                             ta2_zz_z_xxxz_1,  \
                             ta2_zz_z_xxyz_0,  \
                             ta2_zz_z_xxyz_1,  \
                             ta2_zz_z_xxz_0,   \
                             ta2_zz_z_xxz_1,   \
                             ta2_zz_z_xxzz_0,  \
                             ta2_zz_z_xxzz_1,  \
                             ta2_zz_z_xyyz_0,  \
                             ta2_zz_z_xyyz_1,  \
                             ta2_zz_z_xyz_0,   \
                             ta2_zz_z_xyz_1,   \
                             ta2_zz_z_xyzz_0,  \
                             ta2_zz_z_xyzz_1,  \
                             ta2_zz_z_xzz_0,   \
                             ta2_zz_z_xzz_1,   \
                             ta2_zz_z_xzzz_0,  \
                             ta2_zz_z_xzzz_1,  \
                             ta2_zz_z_yyyz_0,  \
                             ta2_zz_z_yyyz_1,  \
                             ta2_zz_z_yyz_0,   \
                             ta2_zz_z_yyz_1,   \
                             ta2_zz_z_yyzz_0,  \
                             ta2_zz_z_yyzz_1,  \
                             ta2_zz_z_yzz_0,   \
                             ta2_zz_z_yzz_1,   \
                             ta2_zz_z_yzzz_0,  \
                             ta2_zz_z_yzzz_1,  \
                             ta2_zz_z_zzz_0,   \
                             ta2_zz_z_zzz_1,   \
                             ta2_zz_z_zzzz_0,  \
                             ta2_zz_z_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yz_xxxx_0[i] = ta2_zz_z_xxxx_0[i] * pa_y[i] - ta2_zz_z_xxxx_1[i] * pc_y[i];

        ta2_zz_yz_xxxy_0[i] = 2.0 * ta1_z_y_xxxy_1[i] + ta2_zz_y_xxxy_0[i] * pa_z[i] - ta2_zz_y_xxxy_1[i] * pc_z[i];

        ta2_zz_yz_xxxz_0[i] = ta2_zz_z_xxxz_0[i] * pa_y[i] - ta2_zz_z_xxxz_1[i] * pc_y[i];

        ta2_zz_yz_xxyy_0[i] = 2.0 * ta1_z_y_xxyy_1[i] + ta2_zz_y_xxyy_0[i] * pa_z[i] - ta2_zz_y_xxyy_1[i] * pc_z[i];

        ta2_zz_yz_xxyz_0[i] = ta2_zz_z_xxz_0[i] * fe_0 - ta2_zz_z_xxz_1[i] * fe_0 + ta2_zz_z_xxyz_0[i] * pa_y[i] - ta2_zz_z_xxyz_1[i] * pc_y[i];

        ta2_zz_yz_xxzz_0[i] = ta2_zz_z_xxzz_0[i] * pa_y[i] - ta2_zz_z_xxzz_1[i] * pc_y[i];

        ta2_zz_yz_xyyy_0[i] = 2.0 * ta1_z_y_xyyy_1[i] + ta2_zz_y_xyyy_0[i] * pa_z[i] - ta2_zz_y_xyyy_1[i] * pc_z[i];

        ta2_zz_yz_xyyz_0[i] =
            2.0 * ta2_zz_z_xyz_0[i] * fe_0 - 2.0 * ta2_zz_z_xyz_1[i] * fe_0 + ta2_zz_z_xyyz_0[i] * pa_y[i] - ta2_zz_z_xyyz_1[i] * pc_y[i];

        ta2_zz_yz_xyzz_0[i] = ta2_zz_z_xzz_0[i] * fe_0 - ta2_zz_z_xzz_1[i] * fe_0 + ta2_zz_z_xyzz_0[i] * pa_y[i] - ta2_zz_z_xyzz_1[i] * pc_y[i];

        ta2_zz_yz_xzzz_0[i] = ta2_zz_z_xzzz_0[i] * pa_y[i] - ta2_zz_z_xzzz_1[i] * pc_y[i];

        ta2_zz_yz_yyyy_0[i] = 2.0 * ta1_z_y_yyyy_1[i] + ta2_zz_y_yyyy_0[i] * pa_z[i] - ta2_zz_y_yyyy_1[i] * pc_z[i];

        ta2_zz_yz_yyyz_0[i] =
            3.0 * ta2_zz_z_yyz_0[i] * fe_0 - 3.0 * ta2_zz_z_yyz_1[i] * fe_0 + ta2_zz_z_yyyz_0[i] * pa_y[i] - ta2_zz_z_yyyz_1[i] * pc_y[i];

        ta2_zz_yz_yyzz_0[i] =
            2.0 * ta2_zz_z_yzz_0[i] * fe_0 - 2.0 * ta2_zz_z_yzz_1[i] * fe_0 + ta2_zz_z_yyzz_0[i] * pa_y[i] - ta2_zz_z_yyzz_1[i] * pc_y[i];

        ta2_zz_yz_yzzz_0[i] = ta2_zz_z_zzz_0[i] * fe_0 - ta2_zz_z_zzz_1[i] * fe_0 + ta2_zz_z_yzzz_0[i] * pa_y[i] - ta2_zz_z_yzzz_1[i] * pc_y[i];

        ta2_zz_yz_zzzz_0[i] = ta2_zz_z_zzzz_0[i] * pa_y[i] - ta2_zz_z_zzzz_1[i] * pc_y[i];
    }

    // Set up 525-540 components of targeted buffer : DG

    auto ta2_zz_zz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 525);

    auto ta2_zz_zz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 526);

    auto ta2_zz_zz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 527);

    auto ta2_zz_zz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 528);

    auto ta2_zz_zz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 529);

    auto ta2_zz_zz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 530);

    auto ta2_zz_zz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 531);

    auto ta2_zz_zz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 532);

    auto ta2_zz_zz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 533);

    auto ta2_zz_zz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 534);

    auto ta2_zz_zz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 535);

    auto ta2_zz_zz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 536);

    auto ta2_zz_zz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 537);

    auto ta2_zz_zz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 538);

    auto ta2_zz_zz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 539);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta1_z_z_xxxx_1,   \
                             ta1_z_z_xxxy_1,   \
                             ta1_z_z_xxxz_1,   \
                             ta1_z_z_xxyy_1,   \
                             ta1_z_z_xxyz_1,   \
                             ta1_z_z_xxzz_1,   \
                             ta1_z_z_xyyy_1,   \
                             ta1_z_z_xyyz_1,   \
                             ta1_z_z_xyzz_1,   \
                             ta1_z_z_xzzz_1,   \
                             ta1_z_z_yyyy_1,   \
                             ta1_z_z_yyyz_1,   \
                             ta1_z_z_yyzz_1,   \
                             ta1_z_z_yzzz_1,   \
                             ta1_z_z_zzzz_1,   \
                             ta2_zz_0_xxxx_0,  \
                             ta2_zz_0_xxxx_1,  \
                             ta2_zz_0_xxxy_0,  \
                             ta2_zz_0_xxxy_1,  \
                             ta2_zz_0_xxxz_0,  \
                             ta2_zz_0_xxxz_1,  \
                             ta2_zz_0_xxyy_0,  \
                             ta2_zz_0_xxyy_1,  \
                             ta2_zz_0_xxyz_0,  \
                             ta2_zz_0_xxyz_1,  \
                             ta2_zz_0_xxzz_0,  \
                             ta2_zz_0_xxzz_1,  \
                             ta2_zz_0_xyyy_0,  \
                             ta2_zz_0_xyyy_1,  \
                             ta2_zz_0_xyyz_0,  \
                             ta2_zz_0_xyyz_1,  \
                             ta2_zz_0_xyzz_0,  \
                             ta2_zz_0_xyzz_1,  \
                             ta2_zz_0_xzzz_0,  \
                             ta2_zz_0_xzzz_1,  \
                             ta2_zz_0_yyyy_0,  \
                             ta2_zz_0_yyyy_1,  \
                             ta2_zz_0_yyyz_0,  \
                             ta2_zz_0_yyyz_1,  \
                             ta2_zz_0_yyzz_0,  \
                             ta2_zz_0_yyzz_1,  \
                             ta2_zz_0_yzzz_0,  \
                             ta2_zz_0_yzzz_1,  \
                             ta2_zz_0_zzzz_0,  \
                             ta2_zz_0_zzzz_1,  \
                             ta2_zz_z_xxx_0,   \
                             ta2_zz_z_xxx_1,   \
                             ta2_zz_z_xxxx_0,  \
                             ta2_zz_z_xxxx_1,  \
                             ta2_zz_z_xxxy_0,  \
                             ta2_zz_z_xxxy_1,  \
                             ta2_zz_z_xxxz_0,  \
                             ta2_zz_z_xxxz_1,  \
                             ta2_zz_z_xxy_0,   \
                             ta2_zz_z_xxy_1,   \
                             ta2_zz_z_xxyy_0,  \
                             ta2_zz_z_xxyy_1,  \
                             ta2_zz_z_xxyz_0,  \
                             ta2_zz_z_xxyz_1,  \
                             ta2_zz_z_xxz_0,   \
                             ta2_zz_z_xxz_1,   \
                             ta2_zz_z_xxzz_0,  \
                             ta2_zz_z_xxzz_1,  \
                             ta2_zz_z_xyy_0,   \
                             ta2_zz_z_xyy_1,   \
                             ta2_zz_z_xyyy_0,  \
                             ta2_zz_z_xyyy_1,  \
                             ta2_zz_z_xyyz_0,  \
                             ta2_zz_z_xyyz_1,  \
                             ta2_zz_z_xyz_0,   \
                             ta2_zz_z_xyz_1,   \
                             ta2_zz_z_xyzz_0,  \
                             ta2_zz_z_xyzz_1,  \
                             ta2_zz_z_xzz_0,   \
                             ta2_zz_z_xzz_1,   \
                             ta2_zz_z_xzzz_0,  \
                             ta2_zz_z_xzzz_1,  \
                             ta2_zz_z_yyy_0,   \
                             ta2_zz_z_yyy_1,   \
                             ta2_zz_z_yyyy_0,  \
                             ta2_zz_z_yyyy_1,  \
                             ta2_zz_z_yyyz_0,  \
                             ta2_zz_z_yyyz_1,  \
                             ta2_zz_z_yyz_0,   \
                             ta2_zz_z_yyz_1,   \
                             ta2_zz_z_yyzz_0,  \
                             ta2_zz_z_yyzz_1,  \
                             ta2_zz_z_yzz_0,   \
                             ta2_zz_z_yzz_1,   \
                             ta2_zz_z_yzzz_0,  \
                             ta2_zz_z_yzzz_1,  \
                             ta2_zz_z_zzz_0,   \
                             ta2_zz_z_zzz_1,   \
                             ta2_zz_z_zzzz_0,  \
                             ta2_zz_z_zzzz_1,  \
                             ta2_zz_zz_xxxx_0, \
                             ta2_zz_zz_xxxy_0, \
                             ta2_zz_zz_xxxz_0, \
                             ta2_zz_zz_xxyy_0, \
                             ta2_zz_zz_xxyz_0, \
                             ta2_zz_zz_xxzz_0, \
                             ta2_zz_zz_xyyy_0, \
                             ta2_zz_zz_xyyz_0, \
                             ta2_zz_zz_xyzz_0, \
                             ta2_zz_zz_xzzz_0, \
                             ta2_zz_zz_yyyy_0, \
                             ta2_zz_zz_yyyz_0, \
                             ta2_zz_zz_yyzz_0, \
                             ta2_zz_zz_yzzz_0, \
                             ta2_zz_zz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_zz_xxxx_0[i] = ta2_zz_0_xxxx_0[i] * fe_0 - ta2_zz_0_xxxx_1[i] * fe_0 + 2.0 * ta1_z_z_xxxx_1[i] + ta2_zz_z_xxxx_0[i] * pa_z[i] -
                              ta2_zz_z_xxxx_1[i] * pc_z[i];

        ta2_zz_zz_xxxy_0[i] = ta2_zz_0_xxxy_0[i] * fe_0 - ta2_zz_0_xxxy_1[i] * fe_0 + 2.0 * ta1_z_z_xxxy_1[i] + ta2_zz_z_xxxy_0[i] * pa_z[i] -
                              ta2_zz_z_xxxy_1[i] * pc_z[i];

        ta2_zz_zz_xxxz_0[i] = ta2_zz_0_xxxz_0[i] * fe_0 - ta2_zz_0_xxxz_1[i] * fe_0 + ta2_zz_z_xxx_0[i] * fe_0 - ta2_zz_z_xxx_1[i] * fe_0 +
                              2.0 * ta1_z_z_xxxz_1[i] + ta2_zz_z_xxxz_0[i] * pa_z[i] - ta2_zz_z_xxxz_1[i] * pc_z[i];

        ta2_zz_zz_xxyy_0[i] = ta2_zz_0_xxyy_0[i] * fe_0 - ta2_zz_0_xxyy_1[i] * fe_0 + 2.0 * ta1_z_z_xxyy_1[i] + ta2_zz_z_xxyy_0[i] * pa_z[i] -
                              ta2_zz_z_xxyy_1[i] * pc_z[i];

        ta2_zz_zz_xxyz_0[i] = ta2_zz_0_xxyz_0[i] * fe_0 - ta2_zz_0_xxyz_1[i] * fe_0 + ta2_zz_z_xxy_0[i] * fe_0 - ta2_zz_z_xxy_1[i] * fe_0 +
                              2.0 * ta1_z_z_xxyz_1[i] + ta2_zz_z_xxyz_0[i] * pa_z[i] - ta2_zz_z_xxyz_1[i] * pc_z[i];

        ta2_zz_zz_xxzz_0[i] = ta2_zz_0_xxzz_0[i] * fe_0 - ta2_zz_0_xxzz_1[i] * fe_0 + 2.0 * ta2_zz_z_xxz_0[i] * fe_0 -
                              2.0 * ta2_zz_z_xxz_1[i] * fe_0 + 2.0 * ta1_z_z_xxzz_1[i] + ta2_zz_z_xxzz_0[i] * pa_z[i] - ta2_zz_z_xxzz_1[i] * pc_z[i];

        ta2_zz_zz_xyyy_0[i] = ta2_zz_0_xyyy_0[i] * fe_0 - ta2_zz_0_xyyy_1[i] * fe_0 + 2.0 * ta1_z_z_xyyy_1[i] + ta2_zz_z_xyyy_0[i] * pa_z[i] -
                              ta2_zz_z_xyyy_1[i] * pc_z[i];

        ta2_zz_zz_xyyz_0[i] = ta2_zz_0_xyyz_0[i] * fe_0 - ta2_zz_0_xyyz_1[i] * fe_0 + ta2_zz_z_xyy_0[i] * fe_0 - ta2_zz_z_xyy_1[i] * fe_0 +
                              2.0 * ta1_z_z_xyyz_1[i] + ta2_zz_z_xyyz_0[i] * pa_z[i] - ta2_zz_z_xyyz_1[i] * pc_z[i];

        ta2_zz_zz_xyzz_0[i] = ta2_zz_0_xyzz_0[i] * fe_0 - ta2_zz_0_xyzz_1[i] * fe_0 + 2.0 * ta2_zz_z_xyz_0[i] * fe_0 -
                              2.0 * ta2_zz_z_xyz_1[i] * fe_0 + 2.0 * ta1_z_z_xyzz_1[i] + ta2_zz_z_xyzz_0[i] * pa_z[i] - ta2_zz_z_xyzz_1[i] * pc_z[i];

        ta2_zz_zz_xzzz_0[i] = ta2_zz_0_xzzz_0[i] * fe_0 - ta2_zz_0_xzzz_1[i] * fe_0 + 3.0 * ta2_zz_z_xzz_0[i] * fe_0 -
                              3.0 * ta2_zz_z_xzz_1[i] * fe_0 + 2.0 * ta1_z_z_xzzz_1[i] + ta2_zz_z_xzzz_0[i] * pa_z[i] - ta2_zz_z_xzzz_1[i] * pc_z[i];

        ta2_zz_zz_yyyy_0[i] = ta2_zz_0_yyyy_0[i] * fe_0 - ta2_zz_0_yyyy_1[i] * fe_0 + 2.0 * ta1_z_z_yyyy_1[i] + ta2_zz_z_yyyy_0[i] * pa_z[i] -
                              ta2_zz_z_yyyy_1[i] * pc_z[i];

        ta2_zz_zz_yyyz_0[i] = ta2_zz_0_yyyz_0[i] * fe_0 - ta2_zz_0_yyyz_1[i] * fe_0 + ta2_zz_z_yyy_0[i] * fe_0 - ta2_zz_z_yyy_1[i] * fe_0 +
                              2.0 * ta1_z_z_yyyz_1[i] + ta2_zz_z_yyyz_0[i] * pa_z[i] - ta2_zz_z_yyyz_1[i] * pc_z[i];

        ta2_zz_zz_yyzz_0[i] = ta2_zz_0_yyzz_0[i] * fe_0 - ta2_zz_0_yyzz_1[i] * fe_0 + 2.0 * ta2_zz_z_yyz_0[i] * fe_0 -
                              2.0 * ta2_zz_z_yyz_1[i] * fe_0 + 2.0 * ta1_z_z_yyzz_1[i] + ta2_zz_z_yyzz_0[i] * pa_z[i] - ta2_zz_z_yyzz_1[i] * pc_z[i];

        ta2_zz_zz_yzzz_0[i] = ta2_zz_0_yzzz_0[i] * fe_0 - ta2_zz_0_yzzz_1[i] * fe_0 + 3.0 * ta2_zz_z_yzz_0[i] * fe_0 -
                              3.0 * ta2_zz_z_yzz_1[i] * fe_0 + 2.0 * ta1_z_z_yzzz_1[i] + ta2_zz_z_yzzz_0[i] * pa_z[i] - ta2_zz_z_yzzz_1[i] * pc_z[i];

        ta2_zz_zz_zzzz_0[i] = ta2_zz_0_zzzz_0[i] * fe_0 - ta2_zz_0_zzzz_1[i] * fe_0 + 4.0 * ta2_zz_z_zzz_0[i] * fe_0 -
                              4.0 * ta2_zz_z_zzz_1[i] * fe_0 + 2.0 * ta1_z_z_zzzz_1[i] + ta2_zz_z_zzzz_0[i] * pa_z[i] - ta2_zz_z_zzzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
