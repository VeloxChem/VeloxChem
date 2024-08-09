#include "NuclearPotentialGeom010PrimRecFG.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_fg(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_fg,
                                        const size_t              idx_npot_geom_010_0_pg,
                                        const size_t              idx_npot_geom_010_1_pg,
                                        const size_t              idx_npot_geom_010_0_df,
                                        const size_t              idx_npot_geom_010_1_df,
                                        const size_t              idx_npot_1_dg,
                                        const size_t              idx_npot_geom_010_0_dg,
                                        const size_t              idx_npot_geom_010_1_dg,
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

    // Set up components of auxiliary buffer : PG

    auto ta1_x_x_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_pg);

    auto ta1_x_x_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 1);

    auto ta1_x_x_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 2);

    auto ta1_x_x_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 3);

    auto ta1_x_x_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 4);

    auto ta1_x_x_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 5);

    auto ta1_x_x_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 6);

    auto ta1_x_x_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 7);

    auto ta1_x_x_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 8);

    auto ta1_x_x_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 9);

    auto ta1_x_x_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 10);

    auto ta1_x_x_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 11);

    auto ta1_x_x_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 12);

    auto ta1_x_x_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 13);

    auto ta1_x_x_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 14);

    auto ta1_x_y_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_pg + 15);

    auto ta1_x_y_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 16);

    auto ta1_x_y_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 17);

    auto ta1_x_y_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 18);

    auto ta1_x_y_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 19);

    auto ta1_x_y_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 20);

    auto ta1_x_y_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 21);

    auto ta1_x_y_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 22);

    auto ta1_x_y_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 23);

    auto ta1_x_y_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 24);

    auto ta1_x_y_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 25);

    auto ta1_x_y_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 26);

    auto ta1_x_y_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 27);

    auto ta1_x_y_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 28);

    auto ta1_x_y_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 29);

    auto ta1_x_z_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_pg + 30);

    auto ta1_x_z_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 31);

    auto ta1_x_z_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 32);

    auto ta1_x_z_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 33);

    auto ta1_x_z_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 34);

    auto ta1_x_z_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 35);

    auto ta1_x_z_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 36);

    auto ta1_x_z_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 37);

    auto ta1_x_z_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 38);

    auto ta1_x_z_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 39);

    auto ta1_x_z_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 40);

    auto ta1_x_z_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 41);

    auto ta1_x_z_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 42);

    auto ta1_x_z_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 43);

    auto ta1_x_z_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 44);

    auto ta1_y_x_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_pg + 45);

    auto ta1_y_x_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 46);

    auto ta1_y_x_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 47);

    auto ta1_y_x_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 48);

    auto ta1_y_x_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 49);

    auto ta1_y_x_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 50);

    auto ta1_y_x_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 51);

    auto ta1_y_x_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 52);

    auto ta1_y_x_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 53);

    auto ta1_y_x_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 54);

    auto ta1_y_x_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 55);

    auto ta1_y_x_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 56);

    auto ta1_y_x_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 57);

    auto ta1_y_x_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 58);

    auto ta1_y_x_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 59);

    auto ta1_y_y_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_pg + 60);

    auto ta1_y_y_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 61);

    auto ta1_y_y_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 62);

    auto ta1_y_y_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 63);

    auto ta1_y_y_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 64);

    auto ta1_y_y_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 65);

    auto ta1_y_y_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 66);

    auto ta1_y_y_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 67);

    auto ta1_y_y_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 68);

    auto ta1_y_y_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 69);

    auto ta1_y_y_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 70);

    auto ta1_y_y_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 71);

    auto ta1_y_y_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 72);

    auto ta1_y_y_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 73);

    auto ta1_y_y_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 74);

    auto ta1_y_z_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_pg + 75);

    auto ta1_y_z_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 76);

    auto ta1_y_z_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 77);

    auto ta1_y_z_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 78);

    auto ta1_y_z_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 79);

    auto ta1_y_z_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 80);

    auto ta1_y_z_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 81);

    auto ta1_y_z_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 82);

    auto ta1_y_z_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 83);

    auto ta1_y_z_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 84);

    auto ta1_y_z_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 85);

    auto ta1_y_z_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 86);

    auto ta1_y_z_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 87);

    auto ta1_y_z_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 88);

    auto ta1_y_z_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 89);

    auto ta1_z_x_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_pg + 90);

    auto ta1_z_x_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 91);

    auto ta1_z_x_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 92);

    auto ta1_z_x_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 93);

    auto ta1_z_x_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 94);

    auto ta1_z_x_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 95);

    auto ta1_z_x_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 96);

    auto ta1_z_x_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 97);

    auto ta1_z_x_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 98);

    auto ta1_z_x_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 99);

    auto ta1_z_x_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 100);

    auto ta1_z_x_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 101);

    auto ta1_z_x_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 102);

    auto ta1_z_x_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 103);

    auto ta1_z_x_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 104);

    auto ta1_z_y_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_pg + 105);

    auto ta1_z_y_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 106);

    auto ta1_z_y_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 107);

    auto ta1_z_y_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 108);

    auto ta1_z_y_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 109);

    auto ta1_z_y_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 110);

    auto ta1_z_y_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 111);

    auto ta1_z_y_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 112);

    auto ta1_z_y_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 113);

    auto ta1_z_y_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 114);

    auto ta1_z_y_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 115);

    auto ta1_z_y_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 116);

    auto ta1_z_y_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 117);

    auto ta1_z_y_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 118);

    auto ta1_z_y_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 119);

    auto ta1_z_z_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_pg + 120);

    auto ta1_z_z_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 121);

    auto ta1_z_z_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 122);

    auto ta1_z_z_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 123);

    auto ta1_z_z_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 124);

    auto ta1_z_z_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 125);

    auto ta1_z_z_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 126);

    auto ta1_z_z_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 127);

    auto ta1_z_z_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 128);

    auto ta1_z_z_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 129);

    auto ta1_z_z_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 130);

    auto ta1_z_z_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 131);

    auto ta1_z_z_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 132);

    auto ta1_z_z_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 133);

    auto ta1_z_z_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 134);

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

    // Set up components of auxiliary buffer : DF

    auto ta1_x_xx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df);

    auto ta1_x_xx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 1);

    auto ta1_x_xx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 2);

    auto ta1_x_xx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 3);

    auto ta1_x_xx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 4);

    auto ta1_x_xx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 5);

    auto ta1_x_xx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 6);

    auto ta1_x_xx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 7);

    auto ta1_x_xx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 8);

    auto ta1_x_xx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 9);

    auto ta1_x_xz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 22);

    auto ta1_x_xz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 24);

    auto ta1_x_xz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 25);

    auto ta1_x_yy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 30);

    auto ta1_x_yy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 31);

    auto ta1_x_yy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 32);

    auto ta1_x_yy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 33);

    auto ta1_x_yy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 34);

    auto ta1_x_yy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 35);

    auto ta1_x_yy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 36);

    auto ta1_x_yy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 37);

    auto ta1_x_yy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 38);

    auto ta1_x_yy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 39);

    auto ta1_x_zz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 50);

    auto ta1_x_zz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 51);

    auto ta1_x_zz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 52);

    auto ta1_x_zz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 53);

    auto ta1_x_zz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 54);

    auto ta1_x_zz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 55);

    auto ta1_x_zz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 56);

    auto ta1_x_zz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 57);

    auto ta1_x_zz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 58);

    auto ta1_x_zz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 59);

    auto ta1_y_xx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 60);

    auto ta1_y_xx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 61);

    auto ta1_y_xx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 62);

    auto ta1_y_xx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 63);

    auto ta1_y_xx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 64);

    auto ta1_y_xx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 65);

    auto ta1_y_xx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 66);

    auto ta1_y_xx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 67);

    auto ta1_y_xx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 68);

    auto ta1_y_xx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 69);

    auto ta1_y_yy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 90);

    auto ta1_y_yy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 91);

    auto ta1_y_yy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 92);

    auto ta1_y_yy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 93);

    auto ta1_y_yy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 94);

    auto ta1_y_yy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 95);

    auto ta1_y_yy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 96);

    auto ta1_y_yy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 97);

    auto ta1_y_yy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 98);

    auto ta1_y_yy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 99);

    auto ta1_y_yz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 104);

    auto ta1_y_yz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 107);

    auto ta1_y_yz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 108);

    auto ta1_y_zz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 110);

    auto ta1_y_zz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 111);

    auto ta1_y_zz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 112);

    auto ta1_y_zz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 113);

    auto ta1_y_zz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 114);

    auto ta1_y_zz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 115);

    auto ta1_y_zz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 116);

    auto ta1_y_zz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 117);

    auto ta1_y_zz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 118);

    auto ta1_y_zz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 119);

    auto ta1_z_xx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 120);

    auto ta1_z_xx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 121);

    auto ta1_z_xx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 122);

    auto ta1_z_xx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 123);

    auto ta1_z_xx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 124);

    auto ta1_z_xx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 125);

    auto ta1_z_xx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 126);

    auto ta1_z_xx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 127);

    auto ta1_z_xx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 128);

    auto ta1_z_xx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 129);

    auto ta1_z_yy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 150);

    auto ta1_z_yy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 151);

    auto ta1_z_yy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 152);

    auto ta1_z_yy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 153);

    auto ta1_z_yy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 154);

    auto ta1_z_yy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 155);

    auto ta1_z_yy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 156);

    auto ta1_z_yy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 157);

    auto ta1_z_yy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 158);

    auto ta1_z_yy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 159);

    auto ta1_z_yz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 164);

    auto ta1_z_yz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 167);

    auto ta1_z_yz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 168);

    auto ta1_z_zz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 170);

    auto ta1_z_zz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 171);

    auto ta1_z_zz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 172);

    auto ta1_z_zz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 173);

    auto ta1_z_zz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 174);

    auto ta1_z_zz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 175);

    auto ta1_z_zz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 176);

    auto ta1_z_zz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 177);

    auto ta1_z_zz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 178);

    auto ta1_z_zz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 179);

    // Set up components of auxiliary buffer : DF

    auto ta1_x_xx_xxx_1 = pbuffer.data(idx_npot_geom_010_1_df);

    auto ta1_x_xx_xxy_1 = pbuffer.data(idx_npot_geom_010_1_df + 1);

    auto ta1_x_xx_xxz_1 = pbuffer.data(idx_npot_geom_010_1_df + 2);

    auto ta1_x_xx_xyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 3);

    auto ta1_x_xx_xyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 4);

    auto ta1_x_xx_xzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 5);

    auto ta1_x_xx_yyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 6);

    auto ta1_x_xx_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 7);

    auto ta1_x_xx_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 8);

    auto ta1_x_xx_zzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 9);

    auto ta1_x_xz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_df + 22);

    auto ta1_x_xz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 24);

    auto ta1_x_xz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 25);

    auto ta1_x_yy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_df + 30);

    auto ta1_x_yy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_df + 31);

    auto ta1_x_yy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_df + 32);

    auto ta1_x_yy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 33);

    auto ta1_x_yy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 34);

    auto ta1_x_yy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 35);

    auto ta1_x_yy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 36);

    auto ta1_x_yy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 37);

    auto ta1_x_yy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 38);

    auto ta1_x_yy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 39);

    auto ta1_x_zz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_df + 50);

    auto ta1_x_zz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_df + 51);

    auto ta1_x_zz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_df + 52);

    auto ta1_x_zz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 53);

    auto ta1_x_zz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 54);

    auto ta1_x_zz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 55);

    auto ta1_x_zz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 56);

    auto ta1_x_zz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 57);

    auto ta1_x_zz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 58);

    auto ta1_x_zz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 59);

    auto ta1_y_xx_xxx_1 = pbuffer.data(idx_npot_geom_010_1_df + 60);

    auto ta1_y_xx_xxy_1 = pbuffer.data(idx_npot_geom_010_1_df + 61);

    auto ta1_y_xx_xxz_1 = pbuffer.data(idx_npot_geom_010_1_df + 62);

    auto ta1_y_xx_xyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 63);

    auto ta1_y_xx_xyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 64);

    auto ta1_y_xx_xzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 65);

    auto ta1_y_xx_yyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 66);

    auto ta1_y_xx_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 67);

    auto ta1_y_xx_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 68);

    auto ta1_y_xx_zzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 69);

    auto ta1_y_yy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_df + 90);

    auto ta1_y_yy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_df + 91);

    auto ta1_y_yy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_df + 92);

    auto ta1_y_yy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 93);

    auto ta1_y_yy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 94);

    auto ta1_y_yy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 95);

    auto ta1_y_yy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 96);

    auto ta1_y_yy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 97);

    auto ta1_y_yy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 98);

    auto ta1_y_yy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 99);

    auto ta1_y_yz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 104);

    auto ta1_y_yz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 107);

    auto ta1_y_yz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 108);

    auto ta1_y_zz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_df + 110);

    auto ta1_y_zz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_df + 111);

    auto ta1_y_zz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_df + 112);

    auto ta1_y_zz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 113);

    auto ta1_y_zz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 114);

    auto ta1_y_zz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 115);

    auto ta1_y_zz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 116);

    auto ta1_y_zz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 117);

    auto ta1_y_zz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 118);

    auto ta1_y_zz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 119);

    auto ta1_z_xx_xxx_1 = pbuffer.data(idx_npot_geom_010_1_df + 120);

    auto ta1_z_xx_xxy_1 = pbuffer.data(idx_npot_geom_010_1_df + 121);

    auto ta1_z_xx_xxz_1 = pbuffer.data(idx_npot_geom_010_1_df + 122);

    auto ta1_z_xx_xyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 123);

    auto ta1_z_xx_xyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 124);

    auto ta1_z_xx_xzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 125);

    auto ta1_z_xx_yyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 126);

    auto ta1_z_xx_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 127);

    auto ta1_z_xx_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 128);

    auto ta1_z_xx_zzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 129);

    auto ta1_z_yy_xxx_1 = pbuffer.data(idx_npot_geom_010_1_df + 150);

    auto ta1_z_yy_xxy_1 = pbuffer.data(idx_npot_geom_010_1_df + 151);

    auto ta1_z_yy_xxz_1 = pbuffer.data(idx_npot_geom_010_1_df + 152);

    auto ta1_z_yy_xyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 153);

    auto ta1_z_yy_xyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 154);

    auto ta1_z_yy_xzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 155);

    auto ta1_z_yy_yyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 156);

    auto ta1_z_yy_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 157);

    auto ta1_z_yy_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 158);

    auto ta1_z_yy_zzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 159);

    auto ta1_z_yz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 164);

    auto ta1_z_yz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 167);

    auto ta1_z_yz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 168);

    auto ta1_z_zz_xxx_1 = pbuffer.data(idx_npot_geom_010_1_df + 170);

    auto ta1_z_zz_xxy_1 = pbuffer.data(idx_npot_geom_010_1_df + 171);

    auto ta1_z_zz_xxz_1 = pbuffer.data(idx_npot_geom_010_1_df + 172);

    auto ta1_z_zz_xyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 173);

    auto ta1_z_zz_xyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 174);

    auto ta1_z_zz_xzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 175);

    auto ta1_z_zz_yyy_1 = pbuffer.data(idx_npot_geom_010_1_df + 176);

    auto ta1_z_zz_yyz_1 = pbuffer.data(idx_npot_geom_010_1_df + 177);

    auto ta1_z_zz_yzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 178);

    auto ta1_z_zz_zzz_1 = pbuffer.data(idx_npot_geom_010_1_df + 179);

    // Set up components of auxiliary buffer : DG

    auto ta_xx_xxxx_1 = pbuffer.data(idx_npot_1_dg);

    auto ta_xx_xxxy_1 = pbuffer.data(idx_npot_1_dg + 1);

    auto ta_xx_xxxz_1 = pbuffer.data(idx_npot_1_dg + 2);

    auto ta_xx_xxyy_1 = pbuffer.data(idx_npot_1_dg + 3);

    auto ta_xx_xxyz_1 = pbuffer.data(idx_npot_1_dg + 4);

    auto ta_xx_xxzz_1 = pbuffer.data(idx_npot_1_dg + 5);

    auto ta_xx_xyyy_1 = pbuffer.data(idx_npot_1_dg + 6);

    auto ta_xx_xyyz_1 = pbuffer.data(idx_npot_1_dg + 7);

    auto ta_xx_xyzz_1 = pbuffer.data(idx_npot_1_dg + 8);

    auto ta_xx_xzzz_1 = pbuffer.data(idx_npot_1_dg + 9);

    auto ta_xx_yyyy_1 = pbuffer.data(idx_npot_1_dg + 10);

    auto ta_xx_yyyz_1 = pbuffer.data(idx_npot_1_dg + 11);

    auto ta_xx_yyzz_1 = pbuffer.data(idx_npot_1_dg + 12);

    auto ta_xx_yzzz_1 = pbuffer.data(idx_npot_1_dg + 13);

    auto ta_xx_zzzz_1 = pbuffer.data(idx_npot_1_dg + 14);

    auto ta_xy_xxxy_1 = pbuffer.data(idx_npot_1_dg + 16);

    auto ta_xy_xxyy_1 = pbuffer.data(idx_npot_1_dg + 18);

    auto ta_xy_xyyy_1 = pbuffer.data(idx_npot_1_dg + 21);

    auto ta_xz_xxxz_1 = pbuffer.data(idx_npot_1_dg + 32);

    auto ta_xz_xxzz_1 = pbuffer.data(idx_npot_1_dg + 35);

    auto ta_xz_xzzz_1 = pbuffer.data(idx_npot_1_dg + 39);

    auto ta_yy_xxxx_1 = pbuffer.data(idx_npot_1_dg + 45);

    auto ta_yy_xxxy_1 = pbuffer.data(idx_npot_1_dg + 46);

    auto ta_yy_xxxz_1 = pbuffer.data(idx_npot_1_dg + 47);

    auto ta_yy_xxyy_1 = pbuffer.data(idx_npot_1_dg + 48);

    auto ta_yy_xxyz_1 = pbuffer.data(idx_npot_1_dg + 49);

    auto ta_yy_xxzz_1 = pbuffer.data(idx_npot_1_dg + 50);

    auto ta_yy_xyyy_1 = pbuffer.data(idx_npot_1_dg + 51);

    auto ta_yy_xyyz_1 = pbuffer.data(idx_npot_1_dg + 52);

    auto ta_yy_xyzz_1 = pbuffer.data(idx_npot_1_dg + 53);

    auto ta_yy_xzzz_1 = pbuffer.data(idx_npot_1_dg + 54);

    auto ta_yy_yyyy_1 = pbuffer.data(idx_npot_1_dg + 55);

    auto ta_yy_yyyz_1 = pbuffer.data(idx_npot_1_dg + 56);

    auto ta_yy_yyzz_1 = pbuffer.data(idx_npot_1_dg + 57);

    auto ta_yy_yzzz_1 = pbuffer.data(idx_npot_1_dg + 58);

    auto ta_yy_zzzz_1 = pbuffer.data(idx_npot_1_dg + 59);

    auto ta_yz_yyyz_1 = pbuffer.data(idx_npot_1_dg + 71);

    auto ta_yz_yyzz_1 = pbuffer.data(idx_npot_1_dg + 72);

    auto ta_yz_yzzz_1 = pbuffer.data(idx_npot_1_dg + 73);

    auto ta_zz_xxxx_1 = pbuffer.data(idx_npot_1_dg + 75);

    auto ta_zz_xxxy_1 = pbuffer.data(idx_npot_1_dg + 76);

    auto ta_zz_xxxz_1 = pbuffer.data(idx_npot_1_dg + 77);

    auto ta_zz_xxyy_1 = pbuffer.data(idx_npot_1_dg + 78);

    auto ta_zz_xxyz_1 = pbuffer.data(idx_npot_1_dg + 79);

    auto ta_zz_xxzz_1 = pbuffer.data(idx_npot_1_dg + 80);

    auto ta_zz_xyyy_1 = pbuffer.data(idx_npot_1_dg + 81);

    auto ta_zz_xyyz_1 = pbuffer.data(idx_npot_1_dg + 82);

    auto ta_zz_xyzz_1 = pbuffer.data(idx_npot_1_dg + 83);

    auto ta_zz_xzzz_1 = pbuffer.data(idx_npot_1_dg + 84);

    auto ta_zz_yyyy_1 = pbuffer.data(idx_npot_1_dg + 85);

    auto ta_zz_yyyz_1 = pbuffer.data(idx_npot_1_dg + 86);

    auto ta_zz_yyzz_1 = pbuffer.data(idx_npot_1_dg + 87);

    auto ta_zz_yzzz_1 = pbuffer.data(idx_npot_1_dg + 88);

    auto ta_zz_zzzz_1 = pbuffer.data(idx_npot_1_dg + 89);

    // Set up components of auxiliary buffer : DG

    auto ta1_x_xx_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg);

    auto ta1_x_xx_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 1);

    auto ta1_x_xx_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 2);

    auto ta1_x_xx_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 3);

    auto ta1_x_xx_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 4);

    auto ta1_x_xx_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 5);

    auto ta1_x_xx_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 6);

    auto ta1_x_xx_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 7);

    auto ta1_x_xx_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 8);

    auto ta1_x_xx_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 9);

    auto ta1_x_xx_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 10);

    auto ta1_x_xx_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 11);

    auto ta1_x_xx_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 12);

    auto ta1_x_xx_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 13);

    auto ta1_x_xx_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 14);

    auto ta1_x_xy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 15);

    auto ta1_x_xy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 16);

    auto ta1_x_xy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 17);

    auto ta1_x_xy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 18);

    auto ta1_x_xy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 20);

    auto ta1_x_xy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 21);

    auto ta1_x_xy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 24);

    auto ta1_x_xy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 25);

    auto ta1_x_xz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 30);

    auto ta1_x_xz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 31);

    auto ta1_x_xz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 32);

    auto ta1_x_xz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 33);

    auto ta1_x_xz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 34);

    auto ta1_x_xz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 35);

    auto ta1_x_xz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 36);

    auto ta1_x_xz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 37);

    auto ta1_x_xz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 38);

    auto ta1_x_xz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 39);

    auto ta1_x_xz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 44);

    auto ta1_x_yy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 45);

    auto ta1_x_yy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 46);

    auto ta1_x_yy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 47);

    auto ta1_x_yy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 48);

    auto ta1_x_yy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 49);

    auto ta1_x_yy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 50);

    auto ta1_x_yy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 51);

    auto ta1_x_yy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 52);

    auto ta1_x_yy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 53);

    auto ta1_x_yy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 54);

    auto ta1_x_yy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 55);

    auto ta1_x_yy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 56);

    auto ta1_x_yy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 57);

    auto ta1_x_yy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 58);

    auto ta1_x_yy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 59);

    auto ta1_x_yz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 62);

    auto ta1_x_yz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 65);

    auto ta1_x_yz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 69);

    auto ta1_x_yz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 71);

    auto ta1_x_yz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 72);

    auto ta1_x_yz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 73);

    auto ta1_x_yz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 74);

    auto ta1_x_zz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 75);

    auto ta1_x_zz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 76);

    auto ta1_x_zz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 77);

    auto ta1_x_zz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 78);

    auto ta1_x_zz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 79);

    auto ta1_x_zz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 80);

    auto ta1_x_zz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 81);

    auto ta1_x_zz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 82);

    auto ta1_x_zz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 83);

    auto ta1_x_zz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 84);

    auto ta1_x_zz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 85);

    auto ta1_x_zz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 86);

    auto ta1_x_zz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 87);

    auto ta1_x_zz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 88);

    auto ta1_x_zz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 89);

    auto ta1_y_xx_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 90);

    auto ta1_y_xx_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 91);

    auto ta1_y_xx_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 92);

    auto ta1_y_xx_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 93);

    auto ta1_y_xx_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 94);

    auto ta1_y_xx_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 95);

    auto ta1_y_xx_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 96);

    auto ta1_y_xx_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 97);

    auto ta1_y_xx_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 98);

    auto ta1_y_xx_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 99);

    auto ta1_y_xx_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 100);

    auto ta1_y_xx_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 101);

    auto ta1_y_xx_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 102);

    auto ta1_y_xx_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 103);

    auto ta1_y_xx_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 104);

    auto ta1_y_xy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 105);

    auto ta1_y_xy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 106);

    auto ta1_y_xy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 108);

    auto ta1_y_xy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 111);

    auto ta1_y_xy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 115);

    auto ta1_y_xy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 116);

    auto ta1_y_xy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 117);

    auto ta1_y_xy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 118);

    auto ta1_y_xz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 122);

    auto ta1_y_xz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 125);

    auto ta1_y_xz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 129);

    auto ta1_y_xz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 131);

    auto ta1_y_xz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 132);

    auto ta1_y_xz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 133);

    auto ta1_y_xz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 134);

    auto ta1_y_yy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 135);

    auto ta1_y_yy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 136);

    auto ta1_y_yy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 137);

    auto ta1_y_yy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 138);

    auto ta1_y_yy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 139);

    auto ta1_y_yy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 140);

    auto ta1_y_yy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 141);

    auto ta1_y_yy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 142);

    auto ta1_y_yy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 143);

    auto ta1_y_yy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 144);

    auto ta1_y_yy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 145);

    auto ta1_y_yy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 146);

    auto ta1_y_yy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 147);

    auto ta1_y_yy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 148);

    auto ta1_y_yy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 149);

    auto ta1_y_yz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 151);

    auto ta1_y_yz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 153);

    auto ta1_y_yz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 154);

    auto ta1_y_yz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 156);

    auto ta1_y_yz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 157);

    auto ta1_y_yz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 158);

    auto ta1_y_yz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 160);

    auto ta1_y_yz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 161);

    auto ta1_y_yz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 162);

    auto ta1_y_yz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 163);

    auto ta1_y_yz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 164);

    auto ta1_y_zz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 165);

    auto ta1_y_zz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 166);

    auto ta1_y_zz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 167);

    auto ta1_y_zz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 168);

    auto ta1_y_zz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 169);

    auto ta1_y_zz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 170);

    auto ta1_y_zz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 171);

    auto ta1_y_zz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 172);

    auto ta1_y_zz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 173);

    auto ta1_y_zz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 174);

    auto ta1_y_zz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 175);

    auto ta1_y_zz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 176);

    auto ta1_y_zz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 177);

    auto ta1_y_zz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 178);

    auto ta1_y_zz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 179);

    auto ta1_z_xx_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 180);

    auto ta1_z_xx_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 181);

    auto ta1_z_xx_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 182);

    auto ta1_z_xx_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 183);

    auto ta1_z_xx_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 184);

    auto ta1_z_xx_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 185);

    auto ta1_z_xx_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 186);

    auto ta1_z_xx_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 187);

    auto ta1_z_xx_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 188);

    auto ta1_z_xx_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 189);

    auto ta1_z_xx_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 190);

    auto ta1_z_xx_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 191);

    auto ta1_z_xx_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 192);

    auto ta1_z_xx_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 193);

    auto ta1_z_xx_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 194);

    auto ta1_z_xy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 196);

    auto ta1_z_xy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 198);

    auto ta1_z_xy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 201);

    auto ta1_z_xy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 205);

    auto ta1_z_xy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 206);

    auto ta1_z_xy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 207);

    auto ta1_z_xy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 208);

    auto ta1_z_xz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 210);

    auto ta1_z_xz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 212);

    auto ta1_z_xz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 215);

    auto ta1_z_xz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 219);

    auto ta1_z_xz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 221);

    auto ta1_z_xz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 222);

    auto ta1_z_xz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 223);

    auto ta1_z_xz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 224);

    auto ta1_z_yy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 225);

    auto ta1_z_yy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 226);

    auto ta1_z_yy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 227);

    auto ta1_z_yy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 228);

    auto ta1_z_yy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 229);

    auto ta1_z_yy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 230);

    auto ta1_z_yy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 231);

    auto ta1_z_yy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 232);

    auto ta1_z_yy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 233);

    auto ta1_z_yy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 234);

    auto ta1_z_yy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 235);

    auto ta1_z_yy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 236);

    auto ta1_z_yy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 237);

    auto ta1_z_yy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 238);

    auto ta1_z_yy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 239);

    auto ta1_z_yz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 242);

    auto ta1_z_yz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 244);

    auto ta1_z_yz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 245);

    auto ta1_z_yz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 247);

    auto ta1_z_yz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 248);

    auto ta1_z_yz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 249);

    auto ta1_z_yz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 250);

    auto ta1_z_yz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 251);

    auto ta1_z_yz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 252);

    auto ta1_z_yz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 253);

    auto ta1_z_yz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 254);

    auto ta1_z_zz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_dg + 255);

    auto ta1_z_zz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 256);

    auto ta1_z_zz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 257);

    auto ta1_z_zz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 258);

    auto ta1_z_zz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 259);

    auto ta1_z_zz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 260);

    auto ta1_z_zz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 261);

    auto ta1_z_zz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 262);

    auto ta1_z_zz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 263);

    auto ta1_z_zz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 264);

    auto ta1_z_zz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_dg + 265);

    auto ta1_z_zz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 266);

    auto ta1_z_zz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 267);

    auto ta1_z_zz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 268);

    auto ta1_z_zz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_dg + 269);

    // Set up components of auxiliary buffer : DG

    auto ta1_x_xx_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg);

    auto ta1_x_xx_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 1);

    auto ta1_x_xx_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 2);

    auto ta1_x_xx_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 3);

    auto ta1_x_xx_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 4);

    auto ta1_x_xx_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 5);

    auto ta1_x_xx_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 6);

    auto ta1_x_xx_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 7);

    auto ta1_x_xx_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 8);

    auto ta1_x_xx_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 9);

    auto ta1_x_xx_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 10);

    auto ta1_x_xx_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 11);

    auto ta1_x_xx_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 12);

    auto ta1_x_xx_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 13);

    auto ta1_x_xx_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 14);

    auto ta1_x_xy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 15);

    auto ta1_x_xy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 16);

    auto ta1_x_xy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 17);

    auto ta1_x_xy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 18);

    auto ta1_x_xy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 20);

    auto ta1_x_xy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 21);

    auto ta1_x_xy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 24);

    auto ta1_x_xy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 25);

    auto ta1_x_xz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 30);

    auto ta1_x_xz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 31);

    auto ta1_x_xz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 32);

    auto ta1_x_xz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 33);

    auto ta1_x_xz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 34);

    auto ta1_x_xz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 35);

    auto ta1_x_xz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 36);

    auto ta1_x_xz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 37);

    auto ta1_x_xz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 38);

    auto ta1_x_xz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 39);

    auto ta1_x_xz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 44);

    auto ta1_x_yy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 45);

    auto ta1_x_yy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 46);

    auto ta1_x_yy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 47);

    auto ta1_x_yy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 48);

    auto ta1_x_yy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 49);

    auto ta1_x_yy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 50);

    auto ta1_x_yy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 51);

    auto ta1_x_yy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 52);

    auto ta1_x_yy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 53);

    auto ta1_x_yy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 54);

    auto ta1_x_yy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 55);

    auto ta1_x_yy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 56);

    auto ta1_x_yy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 57);

    auto ta1_x_yy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 58);

    auto ta1_x_yy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 59);

    auto ta1_x_yz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 62);

    auto ta1_x_yz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 65);

    auto ta1_x_yz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 69);

    auto ta1_x_yz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 71);

    auto ta1_x_yz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 72);

    auto ta1_x_yz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 73);

    auto ta1_x_yz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 74);

    auto ta1_x_zz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 75);

    auto ta1_x_zz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 76);

    auto ta1_x_zz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 77);

    auto ta1_x_zz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 78);

    auto ta1_x_zz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 79);

    auto ta1_x_zz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 80);

    auto ta1_x_zz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 81);

    auto ta1_x_zz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 82);

    auto ta1_x_zz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 83);

    auto ta1_x_zz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 84);

    auto ta1_x_zz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 85);

    auto ta1_x_zz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 86);

    auto ta1_x_zz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 87);

    auto ta1_x_zz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 88);

    auto ta1_x_zz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 89);

    auto ta1_y_xx_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 90);

    auto ta1_y_xx_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 91);

    auto ta1_y_xx_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 92);

    auto ta1_y_xx_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 93);

    auto ta1_y_xx_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 94);

    auto ta1_y_xx_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 95);

    auto ta1_y_xx_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 96);

    auto ta1_y_xx_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 97);

    auto ta1_y_xx_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 98);

    auto ta1_y_xx_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 99);

    auto ta1_y_xx_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 100);

    auto ta1_y_xx_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 101);

    auto ta1_y_xx_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 102);

    auto ta1_y_xx_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 103);

    auto ta1_y_xx_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 104);

    auto ta1_y_xy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 105);

    auto ta1_y_xy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 106);

    auto ta1_y_xy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 108);

    auto ta1_y_xy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 111);

    auto ta1_y_xy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 115);

    auto ta1_y_xy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 116);

    auto ta1_y_xy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 117);

    auto ta1_y_xy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 118);

    auto ta1_y_xz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 122);

    auto ta1_y_xz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 125);

    auto ta1_y_xz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 129);

    auto ta1_y_xz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 131);

    auto ta1_y_xz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 132);

    auto ta1_y_xz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 133);

    auto ta1_y_xz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 134);

    auto ta1_y_yy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 135);

    auto ta1_y_yy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 136);

    auto ta1_y_yy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 137);

    auto ta1_y_yy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 138);

    auto ta1_y_yy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 139);

    auto ta1_y_yy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 140);

    auto ta1_y_yy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 141);

    auto ta1_y_yy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 142);

    auto ta1_y_yy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 143);

    auto ta1_y_yy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 144);

    auto ta1_y_yy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 145);

    auto ta1_y_yy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 146);

    auto ta1_y_yy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 147);

    auto ta1_y_yy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 148);

    auto ta1_y_yy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 149);

    auto ta1_y_yz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 151);

    auto ta1_y_yz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 153);

    auto ta1_y_yz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 154);

    auto ta1_y_yz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 156);

    auto ta1_y_yz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 157);

    auto ta1_y_yz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 158);

    auto ta1_y_yz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 160);

    auto ta1_y_yz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 161);

    auto ta1_y_yz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 162);

    auto ta1_y_yz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 163);

    auto ta1_y_yz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 164);

    auto ta1_y_zz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 165);

    auto ta1_y_zz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 166);

    auto ta1_y_zz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 167);

    auto ta1_y_zz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 168);

    auto ta1_y_zz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 169);

    auto ta1_y_zz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 170);

    auto ta1_y_zz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 171);

    auto ta1_y_zz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 172);

    auto ta1_y_zz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 173);

    auto ta1_y_zz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 174);

    auto ta1_y_zz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 175);

    auto ta1_y_zz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 176);

    auto ta1_y_zz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 177);

    auto ta1_y_zz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 178);

    auto ta1_y_zz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 179);

    auto ta1_z_xx_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 180);

    auto ta1_z_xx_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 181);

    auto ta1_z_xx_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 182);

    auto ta1_z_xx_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 183);

    auto ta1_z_xx_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 184);

    auto ta1_z_xx_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 185);

    auto ta1_z_xx_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 186);

    auto ta1_z_xx_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 187);

    auto ta1_z_xx_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 188);

    auto ta1_z_xx_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 189);

    auto ta1_z_xx_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 190);

    auto ta1_z_xx_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 191);

    auto ta1_z_xx_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 192);

    auto ta1_z_xx_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 193);

    auto ta1_z_xx_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 194);

    auto ta1_z_xy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 196);

    auto ta1_z_xy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 198);

    auto ta1_z_xy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 201);

    auto ta1_z_xy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 205);

    auto ta1_z_xy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 206);

    auto ta1_z_xy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 207);

    auto ta1_z_xy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 208);

    auto ta1_z_xz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 210);

    auto ta1_z_xz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 212);

    auto ta1_z_xz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 215);

    auto ta1_z_xz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 219);

    auto ta1_z_xz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 221);

    auto ta1_z_xz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 222);

    auto ta1_z_xz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 223);

    auto ta1_z_xz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 224);

    auto ta1_z_yy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 225);

    auto ta1_z_yy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 226);

    auto ta1_z_yy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 227);

    auto ta1_z_yy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 228);

    auto ta1_z_yy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 229);

    auto ta1_z_yy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 230);

    auto ta1_z_yy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 231);

    auto ta1_z_yy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 232);

    auto ta1_z_yy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 233);

    auto ta1_z_yy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 234);

    auto ta1_z_yy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 235);

    auto ta1_z_yy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 236);

    auto ta1_z_yy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 237);

    auto ta1_z_yy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 238);

    auto ta1_z_yy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 239);

    auto ta1_z_yz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 242);

    auto ta1_z_yz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 244);

    auto ta1_z_yz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 245);

    auto ta1_z_yz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 247);

    auto ta1_z_yz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 248);

    auto ta1_z_yz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 249);

    auto ta1_z_yz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 250);

    auto ta1_z_yz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 251);

    auto ta1_z_yz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 252);

    auto ta1_z_yz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 253);

    auto ta1_z_yz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 254);

    auto ta1_z_zz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_dg + 255);

    auto ta1_z_zz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 256);

    auto ta1_z_zz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 257);

    auto ta1_z_zz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 258);

    auto ta1_z_zz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 259);

    auto ta1_z_zz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 260);

    auto ta1_z_zz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 261);

    auto ta1_z_zz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 262);

    auto ta1_z_zz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 263);

    auto ta1_z_zz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 264);

    auto ta1_z_zz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_dg + 265);

    auto ta1_z_zz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 266);

    auto ta1_z_zz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 267);

    auto ta1_z_zz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 268);

    auto ta1_z_zz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_dg + 269);

    // Set up 0-15 components of targeted buffer : FG

    auto ta1_x_xxx_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg);

    auto ta1_x_xxx_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 1);

    auto ta1_x_xxx_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 2);

    auto ta1_x_xxx_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 3);

    auto ta1_x_xxx_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 4);

    auto ta1_x_xxx_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 5);

    auto ta1_x_xxx_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 6);

    auto ta1_x_xxx_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 7);

    auto ta1_x_xxx_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 8);

    auto ta1_x_xxx_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 9);

    auto ta1_x_xxx_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 10);

    auto ta1_x_xxx_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 11);

    auto ta1_x_xxx_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 12);

    auto ta1_x_xxx_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 13);

    auto ta1_x_xxx_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 14);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_x_x_xxxx_0,   \
                             ta1_x_x_xxxx_1,   \
                             ta1_x_x_xxxy_0,   \
                             ta1_x_x_xxxy_1,   \
                             ta1_x_x_xxxz_0,   \
                             ta1_x_x_xxxz_1,   \
                             ta1_x_x_xxyy_0,   \
                             ta1_x_x_xxyy_1,   \
                             ta1_x_x_xxyz_0,   \
                             ta1_x_x_xxyz_1,   \
                             ta1_x_x_xxzz_0,   \
                             ta1_x_x_xxzz_1,   \
                             ta1_x_x_xyyy_0,   \
                             ta1_x_x_xyyy_1,   \
                             ta1_x_x_xyyz_0,   \
                             ta1_x_x_xyyz_1,   \
                             ta1_x_x_xyzz_0,   \
                             ta1_x_x_xyzz_1,   \
                             ta1_x_x_xzzz_0,   \
                             ta1_x_x_xzzz_1,   \
                             ta1_x_x_yyyy_0,   \
                             ta1_x_x_yyyy_1,   \
                             ta1_x_x_yyyz_0,   \
                             ta1_x_x_yyyz_1,   \
                             ta1_x_x_yyzz_0,   \
                             ta1_x_x_yyzz_1,   \
                             ta1_x_x_yzzz_0,   \
                             ta1_x_x_yzzz_1,   \
                             ta1_x_x_zzzz_0,   \
                             ta1_x_x_zzzz_1,   \
                             ta1_x_xx_xxx_0,   \
                             ta1_x_xx_xxx_1,   \
                             ta1_x_xx_xxxx_0,  \
                             ta1_x_xx_xxxx_1,  \
                             ta1_x_xx_xxxy_0,  \
                             ta1_x_xx_xxxy_1,  \
                             ta1_x_xx_xxxz_0,  \
                             ta1_x_xx_xxxz_1,  \
                             ta1_x_xx_xxy_0,   \
                             ta1_x_xx_xxy_1,   \
                             ta1_x_xx_xxyy_0,  \
                             ta1_x_xx_xxyy_1,  \
                             ta1_x_xx_xxyz_0,  \
                             ta1_x_xx_xxyz_1,  \
                             ta1_x_xx_xxz_0,   \
                             ta1_x_xx_xxz_1,   \
                             ta1_x_xx_xxzz_0,  \
                             ta1_x_xx_xxzz_1,  \
                             ta1_x_xx_xyy_0,   \
                             ta1_x_xx_xyy_1,   \
                             ta1_x_xx_xyyy_0,  \
                             ta1_x_xx_xyyy_1,  \
                             ta1_x_xx_xyyz_0,  \
                             ta1_x_xx_xyyz_1,  \
                             ta1_x_xx_xyz_0,   \
                             ta1_x_xx_xyz_1,   \
                             ta1_x_xx_xyzz_0,  \
                             ta1_x_xx_xyzz_1,  \
                             ta1_x_xx_xzz_0,   \
                             ta1_x_xx_xzz_1,   \
                             ta1_x_xx_xzzz_0,  \
                             ta1_x_xx_xzzz_1,  \
                             ta1_x_xx_yyy_0,   \
                             ta1_x_xx_yyy_1,   \
                             ta1_x_xx_yyyy_0,  \
                             ta1_x_xx_yyyy_1,  \
                             ta1_x_xx_yyyz_0,  \
                             ta1_x_xx_yyyz_1,  \
                             ta1_x_xx_yyz_0,   \
                             ta1_x_xx_yyz_1,   \
                             ta1_x_xx_yyzz_0,  \
                             ta1_x_xx_yyzz_1,  \
                             ta1_x_xx_yzz_0,   \
                             ta1_x_xx_yzz_1,   \
                             ta1_x_xx_yzzz_0,  \
                             ta1_x_xx_yzzz_1,  \
                             ta1_x_xx_zzz_0,   \
                             ta1_x_xx_zzz_1,   \
                             ta1_x_xx_zzzz_0,  \
                             ta1_x_xx_zzzz_1,  \
                             ta1_x_xxx_xxxx_0, \
                             ta1_x_xxx_xxxy_0, \
                             ta1_x_xxx_xxxz_0, \
                             ta1_x_xxx_xxyy_0, \
                             ta1_x_xxx_xxyz_0, \
                             ta1_x_xxx_xxzz_0, \
                             ta1_x_xxx_xyyy_0, \
                             ta1_x_xxx_xyyz_0, \
                             ta1_x_xxx_xyzz_0, \
                             ta1_x_xxx_xzzz_0, \
                             ta1_x_xxx_yyyy_0, \
                             ta1_x_xxx_yyyz_0, \
                             ta1_x_xxx_yyzz_0, \
                             ta1_x_xxx_yzzz_0, \
                             ta1_x_xxx_zzzz_0, \
                             ta_xx_xxxx_1,     \
                             ta_xx_xxxy_1,     \
                             ta_xx_xxxz_1,     \
                             ta_xx_xxyy_1,     \
                             ta_xx_xxyz_1,     \
                             ta_xx_xxzz_1,     \
                             ta_xx_xyyy_1,     \
                             ta_xx_xyyz_1,     \
                             ta_xx_xyzz_1,     \
                             ta_xx_xzzz_1,     \
                             ta_xx_yyyy_1,     \
                             ta_xx_yyyz_1,     \
                             ta_xx_yyzz_1,     \
                             ta_xx_yzzz_1,     \
                             ta_xx_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxx_xxxx_0[i] = 2.0 * ta1_x_x_xxxx_0[i] * fe_0 - 2.0 * ta1_x_x_xxxx_1[i] * fe_0 +
                              4.0 * ta1_x_xx_xxx_0[i] * fe_0 - 4.0 * ta1_x_xx_xxx_1[i] * fe_0 + ta_xx_xxxx_1[i] +
                              ta1_x_xx_xxxx_0[i] * pa_x[i] - ta1_x_xx_xxxx_1[i] * pc_x[i];

        ta1_x_xxx_xxxy_0[i] = 2.0 * ta1_x_x_xxxy_0[i] * fe_0 - 2.0 * ta1_x_x_xxxy_1[i] * fe_0 +
                              3.0 * ta1_x_xx_xxy_0[i] * fe_0 - 3.0 * ta1_x_xx_xxy_1[i] * fe_0 + ta_xx_xxxy_1[i] +
                              ta1_x_xx_xxxy_0[i] * pa_x[i] - ta1_x_xx_xxxy_1[i] * pc_x[i];

        ta1_x_xxx_xxxz_0[i] = 2.0 * ta1_x_x_xxxz_0[i] * fe_0 - 2.0 * ta1_x_x_xxxz_1[i] * fe_0 +
                              3.0 * ta1_x_xx_xxz_0[i] * fe_0 - 3.0 * ta1_x_xx_xxz_1[i] * fe_0 + ta_xx_xxxz_1[i] +
                              ta1_x_xx_xxxz_0[i] * pa_x[i] - ta1_x_xx_xxxz_1[i] * pc_x[i];

        ta1_x_xxx_xxyy_0[i] = 2.0 * ta1_x_x_xxyy_0[i] * fe_0 - 2.0 * ta1_x_x_xxyy_1[i] * fe_0 +
                              2.0 * ta1_x_xx_xyy_0[i] * fe_0 - 2.0 * ta1_x_xx_xyy_1[i] * fe_0 + ta_xx_xxyy_1[i] +
                              ta1_x_xx_xxyy_0[i] * pa_x[i] - ta1_x_xx_xxyy_1[i] * pc_x[i];

        ta1_x_xxx_xxyz_0[i] = 2.0 * ta1_x_x_xxyz_0[i] * fe_0 - 2.0 * ta1_x_x_xxyz_1[i] * fe_0 +
                              2.0 * ta1_x_xx_xyz_0[i] * fe_0 - 2.0 * ta1_x_xx_xyz_1[i] * fe_0 + ta_xx_xxyz_1[i] +
                              ta1_x_xx_xxyz_0[i] * pa_x[i] - ta1_x_xx_xxyz_1[i] * pc_x[i];

        ta1_x_xxx_xxzz_0[i] = 2.0 * ta1_x_x_xxzz_0[i] * fe_0 - 2.0 * ta1_x_x_xxzz_1[i] * fe_0 +
                              2.0 * ta1_x_xx_xzz_0[i] * fe_0 - 2.0 * ta1_x_xx_xzz_1[i] * fe_0 + ta_xx_xxzz_1[i] +
                              ta1_x_xx_xxzz_0[i] * pa_x[i] - ta1_x_xx_xxzz_1[i] * pc_x[i];

        ta1_x_xxx_xyyy_0[i] = 2.0 * ta1_x_x_xyyy_0[i] * fe_0 - 2.0 * ta1_x_x_xyyy_1[i] * fe_0 +
                              ta1_x_xx_yyy_0[i] * fe_0 - ta1_x_xx_yyy_1[i] * fe_0 + ta_xx_xyyy_1[i] +
                              ta1_x_xx_xyyy_0[i] * pa_x[i] - ta1_x_xx_xyyy_1[i] * pc_x[i];

        ta1_x_xxx_xyyz_0[i] = 2.0 * ta1_x_x_xyyz_0[i] * fe_0 - 2.0 * ta1_x_x_xyyz_1[i] * fe_0 +
                              ta1_x_xx_yyz_0[i] * fe_0 - ta1_x_xx_yyz_1[i] * fe_0 + ta_xx_xyyz_1[i] +
                              ta1_x_xx_xyyz_0[i] * pa_x[i] - ta1_x_xx_xyyz_1[i] * pc_x[i];

        ta1_x_xxx_xyzz_0[i] = 2.0 * ta1_x_x_xyzz_0[i] * fe_0 - 2.0 * ta1_x_x_xyzz_1[i] * fe_0 +
                              ta1_x_xx_yzz_0[i] * fe_0 - ta1_x_xx_yzz_1[i] * fe_0 + ta_xx_xyzz_1[i] +
                              ta1_x_xx_xyzz_0[i] * pa_x[i] - ta1_x_xx_xyzz_1[i] * pc_x[i];

        ta1_x_xxx_xzzz_0[i] = 2.0 * ta1_x_x_xzzz_0[i] * fe_0 - 2.0 * ta1_x_x_xzzz_1[i] * fe_0 +
                              ta1_x_xx_zzz_0[i] * fe_0 - ta1_x_xx_zzz_1[i] * fe_0 + ta_xx_xzzz_1[i] +
                              ta1_x_xx_xzzz_0[i] * pa_x[i] - ta1_x_xx_xzzz_1[i] * pc_x[i];

        ta1_x_xxx_yyyy_0[i] = 2.0 * ta1_x_x_yyyy_0[i] * fe_0 - 2.0 * ta1_x_x_yyyy_1[i] * fe_0 + ta_xx_yyyy_1[i] +
                              ta1_x_xx_yyyy_0[i] * pa_x[i] - ta1_x_xx_yyyy_1[i] * pc_x[i];

        ta1_x_xxx_yyyz_0[i] = 2.0 * ta1_x_x_yyyz_0[i] * fe_0 - 2.0 * ta1_x_x_yyyz_1[i] * fe_0 + ta_xx_yyyz_1[i] +
                              ta1_x_xx_yyyz_0[i] * pa_x[i] - ta1_x_xx_yyyz_1[i] * pc_x[i];

        ta1_x_xxx_yyzz_0[i] = 2.0 * ta1_x_x_yyzz_0[i] * fe_0 - 2.0 * ta1_x_x_yyzz_1[i] * fe_0 + ta_xx_yyzz_1[i] +
                              ta1_x_xx_yyzz_0[i] * pa_x[i] - ta1_x_xx_yyzz_1[i] * pc_x[i];

        ta1_x_xxx_yzzz_0[i] = 2.0 * ta1_x_x_yzzz_0[i] * fe_0 - 2.0 * ta1_x_x_yzzz_1[i] * fe_0 + ta_xx_yzzz_1[i] +
                              ta1_x_xx_yzzz_0[i] * pa_x[i] - ta1_x_xx_yzzz_1[i] * pc_x[i];

        ta1_x_xxx_zzzz_0[i] = 2.0 * ta1_x_x_zzzz_0[i] * fe_0 - 2.0 * ta1_x_x_zzzz_1[i] * fe_0 + ta_xx_zzzz_1[i] +
                              ta1_x_xx_zzzz_0[i] * pa_x[i] - ta1_x_xx_zzzz_1[i] * pc_x[i];
    }

    // Set up 15-30 components of targeted buffer : FG

    auto ta1_x_xxy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 15);

    auto ta1_x_xxy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 16);

    auto ta1_x_xxy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 17);

    auto ta1_x_xxy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 18);

    auto ta1_x_xxy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 19);

    auto ta1_x_xxy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 20);

    auto ta1_x_xxy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 21);

    auto ta1_x_xxy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 22);

    auto ta1_x_xxy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 23);

    auto ta1_x_xxy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 24);

    auto ta1_x_xxy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 25);

    auto ta1_x_xxy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 26);

    auto ta1_x_xxy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 27);

    auto ta1_x_xxy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 28);

    auto ta1_x_xxy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 29);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_x_xx_xxx_0,   \
                             ta1_x_xx_xxx_1,   \
                             ta1_x_xx_xxxx_0,  \
                             ta1_x_xx_xxxx_1,  \
                             ta1_x_xx_xxxy_0,  \
                             ta1_x_xx_xxxy_1,  \
                             ta1_x_xx_xxxz_0,  \
                             ta1_x_xx_xxxz_1,  \
                             ta1_x_xx_xxy_0,   \
                             ta1_x_xx_xxy_1,   \
                             ta1_x_xx_xxyy_0,  \
                             ta1_x_xx_xxyy_1,  \
                             ta1_x_xx_xxyz_0,  \
                             ta1_x_xx_xxyz_1,  \
                             ta1_x_xx_xxz_0,   \
                             ta1_x_xx_xxz_1,   \
                             ta1_x_xx_xxzz_0,  \
                             ta1_x_xx_xxzz_1,  \
                             ta1_x_xx_xyy_0,   \
                             ta1_x_xx_xyy_1,   \
                             ta1_x_xx_xyyy_0,  \
                             ta1_x_xx_xyyy_1,  \
                             ta1_x_xx_xyyz_0,  \
                             ta1_x_xx_xyyz_1,  \
                             ta1_x_xx_xyz_0,   \
                             ta1_x_xx_xyz_1,   \
                             ta1_x_xx_xyzz_0,  \
                             ta1_x_xx_xyzz_1,  \
                             ta1_x_xx_xzz_0,   \
                             ta1_x_xx_xzz_1,   \
                             ta1_x_xx_xzzz_0,  \
                             ta1_x_xx_xzzz_1,  \
                             ta1_x_xx_yyy_0,   \
                             ta1_x_xx_yyy_1,   \
                             ta1_x_xx_yyyy_0,  \
                             ta1_x_xx_yyyy_1,  \
                             ta1_x_xx_yyyz_0,  \
                             ta1_x_xx_yyyz_1,  \
                             ta1_x_xx_yyz_0,   \
                             ta1_x_xx_yyz_1,   \
                             ta1_x_xx_yyzz_0,  \
                             ta1_x_xx_yyzz_1,  \
                             ta1_x_xx_yzz_0,   \
                             ta1_x_xx_yzz_1,   \
                             ta1_x_xx_yzzz_0,  \
                             ta1_x_xx_yzzz_1,  \
                             ta1_x_xx_zzz_0,   \
                             ta1_x_xx_zzz_1,   \
                             ta1_x_xx_zzzz_0,  \
                             ta1_x_xx_zzzz_1,  \
                             ta1_x_xxy_xxxx_0, \
                             ta1_x_xxy_xxxy_0, \
                             ta1_x_xxy_xxxz_0, \
                             ta1_x_xxy_xxyy_0, \
                             ta1_x_xxy_xxyz_0, \
                             ta1_x_xxy_xxzz_0, \
                             ta1_x_xxy_xyyy_0, \
                             ta1_x_xxy_xyyz_0, \
                             ta1_x_xxy_xyzz_0, \
                             ta1_x_xxy_xzzz_0, \
                             ta1_x_xxy_yyyy_0, \
                             ta1_x_xxy_yyyz_0, \
                             ta1_x_xxy_yyzz_0, \
                             ta1_x_xxy_yzzz_0, \
                             ta1_x_xxy_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxy_xxxx_0[i] = ta1_x_xx_xxxx_0[i] * pa_y[i] - ta1_x_xx_xxxx_1[i] * pc_y[i];

        ta1_x_xxy_xxxy_0[i] = ta1_x_xx_xxx_0[i] * fe_0 - ta1_x_xx_xxx_1[i] * fe_0 + ta1_x_xx_xxxy_0[i] * pa_y[i] -
                              ta1_x_xx_xxxy_1[i] * pc_y[i];

        ta1_x_xxy_xxxz_0[i] = ta1_x_xx_xxxz_0[i] * pa_y[i] - ta1_x_xx_xxxz_1[i] * pc_y[i];

        ta1_x_xxy_xxyy_0[i] = 2.0 * ta1_x_xx_xxy_0[i] * fe_0 - 2.0 * ta1_x_xx_xxy_1[i] * fe_0 +
                              ta1_x_xx_xxyy_0[i] * pa_y[i] - ta1_x_xx_xxyy_1[i] * pc_y[i];

        ta1_x_xxy_xxyz_0[i] = ta1_x_xx_xxz_0[i] * fe_0 - ta1_x_xx_xxz_1[i] * fe_0 + ta1_x_xx_xxyz_0[i] * pa_y[i] -
                              ta1_x_xx_xxyz_1[i] * pc_y[i];

        ta1_x_xxy_xxzz_0[i] = ta1_x_xx_xxzz_0[i] * pa_y[i] - ta1_x_xx_xxzz_1[i] * pc_y[i];

        ta1_x_xxy_xyyy_0[i] = 3.0 * ta1_x_xx_xyy_0[i] * fe_0 - 3.0 * ta1_x_xx_xyy_1[i] * fe_0 +
                              ta1_x_xx_xyyy_0[i] * pa_y[i] - ta1_x_xx_xyyy_1[i] * pc_y[i];

        ta1_x_xxy_xyyz_0[i] = 2.0 * ta1_x_xx_xyz_0[i] * fe_0 - 2.0 * ta1_x_xx_xyz_1[i] * fe_0 +
                              ta1_x_xx_xyyz_0[i] * pa_y[i] - ta1_x_xx_xyyz_1[i] * pc_y[i];

        ta1_x_xxy_xyzz_0[i] = ta1_x_xx_xzz_0[i] * fe_0 - ta1_x_xx_xzz_1[i] * fe_0 + ta1_x_xx_xyzz_0[i] * pa_y[i] -
                              ta1_x_xx_xyzz_1[i] * pc_y[i];

        ta1_x_xxy_xzzz_0[i] = ta1_x_xx_xzzz_0[i] * pa_y[i] - ta1_x_xx_xzzz_1[i] * pc_y[i];

        ta1_x_xxy_yyyy_0[i] = 4.0 * ta1_x_xx_yyy_0[i] * fe_0 - 4.0 * ta1_x_xx_yyy_1[i] * fe_0 +
                              ta1_x_xx_yyyy_0[i] * pa_y[i] - ta1_x_xx_yyyy_1[i] * pc_y[i];

        ta1_x_xxy_yyyz_0[i] = 3.0 * ta1_x_xx_yyz_0[i] * fe_0 - 3.0 * ta1_x_xx_yyz_1[i] * fe_0 +
                              ta1_x_xx_yyyz_0[i] * pa_y[i] - ta1_x_xx_yyyz_1[i] * pc_y[i];

        ta1_x_xxy_yyzz_0[i] = 2.0 * ta1_x_xx_yzz_0[i] * fe_0 - 2.0 * ta1_x_xx_yzz_1[i] * fe_0 +
                              ta1_x_xx_yyzz_0[i] * pa_y[i] - ta1_x_xx_yyzz_1[i] * pc_y[i];

        ta1_x_xxy_yzzz_0[i] = ta1_x_xx_zzz_0[i] * fe_0 - ta1_x_xx_zzz_1[i] * fe_0 + ta1_x_xx_yzzz_0[i] * pa_y[i] -
                              ta1_x_xx_yzzz_1[i] * pc_y[i];

        ta1_x_xxy_zzzz_0[i] = ta1_x_xx_zzzz_0[i] * pa_y[i] - ta1_x_xx_zzzz_1[i] * pc_y[i];
    }

    // Set up 30-45 components of targeted buffer : FG

    auto ta1_x_xxz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 30);

    auto ta1_x_xxz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 31);

    auto ta1_x_xxz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 32);

    auto ta1_x_xxz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 33);

    auto ta1_x_xxz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 34);

    auto ta1_x_xxz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 35);

    auto ta1_x_xxz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 36);

    auto ta1_x_xxz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 37);

    auto ta1_x_xxz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 38);

    auto ta1_x_xxz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 39);

    auto ta1_x_xxz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 40);

    auto ta1_x_xxz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 41);

    auto ta1_x_xxz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 42);

    auto ta1_x_xxz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 43);

    auto ta1_x_xxz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 44);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta1_x_xx_xxx_0,   \
                             ta1_x_xx_xxx_1,   \
                             ta1_x_xx_xxxx_0,  \
                             ta1_x_xx_xxxx_1,  \
                             ta1_x_xx_xxxy_0,  \
                             ta1_x_xx_xxxy_1,  \
                             ta1_x_xx_xxxz_0,  \
                             ta1_x_xx_xxxz_1,  \
                             ta1_x_xx_xxy_0,   \
                             ta1_x_xx_xxy_1,   \
                             ta1_x_xx_xxyy_0,  \
                             ta1_x_xx_xxyy_1,  \
                             ta1_x_xx_xxyz_0,  \
                             ta1_x_xx_xxyz_1,  \
                             ta1_x_xx_xxz_0,   \
                             ta1_x_xx_xxz_1,   \
                             ta1_x_xx_xxzz_0,  \
                             ta1_x_xx_xxzz_1,  \
                             ta1_x_xx_xyy_0,   \
                             ta1_x_xx_xyy_1,   \
                             ta1_x_xx_xyyy_0,  \
                             ta1_x_xx_xyyy_1,  \
                             ta1_x_xx_xyyz_0,  \
                             ta1_x_xx_xyyz_1,  \
                             ta1_x_xx_xyz_0,   \
                             ta1_x_xx_xyz_1,   \
                             ta1_x_xx_xyzz_0,  \
                             ta1_x_xx_xyzz_1,  \
                             ta1_x_xx_xzz_0,   \
                             ta1_x_xx_xzz_1,   \
                             ta1_x_xx_xzzz_0,  \
                             ta1_x_xx_xzzz_1,  \
                             ta1_x_xx_yyy_0,   \
                             ta1_x_xx_yyy_1,   \
                             ta1_x_xx_yyyy_0,  \
                             ta1_x_xx_yyyy_1,  \
                             ta1_x_xx_yyyz_0,  \
                             ta1_x_xx_yyyz_1,  \
                             ta1_x_xx_yyz_0,   \
                             ta1_x_xx_yyz_1,   \
                             ta1_x_xx_yyzz_0,  \
                             ta1_x_xx_yyzz_1,  \
                             ta1_x_xx_yzz_0,   \
                             ta1_x_xx_yzz_1,   \
                             ta1_x_xx_yzzz_0,  \
                             ta1_x_xx_yzzz_1,  \
                             ta1_x_xx_zzz_0,   \
                             ta1_x_xx_zzz_1,   \
                             ta1_x_xx_zzzz_0,  \
                             ta1_x_xx_zzzz_1,  \
                             ta1_x_xxz_xxxx_0, \
                             ta1_x_xxz_xxxy_0, \
                             ta1_x_xxz_xxxz_0, \
                             ta1_x_xxz_xxyy_0, \
                             ta1_x_xxz_xxyz_0, \
                             ta1_x_xxz_xxzz_0, \
                             ta1_x_xxz_xyyy_0, \
                             ta1_x_xxz_xyyz_0, \
                             ta1_x_xxz_xyzz_0, \
                             ta1_x_xxz_xzzz_0, \
                             ta1_x_xxz_yyyy_0, \
                             ta1_x_xxz_yyyz_0, \
                             ta1_x_xxz_yyzz_0, \
                             ta1_x_xxz_yzzz_0, \
                             ta1_x_xxz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxz_xxxx_0[i] = ta1_x_xx_xxxx_0[i] * pa_z[i] - ta1_x_xx_xxxx_1[i] * pc_z[i];

        ta1_x_xxz_xxxy_0[i] = ta1_x_xx_xxxy_0[i] * pa_z[i] - ta1_x_xx_xxxy_1[i] * pc_z[i];

        ta1_x_xxz_xxxz_0[i] = ta1_x_xx_xxx_0[i] * fe_0 - ta1_x_xx_xxx_1[i] * fe_0 + ta1_x_xx_xxxz_0[i] * pa_z[i] -
                              ta1_x_xx_xxxz_1[i] * pc_z[i];

        ta1_x_xxz_xxyy_0[i] = ta1_x_xx_xxyy_0[i] * pa_z[i] - ta1_x_xx_xxyy_1[i] * pc_z[i];

        ta1_x_xxz_xxyz_0[i] = ta1_x_xx_xxy_0[i] * fe_0 - ta1_x_xx_xxy_1[i] * fe_0 + ta1_x_xx_xxyz_0[i] * pa_z[i] -
                              ta1_x_xx_xxyz_1[i] * pc_z[i];

        ta1_x_xxz_xxzz_0[i] = 2.0 * ta1_x_xx_xxz_0[i] * fe_0 - 2.0 * ta1_x_xx_xxz_1[i] * fe_0 +
                              ta1_x_xx_xxzz_0[i] * pa_z[i] - ta1_x_xx_xxzz_1[i] * pc_z[i];

        ta1_x_xxz_xyyy_0[i] = ta1_x_xx_xyyy_0[i] * pa_z[i] - ta1_x_xx_xyyy_1[i] * pc_z[i];

        ta1_x_xxz_xyyz_0[i] = ta1_x_xx_xyy_0[i] * fe_0 - ta1_x_xx_xyy_1[i] * fe_0 + ta1_x_xx_xyyz_0[i] * pa_z[i] -
                              ta1_x_xx_xyyz_1[i] * pc_z[i];

        ta1_x_xxz_xyzz_0[i] = 2.0 * ta1_x_xx_xyz_0[i] * fe_0 - 2.0 * ta1_x_xx_xyz_1[i] * fe_0 +
                              ta1_x_xx_xyzz_0[i] * pa_z[i] - ta1_x_xx_xyzz_1[i] * pc_z[i];

        ta1_x_xxz_xzzz_0[i] = 3.0 * ta1_x_xx_xzz_0[i] * fe_0 - 3.0 * ta1_x_xx_xzz_1[i] * fe_0 +
                              ta1_x_xx_xzzz_0[i] * pa_z[i] - ta1_x_xx_xzzz_1[i] * pc_z[i];

        ta1_x_xxz_yyyy_0[i] = ta1_x_xx_yyyy_0[i] * pa_z[i] - ta1_x_xx_yyyy_1[i] * pc_z[i];

        ta1_x_xxz_yyyz_0[i] = ta1_x_xx_yyy_0[i] * fe_0 - ta1_x_xx_yyy_1[i] * fe_0 + ta1_x_xx_yyyz_0[i] * pa_z[i] -
                              ta1_x_xx_yyyz_1[i] * pc_z[i];

        ta1_x_xxz_yyzz_0[i] = 2.0 * ta1_x_xx_yyz_0[i] * fe_0 - 2.0 * ta1_x_xx_yyz_1[i] * fe_0 +
                              ta1_x_xx_yyzz_0[i] * pa_z[i] - ta1_x_xx_yyzz_1[i] * pc_z[i];

        ta1_x_xxz_yzzz_0[i] = 3.0 * ta1_x_xx_yzz_0[i] * fe_0 - 3.0 * ta1_x_xx_yzz_1[i] * fe_0 +
                              ta1_x_xx_yzzz_0[i] * pa_z[i] - ta1_x_xx_yzzz_1[i] * pc_z[i];

        ta1_x_xxz_zzzz_0[i] = 4.0 * ta1_x_xx_zzz_0[i] * fe_0 - 4.0 * ta1_x_xx_zzz_1[i] * fe_0 +
                              ta1_x_xx_zzzz_0[i] * pa_z[i] - ta1_x_xx_zzzz_1[i] * pc_z[i];
    }

    // Set up 45-60 components of targeted buffer : FG

    auto ta1_x_xyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 45);

    auto ta1_x_xyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 46);

    auto ta1_x_xyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 47);

    auto ta1_x_xyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 48);

    auto ta1_x_xyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 49);

    auto ta1_x_xyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 50);

    auto ta1_x_xyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 51);

    auto ta1_x_xyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 52);

    auto ta1_x_xyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 53);

    auto ta1_x_xyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 54);

    auto ta1_x_xyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 55);

    auto ta1_x_xyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 56);

    auto ta1_x_xyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 57);

    auto ta1_x_xyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 58);

    auto ta1_x_xyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 59);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_x_x_xxxx_0,   \
                             ta1_x_x_xxxx_1,   \
                             ta1_x_x_xxxz_0,   \
                             ta1_x_x_xxxz_1,   \
                             ta1_x_x_xxzz_0,   \
                             ta1_x_x_xxzz_1,   \
                             ta1_x_x_xzzz_0,   \
                             ta1_x_x_xzzz_1,   \
                             ta1_x_xy_xxxx_0,  \
                             ta1_x_xy_xxxx_1,  \
                             ta1_x_xy_xxxz_0,  \
                             ta1_x_xy_xxxz_1,  \
                             ta1_x_xy_xxzz_0,  \
                             ta1_x_xy_xxzz_1,  \
                             ta1_x_xy_xzzz_0,  \
                             ta1_x_xy_xzzz_1,  \
                             ta1_x_xyy_xxxx_0, \
                             ta1_x_xyy_xxxy_0, \
                             ta1_x_xyy_xxxz_0, \
                             ta1_x_xyy_xxyy_0, \
                             ta1_x_xyy_xxyz_0, \
                             ta1_x_xyy_xxzz_0, \
                             ta1_x_xyy_xyyy_0, \
                             ta1_x_xyy_xyyz_0, \
                             ta1_x_xyy_xyzz_0, \
                             ta1_x_xyy_xzzz_0, \
                             ta1_x_xyy_yyyy_0, \
                             ta1_x_xyy_yyyz_0, \
                             ta1_x_xyy_yyzz_0, \
                             ta1_x_xyy_yzzz_0, \
                             ta1_x_xyy_zzzz_0, \
                             ta1_x_yy_xxxy_0,  \
                             ta1_x_yy_xxxy_1,  \
                             ta1_x_yy_xxy_0,   \
                             ta1_x_yy_xxy_1,   \
                             ta1_x_yy_xxyy_0,  \
                             ta1_x_yy_xxyy_1,  \
                             ta1_x_yy_xxyz_0,  \
                             ta1_x_yy_xxyz_1,  \
                             ta1_x_yy_xyy_0,   \
                             ta1_x_yy_xyy_1,   \
                             ta1_x_yy_xyyy_0,  \
                             ta1_x_yy_xyyy_1,  \
                             ta1_x_yy_xyyz_0,  \
                             ta1_x_yy_xyyz_1,  \
                             ta1_x_yy_xyz_0,   \
                             ta1_x_yy_xyz_1,   \
                             ta1_x_yy_xyzz_0,  \
                             ta1_x_yy_xyzz_1,  \
                             ta1_x_yy_yyy_0,   \
                             ta1_x_yy_yyy_1,   \
                             ta1_x_yy_yyyy_0,  \
                             ta1_x_yy_yyyy_1,  \
                             ta1_x_yy_yyyz_0,  \
                             ta1_x_yy_yyyz_1,  \
                             ta1_x_yy_yyz_0,   \
                             ta1_x_yy_yyz_1,   \
                             ta1_x_yy_yyzz_0,  \
                             ta1_x_yy_yyzz_1,  \
                             ta1_x_yy_yzz_0,   \
                             ta1_x_yy_yzz_1,   \
                             ta1_x_yy_yzzz_0,  \
                             ta1_x_yy_yzzz_1,  \
                             ta1_x_yy_zzzz_0,  \
                             ta1_x_yy_zzzz_1,  \
                             ta_yy_xxxy_1,     \
                             ta_yy_xxyy_1,     \
                             ta_yy_xxyz_1,     \
                             ta_yy_xyyy_1,     \
                             ta_yy_xyyz_1,     \
                             ta_yy_xyzz_1,     \
                             ta_yy_yyyy_1,     \
                             ta_yy_yyyz_1,     \
                             ta_yy_yyzz_1,     \
                             ta_yy_yzzz_1,     \
                             ta_yy_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyy_xxxx_0[i] = ta1_x_x_xxxx_0[i] * fe_0 - ta1_x_x_xxxx_1[i] * fe_0 + ta1_x_xy_xxxx_0[i] * pa_y[i] -
                              ta1_x_xy_xxxx_1[i] * pc_y[i];

        ta1_x_xyy_xxxy_0[i] = 3.0 * ta1_x_yy_xxy_0[i] * fe_0 - 3.0 * ta1_x_yy_xxy_1[i] * fe_0 + ta_yy_xxxy_1[i] +
                              ta1_x_yy_xxxy_0[i] * pa_x[i] - ta1_x_yy_xxxy_1[i] * pc_x[i];

        ta1_x_xyy_xxxz_0[i] = ta1_x_x_xxxz_0[i] * fe_0 - ta1_x_x_xxxz_1[i] * fe_0 + ta1_x_xy_xxxz_0[i] * pa_y[i] -
                              ta1_x_xy_xxxz_1[i] * pc_y[i];

        ta1_x_xyy_xxyy_0[i] = 2.0 * ta1_x_yy_xyy_0[i] * fe_0 - 2.0 * ta1_x_yy_xyy_1[i] * fe_0 + ta_yy_xxyy_1[i] +
                              ta1_x_yy_xxyy_0[i] * pa_x[i] - ta1_x_yy_xxyy_1[i] * pc_x[i];

        ta1_x_xyy_xxyz_0[i] = 2.0 * ta1_x_yy_xyz_0[i] * fe_0 - 2.0 * ta1_x_yy_xyz_1[i] * fe_0 + ta_yy_xxyz_1[i] +
                              ta1_x_yy_xxyz_0[i] * pa_x[i] - ta1_x_yy_xxyz_1[i] * pc_x[i];

        ta1_x_xyy_xxzz_0[i] = ta1_x_x_xxzz_0[i] * fe_0 - ta1_x_x_xxzz_1[i] * fe_0 + ta1_x_xy_xxzz_0[i] * pa_y[i] -
                              ta1_x_xy_xxzz_1[i] * pc_y[i];

        ta1_x_xyy_xyyy_0[i] = ta1_x_yy_yyy_0[i] * fe_0 - ta1_x_yy_yyy_1[i] * fe_0 + ta_yy_xyyy_1[i] +
                              ta1_x_yy_xyyy_0[i] * pa_x[i] - ta1_x_yy_xyyy_1[i] * pc_x[i];

        ta1_x_xyy_xyyz_0[i] = ta1_x_yy_yyz_0[i] * fe_0 - ta1_x_yy_yyz_1[i] * fe_0 + ta_yy_xyyz_1[i] +
                              ta1_x_yy_xyyz_0[i] * pa_x[i] - ta1_x_yy_xyyz_1[i] * pc_x[i];

        ta1_x_xyy_xyzz_0[i] = ta1_x_yy_yzz_0[i] * fe_0 - ta1_x_yy_yzz_1[i] * fe_0 + ta_yy_xyzz_1[i] +
                              ta1_x_yy_xyzz_0[i] * pa_x[i] - ta1_x_yy_xyzz_1[i] * pc_x[i];

        ta1_x_xyy_xzzz_0[i] = ta1_x_x_xzzz_0[i] * fe_0 - ta1_x_x_xzzz_1[i] * fe_0 + ta1_x_xy_xzzz_0[i] * pa_y[i] -
                              ta1_x_xy_xzzz_1[i] * pc_y[i];

        ta1_x_xyy_yyyy_0[i] = ta_yy_yyyy_1[i] + ta1_x_yy_yyyy_0[i] * pa_x[i] - ta1_x_yy_yyyy_1[i] * pc_x[i];

        ta1_x_xyy_yyyz_0[i] = ta_yy_yyyz_1[i] + ta1_x_yy_yyyz_0[i] * pa_x[i] - ta1_x_yy_yyyz_1[i] * pc_x[i];

        ta1_x_xyy_yyzz_0[i] = ta_yy_yyzz_1[i] + ta1_x_yy_yyzz_0[i] * pa_x[i] - ta1_x_yy_yyzz_1[i] * pc_x[i];

        ta1_x_xyy_yzzz_0[i] = ta_yy_yzzz_1[i] + ta1_x_yy_yzzz_0[i] * pa_x[i] - ta1_x_yy_yzzz_1[i] * pc_x[i];

        ta1_x_xyy_zzzz_0[i] = ta_yy_zzzz_1[i] + ta1_x_yy_zzzz_0[i] * pa_x[i] - ta1_x_yy_zzzz_1[i] * pc_x[i];
    }

    // Set up 60-75 components of targeted buffer : FG

    auto ta1_x_xyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 60);

    auto ta1_x_xyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 61);

    auto ta1_x_xyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 62);

    auto ta1_x_xyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 63);

    auto ta1_x_xyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 64);

    auto ta1_x_xyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 65);

    auto ta1_x_xyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 66);

    auto ta1_x_xyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 67);

    auto ta1_x_xyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 68);

    auto ta1_x_xyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 69);

    auto ta1_x_xyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 70);

    auto ta1_x_xyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 71);

    auto ta1_x_xyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 72);

    auto ta1_x_xyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 73);

    auto ta1_x_xyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 74);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             pc_x,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_x_xy_xxxy_0,  \
                             ta1_x_xy_xxxy_1,  \
                             ta1_x_xy_xxyy_0,  \
                             ta1_x_xy_xxyy_1,  \
                             ta1_x_xy_xyyy_0,  \
                             ta1_x_xy_xyyy_1,  \
                             ta1_x_xy_yyyy_0,  \
                             ta1_x_xy_yyyy_1,  \
                             ta1_x_xyz_xxxx_0, \
                             ta1_x_xyz_xxxy_0, \
                             ta1_x_xyz_xxxz_0, \
                             ta1_x_xyz_xxyy_0, \
                             ta1_x_xyz_xxyz_0, \
                             ta1_x_xyz_xxzz_0, \
                             ta1_x_xyz_xyyy_0, \
                             ta1_x_xyz_xyyz_0, \
                             ta1_x_xyz_xyzz_0, \
                             ta1_x_xyz_xzzz_0, \
                             ta1_x_xyz_yyyy_0, \
                             ta1_x_xyz_yyyz_0, \
                             ta1_x_xyz_yyzz_0, \
                             ta1_x_xyz_yzzz_0, \
                             ta1_x_xyz_zzzz_0, \
                             ta1_x_xz_xxxx_0,  \
                             ta1_x_xz_xxxx_1,  \
                             ta1_x_xz_xxxz_0,  \
                             ta1_x_xz_xxxz_1,  \
                             ta1_x_xz_xxyz_0,  \
                             ta1_x_xz_xxyz_1,  \
                             ta1_x_xz_xxz_0,   \
                             ta1_x_xz_xxz_1,   \
                             ta1_x_xz_xxzz_0,  \
                             ta1_x_xz_xxzz_1,  \
                             ta1_x_xz_xyyz_0,  \
                             ta1_x_xz_xyyz_1,  \
                             ta1_x_xz_xyz_0,   \
                             ta1_x_xz_xyz_1,   \
                             ta1_x_xz_xyzz_0,  \
                             ta1_x_xz_xyzz_1,  \
                             ta1_x_xz_xzz_0,   \
                             ta1_x_xz_xzz_1,   \
                             ta1_x_xz_xzzz_0,  \
                             ta1_x_xz_xzzz_1,  \
                             ta1_x_xz_zzzz_0,  \
                             ta1_x_xz_zzzz_1,  \
                             ta1_x_yz_yyyz_0,  \
                             ta1_x_yz_yyyz_1,  \
                             ta1_x_yz_yyzz_0,  \
                             ta1_x_yz_yyzz_1,  \
                             ta1_x_yz_yzzz_0,  \
                             ta1_x_yz_yzzz_1,  \
                             ta_yz_yyyz_1,     \
                             ta_yz_yyzz_1,     \
                             ta_yz_yzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyz_xxxx_0[i] = ta1_x_xz_xxxx_0[i] * pa_y[i] - ta1_x_xz_xxxx_1[i] * pc_y[i];

        ta1_x_xyz_xxxy_0[i] = ta1_x_xy_xxxy_0[i] * pa_z[i] - ta1_x_xy_xxxy_1[i] * pc_z[i];

        ta1_x_xyz_xxxz_0[i] = ta1_x_xz_xxxz_0[i] * pa_y[i] - ta1_x_xz_xxxz_1[i] * pc_y[i];

        ta1_x_xyz_xxyy_0[i] = ta1_x_xy_xxyy_0[i] * pa_z[i] - ta1_x_xy_xxyy_1[i] * pc_z[i];

        ta1_x_xyz_xxyz_0[i] = ta1_x_xz_xxz_0[i] * fe_0 - ta1_x_xz_xxz_1[i] * fe_0 + ta1_x_xz_xxyz_0[i] * pa_y[i] -
                              ta1_x_xz_xxyz_1[i] * pc_y[i];

        ta1_x_xyz_xxzz_0[i] = ta1_x_xz_xxzz_0[i] * pa_y[i] - ta1_x_xz_xxzz_1[i] * pc_y[i];

        ta1_x_xyz_xyyy_0[i] = ta1_x_xy_xyyy_0[i] * pa_z[i] - ta1_x_xy_xyyy_1[i] * pc_z[i];

        ta1_x_xyz_xyyz_0[i] = 2.0 * ta1_x_xz_xyz_0[i] * fe_0 - 2.0 * ta1_x_xz_xyz_1[i] * fe_0 +
                              ta1_x_xz_xyyz_0[i] * pa_y[i] - ta1_x_xz_xyyz_1[i] * pc_y[i];

        ta1_x_xyz_xyzz_0[i] = ta1_x_xz_xzz_0[i] * fe_0 - ta1_x_xz_xzz_1[i] * fe_0 + ta1_x_xz_xyzz_0[i] * pa_y[i] -
                              ta1_x_xz_xyzz_1[i] * pc_y[i];

        ta1_x_xyz_xzzz_0[i] = ta1_x_xz_xzzz_0[i] * pa_y[i] - ta1_x_xz_xzzz_1[i] * pc_y[i];

        ta1_x_xyz_yyyy_0[i] = ta1_x_xy_yyyy_0[i] * pa_z[i] - ta1_x_xy_yyyy_1[i] * pc_z[i];

        ta1_x_xyz_yyyz_0[i] = ta_yz_yyyz_1[i] + ta1_x_yz_yyyz_0[i] * pa_x[i] - ta1_x_yz_yyyz_1[i] * pc_x[i];

        ta1_x_xyz_yyzz_0[i] = ta_yz_yyzz_1[i] + ta1_x_yz_yyzz_0[i] * pa_x[i] - ta1_x_yz_yyzz_1[i] * pc_x[i];

        ta1_x_xyz_yzzz_0[i] = ta_yz_yzzz_1[i] + ta1_x_yz_yzzz_0[i] * pa_x[i] - ta1_x_yz_yzzz_1[i] * pc_x[i];

        ta1_x_xyz_zzzz_0[i] = ta1_x_xz_zzzz_0[i] * pa_y[i] - ta1_x_xz_zzzz_1[i] * pc_y[i];
    }

    // Set up 75-90 components of targeted buffer : FG

    auto ta1_x_xzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 75);

    auto ta1_x_xzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 76);

    auto ta1_x_xzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 77);

    auto ta1_x_xzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 78);

    auto ta1_x_xzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 79);

    auto ta1_x_xzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 80);

    auto ta1_x_xzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 81);

    auto ta1_x_xzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 82);

    auto ta1_x_xzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 83);

    auto ta1_x_xzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 84);

    auto ta1_x_xzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 85);

    auto ta1_x_xzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 86);

    auto ta1_x_xzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 87);

    auto ta1_x_xzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 88);

    auto ta1_x_xzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 89);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_x_x_xxxx_0,   \
                             ta1_x_x_xxxx_1,   \
                             ta1_x_x_xxxy_0,   \
                             ta1_x_x_xxxy_1,   \
                             ta1_x_x_xxyy_0,   \
                             ta1_x_x_xxyy_1,   \
                             ta1_x_x_xyyy_0,   \
                             ta1_x_x_xyyy_1,   \
                             ta1_x_xz_xxxx_0,  \
                             ta1_x_xz_xxxx_1,  \
                             ta1_x_xz_xxxy_0,  \
                             ta1_x_xz_xxxy_1,  \
                             ta1_x_xz_xxyy_0,  \
                             ta1_x_xz_xxyy_1,  \
                             ta1_x_xz_xyyy_0,  \
                             ta1_x_xz_xyyy_1,  \
                             ta1_x_xzz_xxxx_0, \
                             ta1_x_xzz_xxxy_0, \
                             ta1_x_xzz_xxxz_0, \
                             ta1_x_xzz_xxyy_0, \
                             ta1_x_xzz_xxyz_0, \
                             ta1_x_xzz_xxzz_0, \
                             ta1_x_xzz_xyyy_0, \
                             ta1_x_xzz_xyyz_0, \
                             ta1_x_xzz_xyzz_0, \
                             ta1_x_xzz_xzzz_0, \
                             ta1_x_xzz_yyyy_0, \
                             ta1_x_xzz_yyyz_0, \
                             ta1_x_xzz_yyzz_0, \
                             ta1_x_xzz_yzzz_0, \
                             ta1_x_xzz_zzzz_0, \
                             ta1_x_zz_xxxz_0,  \
                             ta1_x_zz_xxxz_1,  \
                             ta1_x_zz_xxyz_0,  \
                             ta1_x_zz_xxyz_1,  \
                             ta1_x_zz_xxz_0,   \
                             ta1_x_zz_xxz_1,   \
                             ta1_x_zz_xxzz_0,  \
                             ta1_x_zz_xxzz_1,  \
                             ta1_x_zz_xyyz_0,  \
                             ta1_x_zz_xyyz_1,  \
                             ta1_x_zz_xyz_0,   \
                             ta1_x_zz_xyz_1,   \
                             ta1_x_zz_xyzz_0,  \
                             ta1_x_zz_xyzz_1,  \
                             ta1_x_zz_xzz_0,   \
                             ta1_x_zz_xzz_1,   \
                             ta1_x_zz_xzzz_0,  \
                             ta1_x_zz_xzzz_1,  \
                             ta1_x_zz_yyyy_0,  \
                             ta1_x_zz_yyyy_1,  \
                             ta1_x_zz_yyyz_0,  \
                             ta1_x_zz_yyyz_1,  \
                             ta1_x_zz_yyz_0,   \
                             ta1_x_zz_yyz_1,   \
                             ta1_x_zz_yyzz_0,  \
                             ta1_x_zz_yyzz_1,  \
                             ta1_x_zz_yzz_0,   \
                             ta1_x_zz_yzz_1,   \
                             ta1_x_zz_yzzz_0,  \
                             ta1_x_zz_yzzz_1,  \
                             ta1_x_zz_zzz_0,   \
                             ta1_x_zz_zzz_1,   \
                             ta1_x_zz_zzzz_0,  \
                             ta1_x_zz_zzzz_1,  \
                             ta_zz_xxxz_1,     \
                             ta_zz_xxyz_1,     \
                             ta_zz_xxzz_1,     \
                             ta_zz_xyyz_1,     \
                             ta_zz_xyzz_1,     \
                             ta_zz_xzzz_1,     \
                             ta_zz_yyyy_1,     \
                             ta_zz_yyyz_1,     \
                             ta_zz_yyzz_1,     \
                             ta_zz_yzzz_1,     \
                             ta_zz_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xzz_xxxx_0[i] = ta1_x_x_xxxx_0[i] * fe_0 - ta1_x_x_xxxx_1[i] * fe_0 + ta1_x_xz_xxxx_0[i] * pa_z[i] -
                              ta1_x_xz_xxxx_1[i] * pc_z[i];

        ta1_x_xzz_xxxy_0[i] = ta1_x_x_xxxy_0[i] * fe_0 - ta1_x_x_xxxy_1[i] * fe_0 + ta1_x_xz_xxxy_0[i] * pa_z[i] -
                              ta1_x_xz_xxxy_1[i] * pc_z[i];

        ta1_x_xzz_xxxz_0[i] = 3.0 * ta1_x_zz_xxz_0[i] * fe_0 - 3.0 * ta1_x_zz_xxz_1[i] * fe_0 + ta_zz_xxxz_1[i] +
                              ta1_x_zz_xxxz_0[i] * pa_x[i] - ta1_x_zz_xxxz_1[i] * pc_x[i];

        ta1_x_xzz_xxyy_0[i] = ta1_x_x_xxyy_0[i] * fe_0 - ta1_x_x_xxyy_1[i] * fe_0 + ta1_x_xz_xxyy_0[i] * pa_z[i] -
                              ta1_x_xz_xxyy_1[i] * pc_z[i];

        ta1_x_xzz_xxyz_0[i] = 2.0 * ta1_x_zz_xyz_0[i] * fe_0 - 2.0 * ta1_x_zz_xyz_1[i] * fe_0 + ta_zz_xxyz_1[i] +
                              ta1_x_zz_xxyz_0[i] * pa_x[i] - ta1_x_zz_xxyz_1[i] * pc_x[i];

        ta1_x_xzz_xxzz_0[i] = 2.0 * ta1_x_zz_xzz_0[i] * fe_0 - 2.0 * ta1_x_zz_xzz_1[i] * fe_0 + ta_zz_xxzz_1[i] +
                              ta1_x_zz_xxzz_0[i] * pa_x[i] - ta1_x_zz_xxzz_1[i] * pc_x[i];

        ta1_x_xzz_xyyy_0[i] = ta1_x_x_xyyy_0[i] * fe_0 - ta1_x_x_xyyy_1[i] * fe_0 + ta1_x_xz_xyyy_0[i] * pa_z[i] -
                              ta1_x_xz_xyyy_1[i] * pc_z[i];

        ta1_x_xzz_xyyz_0[i] = ta1_x_zz_yyz_0[i] * fe_0 - ta1_x_zz_yyz_1[i] * fe_0 + ta_zz_xyyz_1[i] +
                              ta1_x_zz_xyyz_0[i] * pa_x[i] - ta1_x_zz_xyyz_1[i] * pc_x[i];

        ta1_x_xzz_xyzz_0[i] = ta1_x_zz_yzz_0[i] * fe_0 - ta1_x_zz_yzz_1[i] * fe_0 + ta_zz_xyzz_1[i] +
                              ta1_x_zz_xyzz_0[i] * pa_x[i] - ta1_x_zz_xyzz_1[i] * pc_x[i];

        ta1_x_xzz_xzzz_0[i] = ta1_x_zz_zzz_0[i] * fe_0 - ta1_x_zz_zzz_1[i] * fe_0 + ta_zz_xzzz_1[i] +
                              ta1_x_zz_xzzz_0[i] * pa_x[i] - ta1_x_zz_xzzz_1[i] * pc_x[i];

        ta1_x_xzz_yyyy_0[i] = ta_zz_yyyy_1[i] + ta1_x_zz_yyyy_0[i] * pa_x[i] - ta1_x_zz_yyyy_1[i] * pc_x[i];

        ta1_x_xzz_yyyz_0[i] = ta_zz_yyyz_1[i] + ta1_x_zz_yyyz_0[i] * pa_x[i] - ta1_x_zz_yyyz_1[i] * pc_x[i];

        ta1_x_xzz_yyzz_0[i] = ta_zz_yyzz_1[i] + ta1_x_zz_yyzz_0[i] * pa_x[i] - ta1_x_zz_yyzz_1[i] * pc_x[i];

        ta1_x_xzz_yzzz_0[i] = ta_zz_yzzz_1[i] + ta1_x_zz_yzzz_0[i] * pa_x[i] - ta1_x_zz_yzzz_1[i] * pc_x[i];

        ta1_x_xzz_zzzz_0[i] = ta_zz_zzzz_1[i] + ta1_x_zz_zzzz_0[i] * pa_x[i] - ta1_x_zz_zzzz_1[i] * pc_x[i];
    }

    // Set up 90-105 components of targeted buffer : FG

    auto ta1_x_yyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 90);

    auto ta1_x_yyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 91);

    auto ta1_x_yyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 92);

    auto ta1_x_yyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 93);

    auto ta1_x_yyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 94);

    auto ta1_x_yyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 95);

    auto ta1_x_yyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 96);

    auto ta1_x_yyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 97);

    auto ta1_x_yyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 98);

    auto ta1_x_yyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 99);

    auto ta1_x_yyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 100);

    auto ta1_x_yyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 101);

    auto ta1_x_yyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 102);

    auto ta1_x_yyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 103);

    auto ta1_x_yyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 104);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_x_y_xxxx_0,   \
                             ta1_x_y_xxxx_1,   \
                             ta1_x_y_xxxy_0,   \
                             ta1_x_y_xxxy_1,   \
                             ta1_x_y_xxxz_0,   \
                             ta1_x_y_xxxz_1,   \
                             ta1_x_y_xxyy_0,   \
                             ta1_x_y_xxyy_1,   \
                             ta1_x_y_xxyz_0,   \
                             ta1_x_y_xxyz_1,   \
                             ta1_x_y_xxzz_0,   \
                             ta1_x_y_xxzz_1,   \
                             ta1_x_y_xyyy_0,   \
                             ta1_x_y_xyyy_1,   \
                             ta1_x_y_xyyz_0,   \
                             ta1_x_y_xyyz_1,   \
                             ta1_x_y_xyzz_0,   \
                             ta1_x_y_xyzz_1,   \
                             ta1_x_y_xzzz_0,   \
                             ta1_x_y_xzzz_1,   \
                             ta1_x_y_yyyy_0,   \
                             ta1_x_y_yyyy_1,   \
                             ta1_x_y_yyyz_0,   \
                             ta1_x_y_yyyz_1,   \
                             ta1_x_y_yyzz_0,   \
                             ta1_x_y_yyzz_1,   \
                             ta1_x_y_yzzz_0,   \
                             ta1_x_y_yzzz_1,   \
                             ta1_x_y_zzzz_0,   \
                             ta1_x_y_zzzz_1,   \
                             ta1_x_yy_xxx_0,   \
                             ta1_x_yy_xxx_1,   \
                             ta1_x_yy_xxxx_0,  \
                             ta1_x_yy_xxxx_1,  \
                             ta1_x_yy_xxxy_0,  \
                             ta1_x_yy_xxxy_1,  \
                             ta1_x_yy_xxxz_0,  \
                             ta1_x_yy_xxxz_1,  \
                             ta1_x_yy_xxy_0,   \
                             ta1_x_yy_xxy_1,   \
                             ta1_x_yy_xxyy_0,  \
                             ta1_x_yy_xxyy_1,  \
                             ta1_x_yy_xxyz_0,  \
                             ta1_x_yy_xxyz_1,  \
                             ta1_x_yy_xxz_0,   \
                             ta1_x_yy_xxz_1,   \
                             ta1_x_yy_xxzz_0,  \
                             ta1_x_yy_xxzz_1,  \
                             ta1_x_yy_xyy_0,   \
                             ta1_x_yy_xyy_1,   \
                             ta1_x_yy_xyyy_0,  \
                             ta1_x_yy_xyyy_1,  \
                             ta1_x_yy_xyyz_0,  \
                             ta1_x_yy_xyyz_1,  \
                             ta1_x_yy_xyz_0,   \
                             ta1_x_yy_xyz_1,   \
                             ta1_x_yy_xyzz_0,  \
                             ta1_x_yy_xyzz_1,  \
                             ta1_x_yy_xzz_0,   \
                             ta1_x_yy_xzz_1,   \
                             ta1_x_yy_xzzz_0,  \
                             ta1_x_yy_xzzz_1,  \
                             ta1_x_yy_yyy_0,   \
                             ta1_x_yy_yyy_1,   \
                             ta1_x_yy_yyyy_0,  \
                             ta1_x_yy_yyyy_1,  \
                             ta1_x_yy_yyyz_0,  \
                             ta1_x_yy_yyyz_1,  \
                             ta1_x_yy_yyz_0,   \
                             ta1_x_yy_yyz_1,   \
                             ta1_x_yy_yyzz_0,  \
                             ta1_x_yy_yyzz_1,  \
                             ta1_x_yy_yzz_0,   \
                             ta1_x_yy_yzz_1,   \
                             ta1_x_yy_yzzz_0,  \
                             ta1_x_yy_yzzz_1,  \
                             ta1_x_yy_zzz_0,   \
                             ta1_x_yy_zzz_1,   \
                             ta1_x_yy_zzzz_0,  \
                             ta1_x_yy_zzzz_1,  \
                             ta1_x_yyy_xxxx_0, \
                             ta1_x_yyy_xxxy_0, \
                             ta1_x_yyy_xxxz_0, \
                             ta1_x_yyy_xxyy_0, \
                             ta1_x_yyy_xxyz_0, \
                             ta1_x_yyy_xxzz_0, \
                             ta1_x_yyy_xyyy_0, \
                             ta1_x_yyy_xyyz_0, \
                             ta1_x_yyy_xyzz_0, \
                             ta1_x_yyy_xzzz_0, \
                             ta1_x_yyy_yyyy_0, \
                             ta1_x_yyy_yyyz_0, \
                             ta1_x_yyy_yyzz_0, \
                             ta1_x_yyy_yzzz_0, \
                             ta1_x_yyy_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyy_xxxx_0[i] = 2.0 * ta1_x_y_xxxx_0[i] * fe_0 - 2.0 * ta1_x_y_xxxx_1[i] * fe_0 +
                              ta1_x_yy_xxxx_0[i] * pa_y[i] - ta1_x_yy_xxxx_1[i] * pc_y[i];

        ta1_x_yyy_xxxy_0[i] = 2.0 * ta1_x_y_xxxy_0[i] * fe_0 - 2.0 * ta1_x_y_xxxy_1[i] * fe_0 +
                              ta1_x_yy_xxx_0[i] * fe_0 - ta1_x_yy_xxx_1[i] * fe_0 + ta1_x_yy_xxxy_0[i] * pa_y[i] -
                              ta1_x_yy_xxxy_1[i] * pc_y[i];

        ta1_x_yyy_xxxz_0[i] = 2.0 * ta1_x_y_xxxz_0[i] * fe_0 - 2.0 * ta1_x_y_xxxz_1[i] * fe_0 +
                              ta1_x_yy_xxxz_0[i] * pa_y[i] - ta1_x_yy_xxxz_1[i] * pc_y[i];

        ta1_x_yyy_xxyy_0[i] = 2.0 * ta1_x_y_xxyy_0[i] * fe_0 - 2.0 * ta1_x_y_xxyy_1[i] * fe_0 +
                              2.0 * ta1_x_yy_xxy_0[i] * fe_0 - 2.0 * ta1_x_yy_xxy_1[i] * fe_0 +
                              ta1_x_yy_xxyy_0[i] * pa_y[i] - ta1_x_yy_xxyy_1[i] * pc_y[i];

        ta1_x_yyy_xxyz_0[i] = 2.0 * ta1_x_y_xxyz_0[i] * fe_0 - 2.0 * ta1_x_y_xxyz_1[i] * fe_0 +
                              ta1_x_yy_xxz_0[i] * fe_0 - ta1_x_yy_xxz_1[i] * fe_0 + ta1_x_yy_xxyz_0[i] * pa_y[i] -
                              ta1_x_yy_xxyz_1[i] * pc_y[i];

        ta1_x_yyy_xxzz_0[i] = 2.0 * ta1_x_y_xxzz_0[i] * fe_0 - 2.0 * ta1_x_y_xxzz_1[i] * fe_0 +
                              ta1_x_yy_xxzz_0[i] * pa_y[i] - ta1_x_yy_xxzz_1[i] * pc_y[i];

        ta1_x_yyy_xyyy_0[i] = 2.0 * ta1_x_y_xyyy_0[i] * fe_0 - 2.0 * ta1_x_y_xyyy_1[i] * fe_0 +
                              3.0 * ta1_x_yy_xyy_0[i] * fe_0 - 3.0 * ta1_x_yy_xyy_1[i] * fe_0 +
                              ta1_x_yy_xyyy_0[i] * pa_y[i] - ta1_x_yy_xyyy_1[i] * pc_y[i];

        ta1_x_yyy_xyyz_0[i] = 2.0 * ta1_x_y_xyyz_0[i] * fe_0 - 2.0 * ta1_x_y_xyyz_1[i] * fe_0 +
                              2.0 * ta1_x_yy_xyz_0[i] * fe_0 - 2.0 * ta1_x_yy_xyz_1[i] * fe_0 +
                              ta1_x_yy_xyyz_0[i] * pa_y[i] - ta1_x_yy_xyyz_1[i] * pc_y[i];

        ta1_x_yyy_xyzz_0[i] = 2.0 * ta1_x_y_xyzz_0[i] * fe_0 - 2.0 * ta1_x_y_xyzz_1[i] * fe_0 +
                              ta1_x_yy_xzz_0[i] * fe_0 - ta1_x_yy_xzz_1[i] * fe_0 + ta1_x_yy_xyzz_0[i] * pa_y[i] -
                              ta1_x_yy_xyzz_1[i] * pc_y[i];

        ta1_x_yyy_xzzz_0[i] = 2.0 * ta1_x_y_xzzz_0[i] * fe_0 - 2.0 * ta1_x_y_xzzz_1[i] * fe_0 +
                              ta1_x_yy_xzzz_0[i] * pa_y[i] - ta1_x_yy_xzzz_1[i] * pc_y[i];

        ta1_x_yyy_yyyy_0[i] = 2.0 * ta1_x_y_yyyy_0[i] * fe_0 - 2.0 * ta1_x_y_yyyy_1[i] * fe_0 +
                              4.0 * ta1_x_yy_yyy_0[i] * fe_0 - 4.0 * ta1_x_yy_yyy_1[i] * fe_0 +
                              ta1_x_yy_yyyy_0[i] * pa_y[i] - ta1_x_yy_yyyy_1[i] * pc_y[i];

        ta1_x_yyy_yyyz_0[i] = 2.0 * ta1_x_y_yyyz_0[i] * fe_0 - 2.0 * ta1_x_y_yyyz_1[i] * fe_0 +
                              3.0 * ta1_x_yy_yyz_0[i] * fe_0 - 3.0 * ta1_x_yy_yyz_1[i] * fe_0 +
                              ta1_x_yy_yyyz_0[i] * pa_y[i] - ta1_x_yy_yyyz_1[i] * pc_y[i];

        ta1_x_yyy_yyzz_0[i] = 2.0 * ta1_x_y_yyzz_0[i] * fe_0 - 2.0 * ta1_x_y_yyzz_1[i] * fe_0 +
                              2.0 * ta1_x_yy_yzz_0[i] * fe_0 - 2.0 * ta1_x_yy_yzz_1[i] * fe_0 +
                              ta1_x_yy_yyzz_0[i] * pa_y[i] - ta1_x_yy_yyzz_1[i] * pc_y[i];

        ta1_x_yyy_yzzz_0[i] = 2.0 * ta1_x_y_yzzz_0[i] * fe_0 - 2.0 * ta1_x_y_yzzz_1[i] * fe_0 +
                              ta1_x_yy_zzz_0[i] * fe_0 - ta1_x_yy_zzz_1[i] * fe_0 + ta1_x_yy_yzzz_0[i] * pa_y[i] -
                              ta1_x_yy_yzzz_1[i] * pc_y[i];

        ta1_x_yyy_zzzz_0[i] = 2.0 * ta1_x_y_zzzz_0[i] * fe_0 - 2.0 * ta1_x_y_zzzz_1[i] * fe_0 +
                              ta1_x_yy_zzzz_0[i] * pa_y[i] - ta1_x_yy_zzzz_1[i] * pc_y[i];
    }

    // Set up 105-120 components of targeted buffer : FG

    auto ta1_x_yyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 105);

    auto ta1_x_yyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 106);

    auto ta1_x_yyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 107);

    auto ta1_x_yyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 108);

    auto ta1_x_yyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 109);

    auto ta1_x_yyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 110);

    auto ta1_x_yyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 111);

    auto ta1_x_yyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 112);

    auto ta1_x_yyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 113);

    auto ta1_x_yyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 114);

    auto ta1_x_yyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 115);

    auto ta1_x_yyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 116);

    auto ta1_x_yyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 117);

    auto ta1_x_yyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 118);

    auto ta1_x_yyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 119);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_x_yy_xxxx_0,  \
                             ta1_x_yy_xxxx_1,  \
                             ta1_x_yy_xxxy_0,  \
                             ta1_x_yy_xxxy_1,  \
                             ta1_x_yy_xxy_0,   \
                             ta1_x_yy_xxy_1,   \
                             ta1_x_yy_xxyy_0,  \
                             ta1_x_yy_xxyy_1,  \
                             ta1_x_yy_xxyz_0,  \
                             ta1_x_yy_xxyz_1,  \
                             ta1_x_yy_xyy_0,   \
                             ta1_x_yy_xyy_1,   \
                             ta1_x_yy_xyyy_0,  \
                             ta1_x_yy_xyyy_1,  \
                             ta1_x_yy_xyyz_0,  \
                             ta1_x_yy_xyyz_1,  \
                             ta1_x_yy_xyz_0,   \
                             ta1_x_yy_xyz_1,   \
                             ta1_x_yy_xyzz_0,  \
                             ta1_x_yy_xyzz_1,  \
                             ta1_x_yy_yyy_0,   \
                             ta1_x_yy_yyy_1,   \
                             ta1_x_yy_yyyy_0,  \
                             ta1_x_yy_yyyy_1,  \
                             ta1_x_yy_yyyz_0,  \
                             ta1_x_yy_yyyz_1,  \
                             ta1_x_yy_yyz_0,   \
                             ta1_x_yy_yyz_1,   \
                             ta1_x_yy_yyzz_0,  \
                             ta1_x_yy_yyzz_1,  \
                             ta1_x_yy_yzz_0,   \
                             ta1_x_yy_yzz_1,   \
                             ta1_x_yy_yzzz_0,  \
                             ta1_x_yy_yzzz_1,  \
                             ta1_x_yyz_xxxx_0, \
                             ta1_x_yyz_xxxy_0, \
                             ta1_x_yyz_xxxz_0, \
                             ta1_x_yyz_xxyy_0, \
                             ta1_x_yyz_xxyz_0, \
                             ta1_x_yyz_xxzz_0, \
                             ta1_x_yyz_xyyy_0, \
                             ta1_x_yyz_xyyz_0, \
                             ta1_x_yyz_xyzz_0, \
                             ta1_x_yyz_xzzz_0, \
                             ta1_x_yyz_yyyy_0, \
                             ta1_x_yyz_yyyz_0, \
                             ta1_x_yyz_yyzz_0, \
                             ta1_x_yyz_yzzz_0, \
                             ta1_x_yyz_zzzz_0, \
                             ta1_x_yz_xxxz_0,  \
                             ta1_x_yz_xxxz_1,  \
                             ta1_x_yz_xxzz_0,  \
                             ta1_x_yz_xxzz_1,  \
                             ta1_x_yz_xzzz_0,  \
                             ta1_x_yz_xzzz_1,  \
                             ta1_x_yz_zzzz_0,  \
                             ta1_x_yz_zzzz_1,  \
                             ta1_x_z_xxxz_0,   \
                             ta1_x_z_xxxz_1,   \
                             ta1_x_z_xxzz_0,   \
                             ta1_x_z_xxzz_1,   \
                             ta1_x_z_xzzz_0,   \
                             ta1_x_z_xzzz_1,   \
                             ta1_x_z_zzzz_0,   \
                             ta1_x_z_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyz_xxxx_0[i] = ta1_x_yy_xxxx_0[i] * pa_z[i] - ta1_x_yy_xxxx_1[i] * pc_z[i];

        ta1_x_yyz_xxxy_0[i] = ta1_x_yy_xxxy_0[i] * pa_z[i] - ta1_x_yy_xxxy_1[i] * pc_z[i];

        ta1_x_yyz_xxxz_0[i] = ta1_x_z_xxxz_0[i] * fe_0 - ta1_x_z_xxxz_1[i] * fe_0 + ta1_x_yz_xxxz_0[i] * pa_y[i] -
                              ta1_x_yz_xxxz_1[i] * pc_y[i];

        ta1_x_yyz_xxyy_0[i] = ta1_x_yy_xxyy_0[i] * pa_z[i] - ta1_x_yy_xxyy_1[i] * pc_z[i];

        ta1_x_yyz_xxyz_0[i] = ta1_x_yy_xxy_0[i] * fe_0 - ta1_x_yy_xxy_1[i] * fe_0 + ta1_x_yy_xxyz_0[i] * pa_z[i] -
                              ta1_x_yy_xxyz_1[i] * pc_z[i];

        ta1_x_yyz_xxzz_0[i] = ta1_x_z_xxzz_0[i] * fe_0 - ta1_x_z_xxzz_1[i] * fe_0 + ta1_x_yz_xxzz_0[i] * pa_y[i] -
                              ta1_x_yz_xxzz_1[i] * pc_y[i];

        ta1_x_yyz_xyyy_0[i] = ta1_x_yy_xyyy_0[i] * pa_z[i] - ta1_x_yy_xyyy_1[i] * pc_z[i];

        ta1_x_yyz_xyyz_0[i] = ta1_x_yy_xyy_0[i] * fe_0 - ta1_x_yy_xyy_1[i] * fe_0 + ta1_x_yy_xyyz_0[i] * pa_z[i] -
                              ta1_x_yy_xyyz_1[i] * pc_z[i];

        ta1_x_yyz_xyzz_0[i] = 2.0 * ta1_x_yy_xyz_0[i] * fe_0 - 2.0 * ta1_x_yy_xyz_1[i] * fe_0 +
                              ta1_x_yy_xyzz_0[i] * pa_z[i] - ta1_x_yy_xyzz_1[i] * pc_z[i];

        ta1_x_yyz_xzzz_0[i] = ta1_x_z_xzzz_0[i] * fe_0 - ta1_x_z_xzzz_1[i] * fe_0 + ta1_x_yz_xzzz_0[i] * pa_y[i] -
                              ta1_x_yz_xzzz_1[i] * pc_y[i];

        ta1_x_yyz_yyyy_0[i] = ta1_x_yy_yyyy_0[i] * pa_z[i] - ta1_x_yy_yyyy_1[i] * pc_z[i];

        ta1_x_yyz_yyyz_0[i] = ta1_x_yy_yyy_0[i] * fe_0 - ta1_x_yy_yyy_1[i] * fe_0 + ta1_x_yy_yyyz_0[i] * pa_z[i] -
                              ta1_x_yy_yyyz_1[i] * pc_z[i];

        ta1_x_yyz_yyzz_0[i] = 2.0 * ta1_x_yy_yyz_0[i] * fe_0 - 2.0 * ta1_x_yy_yyz_1[i] * fe_0 +
                              ta1_x_yy_yyzz_0[i] * pa_z[i] - ta1_x_yy_yyzz_1[i] * pc_z[i];

        ta1_x_yyz_yzzz_0[i] = 3.0 * ta1_x_yy_yzz_0[i] * fe_0 - 3.0 * ta1_x_yy_yzz_1[i] * fe_0 +
                              ta1_x_yy_yzzz_0[i] * pa_z[i] - ta1_x_yy_yzzz_1[i] * pc_z[i];

        ta1_x_yyz_zzzz_0[i] = ta1_x_z_zzzz_0[i] * fe_0 - ta1_x_z_zzzz_1[i] * fe_0 + ta1_x_yz_zzzz_0[i] * pa_y[i] -
                              ta1_x_yz_zzzz_1[i] * pc_y[i];
    }

    // Set up 120-135 components of targeted buffer : FG

    auto ta1_x_yzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 120);

    auto ta1_x_yzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 121);

    auto ta1_x_yzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 122);

    auto ta1_x_yzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 123);

    auto ta1_x_yzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 124);

    auto ta1_x_yzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 125);

    auto ta1_x_yzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 126);

    auto ta1_x_yzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 127);

    auto ta1_x_yzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 128);

    auto ta1_x_yzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 129);

    auto ta1_x_yzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 130);

    auto ta1_x_yzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 131);

    auto ta1_x_yzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 132);

    auto ta1_x_yzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 133);

    auto ta1_x_yzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 134);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_x_yzz_xxxx_0, \
                             ta1_x_yzz_xxxy_0, \
                             ta1_x_yzz_xxxz_0, \
                             ta1_x_yzz_xxyy_0, \
                             ta1_x_yzz_xxyz_0, \
                             ta1_x_yzz_xxzz_0, \
                             ta1_x_yzz_xyyy_0, \
                             ta1_x_yzz_xyyz_0, \
                             ta1_x_yzz_xyzz_0, \
                             ta1_x_yzz_xzzz_0, \
                             ta1_x_yzz_yyyy_0, \
                             ta1_x_yzz_yyyz_0, \
                             ta1_x_yzz_yyzz_0, \
                             ta1_x_yzz_yzzz_0, \
                             ta1_x_yzz_zzzz_0, \
                             ta1_x_zz_xxx_0,   \
                             ta1_x_zz_xxx_1,   \
                             ta1_x_zz_xxxx_0,  \
                             ta1_x_zz_xxxx_1,  \
                             ta1_x_zz_xxxy_0,  \
                             ta1_x_zz_xxxy_1,  \
                             ta1_x_zz_xxxz_0,  \
                             ta1_x_zz_xxxz_1,  \
                             ta1_x_zz_xxy_0,   \
                             ta1_x_zz_xxy_1,   \
                             ta1_x_zz_xxyy_0,  \
                             ta1_x_zz_xxyy_1,  \
                             ta1_x_zz_xxyz_0,  \
                             ta1_x_zz_xxyz_1,  \
                             ta1_x_zz_xxz_0,   \
                             ta1_x_zz_xxz_1,   \
                             ta1_x_zz_xxzz_0,  \
                             ta1_x_zz_xxzz_1,  \
                             ta1_x_zz_xyy_0,   \
                             ta1_x_zz_xyy_1,   \
                             ta1_x_zz_xyyy_0,  \
                             ta1_x_zz_xyyy_1,  \
                             ta1_x_zz_xyyz_0,  \
                             ta1_x_zz_xyyz_1,  \
                             ta1_x_zz_xyz_0,   \
                             ta1_x_zz_xyz_1,   \
                             ta1_x_zz_xyzz_0,  \
                             ta1_x_zz_xyzz_1,  \
                             ta1_x_zz_xzz_0,   \
                             ta1_x_zz_xzz_1,   \
                             ta1_x_zz_xzzz_0,  \
                             ta1_x_zz_xzzz_1,  \
                             ta1_x_zz_yyy_0,   \
                             ta1_x_zz_yyy_1,   \
                             ta1_x_zz_yyyy_0,  \
                             ta1_x_zz_yyyy_1,  \
                             ta1_x_zz_yyyz_0,  \
                             ta1_x_zz_yyyz_1,  \
                             ta1_x_zz_yyz_0,   \
                             ta1_x_zz_yyz_1,   \
                             ta1_x_zz_yyzz_0,  \
                             ta1_x_zz_yyzz_1,  \
                             ta1_x_zz_yzz_0,   \
                             ta1_x_zz_yzz_1,   \
                             ta1_x_zz_yzzz_0,  \
                             ta1_x_zz_yzzz_1,  \
                             ta1_x_zz_zzz_0,   \
                             ta1_x_zz_zzz_1,   \
                             ta1_x_zz_zzzz_0,  \
                             ta1_x_zz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yzz_xxxx_0[i] = ta1_x_zz_xxxx_0[i] * pa_y[i] - ta1_x_zz_xxxx_1[i] * pc_y[i];

        ta1_x_yzz_xxxy_0[i] = ta1_x_zz_xxx_0[i] * fe_0 - ta1_x_zz_xxx_1[i] * fe_0 + ta1_x_zz_xxxy_0[i] * pa_y[i] -
                              ta1_x_zz_xxxy_1[i] * pc_y[i];

        ta1_x_yzz_xxxz_0[i] = ta1_x_zz_xxxz_0[i] * pa_y[i] - ta1_x_zz_xxxz_1[i] * pc_y[i];

        ta1_x_yzz_xxyy_0[i] = 2.0 * ta1_x_zz_xxy_0[i] * fe_0 - 2.0 * ta1_x_zz_xxy_1[i] * fe_0 +
                              ta1_x_zz_xxyy_0[i] * pa_y[i] - ta1_x_zz_xxyy_1[i] * pc_y[i];

        ta1_x_yzz_xxyz_0[i] = ta1_x_zz_xxz_0[i] * fe_0 - ta1_x_zz_xxz_1[i] * fe_0 + ta1_x_zz_xxyz_0[i] * pa_y[i] -
                              ta1_x_zz_xxyz_1[i] * pc_y[i];

        ta1_x_yzz_xxzz_0[i] = ta1_x_zz_xxzz_0[i] * pa_y[i] - ta1_x_zz_xxzz_1[i] * pc_y[i];

        ta1_x_yzz_xyyy_0[i] = 3.0 * ta1_x_zz_xyy_0[i] * fe_0 - 3.0 * ta1_x_zz_xyy_1[i] * fe_0 +
                              ta1_x_zz_xyyy_0[i] * pa_y[i] - ta1_x_zz_xyyy_1[i] * pc_y[i];

        ta1_x_yzz_xyyz_0[i] = 2.0 * ta1_x_zz_xyz_0[i] * fe_0 - 2.0 * ta1_x_zz_xyz_1[i] * fe_0 +
                              ta1_x_zz_xyyz_0[i] * pa_y[i] - ta1_x_zz_xyyz_1[i] * pc_y[i];

        ta1_x_yzz_xyzz_0[i] = ta1_x_zz_xzz_0[i] * fe_0 - ta1_x_zz_xzz_1[i] * fe_0 + ta1_x_zz_xyzz_0[i] * pa_y[i] -
                              ta1_x_zz_xyzz_1[i] * pc_y[i];

        ta1_x_yzz_xzzz_0[i] = ta1_x_zz_xzzz_0[i] * pa_y[i] - ta1_x_zz_xzzz_1[i] * pc_y[i];

        ta1_x_yzz_yyyy_0[i] = 4.0 * ta1_x_zz_yyy_0[i] * fe_0 - 4.0 * ta1_x_zz_yyy_1[i] * fe_0 +
                              ta1_x_zz_yyyy_0[i] * pa_y[i] - ta1_x_zz_yyyy_1[i] * pc_y[i];

        ta1_x_yzz_yyyz_0[i] = 3.0 * ta1_x_zz_yyz_0[i] * fe_0 - 3.0 * ta1_x_zz_yyz_1[i] * fe_0 +
                              ta1_x_zz_yyyz_0[i] * pa_y[i] - ta1_x_zz_yyyz_1[i] * pc_y[i];

        ta1_x_yzz_yyzz_0[i] = 2.0 * ta1_x_zz_yzz_0[i] * fe_0 - 2.0 * ta1_x_zz_yzz_1[i] * fe_0 +
                              ta1_x_zz_yyzz_0[i] * pa_y[i] - ta1_x_zz_yyzz_1[i] * pc_y[i];

        ta1_x_yzz_yzzz_0[i] = ta1_x_zz_zzz_0[i] * fe_0 - ta1_x_zz_zzz_1[i] * fe_0 + ta1_x_zz_yzzz_0[i] * pa_y[i] -
                              ta1_x_zz_yzzz_1[i] * pc_y[i];

        ta1_x_yzz_zzzz_0[i] = ta1_x_zz_zzzz_0[i] * pa_y[i] - ta1_x_zz_zzzz_1[i] * pc_y[i];
    }

    // Set up 135-150 components of targeted buffer : FG

    auto ta1_x_zzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 135);

    auto ta1_x_zzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 136);

    auto ta1_x_zzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 137);

    auto ta1_x_zzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 138);

    auto ta1_x_zzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 139);

    auto ta1_x_zzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 140);

    auto ta1_x_zzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 141);

    auto ta1_x_zzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 142);

    auto ta1_x_zzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 143);

    auto ta1_x_zzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 144);

    auto ta1_x_zzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 145);

    auto ta1_x_zzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 146);

    auto ta1_x_zzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 147);

    auto ta1_x_zzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 148);

    auto ta1_x_zzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 149);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta1_x_z_xxxx_0,   \
                             ta1_x_z_xxxx_1,   \
                             ta1_x_z_xxxy_0,   \
                             ta1_x_z_xxxy_1,   \
                             ta1_x_z_xxxz_0,   \
                             ta1_x_z_xxxz_1,   \
                             ta1_x_z_xxyy_0,   \
                             ta1_x_z_xxyy_1,   \
                             ta1_x_z_xxyz_0,   \
                             ta1_x_z_xxyz_1,   \
                             ta1_x_z_xxzz_0,   \
                             ta1_x_z_xxzz_1,   \
                             ta1_x_z_xyyy_0,   \
                             ta1_x_z_xyyy_1,   \
                             ta1_x_z_xyyz_0,   \
                             ta1_x_z_xyyz_1,   \
                             ta1_x_z_xyzz_0,   \
                             ta1_x_z_xyzz_1,   \
                             ta1_x_z_xzzz_0,   \
                             ta1_x_z_xzzz_1,   \
                             ta1_x_z_yyyy_0,   \
                             ta1_x_z_yyyy_1,   \
                             ta1_x_z_yyyz_0,   \
                             ta1_x_z_yyyz_1,   \
                             ta1_x_z_yyzz_0,   \
                             ta1_x_z_yyzz_1,   \
                             ta1_x_z_yzzz_0,   \
                             ta1_x_z_yzzz_1,   \
                             ta1_x_z_zzzz_0,   \
                             ta1_x_z_zzzz_1,   \
                             ta1_x_zz_xxx_0,   \
                             ta1_x_zz_xxx_1,   \
                             ta1_x_zz_xxxx_0,  \
                             ta1_x_zz_xxxx_1,  \
                             ta1_x_zz_xxxy_0,  \
                             ta1_x_zz_xxxy_1,  \
                             ta1_x_zz_xxxz_0,  \
                             ta1_x_zz_xxxz_1,  \
                             ta1_x_zz_xxy_0,   \
                             ta1_x_zz_xxy_1,   \
                             ta1_x_zz_xxyy_0,  \
                             ta1_x_zz_xxyy_1,  \
                             ta1_x_zz_xxyz_0,  \
                             ta1_x_zz_xxyz_1,  \
                             ta1_x_zz_xxz_0,   \
                             ta1_x_zz_xxz_1,   \
                             ta1_x_zz_xxzz_0,  \
                             ta1_x_zz_xxzz_1,  \
                             ta1_x_zz_xyy_0,   \
                             ta1_x_zz_xyy_1,   \
                             ta1_x_zz_xyyy_0,  \
                             ta1_x_zz_xyyy_1,  \
                             ta1_x_zz_xyyz_0,  \
                             ta1_x_zz_xyyz_1,  \
                             ta1_x_zz_xyz_0,   \
                             ta1_x_zz_xyz_1,   \
                             ta1_x_zz_xyzz_0,  \
                             ta1_x_zz_xyzz_1,  \
                             ta1_x_zz_xzz_0,   \
                             ta1_x_zz_xzz_1,   \
                             ta1_x_zz_xzzz_0,  \
                             ta1_x_zz_xzzz_1,  \
                             ta1_x_zz_yyy_0,   \
                             ta1_x_zz_yyy_1,   \
                             ta1_x_zz_yyyy_0,  \
                             ta1_x_zz_yyyy_1,  \
                             ta1_x_zz_yyyz_0,  \
                             ta1_x_zz_yyyz_1,  \
                             ta1_x_zz_yyz_0,   \
                             ta1_x_zz_yyz_1,   \
                             ta1_x_zz_yyzz_0,  \
                             ta1_x_zz_yyzz_1,  \
                             ta1_x_zz_yzz_0,   \
                             ta1_x_zz_yzz_1,   \
                             ta1_x_zz_yzzz_0,  \
                             ta1_x_zz_yzzz_1,  \
                             ta1_x_zz_zzz_0,   \
                             ta1_x_zz_zzz_1,   \
                             ta1_x_zz_zzzz_0,  \
                             ta1_x_zz_zzzz_1,  \
                             ta1_x_zzz_xxxx_0, \
                             ta1_x_zzz_xxxy_0, \
                             ta1_x_zzz_xxxz_0, \
                             ta1_x_zzz_xxyy_0, \
                             ta1_x_zzz_xxyz_0, \
                             ta1_x_zzz_xxzz_0, \
                             ta1_x_zzz_xyyy_0, \
                             ta1_x_zzz_xyyz_0, \
                             ta1_x_zzz_xyzz_0, \
                             ta1_x_zzz_xzzz_0, \
                             ta1_x_zzz_yyyy_0, \
                             ta1_x_zzz_yyyz_0, \
                             ta1_x_zzz_yyzz_0, \
                             ta1_x_zzz_yzzz_0, \
                             ta1_x_zzz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_zzz_xxxx_0[i] = 2.0 * ta1_x_z_xxxx_0[i] * fe_0 - 2.0 * ta1_x_z_xxxx_1[i] * fe_0 +
                              ta1_x_zz_xxxx_0[i] * pa_z[i] - ta1_x_zz_xxxx_1[i] * pc_z[i];

        ta1_x_zzz_xxxy_0[i] = 2.0 * ta1_x_z_xxxy_0[i] * fe_0 - 2.0 * ta1_x_z_xxxy_1[i] * fe_0 +
                              ta1_x_zz_xxxy_0[i] * pa_z[i] - ta1_x_zz_xxxy_1[i] * pc_z[i];

        ta1_x_zzz_xxxz_0[i] = 2.0 * ta1_x_z_xxxz_0[i] * fe_0 - 2.0 * ta1_x_z_xxxz_1[i] * fe_0 +
                              ta1_x_zz_xxx_0[i] * fe_0 - ta1_x_zz_xxx_1[i] * fe_0 + ta1_x_zz_xxxz_0[i] * pa_z[i] -
                              ta1_x_zz_xxxz_1[i] * pc_z[i];

        ta1_x_zzz_xxyy_0[i] = 2.0 * ta1_x_z_xxyy_0[i] * fe_0 - 2.0 * ta1_x_z_xxyy_1[i] * fe_0 +
                              ta1_x_zz_xxyy_0[i] * pa_z[i] - ta1_x_zz_xxyy_1[i] * pc_z[i];

        ta1_x_zzz_xxyz_0[i] = 2.0 * ta1_x_z_xxyz_0[i] * fe_0 - 2.0 * ta1_x_z_xxyz_1[i] * fe_0 +
                              ta1_x_zz_xxy_0[i] * fe_0 - ta1_x_zz_xxy_1[i] * fe_0 + ta1_x_zz_xxyz_0[i] * pa_z[i] -
                              ta1_x_zz_xxyz_1[i] * pc_z[i];

        ta1_x_zzz_xxzz_0[i] = 2.0 * ta1_x_z_xxzz_0[i] * fe_0 - 2.0 * ta1_x_z_xxzz_1[i] * fe_0 +
                              2.0 * ta1_x_zz_xxz_0[i] * fe_0 - 2.0 * ta1_x_zz_xxz_1[i] * fe_0 +
                              ta1_x_zz_xxzz_0[i] * pa_z[i] - ta1_x_zz_xxzz_1[i] * pc_z[i];

        ta1_x_zzz_xyyy_0[i] = 2.0 * ta1_x_z_xyyy_0[i] * fe_0 - 2.0 * ta1_x_z_xyyy_1[i] * fe_0 +
                              ta1_x_zz_xyyy_0[i] * pa_z[i] - ta1_x_zz_xyyy_1[i] * pc_z[i];

        ta1_x_zzz_xyyz_0[i] = 2.0 * ta1_x_z_xyyz_0[i] * fe_0 - 2.0 * ta1_x_z_xyyz_1[i] * fe_0 +
                              ta1_x_zz_xyy_0[i] * fe_0 - ta1_x_zz_xyy_1[i] * fe_0 + ta1_x_zz_xyyz_0[i] * pa_z[i] -
                              ta1_x_zz_xyyz_1[i] * pc_z[i];

        ta1_x_zzz_xyzz_0[i] = 2.0 * ta1_x_z_xyzz_0[i] * fe_0 - 2.0 * ta1_x_z_xyzz_1[i] * fe_0 +
                              2.0 * ta1_x_zz_xyz_0[i] * fe_0 - 2.0 * ta1_x_zz_xyz_1[i] * fe_0 +
                              ta1_x_zz_xyzz_0[i] * pa_z[i] - ta1_x_zz_xyzz_1[i] * pc_z[i];

        ta1_x_zzz_xzzz_0[i] = 2.0 * ta1_x_z_xzzz_0[i] * fe_0 - 2.0 * ta1_x_z_xzzz_1[i] * fe_0 +
                              3.0 * ta1_x_zz_xzz_0[i] * fe_0 - 3.0 * ta1_x_zz_xzz_1[i] * fe_0 +
                              ta1_x_zz_xzzz_0[i] * pa_z[i] - ta1_x_zz_xzzz_1[i] * pc_z[i];

        ta1_x_zzz_yyyy_0[i] = 2.0 * ta1_x_z_yyyy_0[i] * fe_0 - 2.0 * ta1_x_z_yyyy_1[i] * fe_0 +
                              ta1_x_zz_yyyy_0[i] * pa_z[i] - ta1_x_zz_yyyy_1[i] * pc_z[i];

        ta1_x_zzz_yyyz_0[i] = 2.0 * ta1_x_z_yyyz_0[i] * fe_0 - 2.0 * ta1_x_z_yyyz_1[i] * fe_0 +
                              ta1_x_zz_yyy_0[i] * fe_0 - ta1_x_zz_yyy_1[i] * fe_0 + ta1_x_zz_yyyz_0[i] * pa_z[i] -
                              ta1_x_zz_yyyz_1[i] * pc_z[i];

        ta1_x_zzz_yyzz_0[i] = 2.0 * ta1_x_z_yyzz_0[i] * fe_0 - 2.0 * ta1_x_z_yyzz_1[i] * fe_0 +
                              2.0 * ta1_x_zz_yyz_0[i] * fe_0 - 2.0 * ta1_x_zz_yyz_1[i] * fe_0 +
                              ta1_x_zz_yyzz_0[i] * pa_z[i] - ta1_x_zz_yyzz_1[i] * pc_z[i];

        ta1_x_zzz_yzzz_0[i] = 2.0 * ta1_x_z_yzzz_0[i] * fe_0 - 2.0 * ta1_x_z_yzzz_1[i] * fe_0 +
                              3.0 * ta1_x_zz_yzz_0[i] * fe_0 - 3.0 * ta1_x_zz_yzz_1[i] * fe_0 +
                              ta1_x_zz_yzzz_0[i] * pa_z[i] - ta1_x_zz_yzzz_1[i] * pc_z[i];

        ta1_x_zzz_zzzz_0[i] = 2.0 * ta1_x_z_zzzz_0[i] * fe_0 - 2.0 * ta1_x_z_zzzz_1[i] * fe_0 +
                              4.0 * ta1_x_zz_zzz_0[i] * fe_0 - 4.0 * ta1_x_zz_zzz_1[i] * fe_0 +
                              ta1_x_zz_zzzz_0[i] * pa_z[i] - ta1_x_zz_zzzz_1[i] * pc_z[i];
    }

    // Set up 150-165 components of targeted buffer : FG

    auto ta1_y_xxx_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 150);

    auto ta1_y_xxx_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 151);

    auto ta1_y_xxx_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 152);

    auto ta1_y_xxx_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 153);

    auto ta1_y_xxx_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 154);

    auto ta1_y_xxx_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 155);

    auto ta1_y_xxx_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 156);

    auto ta1_y_xxx_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 157);

    auto ta1_y_xxx_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 158);

    auto ta1_y_xxx_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 159);

    auto ta1_y_xxx_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 160);

    auto ta1_y_xxx_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 161);

    auto ta1_y_xxx_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 162);

    auto ta1_y_xxx_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 163);

    auto ta1_y_xxx_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 164);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_y_x_xxxx_0,   \
                             ta1_y_x_xxxx_1,   \
                             ta1_y_x_xxxy_0,   \
                             ta1_y_x_xxxy_1,   \
                             ta1_y_x_xxxz_0,   \
                             ta1_y_x_xxxz_1,   \
                             ta1_y_x_xxyy_0,   \
                             ta1_y_x_xxyy_1,   \
                             ta1_y_x_xxyz_0,   \
                             ta1_y_x_xxyz_1,   \
                             ta1_y_x_xxzz_0,   \
                             ta1_y_x_xxzz_1,   \
                             ta1_y_x_xyyy_0,   \
                             ta1_y_x_xyyy_1,   \
                             ta1_y_x_xyyz_0,   \
                             ta1_y_x_xyyz_1,   \
                             ta1_y_x_xyzz_0,   \
                             ta1_y_x_xyzz_1,   \
                             ta1_y_x_xzzz_0,   \
                             ta1_y_x_xzzz_1,   \
                             ta1_y_x_yyyy_0,   \
                             ta1_y_x_yyyy_1,   \
                             ta1_y_x_yyyz_0,   \
                             ta1_y_x_yyyz_1,   \
                             ta1_y_x_yyzz_0,   \
                             ta1_y_x_yyzz_1,   \
                             ta1_y_x_yzzz_0,   \
                             ta1_y_x_yzzz_1,   \
                             ta1_y_x_zzzz_0,   \
                             ta1_y_x_zzzz_1,   \
                             ta1_y_xx_xxx_0,   \
                             ta1_y_xx_xxx_1,   \
                             ta1_y_xx_xxxx_0,  \
                             ta1_y_xx_xxxx_1,  \
                             ta1_y_xx_xxxy_0,  \
                             ta1_y_xx_xxxy_1,  \
                             ta1_y_xx_xxxz_0,  \
                             ta1_y_xx_xxxz_1,  \
                             ta1_y_xx_xxy_0,   \
                             ta1_y_xx_xxy_1,   \
                             ta1_y_xx_xxyy_0,  \
                             ta1_y_xx_xxyy_1,  \
                             ta1_y_xx_xxyz_0,  \
                             ta1_y_xx_xxyz_1,  \
                             ta1_y_xx_xxz_0,   \
                             ta1_y_xx_xxz_1,   \
                             ta1_y_xx_xxzz_0,  \
                             ta1_y_xx_xxzz_1,  \
                             ta1_y_xx_xyy_0,   \
                             ta1_y_xx_xyy_1,   \
                             ta1_y_xx_xyyy_0,  \
                             ta1_y_xx_xyyy_1,  \
                             ta1_y_xx_xyyz_0,  \
                             ta1_y_xx_xyyz_1,  \
                             ta1_y_xx_xyz_0,   \
                             ta1_y_xx_xyz_1,   \
                             ta1_y_xx_xyzz_0,  \
                             ta1_y_xx_xyzz_1,  \
                             ta1_y_xx_xzz_0,   \
                             ta1_y_xx_xzz_1,   \
                             ta1_y_xx_xzzz_0,  \
                             ta1_y_xx_xzzz_1,  \
                             ta1_y_xx_yyy_0,   \
                             ta1_y_xx_yyy_1,   \
                             ta1_y_xx_yyyy_0,  \
                             ta1_y_xx_yyyy_1,  \
                             ta1_y_xx_yyyz_0,  \
                             ta1_y_xx_yyyz_1,  \
                             ta1_y_xx_yyz_0,   \
                             ta1_y_xx_yyz_1,   \
                             ta1_y_xx_yyzz_0,  \
                             ta1_y_xx_yyzz_1,  \
                             ta1_y_xx_yzz_0,   \
                             ta1_y_xx_yzz_1,   \
                             ta1_y_xx_yzzz_0,  \
                             ta1_y_xx_yzzz_1,  \
                             ta1_y_xx_zzz_0,   \
                             ta1_y_xx_zzz_1,   \
                             ta1_y_xx_zzzz_0,  \
                             ta1_y_xx_zzzz_1,  \
                             ta1_y_xxx_xxxx_0, \
                             ta1_y_xxx_xxxy_0, \
                             ta1_y_xxx_xxxz_0, \
                             ta1_y_xxx_xxyy_0, \
                             ta1_y_xxx_xxyz_0, \
                             ta1_y_xxx_xxzz_0, \
                             ta1_y_xxx_xyyy_0, \
                             ta1_y_xxx_xyyz_0, \
                             ta1_y_xxx_xyzz_0, \
                             ta1_y_xxx_xzzz_0, \
                             ta1_y_xxx_yyyy_0, \
                             ta1_y_xxx_yyyz_0, \
                             ta1_y_xxx_yyzz_0, \
                             ta1_y_xxx_yzzz_0, \
                             ta1_y_xxx_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxx_xxxx_0[i] = 2.0 * ta1_y_x_xxxx_0[i] * fe_0 - 2.0 * ta1_y_x_xxxx_1[i] * fe_0 +
                              4.0 * ta1_y_xx_xxx_0[i] * fe_0 - 4.0 * ta1_y_xx_xxx_1[i] * fe_0 +
                              ta1_y_xx_xxxx_0[i] * pa_x[i] - ta1_y_xx_xxxx_1[i] * pc_x[i];

        ta1_y_xxx_xxxy_0[i] = 2.0 * ta1_y_x_xxxy_0[i] * fe_0 - 2.0 * ta1_y_x_xxxy_1[i] * fe_0 +
                              3.0 * ta1_y_xx_xxy_0[i] * fe_0 - 3.0 * ta1_y_xx_xxy_1[i] * fe_0 +
                              ta1_y_xx_xxxy_0[i] * pa_x[i] - ta1_y_xx_xxxy_1[i] * pc_x[i];

        ta1_y_xxx_xxxz_0[i] = 2.0 * ta1_y_x_xxxz_0[i] * fe_0 - 2.0 * ta1_y_x_xxxz_1[i] * fe_0 +
                              3.0 * ta1_y_xx_xxz_0[i] * fe_0 - 3.0 * ta1_y_xx_xxz_1[i] * fe_0 +
                              ta1_y_xx_xxxz_0[i] * pa_x[i] - ta1_y_xx_xxxz_1[i] * pc_x[i];

        ta1_y_xxx_xxyy_0[i] = 2.0 * ta1_y_x_xxyy_0[i] * fe_0 - 2.0 * ta1_y_x_xxyy_1[i] * fe_0 +
                              2.0 * ta1_y_xx_xyy_0[i] * fe_0 - 2.0 * ta1_y_xx_xyy_1[i] * fe_0 +
                              ta1_y_xx_xxyy_0[i] * pa_x[i] - ta1_y_xx_xxyy_1[i] * pc_x[i];

        ta1_y_xxx_xxyz_0[i] = 2.0 * ta1_y_x_xxyz_0[i] * fe_0 - 2.0 * ta1_y_x_xxyz_1[i] * fe_0 +
                              2.0 * ta1_y_xx_xyz_0[i] * fe_0 - 2.0 * ta1_y_xx_xyz_1[i] * fe_0 +
                              ta1_y_xx_xxyz_0[i] * pa_x[i] - ta1_y_xx_xxyz_1[i] * pc_x[i];

        ta1_y_xxx_xxzz_0[i] = 2.0 * ta1_y_x_xxzz_0[i] * fe_0 - 2.0 * ta1_y_x_xxzz_1[i] * fe_0 +
                              2.0 * ta1_y_xx_xzz_0[i] * fe_0 - 2.0 * ta1_y_xx_xzz_1[i] * fe_0 +
                              ta1_y_xx_xxzz_0[i] * pa_x[i] - ta1_y_xx_xxzz_1[i] * pc_x[i];

        ta1_y_xxx_xyyy_0[i] = 2.0 * ta1_y_x_xyyy_0[i] * fe_0 - 2.0 * ta1_y_x_xyyy_1[i] * fe_0 +
                              ta1_y_xx_yyy_0[i] * fe_0 - ta1_y_xx_yyy_1[i] * fe_0 + ta1_y_xx_xyyy_0[i] * pa_x[i] -
                              ta1_y_xx_xyyy_1[i] * pc_x[i];

        ta1_y_xxx_xyyz_0[i] = 2.0 * ta1_y_x_xyyz_0[i] * fe_0 - 2.0 * ta1_y_x_xyyz_1[i] * fe_0 +
                              ta1_y_xx_yyz_0[i] * fe_0 - ta1_y_xx_yyz_1[i] * fe_0 + ta1_y_xx_xyyz_0[i] * pa_x[i] -
                              ta1_y_xx_xyyz_1[i] * pc_x[i];

        ta1_y_xxx_xyzz_0[i] = 2.0 * ta1_y_x_xyzz_0[i] * fe_0 - 2.0 * ta1_y_x_xyzz_1[i] * fe_0 +
                              ta1_y_xx_yzz_0[i] * fe_0 - ta1_y_xx_yzz_1[i] * fe_0 + ta1_y_xx_xyzz_0[i] * pa_x[i] -
                              ta1_y_xx_xyzz_1[i] * pc_x[i];

        ta1_y_xxx_xzzz_0[i] = 2.0 * ta1_y_x_xzzz_0[i] * fe_0 - 2.0 * ta1_y_x_xzzz_1[i] * fe_0 +
                              ta1_y_xx_zzz_0[i] * fe_0 - ta1_y_xx_zzz_1[i] * fe_0 + ta1_y_xx_xzzz_0[i] * pa_x[i] -
                              ta1_y_xx_xzzz_1[i] * pc_x[i];

        ta1_y_xxx_yyyy_0[i] = 2.0 * ta1_y_x_yyyy_0[i] * fe_0 - 2.0 * ta1_y_x_yyyy_1[i] * fe_0 +
                              ta1_y_xx_yyyy_0[i] * pa_x[i] - ta1_y_xx_yyyy_1[i] * pc_x[i];

        ta1_y_xxx_yyyz_0[i] = 2.0 * ta1_y_x_yyyz_0[i] * fe_0 - 2.0 * ta1_y_x_yyyz_1[i] * fe_0 +
                              ta1_y_xx_yyyz_0[i] * pa_x[i] - ta1_y_xx_yyyz_1[i] * pc_x[i];

        ta1_y_xxx_yyzz_0[i] = 2.0 * ta1_y_x_yyzz_0[i] * fe_0 - 2.0 * ta1_y_x_yyzz_1[i] * fe_0 +
                              ta1_y_xx_yyzz_0[i] * pa_x[i] - ta1_y_xx_yyzz_1[i] * pc_x[i];

        ta1_y_xxx_yzzz_0[i] = 2.0 * ta1_y_x_yzzz_0[i] * fe_0 - 2.0 * ta1_y_x_yzzz_1[i] * fe_0 +
                              ta1_y_xx_yzzz_0[i] * pa_x[i] - ta1_y_xx_yzzz_1[i] * pc_x[i];

        ta1_y_xxx_zzzz_0[i] = 2.0 * ta1_y_x_zzzz_0[i] * fe_0 - 2.0 * ta1_y_x_zzzz_1[i] * fe_0 +
                              ta1_y_xx_zzzz_0[i] * pa_x[i] - ta1_y_xx_zzzz_1[i] * pc_x[i];
    }

    // Set up 165-180 components of targeted buffer : FG

    auto ta1_y_xxy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 165);

    auto ta1_y_xxy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 166);

    auto ta1_y_xxy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 167);

    auto ta1_y_xxy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 168);

    auto ta1_y_xxy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 169);

    auto ta1_y_xxy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 170);

    auto ta1_y_xxy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 171);

    auto ta1_y_xxy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 172);

    auto ta1_y_xxy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 173);

    auto ta1_y_xxy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 174);

    auto ta1_y_xxy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 175);

    auto ta1_y_xxy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 176);

    auto ta1_y_xxy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 177);

    auto ta1_y_xxy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 178);

    auto ta1_y_xxy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 179);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_y_xx_xxx_0,   \
                             ta1_y_xx_xxx_1,   \
                             ta1_y_xx_xxxx_0,  \
                             ta1_y_xx_xxxx_1,  \
                             ta1_y_xx_xxxy_0,  \
                             ta1_y_xx_xxxy_1,  \
                             ta1_y_xx_xxxz_0,  \
                             ta1_y_xx_xxxz_1,  \
                             ta1_y_xx_xxy_0,   \
                             ta1_y_xx_xxy_1,   \
                             ta1_y_xx_xxyy_0,  \
                             ta1_y_xx_xxyy_1,  \
                             ta1_y_xx_xxyz_0,  \
                             ta1_y_xx_xxyz_1,  \
                             ta1_y_xx_xxz_0,   \
                             ta1_y_xx_xxz_1,   \
                             ta1_y_xx_xxzz_0,  \
                             ta1_y_xx_xxzz_1,  \
                             ta1_y_xx_xyy_0,   \
                             ta1_y_xx_xyy_1,   \
                             ta1_y_xx_xyyy_0,  \
                             ta1_y_xx_xyyy_1,  \
                             ta1_y_xx_xyyz_0,  \
                             ta1_y_xx_xyyz_1,  \
                             ta1_y_xx_xyz_0,   \
                             ta1_y_xx_xyz_1,   \
                             ta1_y_xx_xyzz_0,  \
                             ta1_y_xx_xyzz_1,  \
                             ta1_y_xx_xzz_0,   \
                             ta1_y_xx_xzz_1,   \
                             ta1_y_xx_xzzz_0,  \
                             ta1_y_xx_xzzz_1,  \
                             ta1_y_xx_zzzz_0,  \
                             ta1_y_xx_zzzz_1,  \
                             ta1_y_xxy_xxxx_0, \
                             ta1_y_xxy_xxxy_0, \
                             ta1_y_xxy_xxxz_0, \
                             ta1_y_xxy_xxyy_0, \
                             ta1_y_xxy_xxyz_0, \
                             ta1_y_xxy_xxzz_0, \
                             ta1_y_xxy_xyyy_0, \
                             ta1_y_xxy_xyyz_0, \
                             ta1_y_xxy_xyzz_0, \
                             ta1_y_xxy_xzzz_0, \
                             ta1_y_xxy_yyyy_0, \
                             ta1_y_xxy_yyyz_0, \
                             ta1_y_xxy_yyzz_0, \
                             ta1_y_xxy_yzzz_0, \
                             ta1_y_xxy_zzzz_0, \
                             ta1_y_xy_yyyy_0,  \
                             ta1_y_xy_yyyy_1,  \
                             ta1_y_xy_yyyz_0,  \
                             ta1_y_xy_yyyz_1,  \
                             ta1_y_xy_yyzz_0,  \
                             ta1_y_xy_yyzz_1,  \
                             ta1_y_xy_yzzz_0,  \
                             ta1_y_xy_yzzz_1,  \
                             ta1_y_y_yyyy_0,   \
                             ta1_y_y_yyyy_1,   \
                             ta1_y_y_yyyz_0,   \
                             ta1_y_y_yyyz_1,   \
                             ta1_y_y_yyzz_0,   \
                             ta1_y_y_yyzz_1,   \
                             ta1_y_y_yzzz_0,   \
                             ta1_y_y_yzzz_1,   \
                             ta_xx_xxxx_1,     \
                             ta_xx_xxxy_1,     \
                             ta_xx_xxxz_1,     \
                             ta_xx_xxyy_1,     \
                             ta_xx_xxyz_1,     \
                             ta_xx_xxzz_1,     \
                             ta_xx_xyyy_1,     \
                             ta_xx_xyyz_1,     \
                             ta_xx_xyzz_1,     \
                             ta_xx_xzzz_1,     \
                             ta_xx_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxy_xxxx_0[i] = ta_xx_xxxx_1[i] + ta1_y_xx_xxxx_0[i] * pa_y[i] - ta1_y_xx_xxxx_1[i] * pc_y[i];

        ta1_y_xxy_xxxy_0[i] = ta1_y_xx_xxx_0[i] * fe_0 - ta1_y_xx_xxx_1[i] * fe_0 + ta_xx_xxxy_1[i] +
                              ta1_y_xx_xxxy_0[i] * pa_y[i] - ta1_y_xx_xxxy_1[i] * pc_y[i];

        ta1_y_xxy_xxxz_0[i] = ta_xx_xxxz_1[i] + ta1_y_xx_xxxz_0[i] * pa_y[i] - ta1_y_xx_xxxz_1[i] * pc_y[i];

        ta1_y_xxy_xxyy_0[i] = 2.0 * ta1_y_xx_xxy_0[i] * fe_0 - 2.0 * ta1_y_xx_xxy_1[i] * fe_0 + ta_xx_xxyy_1[i] +
                              ta1_y_xx_xxyy_0[i] * pa_y[i] - ta1_y_xx_xxyy_1[i] * pc_y[i];

        ta1_y_xxy_xxyz_0[i] = ta1_y_xx_xxz_0[i] * fe_0 - ta1_y_xx_xxz_1[i] * fe_0 + ta_xx_xxyz_1[i] +
                              ta1_y_xx_xxyz_0[i] * pa_y[i] - ta1_y_xx_xxyz_1[i] * pc_y[i];

        ta1_y_xxy_xxzz_0[i] = ta_xx_xxzz_1[i] + ta1_y_xx_xxzz_0[i] * pa_y[i] - ta1_y_xx_xxzz_1[i] * pc_y[i];

        ta1_y_xxy_xyyy_0[i] = 3.0 * ta1_y_xx_xyy_0[i] * fe_0 - 3.0 * ta1_y_xx_xyy_1[i] * fe_0 + ta_xx_xyyy_1[i] +
                              ta1_y_xx_xyyy_0[i] * pa_y[i] - ta1_y_xx_xyyy_1[i] * pc_y[i];

        ta1_y_xxy_xyyz_0[i] = 2.0 * ta1_y_xx_xyz_0[i] * fe_0 - 2.0 * ta1_y_xx_xyz_1[i] * fe_0 + ta_xx_xyyz_1[i] +
                              ta1_y_xx_xyyz_0[i] * pa_y[i] - ta1_y_xx_xyyz_1[i] * pc_y[i];

        ta1_y_xxy_xyzz_0[i] = ta1_y_xx_xzz_0[i] * fe_0 - ta1_y_xx_xzz_1[i] * fe_0 + ta_xx_xyzz_1[i] +
                              ta1_y_xx_xyzz_0[i] * pa_y[i] - ta1_y_xx_xyzz_1[i] * pc_y[i];

        ta1_y_xxy_xzzz_0[i] = ta_xx_xzzz_1[i] + ta1_y_xx_xzzz_0[i] * pa_y[i] - ta1_y_xx_xzzz_1[i] * pc_y[i];

        ta1_y_xxy_yyyy_0[i] = ta1_y_y_yyyy_0[i] * fe_0 - ta1_y_y_yyyy_1[i] * fe_0 + ta1_y_xy_yyyy_0[i] * pa_x[i] -
                              ta1_y_xy_yyyy_1[i] * pc_x[i];

        ta1_y_xxy_yyyz_0[i] = ta1_y_y_yyyz_0[i] * fe_0 - ta1_y_y_yyyz_1[i] * fe_0 + ta1_y_xy_yyyz_0[i] * pa_x[i] -
                              ta1_y_xy_yyyz_1[i] * pc_x[i];

        ta1_y_xxy_yyzz_0[i] = ta1_y_y_yyzz_0[i] * fe_0 - ta1_y_y_yyzz_1[i] * fe_0 + ta1_y_xy_yyzz_0[i] * pa_x[i] -
                              ta1_y_xy_yyzz_1[i] * pc_x[i];

        ta1_y_xxy_yzzz_0[i] = ta1_y_y_yzzz_0[i] * fe_0 - ta1_y_y_yzzz_1[i] * fe_0 + ta1_y_xy_yzzz_0[i] * pa_x[i] -
                              ta1_y_xy_yzzz_1[i] * pc_x[i];

        ta1_y_xxy_zzzz_0[i] = ta_xx_zzzz_1[i] + ta1_y_xx_zzzz_0[i] * pa_y[i] - ta1_y_xx_zzzz_1[i] * pc_y[i];
    }

    // Set up 180-195 components of targeted buffer : FG

    auto ta1_y_xxz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 180);

    auto ta1_y_xxz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 181);

    auto ta1_y_xxz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 182);

    auto ta1_y_xxz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 183);

    auto ta1_y_xxz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 184);

    auto ta1_y_xxz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 185);

    auto ta1_y_xxz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 186);

    auto ta1_y_xxz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 187);

    auto ta1_y_xxz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 188);

    auto ta1_y_xxz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 189);

    auto ta1_y_xxz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 190);

    auto ta1_y_xxz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 191);

    auto ta1_y_xxz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 192);

    auto ta1_y_xxz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 193);

    auto ta1_y_xxz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 194);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_y_xx_xxx_0,   \
                             ta1_y_xx_xxx_1,   \
                             ta1_y_xx_xxxx_0,  \
                             ta1_y_xx_xxxx_1,  \
                             ta1_y_xx_xxxy_0,  \
                             ta1_y_xx_xxxy_1,  \
                             ta1_y_xx_xxxz_0,  \
                             ta1_y_xx_xxxz_1,  \
                             ta1_y_xx_xxy_0,   \
                             ta1_y_xx_xxy_1,   \
                             ta1_y_xx_xxyy_0,  \
                             ta1_y_xx_xxyy_1,  \
                             ta1_y_xx_xxyz_0,  \
                             ta1_y_xx_xxyz_1,  \
                             ta1_y_xx_xxz_0,   \
                             ta1_y_xx_xxz_1,   \
                             ta1_y_xx_xxzz_0,  \
                             ta1_y_xx_xxzz_1,  \
                             ta1_y_xx_xyy_0,   \
                             ta1_y_xx_xyy_1,   \
                             ta1_y_xx_xyyy_0,  \
                             ta1_y_xx_xyyy_1,  \
                             ta1_y_xx_xyyz_0,  \
                             ta1_y_xx_xyyz_1,  \
                             ta1_y_xx_xyz_0,   \
                             ta1_y_xx_xyz_1,   \
                             ta1_y_xx_xyzz_0,  \
                             ta1_y_xx_xyzz_1,  \
                             ta1_y_xx_xzz_0,   \
                             ta1_y_xx_xzz_1,   \
                             ta1_y_xx_xzzz_0,  \
                             ta1_y_xx_xzzz_1,  \
                             ta1_y_xx_yyyy_0,  \
                             ta1_y_xx_yyyy_1,  \
                             ta1_y_xxz_xxxx_0, \
                             ta1_y_xxz_xxxy_0, \
                             ta1_y_xxz_xxxz_0, \
                             ta1_y_xxz_xxyy_0, \
                             ta1_y_xxz_xxyz_0, \
                             ta1_y_xxz_xxzz_0, \
                             ta1_y_xxz_xyyy_0, \
                             ta1_y_xxz_xyyz_0, \
                             ta1_y_xxz_xyzz_0, \
                             ta1_y_xxz_xzzz_0, \
                             ta1_y_xxz_yyyy_0, \
                             ta1_y_xxz_yyyz_0, \
                             ta1_y_xxz_yyzz_0, \
                             ta1_y_xxz_yzzz_0, \
                             ta1_y_xxz_zzzz_0, \
                             ta1_y_xz_yyyz_0,  \
                             ta1_y_xz_yyyz_1,  \
                             ta1_y_xz_yyzz_0,  \
                             ta1_y_xz_yyzz_1,  \
                             ta1_y_xz_yzzz_0,  \
                             ta1_y_xz_yzzz_1,  \
                             ta1_y_xz_zzzz_0,  \
                             ta1_y_xz_zzzz_1,  \
                             ta1_y_z_yyyz_0,   \
                             ta1_y_z_yyyz_1,   \
                             ta1_y_z_yyzz_0,   \
                             ta1_y_z_yyzz_1,   \
                             ta1_y_z_yzzz_0,   \
                             ta1_y_z_yzzz_1,   \
                             ta1_y_z_zzzz_0,   \
                             ta1_y_z_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxz_xxxx_0[i] = ta1_y_xx_xxxx_0[i] * pa_z[i] - ta1_y_xx_xxxx_1[i] * pc_z[i];

        ta1_y_xxz_xxxy_0[i] = ta1_y_xx_xxxy_0[i] * pa_z[i] - ta1_y_xx_xxxy_1[i] * pc_z[i];

        ta1_y_xxz_xxxz_0[i] = ta1_y_xx_xxx_0[i] * fe_0 - ta1_y_xx_xxx_1[i] * fe_0 + ta1_y_xx_xxxz_0[i] * pa_z[i] -
                              ta1_y_xx_xxxz_1[i] * pc_z[i];

        ta1_y_xxz_xxyy_0[i] = ta1_y_xx_xxyy_0[i] * pa_z[i] - ta1_y_xx_xxyy_1[i] * pc_z[i];

        ta1_y_xxz_xxyz_0[i] = ta1_y_xx_xxy_0[i] * fe_0 - ta1_y_xx_xxy_1[i] * fe_0 + ta1_y_xx_xxyz_0[i] * pa_z[i] -
                              ta1_y_xx_xxyz_1[i] * pc_z[i];

        ta1_y_xxz_xxzz_0[i] = 2.0 * ta1_y_xx_xxz_0[i] * fe_0 - 2.0 * ta1_y_xx_xxz_1[i] * fe_0 +
                              ta1_y_xx_xxzz_0[i] * pa_z[i] - ta1_y_xx_xxzz_1[i] * pc_z[i];

        ta1_y_xxz_xyyy_0[i] = ta1_y_xx_xyyy_0[i] * pa_z[i] - ta1_y_xx_xyyy_1[i] * pc_z[i];

        ta1_y_xxz_xyyz_0[i] = ta1_y_xx_xyy_0[i] * fe_0 - ta1_y_xx_xyy_1[i] * fe_0 + ta1_y_xx_xyyz_0[i] * pa_z[i] -
                              ta1_y_xx_xyyz_1[i] * pc_z[i];

        ta1_y_xxz_xyzz_0[i] = 2.0 * ta1_y_xx_xyz_0[i] * fe_0 - 2.0 * ta1_y_xx_xyz_1[i] * fe_0 +
                              ta1_y_xx_xyzz_0[i] * pa_z[i] - ta1_y_xx_xyzz_1[i] * pc_z[i];

        ta1_y_xxz_xzzz_0[i] = 3.0 * ta1_y_xx_xzz_0[i] * fe_0 - 3.0 * ta1_y_xx_xzz_1[i] * fe_0 +
                              ta1_y_xx_xzzz_0[i] * pa_z[i] - ta1_y_xx_xzzz_1[i] * pc_z[i];

        ta1_y_xxz_yyyy_0[i] = ta1_y_xx_yyyy_0[i] * pa_z[i] - ta1_y_xx_yyyy_1[i] * pc_z[i];

        ta1_y_xxz_yyyz_0[i] = ta1_y_z_yyyz_0[i] * fe_0 - ta1_y_z_yyyz_1[i] * fe_0 + ta1_y_xz_yyyz_0[i] * pa_x[i] -
                              ta1_y_xz_yyyz_1[i] * pc_x[i];

        ta1_y_xxz_yyzz_0[i] = ta1_y_z_yyzz_0[i] * fe_0 - ta1_y_z_yyzz_1[i] * fe_0 + ta1_y_xz_yyzz_0[i] * pa_x[i] -
                              ta1_y_xz_yyzz_1[i] * pc_x[i];

        ta1_y_xxz_yzzz_0[i] = ta1_y_z_yzzz_0[i] * fe_0 - ta1_y_z_yzzz_1[i] * fe_0 + ta1_y_xz_yzzz_0[i] * pa_x[i] -
                              ta1_y_xz_yzzz_1[i] * pc_x[i];

        ta1_y_xxz_zzzz_0[i] = ta1_y_z_zzzz_0[i] * fe_0 - ta1_y_z_zzzz_1[i] * fe_0 + ta1_y_xz_zzzz_0[i] * pa_x[i] -
                              ta1_y_xz_zzzz_1[i] * pc_x[i];
    }

    // Set up 195-210 components of targeted buffer : FG

    auto ta1_y_xyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 195);

    auto ta1_y_xyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 196);

    auto ta1_y_xyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 197);

    auto ta1_y_xyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 198);

    auto ta1_y_xyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 199);

    auto ta1_y_xyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 200);

    auto ta1_y_xyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 201);

    auto ta1_y_xyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 202);

    auto ta1_y_xyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 203);

    auto ta1_y_xyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 204);

    auto ta1_y_xyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 205);

    auto ta1_y_xyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 206);

    auto ta1_y_xyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 207);

    auto ta1_y_xyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 208);

    auto ta1_y_xyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 209);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_y_xyy_xxxx_0, \
                             ta1_y_xyy_xxxy_0, \
                             ta1_y_xyy_xxxz_0, \
                             ta1_y_xyy_xxyy_0, \
                             ta1_y_xyy_xxyz_0, \
                             ta1_y_xyy_xxzz_0, \
                             ta1_y_xyy_xyyy_0, \
                             ta1_y_xyy_xyyz_0, \
                             ta1_y_xyy_xyzz_0, \
                             ta1_y_xyy_xzzz_0, \
                             ta1_y_xyy_yyyy_0, \
                             ta1_y_xyy_yyyz_0, \
                             ta1_y_xyy_yyzz_0, \
                             ta1_y_xyy_yzzz_0, \
                             ta1_y_xyy_zzzz_0, \
                             ta1_y_yy_xxx_0,   \
                             ta1_y_yy_xxx_1,   \
                             ta1_y_yy_xxxx_0,  \
                             ta1_y_yy_xxxx_1,  \
                             ta1_y_yy_xxxy_0,  \
                             ta1_y_yy_xxxy_1,  \
                             ta1_y_yy_xxxz_0,  \
                             ta1_y_yy_xxxz_1,  \
                             ta1_y_yy_xxy_0,   \
                             ta1_y_yy_xxy_1,   \
                             ta1_y_yy_xxyy_0,  \
                             ta1_y_yy_xxyy_1,  \
                             ta1_y_yy_xxyz_0,  \
                             ta1_y_yy_xxyz_1,  \
                             ta1_y_yy_xxz_0,   \
                             ta1_y_yy_xxz_1,   \
                             ta1_y_yy_xxzz_0,  \
                             ta1_y_yy_xxzz_1,  \
                             ta1_y_yy_xyy_0,   \
                             ta1_y_yy_xyy_1,   \
                             ta1_y_yy_xyyy_0,  \
                             ta1_y_yy_xyyy_1,  \
                             ta1_y_yy_xyyz_0,  \
                             ta1_y_yy_xyyz_1,  \
                             ta1_y_yy_xyz_0,   \
                             ta1_y_yy_xyz_1,   \
                             ta1_y_yy_xyzz_0,  \
                             ta1_y_yy_xyzz_1,  \
                             ta1_y_yy_xzz_0,   \
                             ta1_y_yy_xzz_1,   \
                             ta1_y_yy_xzzz_0,  \
                             ta1_y_yy_xzzz_1,  \
                             ta1_y_yy_yyy_0,   \
                             ta1_y_yy_yyy_1,   \
                             ta1_y_yy_yyyy_0,  \
                             ta1_y_yy_yyyy_1,  \
                             ta1_y_yy_yyyz_0,  \
                             ta1_y_yy_yyyz_1,  \
                             ta1_y_yy_yyz_0,   \
                             ta1_y_yy_yyz_1,   \
                             ta1_y_yy_yyzz_0,  \
                             ta1_y_yy_yyzz_1,  \
                             ta1_y_yy_yzz_0,   \
                             ta1_y_yy_yzz_1,   \
                             ta1_y_yy_yzzz_0,  \
                             ta1_y_yy_yzzz_1,  \
                             ta1_y_yy_zzz_0,   \
                             ta1_y_yy_zzz_1,   \
                             ta1_y_yy_zzzz_0,  \
                             ta1_y_yy_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyy_xxxx_0[i] = 4.0 * ta1_y_yy_xxx_0[i] * fe_0 - 4.0 * ta1_y_yy_xxx_1[i] * fe_0 +
                              ta1_y_yy_xxxx_0[i] * pa_x[i] - ta1_y_yy_xxxx_1[i] * pc_x[i];

        ta1_y_xyy_xxxy_0[i] = 3.0 * ta1_y_yy_xxy_0[i] * fe_0 - 3.0 * ta1_y_yy_xxy_1[i] * fe_0 +
                              ta1_y_yy_xxxy_0[i] * pa_x[i] - ta1_y_yy_xxxy_1[i] * pc_x[i];

        ta1_y_xyy_xxxz_0[i] = 3.0 * ta1_y_yy_xxz_0[i] * fe_0 - 3.0 * ta1_y_yy_xxz_1[i] * fe_0 +
                              ta1_y_yy_xxxz_0[i] * pa_x[i] - ta1_y_yy_xxxz_1[i] * pc_x[i];

        ta1_y_xyy_xxyy_0[i] = 2.0 * ta1_y_yy_xyy_0[i] * fe_0 - 2.0 * ta1_y_yy_xyy_1[i] * fe_0 +
                              ta1_y_yy_xxyy_0[i] * pa_x[i] - ta1_y_yy_xxyy_1[i] * pc_x[i];

        ta1_y_xyy_xxyz_0[i] = 2.0 * ta1_y_yy_xyz_0[i] * fe_0 - 2.0 * ta1_y_yy_xyz_1[i] * fe_0 +
                              ta1_y_yy_xxyz_0[i] * pa_x[i] - ta1_y_yy_xxyz_1[i] * pc_x[i];

        ta1_y_xyy_xxzz_0[i] = 2.0 * ta1_y_yy_xzz_0[i] * fe_0 - 2.0 * ta1_y_yy_xzz_1[i] * fe_0 +
                              ta1_y_yy_xxzz_0[i] * pa_x[i] - ta1_y_yy_xxzz_1[i] * pc_x[i];

        ta1_y_xyy_xyyy_0[i] = ta1_y_yy_yyy_0[i] * fe_0 - ta1_y_yy_yyy_1[i] * fe_0 + ta1_y_yy_xyyy_0[i] * pa_x[i] -
                              ta1_y_yy_xyyy_1[i] * pc_x[i];

        ta1_y_xyy_xyyz_0[i] = ta1_y_yy_yyz_0[i] * fe_0 - ta1_y_yy_yyz_1[i] * fe_0 + ta1_y_yy_xyyz_0[i] * pa_x[i] -
                              ta1_y_yy_xyyz_1[i] * pc_x[i];

        ta1_y_xyy_xyzz_0[i] = ta1_y_yy_yzz_0[i] * fe_0 - ta1_y_yy_yzz_1[i] * fe_0 + ta1_y_yy_xyzz_0[i] * pa_x[i] -
                              ta1_y_yy_xyzz_1[i] * pc_x[i];

        ta1_y_xyy_xzzz_0[i] = ta1_y_yy_zzz_0[i] * fe_0 - ta1_y_yy_zzz_1[i] * fe_0 + ta1_y_yy_xzzz_0[i] * pa_x[i] -
                              ta1_y_yy_xzzz_1[i] * pc_x[i];

        ta1_y_xyy_yyyy_0[i] = ta1_y_yy_yyyy_0[i] * pa_x[i] - ta1_y_yy_yyyy_1[i] * pc_x[i];

        ta1_y_xyy_yyyz_0[i] = ta1_y_yy_yyyz_0[i] * pa_x[i] - ta1_y_yy_yyyz_1[i] * pc_x[i];

        ta1_y_xyy_yyzz_0[i] = ta1_y_yy_yyzz_0[i] * pa_x[i] - ta1_y_yy_yyzz_1[i] * pc_x[i];

        ta1_y_xyy_yzzz_0[i] = ta1_y_yy_yzzz_0[i] * pa_x[i] - ta1_y_yy_yzzz_1[i] * pc_x[i];

        ta1_y_xyy_zzzz_0[i] = ta1_y_yy_zzzz_0[i] * pa_x[i] - ta1_y_yy_zzzz_1[i] * pc_x[i];
    }

    // Set up 210-225 components of targeted buffer : FG

    auto ta1_y_xyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 210);

    auto ta1_y_xyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 211);

    auto ta1_y_xyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 212);

    auto ta1_y_xyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 213);

    auto ta1_y_xyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 214);

    auto ta1_y_xyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 215);

    auto ta1_y_xyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 216);

    auto ta1_y_xyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 217);

    auto ta1_y_xyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 218);

    auto ta1_y_xyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 219);

    auto ta1_y_xyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 220);

    auto ta1_y_xyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 221);

    auto ta1_y_xyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 222);

    auto ta1_y_xyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 223);

    auto ta1_y_xyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 224);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             pc_x,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_y_xy_xxxx_0,  \
                             ta1_y_xy_xxxx_1,  \
                             ta1_y_xy_xxxy_0,  \
                             ta1_y_xy_xxxy_1,  \
                             ta1_y_xy_xxyy_0,  \
                             ta1_y_xy_xxyy_1,  \
                             ta1_y_xy_xyyy_0,  \
                             ta1_y_xy_xyyy_1,  \
                             ta1_y_xyz_xxxx_0, \
                             ta1_y_xyz_xxxy_0, \
                             ta1_y_xyz_xxxz_0, \
                             ta1_y_xyz_xxyy_0, \
                             ta1_y_xyz_xxyz_0, \
                             ta1_y_xyz_xxzz_0, \
                             ta1_y_xyz_xyyy_0, \
                             ta1_y_xyz_xyyz_0, \
                             ta1_y_xyz_xyzz_0, \
                             ta1_y_xyz_xzzz_0, \
                             ta1_y_xyz_yyyy_0, \
                             ta1_y_xyz_yyyz_0, \
                             ta1_y_xyz_yyzz_0, \
                             ta1_y_xyz_yzzz_0, \
                             ta1_y_xyz_zzzz_0, \
                             ta1_y_xz_xxxz_0,  \
                             ta1_y_xz_xxxz_1,  \
                             ta1_y_xz_xxzz_0,  \
                             ta1_y_xz_xxzz_1,  \
                             ta1_y_xz_xzzz_0,  \
                             ta1_y_xz_xzzz_1,  \
                             ta1_y_yz_xxyz_0,  \
                             ta1_y_yz_xxyz_1,  \
                             ta1_y_yz_xyyz_0,  \
                             ta1_y_yz_xyyz_1,  \
                             ta1_y_yz_xyz_0,   \
                             ta1_y_yz_xyz_1,   \
                             ta1_y_yz_xyzz_0,  \
                             ta1_y_yz_xyzz_1,  \
                             ta1_y_yz_yyyy_0,  \
                             ta1_y_yz_yyyy_1,  \
                             ta1_y_yz_yyyz_0,  \
                             ta1_y_yz_yyyz_1,  \
                             ta1_y_yz_yyz_0,   \
                             ta1_y_yz_yyz_1,   \
                             ta1_y_yz_yyzz_0,  \
                             ta1_y_yz_yyzz_1,  \
                             ta1_y_yz_yzz_0,   \
                             ta1_y_yz_yzz_1,   \
                             ta1_y_yz_yzzz_0,  \
                             ta1_y_yz_yzzz_1,  \
                             ta1_y_yz_zzzz_0,  \
                             ta1_y_yz_zzzz_1,  \
                             ta_xz_xxxz_1,     \
                             ta_xz_xxzz_1,     \
                             ta_xz_xzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyz_xxxx_0[i] = ta1_y_xy_xxxx_0[i] * pa_z[i] - ta1_y_xy_xxxx_1[i] * pc_z[i];

        ta1_y_xyz_xxxy_0[i] = ta1_y_xy_xxxy_0[i] * pa_z[i] - ta1_y_xy_xxxy_1[i] * pc_z[i];

        ta1_y_xyz_xxxz_0[i] = ta_xz_xxxz_1[i] + ta1_y_xz_xxxz_0[i] * pa_y[i] - ta1_y_xz_xxxz_1[i] * pc_y[i];

        ta1_y_xyz_xxyy_0[i] = ta1_y_xy_xxyy_0[i] * pa_z[i] - ta1_y_xy_xxyy_1[i] * pc_z[i];

        ta1_y_xyz_xxyz_0[i] = 2.0 * ta1_y_yz_xyz_0[i] * fe_0 - 2.0 * ta1_y_yz_xyz_1[i] * fe_0 +
                              ta1_y_yz_xxyz_0[i] * pa_x[i] - ta1_y_yz_xxyz_1[i] * pc_x[i];

        ta1_y_xyz_xxzz_0[i] = ta_xz_xxzz_1[i] + ta1_y_xz_xxzz_0[i] * pa_y[i] - ta1_y_xz_xxzz_1[i] * pc_y[i];

        ta1_y_xyz_xyyy_0[i] = ta1_y_xy_xyyy_0[i] * pa_z[i] - ta1_y_xy_xyyy_1[i] * pc_z[i];

        ta1_y_xyz_xyyz_0[i] = ta1_y_yz_yyz_0[i] * fe_0 - ta1_y_yz_yyz_1[i] * fe_0 + ta1_y_yz_xyyz_0[i] * pa_x[i] -
                              ta1_y_yz_xyyz_1[i] * pc_x[i];

        ta1_y_xyz_xyzz_0[i] = ta1_y_yz_yzz_0[i] * fe_0 - ta1_y_yz_yzz_1[i] * fe_0 + ta1_y_yz_xyzz_0[i] * pa_x[i] -
                              ta1_y_yz_xyzz_1[i] * pc_x[i];

        ta1_y_xyz_xzzz_0[i] = ta_xz_xzzz_1[i] + ta1_y_xz_xzzz_0[i] * pa_y[i] - ta1_y_xz_xzzz_1[i] * pc_y[i];

        ta1_y_xyz_yyyy_0[i] = ta1_y_yz_yyyy_0[i] * pa_x[i] - ta1_y_yz_yyyy_1[i] * pc_x[i];

        ta1_y_xyz_yyyz_0[i] = ta1_y_yz_yyyz_0[i] * pa_x[i] - ta1_y_yz_yyyz_1[i] * pc_x[i];

        ta1_y_xyz_yyzz_0[i] = ta1_y_yz_yyzz_0[i] * pa_x[i] - ta1_y_yz_yyzz_1[i] * pc_x[i];

        ta1_y_xyz_yzzz_0[i] = ta1_y_yz_yzzz_0[i] * pa_x[i] - ta1_y_yz_yzzz_1[i] * pc_x[i];

        ta1_y_xyz_zzzz_0[i] = ta1_y_yz_zzzz_0[i] * pa_x[i] - ta1_y_yz_zzzz_1[i] * pc_x[i];
    }

    // Set up 225-240 components of targeted buffer : FG

    auto ta1_y_xzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 225);

    auto ta1_y_xzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 226);

    auto ta1_y_xzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 227);

    auto ta1_y_xzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 228);

    auto ta1_y_xzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 229);

    auto ta1_y_xzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 230);

    auto ta1_y_xzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 231);

    auto ta1_y_xzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 232);

    auto ta1_y_xzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 233);

    auto ta1_y_xzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 234);

    auto ta1_y_xzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 235);

    auto ta1_y_xzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 236);

    auto ta1_y_xzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 237);

    auto ta1_y_xzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 238);

    auto ta1_y_xzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 239);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_y_xzz_xxxx_0, \
                             ta1_y_xzz_xxxy_0, \
                             ta1_y_xzz_xxxz_0, \
                             ta1_y_xzz_xxyy_0, \
                             ta1_y_xzz_xxyz_0, \
                             ta1_y_xzz_xxzz_0, \
                             ta1_y_xzz_xyyy_0, \
                             ta1_y_xzz_xyyz_0, \
                             ta1_y_xzz_xyzz_0, \
                             ta1_y_xzz_xzzz_0, \
                             ta1_y_xzz_yyyy_0, \
                             ta1_y_xzz_yyyz_0, \
                             ta1_y_xzz_yyzz_0, \
                             ta1_y_xzz_yzzz_0, \
                             ta1_y_xzz_zzzz_0, \
                             ta1_y_zz_xxx_0,   \
                             ta1_y_zz_xxx_1,   \
                             ta1_y_zz_xxxx_0,  \
                             ta1_y_zz_xxxx_1,  \
                             ta1_y_zz_xxxy_0,  \
                             ta1_y_zz_xxxy_1,  \
                             ta1_y_zz_xxxz_0,  \
                             ta1_y_zz_xxxz_1,  \
                             ta1_y_zz_xxy_0,   \
                             ta1_y_zz_xxy_1,   \
                             ta1_y_zz_xxyy_0,  \
                             ta1_y_zz_xxyy_1,  \
                             ta1_y_zz_xxyz_0,  \
                             ta1_y_zz_xxyz_1,  \
                             ta1_y_zz_xxz_0,   \
                             ta1_y_zz_xxz_1,   \
                             ta1_y_zz_xxzz_0,  \
                             ta1_y_zz_xxzz_1,  \
                             ta1_y_zz_xyy_0,   \
                             ta1_y_zz_xyy_1,   \
                             ta1_y_zz_xyyy_0,  \
                             ta1_y_zz_xyyy_1,  \
                             ta1_y_zz_xyyz_0,  \
                             ta1_y_zz_xyyz_1,  \
                             ta1_y_zz_xyz_0,   \
                             ta1_y_zz_xyz_1,   \
                             ta1_y_zz_xyzz_0,  \
                             ta1_y_zz_xyzz_1,  \
                             ta1_y_zz_xzz_0,   \
                             ta1_y_zz_xzz_1,   \
                             ta1_y_zz_xzzz_0,  \
                             ta1_y_zz_xzzz_1,  \
                             ta1_y_zz_yyy_0,   \
                             ta1_y_zz_yyy_1,   \
                             ta1_y_zz_yyyy_0,  \
                             ta1_y_zz_yyyy_1,  \
                             ta1_y_zz_yyyz_0,  \
                             ta1_y_zz_yyyz_1,  \
                             ta1_y_zz_yyz_0,   \
                             ta1_y_zz_yyz_1,   \
                             ta1_y_zz_yyzz_0,  \
                             ta1_y_zz_yyzz_1,  \
                             ta1_y_zz_yzz_0,   \
                             ta1_y_zz_yzz_1,   \
                             ta1_y_zz_yzzz_0,  \
                             ta1_y_zz_yzzz_1,  \
                             ta1_y_zz_zzz_0,   \
                             ta1_y_zz_zzz_1,   \
                             ta1_y_zz_zzzz_0,  \
                             ta1_y_zz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xzz_xxxx_0[i] = 4.0 * ta1_y_zz_xxx_0[i] * fe_0 - 4.0 * ta1_y_zz_xxx_1[i] * fe_0 +
                              ta1_y_zz_xxxx_0[i] * pa_x[i] - ta1_y_zz_xxxx_1[i] * pc_x[i];

        ta1_y_xzz_xxxy_0[i] = 3.0 * ta1_y_zz_xxy_0[i] * fe_0 - 3.0 * ta1_y_zz_xxy_1[i] * fe_0 +
                              ta1_y_zz_xxxy_0[i] * pa_x[i] - ta1_y_zz_xxxy_1[i] * pc_x[i];

        ta1_y_xzz_xxxz_0[i] = 3.0 * ta1_y_zz_xxz_0[i] * fe_0 - 3.0 * ta1_y_zz_xxz_1[i] * fe_0 +
                              ta1_y_zz_xxxz_0[i] * pa_x[i] - ta1_y_zz_xxxz_1[i] * pc_x[i];

        ta1_y_xzz_xxyy_0[i] = 2.0 * ta1_y_zz_xyy_0[i] * fe_0 - 2.0 * ta1_y_zz_xyy_1[i] * fe_0 +
                              ta1_y_zz_xxyy_0[i] * pa_x[i] - ta1_y_zz_xxyy_1[i] * pc_x[i];

        ta1_y_xzz_xxyz_0[i] = 2.0 * ta1_y_zz_xyz_0[i] * fe_0 - 2.0 * ta1_y_zz_xyz_1[i] * fe_0 +
                              ta1_y_zz_xxyz_0[i] * pa_x[i] - ta1_y_zz_xxyz_1[i] * pc_x[i];

        ta1_y_xzz_xxzz_0[i] = 2.0 * ta1_y_zz_xzz_0[i] * fe_0 - 2.0 * ta1_y_zz_xzz_1[i] * fe_0 +
                              ta1_y_zz_xxzz_0[i] * pa_x[i] - ta1_y_zz_xxzz_1[i] * pc_x[i];

        ta1_y_xzz_xyyy_0[i] = ta1_y_zz_yyy_0[i] * fe_0 - ta1_y_zz_yyy_1[i] * fe_0 + ta1_y_zz_xyyy_0[i] * pa_x[i] -
                              ta1_y_zz_xyyy_1[i] * pc_x[i];

        ta1_y_xzz_xyyz_0[i] = ta1_y_zz_yyz_0[i] * fe_0 - ta1_y_zz_yyz_1[i] * fe_0 + ta1_y_zz_xyyz_0[i] * pa_x[i] -
                              ta1_y_zz_xyyz_1[i] * pc_x[i];

        ta1_y_xzz_xyzz_0[i] = ta1_y_zz_yzz_0[i] * fe_0 - ta1_y_zz_yzz_1[i] * fe_0 + ta1_y_zz_xyzz_0[i] * pa_x[i] -
                              ta1_y_zz_xyzz_1[i] * pc_x[i];

        ta1_y_xzz_xzzz_0[i] = ta1_y_zz_zzz_0[i] * fe_0 - ta1_y_zz_zzz_1[i] * fe_0 + ta1_y_zz_xzzz_0[i] * pa_x[i] -
                              ta1_y_zz_xzzz_1[i] * pc_x[i];

        ta1_y_xzz_yyyy_0[i] = ta1_y_zz_yyyy_0[i] * pa_x[i] - ta1_y_zz_yyyy_1[i] * pc_x[i];

        ta1_y_xzz_yyyz_0[i] = ta1_y_zz_yyyz_0[i] * pa_x[i] - ta1_y_zz_yyyz_1[i] * pc_x[i];

        ta1_y_xzz_yyzz_0[i] = ta1_y_zz_yyzz_0[i] * pa_x[i] - ta1_y_zz_yyzz_1[i] * pc_x[i];

        ta1_y_xzz_yzzz_0[i] = ta1_y_zz_yzzz_0[i] * pa_x[i] - ta1_y_zz_yzzz_1[i] * pc_x[i];

        ta1_y_xzz_zzzz_0[i] = ta1_y_zz_zzzz_0[i] * pa_x[i] - ta1_y_zz_zzzz_1[i] * pc_x[i];
    }

    // Set up 240-255 components of targeted buffer : FG

    auto ta1_y_yyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 240);

    auto ta1_y_yyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 241);

    auto ta1_y_yyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 242);

    auto ta1_y_yyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 243);

    auto ta1_y_yyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 244);

    auto ta1_y_yyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 245);

    auto ta1_y_yyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 246);

    auto ta1_y_yyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 247);

    auto ta1_y_yyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 248);

    auto ta1_y_yyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 249);

    auto ta1_y_yyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 250);

    auto ta1_y_yyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 251);

    auto ta1_y_yyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 252);

    auto ta1_y_yyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 253);

    auto ta1_y_yyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 254);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_y_y_xxxx_0,   \
                             ta1_y_y_xxxx_1,   \
                             ta1_y_y_xxxy_0,   \
                             ta1_y_y_xxxy_1,   \
                             ta1_y_y_xxxz_0,   \
                             ta1_y_y_xxxz_1,   \
                             ta1_y_y_xxyy_0,   \
                             ta1_y_y_xxyy_1,   \
                             ta1_y_y_xxyz_0,   \
                             ta1_y_y_xxyz_1,   \
                             ta1_y_y_xxzz_0,   \
                             ta1_y_y_xxzz_1,   \
                             ta1_y_y_xyyy_0,   \
                             ta1_y_y_xyyy_1,   \
                             ta1_y_y_xyyz_0,   \
                             ta1_y_y_xyyz_1,   \
                             ta1_y_y_xyzz_0,   \
                             ta1_y_y_xyzz_1,   \
                             ta1_y_y_xzzz_0,   \
                             ta1_y_y_xzzz_1,   \
                             ta1_y_y_yyyy_0,   \
                             ta1_y_y_yyyy_1,   \
                             ta1_y_y_yyyz_0,   \
                             ta1_y_y_yyyz_1,   \
                             ta1_y_y_yyzz_0,   \
                             ta1_y_y_yyzz_1,   \
                             ta1_y_y_yzzz_0,   \
                             ta1_y_y_yzzz_1,   \
                             ta1_y_y_zzzz_0,   \
                             ta1_y_y_zzzz_1,   \
                             ta1_y_yy_xxx_0,   \
                             ta1_y_yy_xxx_1,   \
                             ta1_y_yy_xxxx_0,  \
                             ta1_y_yy_xxxx_1,  \
                             ta1_y_yy_xxxy_0,  \
                             ta1_y_yy_xxxy_1,  \
                             ta1_y_yy_xxxz_0,  \
                             ta1_y_yy_xxxz_1,  \
                             ta1_y_yy_xxy_0,   \
                             ta1_y_yy_xxy_1,   \
                             ta1_y_yy_xxyy_0,  \
                             ta1_y_yy_xxyy_1,  \
                             ta1_y_yy_xxyz_0,  \
                             ta1_y_yy_xxyz_1,  \
                             ta1_y_yy_xxz_0,   \
                             ta1_y_yy_xxz_1,   \
                             ta1_y_yy_xxzz_0,  \
                             ta1_y_yy_xxzz_1,  \
                             ta1_y_yy_xyy_0,   \
                             ta1_y_yy_xyy_1,   \
                             ta1_y_yy_xyyy_0,  \
                             ta1_y_yy_xyyy_1,  \
                             ta1_y_yy_xyyz_0,  \
                             ta1_y_yy_xyyz_1,  \
                             ta1_y_yy_xyz_0,   \
                             ta1_y_yy_xyz_1,   \
                             ta1_y_yy_xyzz_0,  \
                             ta1_y_yy_xyzz_1,  \
                             ta1_y_yy_xzz_0,   \
                             ta1_y_yy_xzz_1,   \
                             ta1_y_yy_xzzz_0,  \
                             ta1_y_yy_xzzz_1,  \
                             ta1_y_yy_yyy_0,   \
                             ta1_y_yy_yyy_1,   \
                             ta1_y_yy_yyyy_0,  \
                             ta1_y_yy_yyyy_1,  \
                             ta1_y_yy_yyyz_0,  \
                             ta1_y_yy_yyyz_1,  \
                             ta1_y_yy_yyz_0,   \
                             ta1_y_yy_yyz_1,   \
                             ta1_y_yy_yyzz_0,  \
                             ta1_y_yy_yyzz_1,  \
                             ta1_y_yy_yzz_0,   \
                             ta1_y_yy_yzz_1,   \
                             ta1_y_yy_yzzz_0,  \
                             ta1_y_yy_yzzz_1,  \
                             ta1_y_yy_zzz_0,   \
                             ta1_y_yy_zzz_1,   \
                             ta1_y_yy_zzzz_0,  \
                             ta1_y_yy_zzzz_1,  \
                             ta1_y_yyy_xxxx_0, \
                             ta1_y_yyy_xxxy_0, \
                             ta1_y_yyy_xxxz_0, \
                             ta1_y_yyy_xxyy_0, \
                             ta1_y_yyy_xxyz_0, \
                             ta1_y_yyy_xxzz_0, \
                             ta1_y_yyy_xyyy_0, \
                             ta1_y_yyy_xyyz_0, \
                             ta1_y_yyy_xyzz_0, \
                             ta1_y_yyy_xzzz_0, \
                             ta1_y_yyy_yyyy_0, \
                             ta1_y_yyy_yyyz_0, \
                             ta1_y_yyy_yyzz_0, \
                             ta1_y_yyy_yzzz_0, \
                             ta1_y_yyy_zzzz_0, \
                             ta_yy_xxxx_1,     \
                             ta_yy_xxxy_1,     \
                             ta_yy_xxxz_1,     \
                             ta_yy_xxyy_1,     \
                             ta_yy_xxyz_1,     \
                             ta_yy_xxzz_1,     \
                             ta_yy_xyyy_1,     \
                             ta_yy_xyyz_1,     \
                             ta_yy_xyzz_1,     \
                             ta_yy_xzzz_1,     \
                             ta_yy_yyyy_1,     \
                             ta_yy_yyyz_1,     \
                             ta_yy_yyzz_1,     \
                             ta_yy_yzzz_1,     \
                             ta_yy_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyy_xxxx_0[i] = 2.0 * ta1_y_y_xxxx_0[i] * fe_0 - 2.0 * ta1_y_y_xxxx_1[i] * fe_0 + ta_yy_xxxx_1[i] +
                              ta1_y_yy_xxxx_0[i] * pa_y[i] - ta1_y_yy_xxxx_1[i] * pc_y[i];

        ta1_y_yyy_xxxy_0[i] = 2.0 * ta1_y_y_xxxy_0[i] * fe_0 - 2.0 * ta1_y_y_xxxy_1[i] * fe_0 +
                              ta1_y_yy_xxx_0[i] * fe_0 - ta1_y_yy_xxx_1[i] * fe_0 + ta_yy_xxxy_1[i] +
                              ta1_y_yy_xxxy_0[i] * pa_y[i] - ta1_y_yy_xxxy_1[i] * pc_y[i];

        ta1_y_yyy_xxxz_0[i] = 2.0 * ta1_y_y_xxxz_0[i] * fe_0 - 2.0 * ta1_y_y_xxxz_1[i] * fe_0 + ta_yy_xxxz_1[i] +
                              ta1_y_yy_xxxz_0[i] * pa_y[i] - ta1_y_yy_xxxz_1[i] * pc_y[i];

        ta1_y_yyy_xxyy_0[i] = 2.0 * ta1_y_y_xxyy_0[i] * fe_0 - 2.0 * ta1_y_y_xxyy_1[i] * fe_0 +
                              2.0 * ta1_y_yy_xxy_0[i] * fe_0 - 2.0 * ta1_y_yy_xxy_1[i] * fe_0 + ta_yy_xxyy_1[i] +
                              ta1_y_yy_xxyy_0[i] * pa_y[i] - ta1_y_yy_xxyy_1[i] * pc_y[i];

        ta1_y_yyy_xxyz_0[i] = 2.0 * ta1_y_y_xxyz_0[i] * fe_0 - 2.0 * ta1_y_y_xxyz_1[i] * fe_0 +
                              ta1_y_yy_xxz_0[i] * fe_0 - ta1_y_yy_xxz_1[i] * fe_0 + ta_yy_xxyz_1[i] +
                              ta1_y_yy_xxyz_0[i] * pa_y[i] - ta1_y_yy_xxyz_1[i] * pc_y[i];

        ta1_y_yyy_xxzz_0[i] = 2.0 * ta1_y_y_xxzz_0[i] * fe_0 - 2.0 * ta1_y_y_xxzz_1[i] * fe_0 + ta_yy_xxzz_1[i] +
                              ta1_y_yy_xxzz_0[i] * pa_y[i] - ta1_y_yy_xxzz_1[i] * pc_y[i];

        ta1_y_yyy_xyyy_0[i] = 2.0 * ta1_y_y_xyyy_0[i] * fe_0 - 2.0 * ta1_y_y_xyyy_1[i] * fe_0 +
                              3.0 * ta1_y_yy_xyy_0[i] * fe_0 - 3.0 * ta1_y_yy_xyy_1[i] * fe_0 + ta_yy_xyyy_1[i] +
                              ta1_y_yy_xyyy_0[i] * pa_y[i] - ta1_y_yy_xyyy_1[i] * pc_y[i];

        ta1_y_yyy_xyyz_0[i] = 2.0 * ta1_y_y_xyyz_0[i] * fe_0 - 2.0 * ta1_y_y_xyyz_1[i] * fe_0 +
                              2.0 * ta1_y_yy_xyz_0[i] * fe_0 - 2.0 * ta1_y_yy_xyz_1[i] * fe_0 + ta_yy_xyyz_1[i] +
                              ta1_y_yy_xyyz_0[i] * pa_y[i] - ta1_y_yy_xyyz_1[i] * pc_y[i];

        ta1_y_yyy_xyzz_0[i] = 2.0 * ta1_y_y_xyzz_0[i] * fe_0 - 2.0 * ta1_y_y_xyzz_1[i] * fe_0 +
                              ta1_y_yy_xzz_0[i] * fe_0 - ta1_y_yy_xzz_1[i] * fe_0 + ta_yy_xyzz_1[i] +
                              ta1_y_yy_xyzz_0[i] * pa_y[i] - ta1_y_yy_xyzz_1[i] * pc_y[i];

        ta1_y_yyy_xzzz_0[i] = 2.0 * ta1_y_y_xzzz_0[i] * fe_0 - 2.0 * ta1_y_y_xzzz_1[i] * fe_0 + ta_yy_xzzz_1[i] +
                              ta1_y_yy_xzzz_0[i] * pa_y[i] - ta1_y_yy_xzzz_1[i] * pc_y[i];

        ta1_y_yyy_yyyy_0[i] = 2.0 * ta1_y_y_yyyy_0[i] * fe_0 - 2.0 * ta1_y_y_yyyy_1[i] * fe_0 +
                              4.0 * ta1_y_yy_yyy_0[i] * fe_0 - 4.0 * ta1_y_yy_yyy_1[i] * fe_0 + ta_yy_yyyy_1[i] +
                              ta1_y_yy_yyyy_0[i] * pa_y[i] - ta1_y_yy_yyyy_1[i] * pc_y[i];

        ta1_y_yyy_yyyz_0[i] = 2.0 * ta1_y_y_yyyz_0[i] * fe_0 - 2.0 * ta1_y_y_yyyz_1[i] * fe_0 +
                              3.0 * ta1_y_yy_yyz_0[i] * fe_0 - 3.0 * ta1_y_yy_yyz_1[i] * fe_0 + ta_yy_yyyz_1[i] +
                              ta1_y_yy_yyyz_0[i] * pa_y[i] - ta1_y_yy_yyyz_1[i] * pc_y[i];

        ta1_y_yyy_yyzz_0[i] = 2.0 * ta1_y_y_yyzz_0[i] * fe_0 - 2.0 * ta1_y_y_yyzz_1[i] * fe_0 +
                              2.0 * ta1_y_yy_yzz_0[i] * fe_0 - 2.0 * ta1_y_yy_yzz_1[i] * fe_0 + ta_yy_yyzz_1[i] +
                              ta1_y_yy_yyzz_0[i] * pa_y[i] - ta1_y_yy_yyzz_1[i] * pc_y[i];

        ta1_y_yyy_yzzz_0[i] = 2.0 * ta1_y_y_yzzz_0[i] * fe_0 - 2.0 * ta1_y_y_yzzz_1[i] * fe_0 +
                              ta1_y_yy_zzz_0[i] * fe_0 - ta1_y_yy_zzz_1[i] * fe_0 + ta_yy_yzzz_1[i] +
                              ta1_y_yy_yzzz_0[i] * pa_y[i] - ta1_y_yy_yzzz_1[i] * pc_y[i];

        ta1_y_yyy_zzzz_0[i] = 2.0 * ta1_y_y_zzzz_0[i] * fe_0 - 2.0 * ta1_y_y_zzzz_1[i] * fe_0 + ta_yy_zzzz_1[i] +
                              ta1_y_yy_zzzz_0[i] * pa_y[i] - ta1_y_yy_zzzz_1[i] * pc_y[i];
    }

    // Set up 255-270 components of targeted buffer : FG

    auto ta1_y_yyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 255);

    auto ta1_y_yyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 256);

    auto ta1_y_yyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 257);

    auto ta1_y_yyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 258);

    auto ta1_y_yyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 259);

    auto ta1_y_yyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 260);

    auto ta1_y_yyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 261);

    auto ta1_y_yyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 262);

    auto ta1_y_yyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 263);

    auto ta1_y_yyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 264);

    auto ta1_y_yyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 265);

    auto ta1_y_yyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 266);

    auto ta1_y_yyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 267);

    auto ta1_y_yyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 268);

    auto ta1_y_yyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 269);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta1_y_yy_xxx_0,   \
                             ta1_y_yy_xxx_1,   \
                             ta1_y_yy_xxxx_0,  \
                             ta1_y_yy_xxxx_1,  \
                             ta1_y_yy_xxxy_0,  \
                             ta1_y_yy_xxxy_1,  \
                             ta1_y_yy_xxxz_0,  \
                             ta1_y_yy_xxxz_1,  \
                             ta1_y_yy_xxy_0,   \
                             ta1_y_yy_xxy_1,   \
                             ta1_y_yy_xxyy_0,  \
                             ta1_y_yy_xxyy_1,  \
                             ta1_y_yy_xxyz_0,  \
                             ta1_y_yy_xxyz_1,  \
                             ta1_y_yy_xxz_0,   \
                             ta1_y_yy_xxz_1,   \
                             ta1_y_yy_xxzz_0,  \
                             ta1_y_yy_xxzz_1,  \
                             ta1_y_yy_xyy_0,   \
                             ta1_y_yy_xyy_1,   \
                             ta1_y_yy_xyyy_0,  \
                             ta1_y_yy_xyyy_1,  \
                             ta1_y_yy_xyyz_0,  \
                             ta1_y_yy_xyyz_1,  \
                             ta1_y_yy_xyz_0,   \
                             ta1_y_yy_xyz_1,   \
                             ta1_y_yy_xyzz_0,  \
                             ta1_y_yy_xyzz_1,  \
                             ta1_y_yy_xzz_0,   \
                             ta1_y_yy_xzz_1,   \
                             ta1_y_yy_xzzz_0,  \
                             ta1_y_yy_xzzz_1,  \
                             ta1_y_yy_yyy_0,   \
                             ta1_y_yy_yyy_1,   \
                             ta1_y_yy_yyyy_0,  \
                             ta1_y_yy_yyyy_1,  \
                             ta1_y_yy_yyyz_0,  \
                             ta1_y_yy_yyyz_1,  \
                             ta1_y_yy_yyz_0,   \
                             ta1_y_yy_yyz_1,   \
                             ta1_y_yy_yyzz_0,  \
                             ta1_y_yy_yyzz_1,  \
                             ta1_y_yy_yzz_0,   \
                             ta1_y_yy_yzz_1,   \
                             ta1_y_yy_yzzz_0,  \
                             ta1_y_yy_yzzz_1,  \
                             ta1_y_yy_zzz_0,   \
                             ta1_y_yy_zzz_1,   \
                             ta1_y_yy_zzzz_0,  \
                             ta1_y_yy_zzzz_1,  \
                             ta1_y_yyz_xxxx_0, \
                             ta1_y_yyz_xxxy_0, \
                             ta1_y_yyz_xxxz_0, \
                             ta1_y_yyz_xxyy_0, \
                             ta1_y_yyz_xxyz_0, \
                             ta1_y_yyz_xxzz_0, \
                             ta1_y_yyz_xyyy_0, \
                             ta1_y_yyz_xyyz_0, \
                             ta1_y_yyz_xyzz_0, \
                             ta1_y_yyz_xzzz_0, \
                             ta1_y_yyz_yyyy_0, \
                             ta1_y_yyz_yyyz_0, \
                             ta1_y_yyz_yyzz_0, \
                             ta1_y_yyz_yzzz_0, \
                             ta1_y_yyz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyz_xxxx_0[i] = ta1_y_yy_xxxx_0[i] * pa_z[i] - ta1_y_yy_xxxx_1[i] * pc_z[i];

        ta1_y_yyz_xxxy_0[i] = ta1_y_yy_xxxy_0[i] * pa_z[i] - ta1_y_yy_xxxy_1[i] * pc_z[i];

        ta1_y_yyz_xxxz_0[i] = ta1_y_yy_xxx_0[i] * fe_0 - ta1_y_yy_xxx_1[i] * fe_0 + ta1_y_yy_xxxz_0[i] * pa_z[i] -
                              ta1_y_yy_xxxz_1[i] * pc_z[i];

        ta1_y_yyz_xxyy_0[i] = ta1_y_yy_xxyy_0[i] * pa_z[i] - ta1_y_yy_xxyy_1[i] * pc_z[i];

        ta1_y_yyz_xxyz_0[i] = ta1_y_yy_xxy_0[i] * fe_0 - ta1_y_yy_xxy_1[i] * fe_0 + ta1_y_yy_xxyz_0[i] * pa_z[i] -
                              ta1_y_yy_xxyz_1[i] * pc_z[i];

        ta1_y_yyz_xxzz_0[i] = 2.0 * ta1_y_yy_xxz_0[i] * fe_0 - 2.0 * ta1_y_yy_xxz_1[i] * fe_0 +
                              ta1_y_yy_xxzz_0[i] * pa_z[i] - ta1_y_yy_xxzz_1[i] * pc_z[i];

        ta1_y_yyz_xyyy_0[i] = ta1_y_yy_xyyy_0[i] * pa_z[i] - ta1_y_yy_xyyy_1[i] * pc_z[i];

        ta1_y_yyz_xyyz_0[i] = ta1_y_yy_xyy_0[i] * fe_0 - ta1_y_yy_xyy_1[i] * fe_0 + ta1_y_yy_xyyz_0[i] * pa_z[i] -
                              ta1_y_yy_xyyz_1[i] * pc_z[i];

        ta1_y_yyz_xyzz_0[i] = 2.0 * ta1_y_yy_xyz_0[i] * fe_0 - 2.0 * ta1_y_yy_xyz_1[i] * fe_0 +
                              ta1_y_yy_xyzz_0[i] * pa_z[i] - ta1_y_yy_xyzz_1[i] * pc_z[i];

        ta1_y_yyz_xzzz_0[i] = 3.0 * ta1_y_yy_xzz_0[i] * fe_0 - 3.0 * ta1_y_yy_xzz_1[i] * fe_0 +
                              ta1_y_yy_xzzz_0[i] * pa_z[i] - ta1_y_yy_xzzz_1[i] * pc_z[i];

        ta1_y_yyz_yyyy_0[i] = ta1_y_yy_yyyy_0[i] * pa_z[i] - ta1_y_yy_yyyy_1[i] * pc_z[i];

        ta1_y_yyz_yyyz_0[i] = ta1_y_yy_yyy_0[i] * fe_0 - ta1_y_yy_yyy_1[i] * fe_0 + ta1_y_yy_yyyz_0[i] * pa_z[i] -
                              ta1_y_yy_yyyz_1[i] * pc_z[i];

        ta1_y_yyz_yyzz_0[i] = 2.0 * ta1_y_yy_yyz_0[i] * fe_0 - 2.0 * ta1_y_yy_yyz_1[i] * fe_0 +
                              ta1_y_yy_yyzz_0[i] * pa_z[i] - ta1_y_yy_yyzz_1[i] * pc_z[i];

        ta1_y_yyz_yzzz_0[i] = 3.0 * ta1_y_yy_yzz_0[i] * fe_0 - 3.0 * ta1_y_yy_yzz_1[i] * fe_0 +
                              ta1_y_yy_yzzz_0[i] * pa_z[i] - ta1_y_yy_yzzz_1[i] * pc_z[i];

        ta1_y_yyz_zzzz_0[i] = 4.0 * ta1_y_yy_zzz_0[i] * fe_0 - 4.0 * ta1_y_yy_zzz_1[i] * fe_0 +
                              ta1_y_yy_zzzz_0[i] * pa_z[i] - ta1_y_yy_zzzz_1[i] * pc_z[i];
    }

    // Set up 270-285 components of targeted buffer : FG

    auto ta1_y_yzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 270);

    auto ta1_y_yzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 271);

    auto ta1_y_yzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 272);

    auto ta1_y_yzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 273);

    auto ta1_y_yzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 274);

    auto ta1_y_yzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 275);

    auto ta1_y_yzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 276);

    auto ta1_y_yzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 277);

    auto ta1_y_yzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 278);

    auto ta1_y_yzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 279);

    auto ta1_y_yzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 280);

    auto ta1_y_yzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 281);

    auto ta1_y_yzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 282);

    auto ta1_y_yzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 283);

    auto ta1_y_yzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 284);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_y_y_xxxy_0,   \
                             ta1_y_y_xxxy_1,   \
                             ta1_y_y_xxyy_0,   \
                             ta1_y_y_xxyy_1,   \
                             ta1_y_y_xyyy_0,   \
                             ta1_y_y_xyyy_1,   \
                             ta1_y_y_yyyy_0,   \
                             ta1_y_y_yyyy_1,   \
                             ta1_y_yz_xxxy_0,  \
                             ta1_y_yz_xxxy_1,  \
                             ta1_y_yz_xxyy_0,  \
                             ta1_y_yz_xxyy_1,  \
                             ta1_y_yz_xyyy_0,  \
                             ta1_y_yz_xyyy_1,  \
                             ta1_y_yz_yyyy_0,  \
                             ta1_y_yz_yyyy_1,  \
                             ta1_y_yzz_xxxx_0, \
                             ta1_y_yzz_xxxy_0, \
                             ta1_y_yzz_xxxz_0, \
                             ta1_y_yzz_xxyy_0, \
                             ta1_y_yzz_xxyz_0, \
                             ta1_y_yzz_xxzz_0, \
                             ta1_y_yzz_xyyy_0, \
                             ta1_y_yzz_xyyz_0, \
                             ta1_y_yzz_xyzz_0, \
                             ta1_y_yzz_xzzz_0, \
                             ta1_y_yzz_yyyy_0, \
                             ta1_y_yzz_yyyz_0, \
                             ta1_y_yzz_yyzz_0, \
                             ta1_y_yzz_yzzz_0, \
                             ta1_y_yzz_zzzz_0, \
                             ta1_y_zz_xxxx_0,  \
                             ta1_y_zz_xxxx_1,  \
                             ta1_y_zz_xxxz_0,  \
                             ta1_y_zz_xxxz_1,  \
                             ta1_y_zz_xxyz_0,  \
                             ta1_y_zz_xxyz_1,  \
                             ta1_y_zz_xxz_0,   \
                             ta1_y_zz_xxz_1,   \
                             ta1_y_zz_xxzz_0,  \
                             ta1_y_zz_xxzz_1,  \
                             ta1_y_zz_xyyz_0,  \
                             ta1_y_zz_xyyz_1,  \
                             ta1_y_zz_xyz_0,   \
                             ta1_y_zz_xyz_1,   \
                             ta1_y_zz_xyzz_0,  \
                             ta1_y_zz_xyzz_1,  \
                             ta1_y_zz_xzz_0,   \
                             ta1_y_zz_xzz_1,   \
                             ta1_y_zz_xzzz_0,  \
                             ta1_y_zz_xzzz_1,  \
                             ta1_y_zz_yyyz_0,  \
                             ta1_y_zz_yyyz_1,  \
                             ta1_y_zz_yyz_0,   \
                             ta1_y_zz_yyz_1,   \
                             ta1_y_zz_yyzz_0,  \
                             ta1_y_zz_yyzz_1,  \
                             ta1_y_zz_yzz_0,   \
                             ta1_y_zz_yzz_1,   \
                             ta1_y_zz_yzzz_0,  \
                             ta1_y_zz_yzzz_1,  \
                             ta1_y_zz_zzz_0,   \
                             ta1_y_zz_zzz_1,   \
                             ta1_y_zz_zzzz_0,  \
                             ta1_y_zz_zzzz_1,  \
                             ta_zz_xxxx_1,     \
                             ta_zz_xxxz_1,     \
                             ta_zz_xxyz_1,     \
                             ta_zz_xxzz_1,     \
                             ta_zz_xyyz_1,     \
                             ta_zz_xyzz_1,     \
                             ta_zz_xzzz_1,     \
                             ta_zz_yyyz_1,     \
                             ta_zz_yyzz_1,     \
                             ta_zz_yzzz_1,     \
                             ta_zz_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yzz_xxxx_0[i] = ta_zz_xxxx_1[i] + ta1_y_zz_xxxx_0[i] * pa_y[i] - ta1_y_zz_xxxx_1[i] * pc_y[i];

        ta1_y_yzz_xxxy_0[i] = ta1_y_y_xxxy_0[i] * fe_0 - ta1_y_y_xxxy_1[i] * fe_0 + ta1_y_yz_xxxy_0[i] * pa_z[i] -
                              ta1_y_yz_xxxy_1[i] * pc_z[i];

        ta1_y_yzz_xxxz_0[i] = ta_zz_xxxz_1[i] + ta1_y_zz_xxxz_0[i] * pa_y[i] - ta1_y_zz_xxxz_1[i] * pc_y[i];

        ta1_y_yzz_xxyy_0[i] = ta1_y_y_xxyy_0[i] * fe_0 - ta1_y_y_xxyy_1[i] * fe_0 + ta1_y_yz_xxyy_0[i] * pa_z[i] -
                              ta1_y_yz_xxyy_1[i] * pc_z[i];

        ta1_y_yzz_xxyz_0[i] = ta1_y_zz_xxz_0[i] * fe_0 - ta1_y_zz_xxz_1[i] * fe_0 + ta_zz_xxyz_1[i] +
                              ta1_y_zz_xxyz_0[i] * pa_y[i] - ta1_y_zz_xxyz_1[i] * pc_y[i];

        ta1_y_yzz_xxzz_0[i] = ta_zz_xxzz_1[i] + ta1_y_zz_xxzz_0[i] * pa_y[i] - ta1_y_zz_xxzz_1[i] * pc_y[i];

        ta1_y_yzz_xyyy_0[i] = ta1_y_y_xyyy_0[i] * fe_0 - ta1_y_y_xyyy_1[i] * fe_0 + ta1_y_yz_xyyy_0[i] * pa_z[i] -
                              ta1_y_yz_xyyy_1[i] * pc_z[i];

        ta1_y_yzz_xyyz_0[i] = 2.0 * ta1_y_zz_xyz_0[i] * fe_0 - 2.0 * ta1_y_zz_xyz_1[i] * fe_0 + ta_zz_xyyz_1[i] +
                              ta1_y_zz_xyyz_0[i] * pa_y[i] - ta1_y_zz_xyyz_1[i] * pc_y[i];

        ta1_y_yzz_xyzz_0[i] = ta1_y_zz_xzz_0[i] * fe_0 - ta1_y_zz_xzz_1[i] * fe_0 + ta_zz_xyzz_1[i] +
                              ta1_y_zz_xyzz_0[i] * pa_y[i] - ta1_y_zz_xyzz_1[i] * pc_y[i];

        ta1_y_yzz_xzzz_0[i] = ta_zz_xzzz_1[i] + ta1_y_zz_xzzz_0[i] * pa_y[i] - ta1_y_zz_xzzz_1[i] * pc_y[i];

        ta1_y_yzz_yyyy_0[i] = ta1_y_y_yyyy_0[i] * fe_0 - ta1_y_y_yyyy_1[i] * fe_0 + ta1_y_yz_yyyy_0[i] * pa_z[i] -
                              ta1_y_yz_yyyy_1[i] * pc_z[i];

        ta1_y_yzz_yyyz_0[i] = 3.0 * ta1_y_zz_yyz_0[i] * fe_0 - 3.0 * ta1_y_zz_yyz_1[i] * fe_0 + ta_zz_yyyz_1[i] +
                              ta1_y_zz_yyyz_0[i] * pa_y[i] - ta1_y_zz_yyyz_1[i] * pc_y[i];

        ta1_y_yzz_yyzz_0[i] = 2.0 * ta1_y_zz_yzz_0[i] * fe_0 - 2.0 * ta1_y_zz_yzz_1[i] * fe_0 + ta_zz_yyzz_1[i] +
                              ta1_y_zz_yyzz_0[i] * pa_y[i] - ta1_y_zz_yyzz_1[i] * pc_y[i];

        ta1_y_yzz_yzzz_0[i] = ta1_y_zz_zzz_0[i] * fe_0 - ta1_y_zz_zzz_1[i] * fe_0 + ta_zz_yzzz_1[i] +
                              ta1_y_zz_yzzz_0[i] * pa_y[i] - ta1_y_zz_yzzz_1[i] * pc_y[i];

        ta1_y_yzz_zzzz_0[i] = ta_zz_zzzz_1[i] + ta1_y_zz_zzzz_0[i] * pa_y[i] - ta1_y_zz_zzzz_1[i] * pc_y[i];
    }

    // Set up 285-300 components of targeted buffer : FG

    auto ta1_y_zzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 285);

    auto ta1_y_zzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 286);

    auto ta1_y_zzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 287);

    auto ta1_y_zzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 288);

    auto ta1_y_zzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 289);

    auto ta1_y_zzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 290);

    auto ta1_y_zzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 291);

    auto ta1_y_zzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 292);

    auto ta1_y_zzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 293);

    auto ta1_y_zzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 294);

    auto ta1_y_zzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 295);

    auto ta1_y_zzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 296);

    auto ta1_y_zzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 297);

    auto ta1_y_zzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 298);

    auto ta1_y_zzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 299);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta1_y_z_xxxx_0,   \
                             ta1_y_z_xxxx_1,   \
                             ta1_y_z_xxxy_0,   \
                             ta1_y_z_xxxy_1,   \
                             ta1_y_z_xxxz_0,   \
                             ta1_y_z_xxxz_1,   \
                             ta1_y_z_xxyy_0,   \
                             ta1_y_z_xxyy_1,   \
                             ta1_y_z_xxyz_0,   \
                             ta1_y_z_xxyz_1,   \
                             ta1_y_z_xxzz_0,   \
                             ta1_y_z_xxzz_1,   \
                             ta1_y_z_xyyy_0,   \
                             ta1_y_z_xyyy_1,   \
                             ta1_y_z_xyyz_0,   \
                             ta1_y_z_xyyz_1,   \
                             ta1_y_z_xyzz_0,   \
                             ta1_y_z_xyzz_1,   \
                             ta1_y_z_xzzz_0,   \
                             ta1_y_z_xzzz_1,   \
                             ta1_y_z_yyyy_0,   \
                             ta1_y_z_yyyy_1,   \
                             ta1_y_z_yyyz_0,   \
                             ta1_y_z_yyyz_1,   \
                             ta1_y_z_yyzz_0,   \
                             ta1_y_z_yyzz_1,   \
                             ta1_y_z_yzzz_0,   \
                             ta1_y_z_yzzz_1,   \
                             ta1_y_z_zzzz_0,   \
                             ta1_y_z_zzzz_1,   \
                             ta1_y_zz_xxx_0,   \
                             ta1_y_zz_xxx_1,   \
                             ta1_y_zz_xxxx_0,  \
                             ta1_y_zz_xxxx_1,  \
                             ta1_y_zz_xxxy_0,  \
                             ta1_y_zz_xxxy_1,  \
                             ta1_y_zz_xxxz_0,  \
                             ta1_y_zz_xxxz_1,  \
                             ta1_y_zz_xxy_0,   \
                             ta1_y_zz_xxy_1,   \
                             ta1_y_zz_xxyy_0,  \
                             ta1_y_zz_xxyy_1,  \
                             ta1_y_zz_xxyz_0,  \
                             ta1_y_zz_xxyz_1,  \
                             ta1_y_zz_xxz_0,   \
                             ta1_y_zz_xxz_1,   \
                             ta1_y_zz_xxzz_0,  \
                             ta1_y_zz_xxzz_1,  \
                             ta1_y_zz_xyy_0,   \
                             ta1_y_zz_xyy_1,   \
                             ta1_y_zz_xyyy_0,  \
                             ta1_y_zz_xyyy_1,  \
                             ta1_y_zz_xyyz_0,  \
                             ta1_y_zz_xyyz_1,  \
                             ta1_y_zz_xyz_0,   \
                             ta1_y_zz_xyz_1,   \
                             ta1_y_zz_xyzz_0,  \
                             ta1_y_zz_xyzz_1,  \
                             ta1_y_zz_xzz_0,   \
                             ta1_y_zz_xzz_1,   \
                             ta1_y_zz_xzzz_0,  \
                             ta1_y_zz_xzzz_1,  \
                             ta1_y_zz_yyy_0,   \
                             ta1_y_zz_yyy_1,   \
                             ta1_y_zz_yyyy_0,  \
                             ta1_y_zz_yyyy_1,  \
                             ta1_y_zz_yyyz_0,  \
                             ta1_y_zz_yyyz_1,  \
                             ta1_y_zz_yyz_0,   \
                             ta1_y_zz_yyz_1,   \
                             ta1_y_zz_yyzz_0,  \
                             ta1_y_zz_yyzz_1,  \
                             ta1_y_zz_yzz_0,   \
                             ta1_y_zz_yzz_1,   \
                             ta1_y_zz_yzzz_0,  \
                             ta1_y_zz_yzzz_1,  \
                             ta1_y_zz_zzz_0,   \
                             ta1_y_zz_zzz_1,   \
                             ta1_y_zz_zzzz_0,  \
                             ta1_y_zz_zzzz_1,  \
                             ta1_y_zzz_xxxx_0, \
                             ta1_y_zzz_xxxy_0, \
                             ta1_y_zzz_xxxz_0, \
                             ta1_y_zzz_xxyy_0, \
                             ta1_y_zzz_xxyz_0, \
                             ta1_y_zzz_xxzz_0, \
                             ta1_y_zzz_xyyy_0, \
                             ta1_y_zzz_xyyz_0, \
                             ta1_y_zzz_xyzz_0, \
                             ta1_y_zzz_xzzz_0, \
                             ta1_y_zzz_yyyy_0, \
                             ta1_y_zzz_yyyz_0, \
                             ta1_y_zzz_yyzz_0, \
                             ta1_y_zzz_yzzz_0, \
                             ta1_y_zzz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_zzz_xxxx_0[i] = 2.0 * ta1_y_z_xxxx_0[i] * fe_0 - 2.0 * ta1_y_z_xxxx_1[i] * fe_0 +
                              ta1_y_zz_xxxx_0[i] * pa_z[i] - ta1_y_zz_xxxx_1[i] * pc_z[i];

        ta1_y_zzz_xxxy_0[i] = 2.0 * ta1_y_z_xxxy_0[i] * fe_0 - 2.0 * ta1_y_z_xxxy_1[i] * fe_0 +
                              ta1_y_zz_xxxy_0[i] * pa_z[i] - ta1_y_zz_xxxy_1[i] * pc_z[i];

        ta1_y_zzz_xxxz_0[i] = 2.0 * ta1_y_z_xxxz_0[i] * fe_0 - 2.0 * ta1_y_z_xxxz_1[i] * fe_0 +
                              ta1_y_zz_xxx_0[i] * fe_0 - ta1_y_zz_xxx_1[i] * fe_0 + ta1_y_zz_xxxz_0[i] * pa_z[i] -
                              ta1_y_zz_xxxz_1[i] * pc_z[i];

        ta1_y_zzz_xxyy_0[i] = 2.0 * ta1_y_z_xxyy_0[i] * fe_0 - 2.0 * ta1_y_z_xxyy_1[i] * fe_0 +
                              ta1_y_zz_xxyy_0[i] * pa_z[i] - ta1_y_zz_xxyy_1[i] * pc_z[i];

        ta1_y_zzz_xxyz_0[i] = 2.0 * ta1_y_z_xxyz_0[i] * fe_0 - 2.0 * ta1_y_z_xxyz_1[i] * fe_0 +
                              ta1_y_zz_xxy_0[i] * fe_0 - ta1_y_zz_xxy_1[i] * fe_0 + ta1_y_zz_xxyz_0[i] * pa_z[i] -
                              ta1_y_zz_xxyz_1[i] * pc_z[i];

        ta1_y_zzz_xxzz_0[i] = 2.0 * ta1_y_z_xxzz_0[i] * fe_0 - 2.0 * ta1_y_z_xxzz_1[i] * fe_0 +
                              2.0 * ta1_y_zz_xxz_0[i] * fe_0 - 2.0 * ta1_y_zz_xxz_1[i] * fe_0 +
                              ta1_y_zz_xxzz_0[i] * pa_z[i] - ta1_y_zz_xxzz_1[i] * pc_z[i];

        ta1_y_zzz_xyyy_0[i] = 2.0 * ta1_y_z_xyyy_0[i] * fe_0 - 2.0 * ta1_y_z_xyyy_1[i] * fe_0 +
                              ta1_y_zz_xyyy_0[i] * pa_z[i] - ta1_y_zz_xyyy_1[i] * pc_z[i];

        ta1_y_zzz_xyyz_0[i] = 2.0 * ta1_y_z_xyyz_0[i] * fe_0 - 2.0 * ta1_y_z_xyyz_1[i] * fe_0 +
                              ta1_y_zz_xyy_0[i] * fe_0 - ta1_y_zz_xyy_1[i] * fe_0 + ta1_y_zz_xyyz_0[i] * pa_z[i] -
                              ta1_y_zz_xyyz_1[i] * pc_z[i];

        ta1_y_zzz_xyzz_0[i] = 2.0 * ta1_y_z_xyzz_0[i] * fe_0 - 2.0 * ta1_y_z_xyzz_1[i] * fe_0 +
                              2.0 * ta1_y_zz_xyz_0[i] * fe_0 - 2.0 * ta1_y_zz_xyz_1[i] * fe_0 +
                              ta1_y_zz_xyzz_0[i] * pa_z[i] - ta1_y_zz_xyzz_1[i] * pc_z[i];

        ta1_y_zzz_xzzz_0[i] = 2.0 * ta1_y_z_xzzz_0[i] * fe_0 - 2.0 * ta1_y_z_xzzz_1[i] * fe_0 +
                              3.0 * ta1_y_zz_xzz_0[i] * fe_0 - 3.0 * ta1_y_zz_xzz_1[i] * fe_0 +
                              ta1_y_zz_xzzz_0[i] * pa_z[i] - ta1_y_zz_xzzz_1[i] * pc_z[i];

        ta1_y_zzz_yyyy_0[i] = 2.0 * ta1_y_z_yyyy_0[i] * fe_0 - 2.0 * ta1_y_z_yyyy_1[i] * fe_0 +
                              ta1_y_zz_yyyy_0[i] * pa_z[i] - ta1_y_zz_yyyy_1[i] * pc_z[i];

        ta1_y_zzz_yyyz_0[i] = 2.0 * ta1_y_z_yyyz_0[i] * fe_0 - 2.0 * ta1_y_z_yyyz_1[i] * fe_0 +
                              ta1_y_zz_yyy_0[i] * fe_0 - ta1_y_zz_yyy_1[i] * fe_0 + ta1_y_zz_yyyz_0[i] * pa_z[i] -
                              ta1_y_zz_yyyz_1[i] * pc_z[i];

        ta1_y_zzz_yyzz_0[i] = 2.0 * ta1_y_z_yyzz_0[i] * fe_0 - 2.0 * ta1_y_z_yyzz_1[i] * fe_0 +
                              2.0 * ta1_y_zz_yyz_0[i] * fe_0 - 2.0 * ta1_y_zz_yyz_1[i] * fe_0 +
                              ta1_y_zz_yyzz_0[i] * pa_z[i] - ta1_y_zz_yyzz_1[i] * pc_z[i];

        ta1_y_zzz_yzzz_0[i] = 2.0 * ta1_y_z_yzzz_0[i] * fe_0 - 2.0 * ta1_y_z_yzzz_1[i] * fe_0 +
                              3.0 * ta1_y_zz_yzz_0[i] * fe_0 - 3.0 * ta1_y_zz_yzz_1[i] * fe_0 +
                              ta1_y_zz_yzzz_0[i] * pa_z[i] - ta1_y_zz_yzzz_1[i] * pc_z[i];

        ta1_y_zzz_zzzz_0[i] = 2.0 * ta1_y_z_zzzz_0[i] * fe_0 - 2.0 * ta1_y_z_zzzz_1[i] * fe_0 +
                              4.0 * ta1_y_zz_zzz_0[i] * fe_0 - 4.0 * ta1_y_zz_zzz_1[i] * fe_0 +
                              ta1_y_zz_zzzz_0[i] * pa_z[i] - ta1_y_zz_zzzz_1[i] * pc_z[i];
    }

    // Set up 300-315 components of targeted buffer : FG

    auto ta1_z_xxx_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 300);

    auto ta1_z_xxx_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 301);

    auto ta1_z_xxx_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 302);

    auto ta1_z_xxx_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 303);

    auto ta1_z_xxx_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 304);

    auto ta1_z_xxx_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 305);

    auto ta1_z_xxx_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 306);

    auto ta1_z_xxx_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 307);

    auto ta1_z_xxx_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 308);

    auto ta1_z_xxx_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 309);

    auto ta1_z_xxx_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 310);

    auto ta1_z_xxx_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 311);

    auto ta1_z_xxx_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 312);

    auto ta1_z_xxx_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 313);

    auto ta1_z_xxx_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 314);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_z_x_xxxx_0,   \
                             ta1_z_x_xxxx_1,   \
                             ta1_z_x_xxxy_0,   \
                             ta1_z_x_xxxy_1,   \
                             ta1_z_x_xxxz_0,   \
                             ta1_z_x_xxxz_1,   \
                             ta1_z_x_xxyy_0,   \
                             ta1_z_x_xxyy_1,   \
                             ta1_z_x_xxyz_0,   \
                             ta1_z_x_xxyz_1,   \
                             ta1_z_x_xxzz_0,   \
                             ta1_z_x_xxzz_1,   \
                             ta1_z_x_xyyy_0,   \
                             ta1_z_x_xyyy_1,   \
                             ta1_z_x_xyyz_0,   \
                             ta1_z_x_xyyz_1,   \
                             ta1_z_x_xyzz_0,   \
                             ta1_z_x_xyzz_1,   \
                             ta1_z_x_xzzz_0,   \
                             ta1_z_x_xzzz_1,   \
                             ta1_z_x_yyyy_0,   \
                             ta1_z_x_yyyy_1,   \
                             ta1_z_x_yyyz_0,   \
                             ta1_z_x_yyyz_1,   \
                             ta1_z_x_yyzz_0,   \
                             ta1_z_x_yyzz_1,   \
                             ta1_z_x_yzzz_0,   \
                             ta1_z_x_yzzz_1,   \
                             ta1_z_x_zzzz_0,   \
                             ta1_z_x_zzzz_1,   \
                             ta1_z_xx_xxx_0,   \
                             ta1_z_xx_xxx_1,   \
                             ta1_z_xx_xxxx_0,  \
                             ta1_z_xx_xxxx_1,  \
                             ta1_z_xx_xxxy_0,  \
                             ta1_z_xx_xxxy_1,  \
                             ta1_z_xx_xxxz_0,  \
                             ta1_z_xx_xxxz_1,  \
                             ta1_z_xx_xxy_0,   \
                             ta1_z_xx_xxy_1,   \
                             ta1_z_xx_xxyy_0,  \
                             ta1_z_xx_xxyy_1,  \
                             ta1_z_xx_xxyz_0,  \
                             ta1_z_xx_xxyz_1,  \
                             ta1_z_xx_xxz_0,   \
                             ta1_z_xx_xxz_1,   \
                             ta1_z_xx_xxzz_0,  \
                             ta1_z_xx_xxzz_1,  \
                             ta1_z_xx_xyy_0,   \
                             ta1_z_xx_xyy_1,   \
                             ta1_z_xx_xyyy_0,  \
                             ta1_z_xx_xyyy_1,  \
                             ta1_z_xx_xyyz_0,  \
                             ta1_z_xx_xyyz_1,  \
                             ta1_z_xx_xyz_0,   \
                             ta1_z_xx_xyz_1,   \
                             ta1_z_xx_xyzz_0,  \
                             ta1_z_xx_xyzz_1,  \
                             ta1_z_xx_xzz_0,   \
                             ta1_z_xx_xzz_1,   \
                             ta1_z_xx_xzzz_0,  \
                             ta1_z_xx_xzzz_1,  \
                             ta1_z_xx_yyy_0,   \
                             ta1_z_xx_yyy_1,   \
                             ta1_z_xx_yyyy_0,  \
                             ta1_z_xx_yyyy_1,  \
                             ta1_z_xx_yyyz_0,  \
                             ta1_z_xx_yyyz_1,  \
                             ta1_z_xx_yyz_0,   \
                             ta1_z_xx_yyz_1,   \
                             ta1_z_xx_yyzz_0,  \
                             ta1_z_xx_yyzz_1,  \
                             ta1_z_xx_yzz_0,   \
                             ta1_z_xx_yzz_1,   \
                             ta1_z_xx_yzzz_0,  \
                             ta1_z_xx_yzzz_1,  \
                             ta1_z_xx_zzz_0,   \
                             ta1_z_xx_zzz_1,   \
                             ta1_z_xx_zzzz_0,  \
                             ta1_z_xx_zzzz_1,  \
                             ta1_z_xxx_xxxx_0, \
                             ta1_z_xxx_xxxy_0, \
                             ta1_z_xxx_xxxz_0, \
                             ta1_z_xxx_xxyy_0, \
                             ta1_z_xxx_xxyz_0, \
                             ta1_z_xxx_xxzz_0, \
                             ta1_z_xxx_xyyy_0, \
                             ta1_z_xxx_xyyz_0, \
                             ta1_z_xxx_xyzz_0, \
                             ta1_z_xxx_xzzz_0, \
                             ta1_z_xxx_yyyy_0, \
                             ta1_z_xxx_yyyz_0, \
                             ta1_z_xxx_yyzz_0, \
                             ta1_z_xxx_yzzz_0, \
                             ta1_z_xxx_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxx_xxxx_0[i] = 2.0 * ta1_z_x_xxxx_0[i] * fe_0 - 2.0 * ta1_z_x_xxxx_1[i] * fe_0 +
                              4.0 * ta1_z_xx_xxx_0[i] * fe_0 - 4.0 * ta1_z_xx_xxx_1[i] * fe_0 +
                              ta1_z_xx_xxxx_0[i] * pa_x[i] - ta1_z_xx_xxxx_1[i] * pc_x[i];

        ta1_z_xxx_xxxy_0[i] = 2.0 * ta1_z_x_xxxy_0[i] * fe_0 - 2.0 * ta1_z_x_xxxy_1[i] * fe_0 +
                              3.0 * ta1_z_xx_xxy_0[i] * fe_0 - 3.0 * ta1_z_xx_xxy_1[i] * fe_0 +
                              ta1_z_xx_xxxy_0[i] * pa_x[i] - ta1_z_xx_xxxy_1[i] * pc_x[i];

        ta1_z_xxx_xxxz_0[i] = 2.0 * ta1_z_x_xxxz_0[i] * fe_0 - 2.0 * ta1_z_x_xxxz_1[i] * fe_0 +
                              3.0 * ta1_z_xx_xxz_0[i] * fe_0 - 3.0 * ta1_z_xx_xxz_1[i] * fe_0 +
                              ta1_z_xx_xxxz_0[i] * pa_x[i] - ta1_z_xx_xxxz_1[i] * pc_x[i];

        ta1_z_xxx_xxyy_0[i] = 2.0 * ta1_z_x_xxyy_0[i] * fe_0 - 2.0 * ta1_z_x_xxyy_1[i] * fe_0 +
                              2.0 * ta1_z_xx_xyy_0[i] * fe_0 - 2.0 * ta1_z_xx_xyy_1[i] * fe_0 +
                              ta1_z_xx_xxyy_0[i] * pa_x[i] - ta1_z_xx_xxyy_1[i] * pc_x[i];

        ta1_z_xxx_xxyz_0[i] = 2.0 * ta1_z_x_xxyz_0[i] * fe_0 - 2.0 * ta1_z_x_xxyz_1[i] * fe_0 +
                              2.0 * ta1_z_xx_xyz_0[i] * fe_0 - 2.0 * ta1_z_xx_xyz_1[i] * fe_0 +
                              ta1_z_xx_xxyz_0[i] * pa_x[i] - ta1_z_xx_xxyz_1[i] * pc_x[i];

        ta1_z_xxx_xxzz_0[i] = 2.0 * ta1_z_x_xxzz_0[i] * fe_0 - 2.0 * ta1_z_x_xxzz_1[i] * fe_0 +
                              2.0 * ta1_z_xx_xzz_0[i] * fe_0 - 2.0 * ta1_z_xx_xzz_1[i] * fe_0 +
                              ta1_z_xx_xxzz_0[i] * pa_x[i] - ta1_z_xx_xxzz_1[i] * pc_x[i];

        ta1_z_xxx_xyyy_0[i] = 2.0 * ta1_z_x_xyyy_0[i] * fe_0 - 2.0 * ta1_z_x_xyyy_1[i] * fe_0 +
                              ta1_z_xx_yyy_0[i] * fe_0 - ta1_z_xx_yyy_1[i] * fe_0 + ta1_z_xx_xyyy_0[i] * pa_x[i] -
                              ta1_z_xx_xyyy_1[i] * pc_x[i];

        ta1_z_xxx_xyyz_0[i] = 2.0 * ta1_z_x_xyyz_0[i] * fe_0 - 2.0 * ta1_z_x_xyyz_1[i] * fe_0 +
                              ta1_z_xx_yyz_0[i] * fe_0 - ta1_z_xx_yyz_1[i] * fe_0 + ta1_z_xx_xyyz_0[i] * pa_x[i] -
                              ta1_z_xx_xyyz_1[i] * pc_x[i];

        ta1_z_xxx_xyzz_0[i] = 2.0 * ta1_z_x_xyzz_0[i] * fe_0 - 2.0 * ta1_z_x_xyzz_1[i] * fe_0 +
                              ta1_z_xx_yzz_0[i] * fe_0 - ta1_z_xx_yzz_1[i] * fe_0 + ta1_z_xx_xyzz_0[i] * pa_x[i] -
                              ta1_z_xx_xyzz_1[i] * pc_x[i];

        ta1_z_xxx_xzzz_0[i] = 2.0 * ta1_z_x_xzzz_0[i] * fe_0 - 2.0 * ta1_z_x_xzzz_1[i] * fe_0 +
                              ta1_z_xx_zzz_0[i] * fe_0 - ta1_z_xx_zzz_1[i] * fe_0 + ta1_z_xx_xzzz_0[i] * pa_x[i] -
                              ta1_z_xx_xzzz_1[i] * pc_x[i];

        ta1_z_xxx_yyyy_0[i] = 2.0 * ta1_z_x_yyyy_0[i] * fe_0 - 2.0 * ta1_z_x_yyyy_1[i] * fe_0 +
                              ta1_z_xx_yyyy_0[i] * pa_x[i] - ta1_z_xx_yyyy_1[i] * pc_x[i];

        ta1_z_xxx_yyyz_0[i] = 2.0 * ta1_z_x_yyyz_0[i] * fe_0 - 2.0 * ta1_z_x_yyyz_1[i] * fe_0 +
                              ta1_z_xx_yyyz_0[i] * pa_x[i] - ta1_z_xx_yyyz_1[i] * pc_x[i];

        ta1_z_xxx_yyzz_0[i] = 2.0 * ta1_z_x_yyzz_0[i] * fe_0 - 2.0 * ta1_z_x_yyzz_1[i] * fe_0 +
                              ta1_z_xx_yyzz_0[i] * pa_x[i] - ta1_z_xx_yyzz_1[i] * pc_x[i];

        ta1_z_xxx_yzzz_0[i] = 2.0 * ta1_z_x_yzzz_0[i] * fe_0 - 2.0 * ta1_z_x_yzzz_1[i] * fe_0 +
                              ta1_z_xx_yzzz_0[i] * pa_x[i] - ta1_z_xx_yzzz_1[i] * pc_x[i];

        ta1_z_xxx_zzzz_0[i] = 2.0 * ta1_z_x_zzzz_0[i] * fe_0 - 2.0 * ta1_z_x_zzzz_1[i] * fe_0 +
                              ta1_z_xx_zzzz_0[i] * pa_x[i] - ta1_z_xx_zzzz_1[i] * pc_x[i];
    }

    // Set up 315-330 components of targeted buffer : FG

    auto ta1_z_xxy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 315);

    auto ta1_z_xxy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 316);

    auto ta1_z_xxy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 317);

    auto ta1_z_xxy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 318);

    auto ta1_z_xxy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 319);

    auto ta1_z_xxy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 320);

    auto ta1_z_xxy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 321);

    auto ta1_z_xxy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 322);

    auto ta1_z_xxy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 323);

    auto ta1_z_xxy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 324);

    auto ta1_z_xxy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 325);

    auto ta1_z_xxy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 326);

    auto ta1_z_xxy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 327);

    auto ta1_z_xxy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 328);

    auto ta1_z_xxy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 329);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_z_xx_xxx_0,   \
                             ta1_z_xx_xxx_1,   \
                             ta1_z_xx_xxxx_0,  \
                             ta1_z_xx_xxxx_1,  \
                             ta1_z_xx_xxxy_0,  \
                             ta1_z_xx_xxxy_1,  \
                             ta1_z_xx_xxxz_0,  \
                             ta1_z_xx_xxxz_1,  \
                             ta1_z_xx_xxy_0,   \
                             ta1_z_xx_xxy_1,   \
                             ta1_z_xx_xxyy_0,  \
                             ta1_z_xx_xxyy_1,  \
                             ta1_z_xx_xxyz_0,  \
                             ta1_z_xx_xxyz_1,  \
                             ta1_z_xx_xxz_0,   \
                             ta1_z_xx_xxz_1,   \
                             ta1_z_xx_xxzz_0,  \
                             ta1_z_xx_xxzz_1,  \
                             ta1_z_xx_xyy_0,   \
                             ta1_z_xx_xyy_1,   \
                             ta1_z_xx_xyyy_0,  \
                             ta1_z_xx_xyyy_1,  \
                             ta1_z_xx_xyyz_0,  \
                             ta1_z_xx_xyyz_1,  \
                             ta1_z_xx_xyz_0,   \
                             ta1_z_xx_xyz_1,   \
                             ta1_z_xx_xyzz_0,  \
                             ta1_z_xx_xyzz_1,  \
                             ta1_z_xx_xzz_0,   \
                             ta1_z_xx_xzz_1,   \
                             ta1_z_xx_xzzz_0,  \
                             ta1_z_xx_xzzz_1,  \
                             ta1_z_xx_zzzz_0,  \
                             ta1_z_xx_zzzz_1,  \
                             ta1_z_xxy_xxxx_0, \
                             ta1_z_xxy_xxxy_0, \
                             ta1_z_xxy_xxxz_0, \
                             ta1_z_xxy_xxyy_0, \
                             ta1_z_xxy_xxyz_0, \
                             ta1_z_xxy_xxzz_0, \
                             ta1_z_xxy_xyyy_0, \
                             ta1_z_xxy_xyyz_0, \
                             ta1_z_xxy_xyzz_0, \
                             ta1_z_xxy_xzzz_0, \
                             ta1_z_xxy_yyyy_0, \
                             ta1_z_xxy_yyyz_0, \
                             ta1_z_xxy_yyzz_0, \
                             ta1_z_xxy_yzzz_0, \
                             ta1_z_xxy_zzzz_0, \
                             ta1_z_xy_yyyy_0,  \
                             ta1_z_xy_yyyy_1,  \
                             ta1_z_xy_yyyz_0,  \
                             ta1_z_xy_yyyz_1,  \
                             ta1_z_xy_yyzz_0,  \
                             ta1_z_xy_yyzz_1,  \
                             ta1_z_xy_yzzz_0,  \
                             ta1_z_xy_yzzz_1,  \
                             ta1_z_y_yyyy_0,   \
                             ta1_z_y_yyyy_1,   \
                             ta1_z_y_yyyz_0,   \
                             ta1_z_y_yyyz_1,   \
                             ta1_z_y_yyzz_0,   \
                             ta1_z_y_yyzz_1,   \
                             ta1_z_y_yzzz_0,   \
                             ta1_z_y_yzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxy_xxxx_0[i] = ta1_z_xx_xxxx_0[i] * pa_y[i] - ta1_z_xx_xxxx_1[i] * pc_y[i];

        ta1_z_xxy_xxxy_0[i] = ta1_z_xx_xxx_0[i] * fe_0 - ta1_z_xx_xxx_1[i] * fe_0 + ta1_z_xx_xxxy_0[i] * pa_y[i] -
                              ta1_z_xx_xxxy_1[i] * pc_y[i];

        ta1_z_xxy_xxxz_0[i] = ta1_z_xx_xxxz_0[i] * pa_y[i] - ta1_z_xx_xxxz_1[i] * pc_y[i];

        ta1_z_xxy_xxyy_0[i] = 2.0 * ta1_z_xx_xxy_0[i] * fe_0 - 2.0 * ta1_z_xx_xxy_1[i] * fe_0 +
                              ta1_z_xx_xxyy_0[i] * pa_y[i] - ta1_z_xx_xxyy_1[i] * pc_y[i];

        ta1_z_xxy_xxyz_0[i] = ta1_z_xx_xxz_0[i] * fe_0 - ta1_z_xx_xxz_1[i] * fe_0 + ta1_z_xx_xxyz_0[i] * pa_y[i] -
                              ta1_z_xx_xxyz_1[i] * pc_y[i];

        ta1_z_xxy_xxzz_0[i] = ta1_z_xx_xxzz_0[i] * pa_y[i] - ta1_z_xx_xxzz_1[i] * pc_y[i];

        ta1_z_xxy_xyyy_0[i] = 3.0 * ta1_z_xx_xyy_0[i] * fe_0 - 3.0 * ta1_z_xx_xyy_1[i] * fe_0 +
                              ta1_z_xx_xyyy_0[i] * pa_y[i] - ta1_z_xx_xyyy_1[i] * pc_y[i];

        ta1_z_xxy_xyyz_0[i] = 2.0 * ta1_z_xx_xyz_0[i] * fe_0 - 2.0 * ta1_z_xx_xyz_1[i] * fe_0 +
                              ta1_z_xx_xyyz_0[i] * pa_y[i] - ta1_z_xx_xyyz_1[i] * pc_y[i];

        ta1_z_xxy_xyzz_0[i] = ta1_z_xx_xzz_0[i] * fe_0 - ta1_z_xx_xzz_1[i] * fe_0 + ta1_z_xx_xyzz_0[i] * pa_y[i] -
                              ta1_z_xx_xyzz_1[i] * pc_y[i];

        ta1_z_xxy_xzzz_0[i] = ta1_z_xx_xzzz_0[i] * pa_y[i] - ta1_z_xx_xzzz_1[i] * pc_y[i];

        ta1_z_xxy_yyyy_0[i] = ta1_z_y_yyyy_0[i] * fe_0 - ta1_z_y_yyyy_1[i] * fe_0 + ta1_z_xy_yyyy_0[i] * pa_x[i] -
                              ta1_z_xy_yyyy_1[i] * pc_x[i];

        ta1_z_xxy_yyyz_0[i] = ta1_z_y_yyyz_0[i] * fe_0 - ta1_z_y_yyyz_1[i] * fe_0 + ta1_z_xy_yyyz_0[i] * pa_x[i] -
                              ta1_z_xy_yyyz_1[i] * pc_x[i];

        ta1_z_xxy_yyzz_0[i] = ta1_z_y_yyzz_0[i] * fe_0 - ta1_z_y_yyzz_1[i] * fe_0 + ta1_z_xy_yyzz_0[i] * pa_x[i] -
                              ta1_z_xy_yyzz_1[i] * pc_x[i];

        ta1_z_xxy_yzzz_0[i] = ta1_z_y_yzzz_0[i] * fe_0 - ta1_z_y_yzzz_1[i] * fe_0 + ta1_z_xy_yzzz_0[i] * pa_x[i] -
                              ta1_z_xy_yzzz_1[i] * pc_x[i];

        ta1_z_xxy_zzzz_0[i] = ta1_z_xx_zzzz_0[i] * pa_y[i] - ta1_z_xx_zzzz_1[i] * pc_y[i];
    }

    // Set up 330-345 components of targeted buffer : FG

    auto ta1_z_xxz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 330);

    auto ta1_z_xxz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 331);

    auto ta1_z_xxz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 332);

    auto ta1_z_xxz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 333);

    auto ta1_z_xxz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 334);

    auto ta1_z_xxz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 335);

    auto ta1_z_xxz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 336);

    auto ta1_z_xxz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 337);

    auto ta1_z_xxz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 338);

    auto ta1_z_xxz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 339);

    auto ta1_z_xxz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 340);

    auto ta1_z_xxz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 341);

    auto ta1_z_xxz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 342);

    auto ta1_z_xxz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 343);

    auto ta1_z_xxz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 344);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_z_xx_xxx_0,   \
                             ta1_z_xx_xxx_1,   \
                             ta1_z_xx_xxxx_0,  \
                             ta1_z_xx_xxxx_1,  \
                             ta1_z_xx_xxxy_0,  \
                             ta1_z_xx_xxxy_1,  \
                             ta1_z_xx_xxxz_0,  \
                             ta1_z_xx_xxxz_1,  \
                             ta1_z_xx_xxy_0,   \
                             ta1_z_xx_xxy_1,   \
                             ta1_z_xx_xxyy_0,  \
                             ta1_z_xx_xxyy_1,  \
                             ta1_z_xx_xxyz_0,  \
                             ta1_z_xx_xxyz_1,  \
                             ta1_z_xx_xxz_0,   \
                             ta1_z_xx_xxz_1,   \
                             ta1_z_xx_xxzz_0,  \
                             ta1_z_xx_xxzz_1,  \
                             ta1_z_xx_xyy_0,   \
                             ta1_z_xx_xyy_1,   \
                             ta1_z_xx_xyyy_0,  \
                             ta1_z_xx_xyyy_1,  \
                             ta1_z_xx_xyyz_0,  \
                             ta1_z_xx_xyyz_1,  \
                             ta1_z_xx_xyz_0,   \
                             ta1_z_xx_xyz_1,   \
                             ta1_z_xx_xyzz_0,  \
                             ta1_z_xx_xyzz_1,  \
                             ta1_z_xx_xzz_0,   \
                             ta1_z_xx_xzz_1,   \
                             ta1_z_xx_xzzz_0,  \
                             ta1_z_xx_xzzz_1,  \
                             ta1_z_xx_yyyy_0,  \
                             ta1_z_xx_yyyy_1,  \
                             ta1_z_xxz_xxxx_0, \
                             ta1_z_xxz_xxxy_0, \
                             ta1_z_xxz_xxxz_0, \
                             ta1_z_xxz_xxyy_0, \
                             ta1_z_xxz_xxyz_0, \
                             ta1_z_xxz_xxzz_0, \
                             ta1_z_xxz_xyyy_0, \
                             ta1_z_xxz_xyyz_0, \
                             ta1_z_xxz_xyzz_0, \
                             ta1_z_xxz_xzzz_0, \
                             ta1_z_xxz_yyyy_0, \
                             ta1_z_xxz_yyyz_0, \
                             ta1_z_xxz_yyzz_0, \
                             ta1_z_xxz_yzzz_0, \
                             ta1_z_xxz_zzzz_0, \
                             ta1_z_xz_yyyz_0,  \
                             ta1_z_xz_yyyz_1,  \
                             ta1_z_xz_yyzz_0,  \
                             ta1_z_xz_yyzz_1,  \
                             ta1_z_xz_yzzz_0,  \
                             ta1_z_xz_yzzz_1,  \
                             ta1_z_xz_zzzz_0,  \
                             ta1_z_xz_zzzz_1,  \
                             ta1_z_z_yyyz_0,   \
                             ta1_z_z_yyyz_1,   \
                             ta1_z_z_yyzz_0,   \
                             ta1_z_z_yyzz_1,   \
                             ta1_z_z_yzzz_0,   \
                             ta1_z_z_yzzz_1,   \
                             ta1_z_z_zzzz_0,   \
                             ta1_z_z_zzzz_1,   \
                             ta_xx_xxxx_1,     \
                             ta_xx_xxxy_1,     \
                             ta_xx_xxxz_1,     \
                             ta_xx_xxyy_1,     \
                             ta_xx_xxyz_1,     \
                             ta_xx_xxzz_1,     \
                             ta_xx_xyyy_1,     \
                             ta_xx_xyyz_1,     \
                             ta_xx_xyzz_1,     \
                             ta_xx_xzzz_1,     \
                             ta_xx_yyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxz_xxxx_0[i] = ta_xx_xxxx_1[i] + ta1_z_xx_xxxx_0[i] * pa_z[i] - ta1_z_xx_xxxx_1[i] * pc_z[i];

        ta1_z_xxz_xxxy_0[i] = ta_xx_xxxy_1[i] + ta1_z_xx_xxxy_0[i] * pa_z[i] - ta1_z_xx_xxxy_1[i] * pc_z[i];

        ta1_z_xxz_xxxz_0[i] = ta1_z_xx_xxx_0[i] * fe_0 - ta1_z_xx_xxx_1[i] * fe_0 + ta_xx_xxxz_1[i] +
                              ta1_z_xx_xxxz_0[i] * pa_z[i] - ta1_z_xx_xxxz_1[i] * pc_z[i];

        ta1_z_xxz_xxyy_0[i] = ta_xx_xxyy_1[i] + ta1_z_xx_xxyy_0[i] * pa_z[i] - ta1_z_xx_xxyy_1[i] * pc_z[i];

        ta1_z_xxz_xxyz_0[i] = ta1_z_xx_xxy_0[i] * fe_0 - ta1_z_xx_xxy_1[i] * fe_0 + ta_xx_xxyz_1[i] +
                              ta1_z_xx_xxyz_0[i] * pa_z[i] - ta1_z_xx_xxyz_1[i] * pc_z[i];

        ta1_z_xxz_xxzz_0[i] = 2.0 * ta1_z_xx_xxz_0[i] * fe_0 - 2.0 * ta1_z_xx_xxz_1[i] * fe_0 + ta_xx_xxzz_1[i] +
                              ta1_z_xx_xxzz_0[i] * pa_z[i] - ta1_z_xx_xxzz_1[i] * pc_z[i];

        ta1_z_xxz_xyyy_0[i] = ta_xx_xyyy_1[i] + ta1_z_xx_xyyy_0[i] * pa_z[i] - ta1_z_xx_xyyy_1[i] * pc_z[i];

        ta1_z_xxz_xyyz_0[i] = ta1_z_xx_xyy_0[i] * fe_0 - ta1_z_xx_xyy_1[i] * fe_0 + ta_xx_xyyz_1[i] +
                              ta1_z_xx_xyyz_0[i] * pa_z[i] - ta1_z_xx_xyyz_1[i] * pc_z[i];

        ta1_z_xxz_xyzz_0[i] = 2.0 * ta1_z_xx_xyz_0[i] * fe_0 - 2.0 * ta1_z_xx_xyz_1[i] * fe_0 + ta_xx_xyzz_1[i] +
                              ta1_z_xx_xyzz_0[i] * pa_z[i] - ta1_z_xx_xyzz_1[i] * pc_z[i];

        ta1_z_xxz_xzzz_0[i] = 3.0 * ta1_z_xx_xzz_0[i] * fe_0 - 3.0 * ta1_z_xx_xzz_1[i] * fe_0 + ta_xx_xzzz_1[i] +
                              ta1_z_xx_xzzz_0[i] * pa_z[i] - ta1_z_xx_xzzz_1[i] * pc_z[i];

        ta1_z_xxz_yyyy_0[i] = ta_xx_yyyy_1[i] + ta1_z_xx_yyyy_0[i] * pa_z[i] - ta1_z_xx_yyyy_1[i] * pc_z[i];

        ta1_z_xxz_yyyz_0[i] = ta1_z_z_yyyz_0[i] * fe_0 - ta1_z_z_yyyz_1[i] * fe_0 + ta1_z_xz_yyyz_0[i] * pa_x[i] -
                              ta1_z_xz_yyyz_1[i] * pc_x[i];

        ta1_z_xxz_yyzz_0[i] = ta1_z_z_yyzz_0[i] * fe_0 - ta1_z_z_yyzz_1[i] * fe_0 + ta1_z_xz_yyzz_0[i] * pa_x[i] -
                              ta1_z_xz_yyzz_1[i] * pc_x[i];

        ta1_z_xxz_yzzz_0[i] = ta1_z_z_yzzz_0[i] * fe_0 - ta1_z_z_yzzz_1[i] * fe_0 + ta1_z_xz_yzzz_0[i] * pa_x[i] -
                              ta1_z_xz_yzzz_1[i] * pc_x[i];

        ta1_z_xxz_zzzz_0[i] = ta1_z_z_zzzz_0[i] * fe_0 - ta1_z_z_zzzz_1[i] * fe_0 + ta1_z_xz_zzzz_0[i] * pa_x[i] -
                              ta1_z_xz_zzzz_1[i] * pc_x[i];
    }

    // Set up 345-360 components of targeted buffer : FG

    auto ta1_z_xyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 345);

    auto ta1_z_xyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 346);

    auto ta1_z_xyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 347);

    auto ta1_z_xyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 348);

    auto ta1_z_xyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 349);

    auto ta1_z_xyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 350);

    auto ta1_z_xyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 351);

    auto ta1_z_xyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 352);

    auto ta1_z_xyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 353);

    auto ta1_z_xyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 354);

    auto ta1_z_xyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 355);

    auto ta1_z_xyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 356);

    auto ta1_z_xyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 357);

    auto ta1_z_xyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 358);

    auto ta1_z_xyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 359);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_z_xyy_xxxx_0, \
                             ta1_z_xyy_xxxy_0, \
                             ta1_z_xyy_xxxz_0, \
                             ta1_z_xyy_xxyy_0, \
                             ta1_z_xyy_xxyz_0, \
                             ta1_z_xyy_xxzz_0, \
                             ta1_z_xyy_xyyy_0, \
                             ta1_z_xyy_xyyz_0, \
                             ta1_z_xyy_xyzz_0, \
                             ta1_z_xyy_xzzz_0, \
                             ta1_z_xyy_yyyy_0, \
                             ta1_z_xyy_yyyz_0, \
                             ta1_z_xyy_yyzz_0, \
                             ta1_z_xyy_yzzz_0, \
                             ta1_z_xyy_zzzz_0, \
                             ta1_z_yy_xxx_0,   \
                             ta1_z_yy_xxx_1,   \
                             ta1_z_yy_xxxx_0,  \
                             ta1_z_yy_xxxx_1,  \
                             ta1_z_yy_xxxy_0,  \
                             ta1_z_yy_xxxy_1,  \
                             ta1_z_yy_xxxz_0,  \
                             ta1_z_yy_xxxz_1,  \
                             ta1_z_yy_xxy_0,   \
                             ta1_z_yy_xxy_1,   \
                             ta1_z_yy_xxyy_0,  \
                             ta1_z_yy_xxyy_1,  \
                             ta1_z_yy_xxyz_0,  \
                             ta1_z_yy_xxyz_1,  \
                             ta1_z_yy_xxz_0,   \
                             ta1_z_yy_xxz_1,   \
                             ta1_z_yy_xxzz_0,  \
                             ta1_z_yy_xxzz_1,  \
                             ta1_z_yy_xyy_0,   \
                             ta1_z_yy_xyy_1,   \
                             ta1_z_yy_xyyy_0,  \
                             ta1_z_yy_xyyy_1,  \
                             ta1_z_yy_xyyz_0,  \
                             ta1_z_yy_xyyz_1,  \
                             ta1_z_yy_xyz_0,   \
                             ta1_z_yy_xyz_1,   \
                             ta1_z_yy_xyzz_0,  \
                             ta1_z_yy_xyzz_1,  \
                             ta1_z_yy_xzz_0,   \
                             ta1_z_yy_xzz_1,   \
                             ta1_z_yy_xzzz_0,  \
                             ta1_z_yy_xzzz_1,  \
                             ta1_z_yy_yyy_0,   \
                             ta1_z_yy_yyy_1,   \
                             ta1_z_yy_yyyy_0,  \
                             ta1_z_yy_yyyy_1,  \
                             ta1_z_yy_yyyz_0,  \
                             ta1_z_yy_yyyz_1,  \
                             ta1_z_yy_yyz_0,   \
                             ta1_z_yy_yyz_1,   \
                             ta1_z_yy_yyzz_0,  \
                             ta1_z_yy_yyzz_1,  \
                             ta1_z_yy_yzz_0,   \
                             ta1_z_yy_yzz_1,   \
                             ta1_z_yy_yzzz_0,  \
                             ta1_z_yy_yzzz_1,  \
                             ta1_z_yy_zzz_0,   \
                             ta1_z_yy_zzz_1,   \
                             ta1_z_yy_zzzz_0,  \
                             ta1_z_yy_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyy_xxxx_0[i] = 4.0 * ta1_z_yy_xxx_0[i] * fe_0 - 4.0 * ta1_z_yy_xxx_1[i] * fe_0 +
                              ta1_z_yy_xxxx_0[i] * pa_x[i] - ta1_z_yy_xxxx_1[i] * pc_x[i];

        ta1_z_xyy_xxxy_0[i] = 3.0 * ta1_z_yy_xxy_0[i] * fe_0 - 3.0 * ta1_z_yy_xxy_1[i] * fe_0 +
                              ta1_z_yy_xxxy_0[i] * pa_x[i] - ta1_z_yy_xxxy_1[i] * pc_x[i];

        ta1_z_xyy_xxxz_0[i] = 3.0 * ta1_z_yy_xxz_0[i] * fe_0 - 3.0 * ta1_z_yy_xxz_1[i] * fe_0 +
                              ta1_z_yy_xxxz_0[i] * pa_x[i] - ta1_z_yy_xxxz_1[i] * pc_x[i];

        ta1_z_xyy_xxyy_0[i] = 2.0 * ta1_z_yy_xyy_0[i] * fe_0 - 2.0 * ta1_z_yy_xyy_1[i] * fe_0 +
                              ta1_z_yy_xxyy_0[i] * pa_x[i] - ta1_z_yy_xxyy_1[i] * pc_x[i];

        ta1_z_xyy_xxyz_0[i] = 2.0 * ta1_z_yy_xyz_0[i] * fe_0 - 2.0 * ta1_z_yy_xyz_1[i] * fe_0 +
                              ta1_z_yy_xxyz_0[i] * pa_x[i] - ta1_z_yy_xxyz_1[i] * pc_x[i];

        ta1_z_xyy_xxzz_0[i] = 2.0 * ta1_z_yy_xzz_0[i] * fe_0 - 2.0 * ta1_z_yy_xzz_1[i] * fe_0 +
                              ta1_z_yy_xxzz_0[i] * pa_x[i] - ta1_z_yy_xxzz_1[i] * pc_x[i];

        ta1_z_xyy_xyyy_0[i] = ta1_z_yy_yyy_0[i] * fe_0 - ta1_z_yy_yyy_1[i] * fe_0 + ta1_z_yy_xyyy_0[i] * pa_x[i] -
                              ta1_z_yy_xyyy_1[i] * pc_x[i];

        ta1_z_xyy_xyyz_0[i] = ta1_z_yy_yyz_0[i] * fe_0 - ta1_z_yy_yyz_1[i] * fe_0 + ta1_z_yy_xyyz_0[i] * pa_x[i] -
                              ta1_z_yy_xyyz_1[i] * pc_x[i];

        ta1_z_xyy_xyzz_0[i] = ta1_z_yy_yzz_0[i] * fe_0 - ta1_z_yy_yzz_1[i] * fe_0 + ta1_z_yy_xyzz_0[i] * pa_x[i] -
                              ta1_z_yy_xyzz_1[i] * pc_x[i];

        ta1_z_xyy_xzzz_0[i] = ta1_z_yy_zzz_0[i] * fe_0 - ta1_z_yy_zzz_1[i] * fe_0 + ta1_z_yy_xzzz_0[i] * pa_x[i] -
                              ta1_z_yy_xzzz_1[i] * pc_x[i];

        ta1_z_xyy_yyyy_0[i] = ta1_z_yy_yyyy_0[i] * pa_x[i] - ta1_z_yy_yyyy_1[i] * pc_x[i];

        ta1_z_xyy_yyyz_0[i] = ta1_z_yy_yyyz_0[i] * pa_x[i] - ta1_z_yy_yyyz_1[i] * pc_x[i];

        ta1_z_xyy_yyzz_0[i] = ta1_z_yy_yyzz_0[i] * pa_x[i] - ta1_z_yy_yyzz_1[i] * pc_x[i];

        ta1_z_xyy_yzzz_0[i] = ta1_z_yy_yzzz_0[i] * pa_x[i] - ta1_z_yy_yzzz_1[i] * pc_x[i];

        ta1_z_xyy_zzzz_0[i] = ta1_z_yy_zzzz_0[i] * pa_x[i] - ta1_z_yy_zzzz_1[i] * pc_x[i];
    }

    // Set up 360-375 components of targeted buffer : FG

    auto ta1_z_xyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 360);

    auto ta1_z_xyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 361);

    auto ta1_z_xyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 362);

    auto ta1_z_xyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 363);

    auto ta1_z_xyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 364);

    auto ta1_z_xyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 365);

    auto ta1_z_xyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 366);

    auto ta1_z_xyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 367);

    auto ta1_z_xyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 368);

    auto ta1_z_xyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 369);

    auto ta1_z_xyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 370);

    auto ta1_z_xyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 371);

    auto ta1_z_xyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 372);

    auto ta1_z_xyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 373);

    auto ta1_z_xyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 374);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             pc_x,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_z_xy_xxxy_0,  \
                             ta1_z_xy_xxxy_1,  \
                             ta1_z_xy_xxyy_0,  \
                             ta1_z_xy_xxyy_1,  \
                             ta1_z_xy_xyyy_0,  \
                             ta1_z_xy_xyyy_1,  \
                             ta1_z_xyz_xxxx_0, \
                             ta1_z_xyz_xxxy_0, \
                             ta1_z_xyz_xxxz_0, \
                             ta1_z_xyz_xxyy_0, \
                             ta1_z_xyz_xxyz_0, \
                             ta1_z_xyz_xxzz_0, \
                             ta1_z_xyz_xyyy_0, \
                             ta1_z_xyz_xyyz_0, \
                             ta1_z_xyz_xyzz_0, \
                             ta1_z_xyz_xzzz_0, \
                             ta1_z_xyz_yyyy_0, \
                             ta1_z_xyz_yyyz_0, \
                             ta1_z_xyz_yyzz_0, \
                             ta1_z_xyz_yzzz_0, \
                             ta1_z_xyz_zzzz_0, \
                             ta1_z_xz_xxxx_0,  \
                             ta1_z_xz_xxxx_1,  \
                             ta1_z_xz_xxxz_0,  \
                             ta1_z_xz_xxxz_1,  \
                             ta1_z_xz_xxzz_0,  \
                             ta1_z_xz_xxzz_1,  \
                             ta1_z_xz_xzzz_0,  \
                             ta1_z_xz_xzzz_1,  \
                             ta1_z_yz_xxyz_0,  \
                             ta1_z_yz_xxyz_1,  \
                             ta1_z_yz_xyyz_0,  \
                             ta1_z_yz_xyyz_1,  \
                             ta1_z_yz_xyz_0,   \
                             ta1_z_yz_xyz_1,   \
                             ta1_z_yz_xyzz_0,  \
                             ta1_z_yz_xyzz_1,  \
                             ta1_z_yz_yyyy_0,  \
                             ta1_z_yz_yyyy_1,  \
                             ta1_z_yz_yyyz_0,  \
                             ta1_z_yz_yyyz_1,  \
                             ta1_z_yz_yyz_0,   \
                             ta1_z_yz_yyz_1,   \
                             ta1_z_yz_yyzz_0,  \
                             ta1_z_yz_yyzz_1,  \
                             ta1_z_yz_yzz_0,   \
                             ta1_z_yz_yzz_1,   \
                             ta1_z_yz_yzzz_0,  \
                             ta1_z_yz_yzzz_1,  \
                             ta1_z_yz_zzzz_0,  \
                             ta1_z_yz_zzzz_1,  \
                             ta_xy_xxxy_1,     \
                             ta_xy_xxyy_1,     \
                             ta_xy_xyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyz_xxxx_0[i] = ta1_z_xz_xxxx_0[i] * pa_y[i] - ta1_z_xz_xxxx_1[i] * pc_y[i];

        ta1_z_xyz_xxxy_0[i] = ta_xy_xxxy_1[i] + ta1_z_xy_xxxy_0[i] * pa_z[i] - ta1_z_xy_xxxy_1[i] * pc_z[i];

        ta1_z_xyz_xxxz_0[i] = ta1_z_xz_xxxz_0[i] * pa_y[i] - ta1_z_xz_xxxz_1[i] * pc_y[i];

        ta1_z_xyz_xxyy_0[i] = ta_xy_xxyy_1[i] + ta1_z_xy_xxyy_0[i] * pa_z[i] - ta1_z_xy_xxyy_1[i] * pc_z[i];

        ta1_z_xyz_xxyz_0[i] = 2.0 * ta1_z_yz_xyz_0[i] * fe_0 - 2.0 * ta1_z_yz_xyz_1[i] * fe_0 +
                              ta1_z_yz_xxyz_0[i] * pa_x[i] - ta1_z_yz_xxyz_1[i] * pc_x[i];

        ta1_z_xyz_xxzz_0[i] = ta1_z_xz_xxzz_0[i] * pa_y[i] - ta1_z_xz_xxzz_1[i] * pc_y[i];

        ta1_z_xyz_xyyy_0[i] = ta_xy_xyyy_1[i] + ta1_z_xy_xyyy_0[i] * pa_z[i] - ta1_z_xy_xyyy_1[i] * pc_z[i];

        ta1_z_xyz_xyyz_0[i] = ta1_z_yz_yyz_0[i] * fe_0 - ta1_z_yz_yyz_1[i] * fe_0 + ta1_z_yz_xyyz_0[i] * pa_x[i] -
                              ta1_z_yz_xyyz_1[i] * pc_x[i];

        ta1_z_xyz_xyzz_0[i] = ta1_z_yz_yzz_0[i] * fe_0 - ta1_z_yz_yzz_1[i] * fe_0 + ta1_z_yz_xyzz_0[i] * pa_x[i] -
                              ta1_z_yz_xyzz_1[i] * pc_x[i];

        ta1_z_xyz_xzzz_0[i] = ta1_z_xz_xzzz_0[i] * pa_y[i] - ta1_z_xz_xzzz_1[i] * pc_y[i];

        ta1_z_xyz_yyyy_0[i] = ta1_z_yz_yyyy_0[i] * pa_x[i] - ta1_z_yz_yyyy_1[i] * pc_x[i];

        ta1_z_xyz_yyyz_0[i] = ta1_z_yz_yyyz_0[i] * pa_x[i] - ta1_z_yz_yyyz_1[i] * pc_x[i];

        ta1_z_xyz_yyzz_0[i] = ta1_z_yz_yyzz_0[i] * pa_x[i] - ta1_z_yz_yyzz_1[i] * pc_x[i];

        ta1_z_xyz_yzzz_0[i] = ta1_z_yz_yzzz_0[i] * pa_x[i] - ta1_z_yz_yzzz_1[i] * pc_x[i];

        ta1_z_xyz_zzzz_0[i] = ta1_z_yz_zzzz_0[i] * pa_x[i] - ta1_z_yz_zzzz_1[i] * pc_x[i];
    }

    // Set up 375-390 components of targeted buffer : FG

    auto ta1_z_xzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 375);

    auto ta1_z_xzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 376);

    auto ta1_z_xzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 377);

    auto ta1_z_xzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 378);

    auto ta1_z_xzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 379);

    auto ta1_z_xzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 380);

    auto ta1_z_xzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 381);

    auto ta1_z_xzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 382);

    auto ta1_z_xzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 383);

    auto ta1_z_xzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 384);

    auto ta1_z_xzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 385);

    auto ta1_z_xzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 386);

    auto ta1_z_xzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 387);

    auto ta1_z_xzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 388);

    auto ta1_z_xzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 389);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_z_xzz_xxxx_0, \
                             ta1_z_xzz_xxxy_0, \
                             ta1_z_xzz_xxxz_0, \
                             ta1_z_xzz_xxyy_0, \
                             ta1_z_xzz_xxyz_0, \
                             ta1_z_xzz_xxzz_0, \
                             ta1_z_xzz_xyyy_0, \
                             ta1_z_xzz_xyyz_0, \
                             ta1_z_xzz_xyzz_0, \
                             ta1_z_xzz_xzzz_0, \
                             ta1_z_xzz_yyyy_0, \
                             ta1_z_xzz_yyyz_0, \
                             ta1_z_xzz_yyzz_0, \
                             ta1_z_xzz_yzzz_0, \
                             ta1_z_xzz_zzzz_0, \
                             ta1_z_zz_xxx_0,   \
                             ta1_z_zz_xxx_1,   \
                             ta1_z_zz_xxxx_0,  \
                             ta1_z_zz_xxxx_1,  \
                             ta1_z_zz_xxxy_0,  \
                             ta1_z_zz_xxxy_1,  \
                             ta1_z_zz_xxxz_0,  \
                             ta1_z_zz_xxxz_1,  \
                             ta1_z_zz_xxy_0,   \
                             ta1_z_zz_xxy_1,   \
                             ta1_z_zz_xxyy_0,  \
                             ta1_z_zz_xxyy_1,  \
                             ta1_z_zz_xxyz_0,  \
                             ta1_z_zz_xxyz_1,  \
                             ta1_z_zz_xxz_0,   \
                             ta1_z_zz_xxz_1,   \
                             ta1_z_zz_xxzz_0,  \
                             ta1_z_zz_xxzz_1,  \
                             ta1_z_zz_xyy_0,   \
                             ta1_z_zz_xyy_1,   \
                             ta1_z_zz_xyyy_0,  \
                             ta1_z_zz_xyyy_1,  \
                             ta1_z_zz_xyyz_0,  \
                             ta1_z_zz_xyyz_1,  \
                             ta1_z_zz_xyz_0,   \
                             ta1_z_zz_xyz_1,   \
                             ta1_z_zz_xyzz_0,  \
                             ta1_z_zz_xyzz_1,  \
                             ta1_z_zz_xzz_0,   \
                             ta1_z_zz_xzz_1,   \
                             ta1_z_zz_xzzz_0,  \
                             ta1_z_zz_xzzz_1,  \
                             ta1_z_zz_yyy_0,   \
                             ta1_z_zz_yyy_1,   \
                             ta1_z_zz_yyyy_0,  \
                             ta1_z_zz_yyyy_1,  \
                             ta1_z_zz_yyyz_0,  \
                             ta1_z_zz_yyyz_1,  \
                             ta1_z_zz_yyz_0,   \
                             ta1_z_zz_yyz_1,   \
                             ta1_z_zz_yyzz_0,  \
                             ta1_z_zz_yyzz_1,  \
                             ta1_z_zz_yzz_0,   \
                             ta1_z_zz_yzz_1,   \
                             ta1_z_zz_yzzz_0,  \
                             ta1_z_zz_yzzz_1,  \
                             ta1_z_zz_zzz_0,   \
                             ta1_z_zz_zzz_1,   \
                             ta1_z_zz_zzzz_0,  \
                             ta1_z_zz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xzz_xxxx_0[i] = 4.0 * ta1_z_zz_xxx_0[i] * fe_0 - 4.0 * ta1_z_zz_xxx_1[i] * fe_0 +
                              ta1_z_zz_xxxx_0[i] * pa_x[i] - ta1_z_zz_xxxx_1[i] * pc_x[i];

        ta1_z_xzz_xxxy_0[i] = 3.0 * ta1_z_zz_xxy_0[i] * fe_0 - 3.0 * ta1_z_zz_xxy_1[i] * fe_0 +
                              ta1_z_zz_xxxy_0[i] * pa_x[i] - ta1_z_zz_xxxy_1[i] * pc_x[i];

        ta1_z_xzz_xxxz_0[i] = 3.0 * ta1_z_zz_xxz_0[i] * fe_0 - 3.0 * ta1_z_zz_xxz_1[i] * fe_0 +
                              ta1_z_zz_xxxz_0[i] * pa_x[i] - ta1_z_zz_xxxz_1[i] * pc_x[i];

        ta1_z_xzz_xxyy_0[i] = 2.0 * ta1_z_zz_xyy_0[i] * fe_0 - 2.0 * ta1_z_zz_xyy_1[i] * fe_0 +
                              ta1_z_zz_xxyy_0[i] * pa_x[i] - ta1_z_zz_xxyy_1[i] * pc_x[i];

        ta1_z_xzz_xxyz_0[i] = 2.0 * ta1_z_zz_xyz_0[i] * fe_0 - 2.0 * ta1_z_zz_xyz_1[i] * fe_0 +
                              ta1_z_zz_xxyz_0[i] * pa_x[i] - ta1_z_zz_xxyz_1[i] * pc_x[i];

        ta1_z_xzz_xxzz_0[i] = 2.0 * ta1_z_zz_xzz_0[i] * fe_0 - 2.0 * ta1_z_zz_xzz_1[i] * fe_0 +
                              ta1_z_zz_xxzz_0[i] * pa_x[i] - ta1_z_zz_xxzz_1[i] * pc_x[i];

        ta1_z_xzz_xyyy_0[i] = ta1_z_zz_yyy_0[i] * fe_0 - ta1_z_zz_yyy_1[i] * fe_0 + ta1_z_zz_xyyy_0[i] * pa_x[i] -
                              ta1_z_zz_xyyy_1[i] * pc_x[i];

        ta1_z_xzz_xyyz_0[i] = ta1_z_zz_yyz_0[i] * fe_0 - ta1_z_zz_yyz_1[i] * fe_0 + ta1_z_zz_xyyz_0[i] * pa_x[i] -
                              ta1_z_zz_xyyz_1[i] * pc_x[i];

        ta1_z_xzz_xyzz_0[i] = ta1_z_zz_yzz_0[i] * fe_0 - ta1_z_zz_yzz_1[i] * fe_0 + ta1_z_zz_xyzz_0[i] * pa_x[i] -
                              ta1_z_zz_xyzz_1[i] * pc_x[i];

        ta1_z_xzz_xzzz_0[i] = ta1_z_zz_zzz_0[i] * fe_0 - ta1_z_zz_zzz_1[i] * fe_0 + ta1_z_zz_xzzz_0[i] * pa_x[i] -
                              ta1_z_zz_xzzz_1[i] * pc_x[i];

        ta1_z_xzz_yyyy_0[i] = ta1_z_zz_yyyy_0[i] * pa_x[i] - ta1_z_zz_yyyy_1[i] * pc_x[i];

        ta1_z_xzz_yyyz_0[i] = ta1_z_zz_yyyz_0[i] * pa_x[i] - ta1_z_zz_yyyz_1[i] * pc_x[i];

        ta1_z_xzz_yyzz_0[i] = ta1_z_zz_yyzz_0[i] * pa_x[i] - ta1_z_zz_yyzz_1[i] * pc_x[i];

        ta1_z_xzz_yzzz_0[i] = ta1_z_zz_yzzz_0[i] * pa_x[i] - ta1_z_zz_yzzz_1[i] * pc_x[i];

        ta1_z_xzz_zzzz_0[i] = ta1_z_zz_zzzz_0[i] * pa_x[i] - ta1_z_zz_zzzz_1[i] * pc_x[i];
    }

    // Set up 390-405 components of targeted buffer : FG

    auto ta1_z_yyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 390);

    auto ta1_z_yyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 391);

    auto ta1_z_yyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 392);

    auto ta1_z_yyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 393);

    auto ta1_z_yyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 394);

    auto ta1_z_yyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 395);

    auto ta1_z_yyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 396);

    auto ta1_z_yyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 397);

    auto ta1_z_yyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 398);

    auto ta1_z_yyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 399);

    auto ta1_z_yyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 400);

    auto ta1_z_yyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 401);

    auto ta1_z_yyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 402);

    auto ta1_z_yyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 403);

    auto ta1_z_yyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 404);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_z_y_xxxx_0,   \
                             ta1_z_y_xxxx_1,   \
                             ta1_z_y_xxxy_0,   \
                             ta1_z_y_xxxy_1,   \
                             ta1_z_y_xxxz_0,   \
                             ta1_z_y_xxxz_1,   \
                             ta1_z_y_xxyy_0,   \
                             ta1_z_y_xxyy_1,   \
                             ta1_z_y_xxyz_0,   \
                             ta1_z_y_xxyz_1,   \
                             ta1_z_y_xxzz_0,   \
                             ta1_z_y_xxzz_1,   \
                             ta1_z_y_xyyy_0,   \
                             ta1_z_y_xyyy_1,   \
                             ta1_z_y_xyyz_0,   \
                             ta1_z_y_xyyz_1,   \
                             ta1_z_y_xyzz_0,   \
                             ta1_z_y_xyzz_1,   \
                             ta1_z_y_xzzz_0,   \
                             ta1_z_y_xzzz_1,   \
                             ta1_z_y_yyyy_0,   \
                             ta1_z_y_yyyy_1,   \
                             ta1_z_y_yyyz_0,   \
                             ta1_z_y_yyyz_1,   \
                             ta1_z_y_yyzz_0,   \
                             ta1_z_y_yyzz_1,   \
                             ta1_z_y_yzzz_0,   \
                             ta1_z_y_yzzz_1,   \
                             ta1_z_y_zzzz_0,   \
                             ta1_z_y_zzzz_1,   \
                             ta1_z_yy_xxx_0,   \
                             ta1_z_yy_xxx_1,   \
                             ta1_z_yy_xxxx_0,  \
                             ta1_z_yy_xxxx_1,  \
                             ta1_z_yy_xxxy_0,  \
                             ta1_z_yy_xxxy_1,  \
                             ta1_z_yy_xxxz_0,  \
                             ta1_z_yy_xxxz_1,  \
                             ta1_z_yy_xxy_0,   \
                             ta1_z_yy_xxy_1,   \
                             ta1_z_yy_xxyy_0,  \
                             ta1_z_yy_xxyy_1,  \
                             ta1_z_yy_xxyz_0,  \
                             ta1_z_yy_xxyz_1,  \
                             ta1_z_yy_xxz_0,   \
                             ta1_z_yy_xxz_1,   \
                             ta1_z_yy_xxzz_0,  \
                             ta1_z_yy_xxzz_1,  \
                             ta1_z_yy_xyy_0,   \
                             ta1_z_yy_xyy_1,   \
                             ta1_z_yy_xyyy_0,  \
                             ta1_z_yy_xyyy_1,  \
                             ta1_z_yy_xyyz_0,  \
                             ta1_z_yy_xyyz_1,  \
                             ta1_z_yy_xyz_0,   \
                             ta1_z_yy_xyz_1,   \
                             ta1_z_yy_xyzz_0,  \
                             ta1_z_yy_xyzz_1,  \
                             ta1_z_yy_xzz_0,   \
                             ta1_z_yy_xzz_1,   \
                             ta1_z_yy_xzzz_0,  \
                             ta1_z_yy_xzzz_1,  \
                             ta1_z_yy_yyy_0,   \
                             ta1_z_yy_yyy_1,   \
                             ta1_z_yy_yyyy_0,  \
                             ta1_z_yy_yyyy_1,  \
                             ta1_z_yy_yyyz_0,  \
                             ta1_z_yy_yyyz_1,  \
                             ta1_z_yy_yyz_0,   \
                             ta1_z_yy_yyz_1,   \
                             ta1_z_yy_yyzz_0,  \
                             ta1_z_yy_yyzz_1,  \
                             ta1_z_yy_yzz_0,   \
                             ta1_z_yy_yzz_1,   \
                             ta1_z_yy_yzzz_0,  \
                             ta1_z_yy_yzzz_1,  \
                             ta1_z_yy_zzz_0,   \
                             ta1_z_yy_zzz_1,   \
                             ta1_z_yy_zzzz_0,  \
                             ta1_z_yy_zzzz_1,  \
                             ta1_z_yyy_xxxx_0, \
                             ta1_z_yyy_xxxy_0, \
                             ta1_z_yyy_xxxz_0, \
                             ta1_z_yyy_xxyy_0, \
                             ta1_z_yyy_xxyz_0, \
                             ta1_z_yyy_xxzz_0, \
                             ta1_z_yyy_xyyy_0, \
                             ta1_z_yyy_xyyz_0, \
                             ta1_z_yyy_xyzz_0, \
                             ta1_z_yyy_xzzz_0, \
                             ta1_z_yyy_yyyy_0, \
                             ta1_z_yyy_yyyz_0, \
                             ta1_z_yyy_yyzz_0, \
                             ta1_z_yyy_yzzz_0, \
                             ta1_z_yyy_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyy_xxxx_0[i] = 2.0 * ta1_z_y_xxxx_0[i] * fe_0 - 2.0 * ta1_z_y_xxxx_1[i] * fe_0 +
                              ta1_z_yy_xxxx_0[i] * pa_y[i] - ta1_z_yy_xxxx_1[i] * pc_y[i];

        ta1_z_yyy_xxxy_0[i] = 2.0 * ta1_z_y_xxxy_0[i] * fe_0 - 2.0 * ta1_z_y_xxxy_1[i] * fe_0 +
                              ta1_z_yy_xxx_0[i] * fe_0 - ta1_z_yy_xxx_1[i] * fe_0 + ta1_z_yy_xxxy_0[i] * pa_y[i] -
                              ta1_z_yy_xxxy_1[i] * pc_y[i];

        ta1_z_yyy_xxxz_0[i] = 2.0 * ta1_z_y_xxxz_0[i] * fe_0 - 2.0 * ta1_z_y_xxxz_1[i] * fe_0 +
                              ta1_z_yy_xxxz_0[i] * pa_y[i] - ta1_z_yy_xxxz_1[i] * pc_y[i];

        ta1_z_yyy_xxyy_0[i] = 2.0 * ta1_z_y_xxyy_0[i] * fe_0 - 2.0 * ta1_z_y_xxyy_1[i] * fe_0 +
                              2.0 * ta1_z_yy_xxy_0[i] * fe_0 - 2.0 * ta1_z_yy_xxy_1[i] * fe_0 +
                              ta1_z_yy_xxyy_0[i] * pa_y[i] - ta1_z_yy_xxyy_1[i] * pc_y[i];

        ta1_z_yyy_xxyz_0[i] = 2.0 * ta1_z_y_xxyz_0[i] * fe_0 - 2.0 * ta1_z_y_xxyz_1[i] * fe_0 +
                              ta1_z_yy_xxz_0[i] * fe_0 - ta1_z_yy_xxz_1[i] * fe_0 + ta1_z_yy_xxyz_0[i] * pa_y[i] -
                              ta1_z_yy_xxyz_1[i] * pc_y[i];

        ta1_z_yyy_xxzz_0[i] = 2.0 * ta1_z_y_xxzz_0[i] * fe_0 - 2.0 * ta1_z_y_xxzz_1[i] * fe_0 +
                              ta1_z_yy_xxzz_0[i] * pa_y[i] - ta1_z_yy_xxzz_1[i] * pc_y[i];

        ta1_z_yyy_xyyy_0[i] = 2.0 * ta1_z_y_xyyy_0[i] * fe_0 - 2.0 * ta1_z_y_xyyy_1[i] * fe_0 +
                              3.0 * ta1_z_yy_xyy_0[i] * fe_0 - 3.0 * ta1_z_yy_xyy_1[i] * fe_0 +
                              ta1_z_yy_xyyy_0[i] * pa_y[i] - ta1_z_yy_xyyy_1[i] * pc_y[i];

        ta1_z_yyy_xyyz_0[i] = 2.0 * ta1_z_y_xyyz_0[i] * fe_0 - 2.0 * ta1_z_y_xyyz_1[i] * fe_0 +
                              2.0 * ta1_z_yy_xyz_0[i] * fe_0 - 2.0 * ta1_z_yy_xyz_1[i] * fe_0 +
                              ta1_z_yy_xyyz_0[i] * pa_y[i] - ta1_z_yy_xyyz_1[i] * pc_y[i];

        ta1_z_yyy_xyzz_0[i] = 2.0 * ta1_z_y_xyzz_0[i] * fe_0 - 2.0 * ta1_z_y_xyzz_1[i] * fe_0 +
                              ta1_z_yy_xzz_0[i] * fe_0 - ta1_z_yy_xzz_1[i] * fe_0 + ta1_z_yy_xyzz_0[i] * pa_y[i] -
                              ta1_z_yy_xyzz_1[i] * pc_y[i];

        ta1_z_yyy_xzzz_0[i] = 2.0 * ta1_z_y_xzzz_0[i] * fe_0 - 2.0 * ta1_z_y_xzzz_1[i] * fe_0 +
                              ta1_z_yy_xzzz_0[i] * pa_y[i] - ta1_z_yy_xzzz_1[i] * pc_y[i];

        ta1_z_yyy_yyyy_0[i] = 2.0 * ta1_z_y_yyyy_0[i] * fe_0 - 2.0 * ta1_z_y_yyyy_1[i] * fe_0 +
                              4.0 * ta1_z_yy_yyy_0[i] * fe_0 - 4.0 * ta1_z_yy_yyy_1[i] * fe_0 +
                              ta1_z_yy_yyyy_0[i] * pa_y[i] - ta1_z_yy_yyyy_1[i] * pc_y[i];

        ta1_z_yyy_yyyz_0[i] = 2.0 * ta1_z_y_yyyz_0[i] * fe_0 - 2.0 * ta1_z_y_yyyz_1[i] * fe_0 +
                              3.0 * ta1_z_yy_yyz_0[i] * fe_0 - 3.0 * ta1_z_yy_yyz_1[i] * fe_0 +
                              ta1_z_yy_yyyz_0[i] * pa_y[i] - ta1_z_yy_yyyz_1[i] * pc_y[i];

        ta1_z_yyy_yyzz_0[i] = 2.0 * ta1_z_y_yyzz_0[i] * fe_0 - 2.0 * ta1_z_y_yyzz_1[i] * fe_0 +
                              2.0 * ta1_z_yy_yzz_0[i] * fe_0 - 2.0 * ta1_z_yy_yzz_1[i] * fe_0 +
                              ta1_z_yy_yyzz_0[i] * pa_y[i] - ta1_z_yy_yyzz_1[i] * pc_y[i];

        ta1_z_yyy_yzzz_0[i] = 2.0 * ta1_z_y_yzzz_0[i] * fe_0 - 2.0 * ta1_z_y_yzzz_1[i] * fe_0 +
                              ta1_z_yy_zzz_0[i] * fe_0 - ta1_z_yy_zzz_1[i] * fe_0 + ta1_z_yy_yzzz_0[i] * pa_y[i] -
                              ta1_z_yy_yzzz_1[i] * pc_y[i];

        ta1_z_yyy_zzzz_0[i] = 2.0 * ta1_z_y_zzzz_0[i] * fe_0 - 2.0 * ta1_z_y_zzzz_1[i] * fe_0 +
                              ta1_z_yy_zzzz_0[i] * pa_y[i] - ta1_z_yy_zzzz_1[i] * pc_y[i];
    }

    // Set up 405-420 components of targeted buffer : FG

    auto ta1_z_yyz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 405);

    auto ta1_z_yyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 406);

    auto ta1_z_yyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 407);

    auto ta1_z_yyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 408);

    auto ta1_z_yyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 409);

    auto ta1_z_yyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 410);

    auto ta1_z_yyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 411);

    auto ta1_z_yyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 412);

    auto ta1_z_yyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 413);

    auto ta1_z_yyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 414);

    auto ta1_z_yyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 415);

    auto ta1_z_yyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 416);

    auto ta1_z_yyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 417);

    auto ta1_z_yyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 418);

    auto ta1_z_yyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 419);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_z_yy_xxxx_0,  \
                             ta1_z_yy_xxxx_1,  \
                             ta1_z_yy_xxxy_0,  \
                             ta1_z_yy_xxxy_1,  \
                             ta1_z_yy_xxy_0,   \
                             ta1_z_yy_xxy_1,   \
                             ta1_z_yy_xxyy_0,  \
                             ta1_z_yy_xxyy_1,  \
                             ta1_z_yy_xxyz_0,  \
                             ta1_z_yy_xxyz_1,  \
                             ta1_z_yy_xyy_0,   \
                             ta1_z_yy_xyy_1,   \
                             ta1_z_yy_xyyy_0,  \
                             ta1_z_yy_xyyy_1,  \
                             ta1_z_yy_xyyz_0,  \
                             ta1_z_yy_xyyz_1,  \
                             ta1_z_yy_xyz_0,   \
                             ta1_z_yy_xyz_1,   \
                             ta1_z_yy_xyzz_0,  \
                             ta1_z_yy_xyzz_1,  \
                             ta1_z_yy_yyy_0,   \
                             ta1_z_yy_yyy_1,   \
                             ta1_z_yy_yyyy_0,  \
                             ta1_z_yy_yyyy_1,  \
                             ta1_z_yy_yyyz_0,  \
                             ta1_z_yy_yyyz_1,  \
                             ta1_z_yy_yyz_0,   \
                             ta1_z_yy_yyz_1,   \
                             ta1_z_yy_yyzz_0,  \
                             ta1_z_yy_yyzz_1,  \
                             ta1_z_yy_yzz_0,   \
                             ta1_z_yy_yzz_1,   \
                             ta1_z_yy_yzzz_0,  \
                             ta1_z_yy_yzzz_1,  \
                             ta1_z_yyz_xxxx_0, \
                             ta1_z_yyz_xxxy_0, \
                             ta1_z_yyz_xxxz_0, \
                             ta1_z_yyz_xxyy_0, \
                             ta1_z_yyz_xxyz_0, \
                             ta1_z_yyz_xxzz_0, \
                             ta1_z_yyz_xyyy_0, \
                             ta1_z_yyz_xyyz_0, \
                             ta1_z_yyz_xyzz_0, \
                             ta1_z_yyz_xzzz_0, \
                             ta1_z_yyz_yyyy_0, \
                             ta1_z_yyz_yyyz_0, \
                             ta1_z_yyz_yyzz_0, \
                             ta1_z_yyz_yzzz_0, \
                             ta1_z_yyz_zzzz_0, \
                             ta1_z_yz_xxxz_0,  \
                             ta1_z_yz_xxxz_1,  \
                             ta1_z_yz_xxzz_0,  \
                             ta1_z_yz_xxzz_1,  \
                             ta1_z_yz_xzzz_0,  \
                             ta1_z_yz_xzzz_1,  \
                             ta1_z_yz_zzzz_0,  \
                             ta1_z_yz_zzzz_1,  \
                             ta1_z_z_xxxz_0,   \
                             ta1_z_z_xxxz_1,   \
                             ta1_z_z_xxzz_0,   \
                             ta1_z_z_xxzz_1,   \
                             ta1_z_z_xzzz_0,   \
                             ta1_z_z_xzzz_1,   \
                             ta1_z_z_zzzz_0,   \
                             ta1_z_z_zzzz_1,   \
                             ta_yy_xxxx_1,     \
                             ta_yy_xxxy_1,     \
                             ta_yy_xxyy_1,     \
                             ta_yy_xxyz_1,     \
                             ta_yy_xyyy_1,     \
                             ta_yy_xyyz_1,     \
                             ta_yy_xyzz_1,     \
                             ta_yy_yyyy_1,     \
                             ta_yy_yyyz_1,     \
                             ta_yy_yyzz_1,     \
                             ta_yy_yzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyz_xxxx_0[i] = ta_yy_xxxx_1[i] + ta1_z_yy_xxxx_0[i] * pa_z[i] - ta1_z_yy_xxxx_1[i] * pc_z[i];

        ta1_z_yyz_xxxy_0[i] = ta_yy_xxxy_1[i] + ta1_z_yy_xxxy_0[i] * pa_z[i] - ta1_z_yy_xxxy_1[i] * pc_z[i];

        ta1_z_yyz_xxxz_0[i] = ta1_z_z_xxxz_0[i] * fe_0 - ta1_z_z_xxxz_1[i] * fe_0 + ta1_z_yz_xxxz_0[i] * pa_y[i] -
                              ta1_z_yz_xxxz_1[i] * pc_y[i];

        ta1_z_yyz_xxyy_0[i] = ta_yy_xxyy_1[i] + ta1_z_yy_xxyy_0[i] * pa_z[i] - ta1_z_yy_xxyy_1[i] * pc_z[i];

        ta1_z_yyz_xxyz_0[i] = ta1_z_yy_xxy_0[i] * fe_0 - ta1_z_yy_xxy_1[i] * fe_0 + ta_yy_xxyz_1[i] +
                              ta1_z_yy_xxyz_0[i] * pa_z[i] - ta1_z_yy_xxyz_1[i] * pc_z[i];

        ta1_z_yyz_xxzz_0[i] = ta1_z_z_xxzz_0[i] * fe_0 - ta1_z_z_xxzz_1[i] * fe_0 + ta1_z_yz_xxzz_0[i] * pa_y[i] -
                              ta1_z_yz_xxzz_1[i] * pc_y[i];

        ta1_z_yyz_xyyy_0[i] = ta_yy_xyyy_1[i] + ta1_z_yy_xyyy_0[i] * pa_z[i] - ta1_z_yy_xyyy_1[i] * pc_z[i];

        ta1_z_yyz_xyyz_0[i] = ta1_z_yy_xyy_0[i] * fe_0 - ta1_z_yy_xyy_1[i] * fe_0 + ta_yy_xyyz_1[i] +
                              ta1_z_yy_xyyz_0[i] * pa_z[i] - ta1_z_yy_xyyz_1[i] * pc_z[i];

        ta1_z_yyz_xyzz_0[i] = 2.0 * ta1_z_yy_xyz_0[i] * fe_0 - 2.0 * ta1_z_yy_xyz_1[i] * fe_0 + ta_yy_xyzz_1[i] +
                              ta1_z_yy_xyzz_0[i] * pa_z[i] - ta1_z_yy_xyzz_1[i] * pc_z[i];

        ta1_z_yyz_xzzz_0[i] = ta1_z_z_xzzz_0[i] * fe_0 - ta1_z_z_xzzz_1[i] * fe_0 + ta1_z_yz_xzzz_0[i] * pa_y[i] -
                              ta1_z_yz_xzzz_1[i] * pc_y[i];

        ta1_z_yyz_yyyy_0[i] = ta_yy_yyyy_1[i] + ta1_z_yy_yyyy_0[i] * pa_z[i] - ta1_z_yy_yyyy_1[i] * pc_z[i];

        ta1_z_yyz_yyyz_0[i] = ta1_z_yy_yyy_0[i] * fe_0 - ta1_z_yy_yyy_1[i] * fe_0 + ta_yy_yyyz_1[i] +
                              ta1_z_yy_yyyz_0[i] * pa_z[i] - ta1_z_yy_yyyz_1[i] * pc_z[i];

        ta1_z_yyz_yyzz_0[i] = 2.0 * ta1_z_yy_yyz_0[i] * fe_0 - 2.0 * ta1_z_yy_yyz_1[i] * fe_0 + ta_yy_yyzz_1[i] +
                              ta1_z_yy_yyzz_0[i] * pa_z[i] - ta1_z_yy_yyzz_1[i] * pc_z[i];

        ta1_z_yyz_yzzz_0[i] = 3.0 * ta1_z_yy_yzz_0[i] * fe_0 - 3.0 * ta1_z_yy_yzz_1[i] * fe_0 + ta_yy_yzzz_1[i] +
                              ta1_z_yy_yzzz_0[i] * pa_z[i] - ta1_z_yy_yzzz_1[i] * pc_z[i];

        ta1_z_yyz_zzzz_0[i] = ta1_z_z_zzzz_0[i] * fe_0 - ta1_z_z_zzzz_1[i] * fe_0 + ta1_z_yz_zzzz_0[i] * pa_y[i] -
                              ta1_z_yz_zzzz_1[i] * pc_y[i];
    }

    // Set up 420-435 components of targeted buffer : FG

    auto ta1_z_yzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 420);

    auto ta1_z_yzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 421);

    auto ta1_z_yzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 422);

    auto ta1_z_yzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 423);

    auto ta1_z_yzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 424);

    auto ta1_z_yzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 425);

    auto ta1_z_yzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 426);

    auto ta1_z_yzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 427);

    auto ta1_z_yzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 428);

    auto ta1_z_yzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 429);

    auto ta1_z_yzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 430);

    auto ta1_z_yzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 431);

    auto ta1_z_yzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 432);

    auto ta1_z_yzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 433);

    auto ta1_z_yzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 434);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_z_yzz_xxxx_0, \
                             ta1_z_yzz_xxxy_0, \
                             ta1_z_yzz_xxxz_0, \
                             ta1_z_yzz_xxyy_0, \
                             ta1_z_yzz_xxyz_0, \
                             ta1_z_yzz_xxzz_0, \
                             ta1_z_yzz_xyyy_0, \
                             ta1_z_yzz_xyyz_0, \
                             ta1_z_yzz_xyzz_0, \
                             ta1_z_yzz_xzzz_0, \
                             ta1_z_yzz_yyyy_0, \
                             ta1_z_yzz_yyyz_0, \
                             ta1_z_yzz_yyzz_0, \
                             ta1_z_yzz_yzzz_0, \
                             ta1_z_yzz_zzzz_0, \
                             ta1_z_zz_xxx_0,   \
                             ta1_z_zz_xxx_1,   \
                             ta1_z_zz_xxxx_0,  \
                             ta1_z_zz_xxxx_1,  \
                             ta1_z_zz_xxxy_0,  \
                             ta1_z_zz_xxxy_1,  \
                             ta1_z_zz_xxxz_0,  \
                             ta1_z_zz_xxxz_1,  \
                             ta1_z_zz_xxy_0,   \
                             ta1_z_zz_xxy_1,   \
                             ta1_z_zz_xxyy_0,  \
                             ta1_z_zz_xxyy_1,  \
                             ta1_z_zz_xxyz_0,  \
                             ta1_z_zz_xxyz_1,  \
                             ta1_z_zz_xxz_0,   \
                             ta1_z_zz_xxz_1,   \
                             ta1_z_zz_xxzz_0,  \
                             ta1_z_zz_xxzz_1,  \
                             ta1_z_zz_xyy_0,   \
                             ta1_z_zz_xyy_1,   \
                             ta1_z_zz_xyyy_0,  \
                             ta1_z_zz_xyyy_1,  \
                             ta1_z_zz_xyyz_0,  \
                             ta1_z_zz_xyyz_1,  \
                             ta1_z_zz_xyz_0,   \
                             ta1_z_zz_xyz_1,   \
                             ta1_z_zz_xyzz_0,  \
                             ta1_z_zz_xyzz_1,  \
                             ta1_z_zz_xzz_0,   \
                             ta1_z_zz_xzz_1,   \
                             ta1_z_zz_xzzz_0,  \
                             ta1_z_zz_xzzz_1,  \
                             ta1_z_zz_yyy_0,   \
                             ta1_z_zz_yyy_1,   \
                             ta1_z_zz_yyyy_0,  \
                             ta1_z_zz_yyyy_1,  \
                             ta1_z_zz_yyyz_0,  \
                             ta1_z_zz_yyyz_1,  \
                             ta1_z_zz_yyz_0,   \
                             ta1_z_zz_yyz_1,   \
                             ta1_z_zz_yyzz_0,  \
                             ta1_z_zz_yyzz_1,  \
                             ta1_z_zz_yzz_0,   \
                             ta1_z_zz_yzz_1,   \
                             ta1_z_zz_yzzz_0,  \
                             ta1_z_zz_yzzz_1,  \
                             ta1_z_zz_zzz_0,   \
                             ta1_z_zz_zzz_1,   \
                             ta1_z_zz_zzzz_0,  \
                             ta1_z_zz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yzz_xxxx_0[i] = ta1_z_zz_xxxx_0[i] * pa_y[i] - ta1_z_zz_xxxx_1[i] * pc_y[i];

        ta1_z_yzz_xxxy_0[i] = ta1_z_zz_xxx_0[i] * fe_0 - ta1_z_zz_xxx_1[i] * fe_0 + ta1_z_zz_xxxy_0[i] * pa_y[i] -
                              ta1_z_zz_xxxy_1[i] * pc_y[i];

        ta1_z_yzz_xxxz_0[i] = ta1_z_zz_xxxz_0[i] * pa_y[i] - ta1_z_zz_xxxz_1[i] * pc_y[i];

        ta1_z_yzz_xxyy_0[i] = 2.0 * ta1_z_zz_xxy_0[i] * fe_0 - 2.0 * ta1_z_zz_xxy_1[i] * fe_0 +
                              ta1_z_zz_xxyy_0[i] * pa_y[i] - ta1_z_zz_xxyy_1[i] * pc_y[i];

        ta1_z_yzz_xxyz_0[i] = ta1_z_zz_xxz_0[i] * fe_0 - ta1_z_zz_xxz_1[i] * fe_0 + ta1_z_zz_xxyz_0[i] * pa_y[i] -
                              ta1_z_zz_xxyz_1[i] * pc_y[i];

        ta1_z_yzz_xxzz_0[i] = ta1_z_zz_xxzz_0[i] * pa_y[i] - ta1_z_zz_xxzz_1[i] * pc_y[i];

        ta1_z_yzz_xyyy_0[i] = 3.0 * ta1_z_zz_xyy_0[i] * fe_0 - 3.0 * ta1_z_zz_xyy_1[i] * fe_0 +
                              ta1_z_zz_xyyy_0[i] * pa_y[i] - ta1_z_zz_xyyy_1[i] * pc_y[i];

        ta1_z_yzz_xyyz_0[i] = 2.0 * ta1_z_zz_xyz_0[i] * fe_0 - 2.0 * ta1_z_zz_xyz_1[i] * fe_0 +
                              ta1_z_zz_xyyz_0[i] * pa_y[i] - ta1_z_zz_xyyz_1[i] * pc_y[i];

        ta1_z_yzz_xyzz_0[i] = ta1_z_zz_xzz_0[i] * fe_0 - ta1_z_zz_xzz_1[i] * fe_0 + ta1_z_zz_xyzz_0[i] * pa_y[i] -
                              ta1_z_zz_xyzz_1[i] * pc_y[i];

        ta1_z_yzz_xzzz_0[i] = ta1_z_zz_xzzz_0[i] * pa_y[i] - ta1_z_zz_xzzz_1[i] * pc_y[i];

        ta1_z_yzz_yyyy_0[i] = 4.0 * ta1_z_zz_yyy_0[i] * fe_0 - 4.0 * ta1_z_zz_yyy_1[i] * fe_0 +
                              ta1_z_zz_yyyy_0[i] * pa_y[i] - ta1_z_zz_yyyy_1[i] * pc_y[i];

        ta1_z_yzz_yyyz_0[i] = 3.0 * ta1_z_zz_yyz_0[i] * fe_0 - 3.0 * ta1_z_zz_yyz_1[i] * fe_0 +
                              ta1_z_zz_yyyz_0[i] * pa_y[i] - ta1_z_zz_yyyz_1[i] * pc_y[i];

        ta1_z_yzz_yyzz_0[i] = 2.0 * ta1_z_zz_yzz_0[i] * fe_0 - 2.0 * ta1_z_zz_yzz_1[i] * fe_0 +
                              ta1_z_zz_yyzz_0[i] * pa_y[i] - ta1_z_zz_yyzz_1[i] * pc_y[i];

        ta1_z_yzz_yzzz_0[i] = ta1_z_zz_zzz_0[i] * fe_0 - ta1_z_zz_zzz_1[i] * fe_0 + ta1_z_zz_yzzz_0[i] * pa_y[i] -
                              ta1_z_zz_yzzz_1[i] * pc_y[i];

        ta1_z_yzz_zzzz_0[i] = ta1_z_zz_zzzz_0[i] * pa_y[i] - ta1_z_zz_zzzz_1[i] * pc_y[i];
    }

    // Set up 435-450 components of targeted buffer : FG

    auto ta1_z_zzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 435);

    auto ta1_z_zzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 436);

    auto ta1_z_zzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 437);

    auto ta1_z_zzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 438);

    auto ta1_z_zzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 439);

    auto ta1_z_zzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 440);

    auto ta1_z_zzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 441);

    auto ta1_z_zzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 442);

    auto ta1_z_zzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 443);

    auto ta1_z_zzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 444);

    auto ta1_z_zzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 445);

    auto ta1_z_zzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 446);

    auto ta1_z_zzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 447);

    auto ta1_z_zzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 448);

    auto ta1_z_zzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 449);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta1_z_z_xxxx_0,   \
                             ta1_z_z_xxxx_1,   \
                             ta1_z_z_xxxy_0,   \
                             ta1_z_z_xxxy_1,   \
                             ta1_z_z_xxxz_0,   \
                             ta1_z_z_xxxz_1,   \
                             ta1_z_z_xxyy_0,   \
                             ta1_z_z_xxyy_1,   \
                             ta1_z_z_xxyz_0,   \
                             ta1_z_z_xxyz_1,   \
                             ta1_z_z_xxzz_0,   \
                             ta1_z_z_xxzz_1,   \
                             ta1_z_z_xyyy_0,   \
                             ta1_z_z_xyyy_1,   \
                             ta1_z_z_xyyz_0,   \
                             ta1_z_z_xyyz_1,   \
                             ta1_z_z_xyzz_0,   \
                             ta1_z_z_xyzz_1,   \
                             ta1_z_z_xzzz_0,   \
                             ta1_z_z_xzzz_1,   \
                             ta1_z_z_yyyy_0,   \
                             ta1_z_z_yyyy_1,   \
                             ta1_z_z_yyyz_0,   \
                             ta1_z_z_yyyz_1,   \
                             ta1_z_z_yyzz_0,   \
                             ta1_z_z_yyzz_1,   \
                             ta1_z_z_yzzz_0,   \
                             ta1_z_z_yzzz_1,   \
                             ta1_z_z_zzzz_0,   \
                             ta1_z_z_zzzz_1,   \
                             ta1_z_zz_xxx_0,   \
                             ta1_z_zz_xxx_1,   \
                             ta1_z_zz_xxxx_0,  \
                             ta1_z_zz_xxxx_1,  \
                             ta1_z_zz_xxxy_0,  \
                             ta1_z_zz_xxxy_1,  \
                             ta1_z_zz_xxxz_0,  \
                             ta1_z_zz_xxxz_1,  \
                             ta1_z_zz_xxy_0,   \
                             ta1_z_zz_xxy_1,   \
                             ta1_z_zz_xxyy_0,  \
                             ta1_z_zz_xxyy_1,  \
                             ta1_z_zz_xxyz_0,  \
                             ta1_z_zz_xxyz_1,  \
                             ta1_z_zz_xxz_0,   \
                             ta1_z_zz_xxz_1,   \
                             ta1_z_zz_xxzz_0,  \
                             ta1_z_zz_xxzz_1,  \
                             ta1_z_zz_xyy_0,   \
                             ta1_z_zz_xyy_1,   \
                             ta1_z_zz_xyyy_0,  \
                             ta1_z_zz_xyyy_1,  \
                             ta1_z_zz_xyyz_0,  \
                             ta1_z_zz_xyyz_1,  \
                             ta1_z_zz_xyz_0,   \
                             ta1_z_zz_xyz_1,   \
                             ta1_z_zz_xyzz_0,  \
                             ta1_z_zz_xyzz_1,  \
                             ta1_z_zz_xzz_0,   \
                             ta1_z_zz_xzz_1,   \
                             ta1_z_zz_xzzz_0,  \
                             ta1_z_zz_xzzz_1,  \
                             ta1_z_zz_yyy_0,   \
                             ta1_z_zz_yyy_1,   \
                             ta1_z_zz_yyyy_0,  \
                             ta1_z_zz_yyyy_1,  \
                             ta1_z_zz_yyyz_0,  \
                             ta1_z_zz_yyyz_1,  \
                             ta1_z_zz_yyz_0,   \
                             ta1_z_zz_yyz_1,   \
                             ta1_z_zz_yyzz_0,  \
                             ta1_z_zz_yyzz_1,  \
                             ta1_z_zz_yzz_0,   \
                             ta1_z_zz_yzz_1,   \
                             ta1_z_zz_yzzz_0,  \
                             ta1_z_zz_yzzz_1,  \
                             ta1_z_zz_zzz_0,   \
                             ta1_z_zz_zzz_1,   \
                             ta1_z_zz_zzzz_0,  \
                             ta1_z_zz_zzzz_1,  \
                             ta1_z_zzz_xxxx_0, \
                             ta1_z_zzz_xxxy_0, \
                             ta1_z_zzz_xxxz_0, \
                             ta1_z_zzz_xxyy_0, \
                             ta1_z_zzz_xxyz_0, \
                             ta1_z_zzz_xxzz_0, \
                             ta1_z_zzz_xyyy_0, \
                             ta1_z_zzz_xyyz_0, \
                             ta1_z_zzz_xyzz_0, \
                             ta1_z_zzz_xzzz_0, \
                             ta1_z_zzz_yyyy_0, \
                             ta1_z_zzz_yyyz_0, \
                             ta1_z_zzz_yyzz_0, \
                             ta1_z_zzz_yzzz_0, \
                             ta1_z_zzz_zzzz_0, \
                             ta_zz_xxxx_1,     \
                             ta_zz_xxxy_1,     \
                             ta_zz_xxxz_1,     \
                             ta_zz_xxyy_1,     \
                             ta_zz_xxyz_1,     \
                             ta_zz_xxzz_1,     \
                             ta_zz_xyyy_1,     \
                             ta_zz_xyyz_1,     \
                             ta_zz_xyzz_1,     \
                             ta_zz_xzzz_1,     \
                             ta_zz_yyyy_1,     \
                             ta_zz_yyyz_1,     \
                             ta_zz_yyzz_1,     \
                             ta_zz_yzzz_1,     \
                             ta_zz_zzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_zzz_xxxx_0[i] = 2.0 * ta1_z_z_xxxx_0[i] * fe_0 - 2.0 * ta1_z_z_xxxx_1[i] * fe_0 + ta_zz_xxxx_1[i] +
                              ta1_z_zz_xxxx_0[i] * pa_z[i] - ta1_z_zz_xxxx_1[i] * pc_z[i];

        ta1_z_zzz_xxxy_0[i] = 2.0 * ta1_z_z_xxxy_0[i] * fe_0 - 2.0 * ta1_z_z_xxxy_1[i] * fe_0 + ta_zz_xxxy_1[i] +
                              ta1_z_zz_xxxy_0[i] * pa_z[i] - ta1_z_zz_xxxy_1[i] * pc_z[i];

        ta1_z_zzz_xxxz_0[i] = 2.0 * ta1_z_z_xxxz_0[i] * fe_0 - 2.0 * ta1_z_z_xxxz_1[i] * fe_0 +
                              ta1_z_zz_xxx_0[i] * fe_0 - ta1_z_zz_xxx_1[i] * fe_0 + ta_zz_xxxz_1[i] +
                              ta1_z_zz_xxxz_0[i] * pa_z[i] - ta1_z_zz_xxxz_1[i] * pc_z[i];

        ta1_z_zzz_xxyy_0[i] = 2.0 * ta1_z_z_xxyy_0[i] * fe_0 - 2.0 * ta1_z_z_xxyy_1[i] * fe_0 + ta_zz_xxyy_1[i] +
                              ta1_z_zz_xxyy_0[i] * pa_z[i] - ta1_z_zz_xxyy_1[i] * pc_z[i];

        ta1_z_zzz_xxyz_0[i] = 2.0 * ta1_z_z_xxyz_0[i] * fe_0 - 2.0 * ta1_z_z_xxyz_1[i] * fe_0 +
                              ta1_z_zz_xxy_0[i] * fe_0 - ta1_z_zz_xxy_1[i] * fe_0 + ta_zz_xxyz_1[i] +
                              ta1_z_zz_xxyz_0[i] * pa_z[i] - ta1_z_zz_xxyz_1[i] * pc_z[i];

        ta1_z_zzz_xxzz_0[i] = 2.0 * ta1_z_z_xxzz_0[i] * fe_0 - 2.0 * ta1_z_z_xxzz_1[i] * fe_0 +
                              2.0 * ta1_z_zz_xxz_0[i] * fe_0 - 2.0 * ta1_z_zz_xxz_1[i] * fe_0 + ta_zz_xxzz_1[i] +
                              ta1_z_zz_xxzz_0[i] * pa_z[i] - ta1_z_zz_xxzz_1[i] * pc_z[i];

        ta1_z_zzz_xyyy_0[i] = 2.0 * ta1_z_z_xyyy_0[i] * fe_0 - 2.0 * ta1_z_z_xyyy_1[i] * fe_0 + ta_zz_xyyy_1[i] +
                              ta1_z_zz_xyyy_0[i] * pa_z[i] - ta1_z_zz_xyyy_1[i] * pc_z[i];

        ta1_z_zzz_xyyz_0[i] = 2.0 * ta1_z_z_xyyz_0[i] * fe_0 - 2.0 * ta1_z_z_xyyz_1[i] * fe_0 +
                              ta1_z_zz_xyy_0[i] * fe_0 - ta1_z_zz_xyy_1[i] * fe_0 + ta_zz_xyyz_1[i] +
                              ta1_z_zz_xyyz_0[i] * pa_z[i] - ta1_z_zz_xyyz_1[i] * pc_z[i];

        ta1_z_zzz_xyzz_0[i] = 2.0 * ta1_z_z_xyzz_0[i] * fe_0 - 2.0 * ta1_z_z_xyzz_1[i] * fe_0 +
                              2.0 * ta1_z_zz_xyz_0[i] * fe_0 - 2.0 * ta1_z_zz_xyz_1[i] * fe_0 + ta_zz_xyzz_1[i] +
                              ta1_z_zz_xyzz_0[i] * pa_z[i] - ta1_z_zz_xyzz_1[i] * pc_z[i];

        ta1_z_zzz_xzzz_0[i] = 2.0 * ta1_z_z_xzzz_0[i] * fe_0 - 2.0 * ta1_z_z_xzzz_1[i] * fe_0 +
                              3.0 * ta1_z_zz_xzz_0[i] * fe_0 - 3.0 * ta1_z_zz_xzz_1[i] * fe_0 + ta_zz_xzzz_1[i] +
                              ta1_z_zz_xzzz_0[i] * pa_z[i] - ta1_z_zz_xzzz_1[i] * pc_z[i];

        ta1_z_zzz_yyyy_0[i] = 2.0 * ta1_z_z_yyyy_0[i] * fe_0 - 2.0 * ta1_z_z_yyyy_1[i] * fe_0 + ta_zz_yyyy_1[i] +
                              ta1_z_zz_yyyy_0[i] * pa_z[i] - ta1_z_zz_yyyy_1[i] * pc_z[i];

        ta1_z_zzz_yyyz_0[i] = 2.0 * ta1_z_z_yyyz_0[i] * fe_0 - 2.0 * ta1_z_z_yyyz_1[i] * fe_0 +
                              ta1_z_zz_yyy_0[i] * fe_0 - ta1_z_zz_yyy_1[i] * fe_0 + ta_zz_yyyz_1[i] +
                              ta1_z_zz_yyyz_0[i] * pa_z[i] - ta1_z_zz_yyyz_1[i] * pc_z[i];

        ta1_z_zzz_yyzz_0[i] = 2.0 * ta1_z_z_yyzz_0[i] * fe_0 - 2.0 * ta1_z_z_yyzz_1[i] * fe_0 +
                              2.0 * ta1_z_zz_yyz_0[i] * fe_0 - 2.0 * ta1_z_zz_yyz_1[i] * fe_0 + ta_zz_yyzz_1[i] +
                              ta1_z_zz_yyzz_0[i] * pa_z[i] - ta1_z_zz_yyzz_1[i] * pc_z[i];

        ta1_z_zzz_yzzz_0[i] = 2.0 * ta1_z_z_yzzz_0[i] * fe_0 - 2.0 * ta1_z_z_yzzz_1[i] * fe_0 +
                              3.0 * ta1_z_zz_yzz_0[i] * fe_0 - 3.0 * ta1_z_zz_yzz_1[i] * fe_0 + ta_zz_yzzz_1[i] +
                              ta1_z_zz_yzzz_0[i] * pa_z[i] - ta1_z_zz_yzzz_1[i] * pc_z[i];

        ta1_z_zzz_zzzz_0[i] = 2.0 * ta1_z_z_zzzz_0[i] * fe_0 - 2.0 * ta1_z_z_zzzz_1[i] * fe_0 +
                              4.0 * ta1_z_zz_zzz_0[i] * fe_0 - 4.0 * ta1_z_zz_zzz_1[i] * fe_0 + ta_zz_zzzz_1[i] +
                              ta1_z_zz_zzzz_0[i] * pa_z[i] - ta1_z_zz_zzzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
