#include "NuclearPotentialGeom010PrimRecIP.hpp"

namespace npotrec { // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_ip(CSimdArray<double>& pbuffer, 
                                        const size_t idx_npot_geom_010_0_ip,
                                        const size_t idx_npot_geom_010_0_gp,
                                        const size_t idx_npot_geom_010_1_gp,
                                        const size_t idx_npot_geom_010_0_hs,
                                        const size_t idx_npot_geom_010_1_hs,
                                        const size_t idx_npot_1_hp,
                                        const size_t idx_npot_geom_010_0_hp,
                                        const size_t idx_npot_geom_010_1_hp,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_rpa,
                                        const size_t idx_rpc,
                                        const double a_exp) -> void
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

    // Set up components of auxiliary buffer : GP

    auto ta1_x_xxxx_x_0 = pbuffer.data(idx_npot_geom_010_0_gp);

    auto ta1_x_xxxx_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 1);

    auto ta1_x_xxxx_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 2);

    auto ta1_x_xxxy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 3);

    auto ta1_x_xxxy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 5);

    auto ta1_x_xxxz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 6);

    auto ta1_x_xxxz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 7);

    auto ta1_x_xxxz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 8);

    auto ta1_x_xxyy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 9);

    auto ta1_x_xxyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 10);

    auto ta1_x_xxyy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 11);

    auto ta1_x_xxyz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 14);

    auto ta1_x_xxzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 15);

    auto ta1_x_xxzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 16);

    auto ta1_x_xxzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 17);

    auto ta1_x_xyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 18);

    auto ta1_x_xyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 19);

    auto ta1_x_xyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 24);

    auto ta1_x_xzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 27);

    auto ta1_x_xzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 29);

    auto ta1_x_yyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 30);

    auto ta1_x_yyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 31);

    auto ta1_x_yyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 32);

    auto ta1_x_yyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 34);

    auto ta1_x_yyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 35);

    auto ta1_x_yyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 36);

    auto ta1_x_yyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 37);

    auto ta1_x_yyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 38);

    auto ta1_x_yzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 39);

    auto ta1_x_yzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 41);

    auto ta1_x_zzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 42);

    auto ta1_x_zzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 43);

    auto ta1_x_zzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 44);

    auto ta1_y_xxxx_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 45);

    auto ta1_y_xxxx_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 46);

    auto ta1_y_xxxx_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 47);

    auto ta1_y_xxxy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 48);

    auto ta1_y_xxxy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 49);

    auto ta1_y_xxxz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 51);

    auto ta1_y_xxxz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 53);

    auto ta1_y_xxyy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 54);

    auto ta1_y_xxyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 55);

    auto ta1_y_xxyy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 56);

    auto ta1_y_xxzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 60);

    auto ta1_y_xxzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 61);

    auto ta1_y_xxzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 62);

    auto ta1_y_xyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 64);

    auto ta1_y_xyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 65);

    auto ta1_y_xyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 68);

    auto ta1_y_xyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 70);

    auto ta1_y_xzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 73);

    auto ta1_y_xzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 74);

    auto ta1_y_yyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 75);

    auto ta1_y_yyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 76);

    auto ta1_y_yyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 77);

    auto ta1_y_yyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 78);

    auto ta1_y_yyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 79);

    auto ta1_y_yyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 80);

    auto ta1_y_yyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 81);

    auto ta1_y_yyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 82);

    auto ta1_y_yyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 83);

    auto ta1_y_yzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 85);

    auto ta1_y_yzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 86);

    auto ta1_y_zzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 87);

    auto ta1_y_zzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 88);

    auto ta1_y_zzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 89);

    auto ta1_z_xxxx_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 90);

    auto ta1_z_xxxx_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 91);

    auto ta1_z_xxxx_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 92);

    auto ta1_z_xxxy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 93);

    auto ta1_z_xxxy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 94);

    auto ta1_z_xxxz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 96);

    auto ta1_z_xxxz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 98);

    auto ta1_z_xxyy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 99);

    auto ta1_z_xxyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 100);

    auto ta1_z_xxyy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 101);

    auto ta1_z_xxzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 105);

    auto ta1_z_xxzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 106);

    auto ta1_z_xxzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 107);

    auto ta1_z_xyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 109);

    auto ta1_z_xyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 110);

    auto ta1_z_xyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 113);

    auto ta1_z_xyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 115);

    auto ta1_z_xzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 118);

    auto ta1_z_xzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 119);

    auto ta1_z_yyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 120);

    auto ta1_z_yyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 121);

    auto ta1_z_yyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 122);

    auto ta1_z_yyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 124);

    auto ta1_z_yyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 125);

    auto ta1_z_yyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 126);

    auto ta1_z_yyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 127);

    auto ta1_z_yyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 128);

    auto ta1_z_yzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 129);

    auto ta1_z_yzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 130);

    auto ta1_z_yzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 131);

    auto ta1_z_zzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 132);

    auto ta1_z_zzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 133);

    auto ta1_z_zzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 134);

    // Set up components of auxiliary buffer : GP

    auto ta1_x_xxxx_x_1 = pbuffer.data(idx_npot_geom_010_1_gp);

    auto ta1_x_xxxx_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 1);

    auto ta1_x_xxxx_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 2);

    auto ta1_x_xxxy_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 3);

    auto ta1_x_xxxy_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 5);

    auto ta1_x_xxxz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 6);

    auto ta1_x_xxxz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 7);

    auto ta1_x_xxxz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 8);

    auto ta1_x_xxyy_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 9);

    auto ta1_x_xxyy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 10);

    auto ta1_x_xxyy_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 11);

    auto ta1_x_xxyz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 14);

    auto ta1_x_xxzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 15);

    auto ta1_x_xxzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 16);

    auto ta1_x_xxzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 17);

    auto ta1_x_xyyy_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 18);

    auto ta1_x_xyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 19);

    auto ta1_x_xyzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 24);

    auto ta1_x_xzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 27);

    auto ta1_x_xzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 29);

    auto ta1_x_yyyy_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 30);

    auto ta1_x_yyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 31);

    auto ta1_x_yyyy_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 32);

    auto ta1_x_yyyz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 34);

    auto ta1_x_yyyz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 35);

    auto ta1_x_yyzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 36);

    auto ta1_x_yyzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 37);

    auto ta1_x_yyzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 38);

    auto ta1_x_yzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 39);

    auto ta1_x_yzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 41);

    auto ta1_x_zzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 42);

    auto ta1_x_zzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 43);

    auto ta1_x_zzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 44);

    auto ta1_y_xxxx_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 45);

    auto ta1_y_xxxx_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 46);

    auto ta1_y_xxxx_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 47);

    auto ta1_y_xxxy_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 48);

    auto ta1_y_xxxy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 49);

    auto ta1_y_xxxz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 51);

    auto ta1_y_xxxz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 53);

    auto ta1_y_xxyy_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 54);

    auto ta1_y_xxyy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 55);

    auto ta1_y_xxyy_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 56);

    auto ta1_y_xxzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 60);

    auto ta1_y_xxzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 61);

    auto ta1_y_xxzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 62);

    auto ta1_y_xyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 64);

    auto ta1_y_xyyy_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 65);

    auto ta1_y_xyyz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 68);

    auto ta1_y_xyzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 70);

    auto ta1_y_xzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 73);

    auto ta1_y_xzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 74);

    auto ta1_y_yyyy_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 75);

    auto ta1_y_yyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 76);

    auto ta1_y_yyyy_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 77);

    auto ta1_y_yyyz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 78);

    auto ta1_y_yyyz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 79);

    auto ta1_y_yyyz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 80);

    auto ta1_y_yyzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 81);

    auto ta1_y_yyzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 82);

    auto ta1_y_yyzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 83);

    auto ta1_y_yzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 85);

    auto ta1_y_yzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 86);

    auto ta1_y_zzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 87);

    auto ta1_y_zzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 88);

    auto ta1_y_zzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 89);

    auto ta1_z_xxxx_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 90);

    auto ta1_z_xxxx_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 91);

    auto ta1_z_xxxx_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 92);

    auto ta1_z_xxxy_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 93);

    auto ta1_z_xxxy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 94);

    auto ta1_z_xxxz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 96);

    auto ta1_z_xxxz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 98);

    auto ta1_z_xxyy_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 99);

    auto ta1_z_xxyy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 100);

    auto ta1_z_xxyy_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 101);

    auto ta1_z_xxzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 105);

    auto ta1_z_xxzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 106);

    auto ta1_z_xxzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 107);

    auto ta1_z_xyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 109);

    auto ta1_z_xyyy_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 110);

    auto ta1_z_xyyz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 113);

    auto ta1_z_xyzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 115);

    auto ta1_z_xzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 118);

    auto ta1_z_xzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 119);

    auto ta1_z_yyyy_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 120);

    auto ta1_z_yyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 121);

    auto ta1_z_yyyy_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 122);

    auto ta1_z_yyyz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 124);

    auto ta1_z_yyyz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 125);

    auto ta1_z_yyzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 126);

    auto ta1_z_yyzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 127);

    auto ta1_z_yyzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 128);

    auto ta1_z_yzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 129);

    auto ta1_z_yzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 130);

    auto ta1_z_yzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 131);

    auto ta1_z_zzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 132);

    auto ta1_z_zzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 133);

    auto ta1_z_zzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 134);

    // Set up components of auxiliary buffer : HS

    auto ta1_x_xxxxx_0_0 = pbuffer.data(idx_npot_geom_010_0_hs);

    auto ta1_x_xxxzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 5);

    auto ta1_x_xxzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 9);

    auto ta1_x_yyyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 15);

    auto ta1_x_zzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 20);

    auto ta1_y_xxxxx_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 21);

    auto ta1_y_yyyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 36);

    auto ta1_y_yyyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 38);

    auto ta1_y_yyzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 39);

    auto ta1_y_zzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 41);

    auto ta1_z_xxxxx_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 42);

    auto ta1_z_yyyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 57);

    auto ta1_z_yyyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 59);

    auto ta1_z_yyzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 60);

    auto ta1_z_zzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 62);

    // Set up components of auxiliary buffer : HS

    auto ta1_x_xxxxx_0_1 = pbuffer.data(idx_npot_geom_010_1_hs);

    auto ta1_x_xxxzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 5);

    auto ta1_x_xxzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 9);

    auto ta1_x_yyyyy_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 15);

    auto ta1_x_zzzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 20);

    auto ta1_y_xxxxx_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 21);

    auto ta1_y_yyyyy_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 36);

    auto ta1_y_yyyzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 38);

    auto ta1_y_yyzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 39);

    auto ta1_y_zzzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 41);

    auto ta1_z_xxxxx_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 42);

    auto ta1_z_yyyyy_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 57);

    auto ta1_z_yyyzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 59);

    auto ta1_z_yyzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 60);

    auto ta1_z_zzzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 62);

    // Set up components of auxiliary buffer : HP

    auto ta_xxxxx_x_1 = pbuffer.data(idx_npot_1_hp);

    auto ta_xxxxx_y_1 = pbuffer.data(idx_npot_1_hp + 1);

    auto ta_xxxxx_z_1 = pbuffer.data(idx_npot_1_hp + 2);

    auto ta_xxxxy_x_1 = pbuffer.data(idx_npot_1_hp + 3);

    auto ta_xxxxy_y_1 = pbuffer.data(idx_npot_1_hp + 4);

    auto ta_xxxxz_x_1 = pbuffer.data(idx_npot_1_hp + 6);

    auto ta_xxxxz_z_1 = pbuffer.data(idx_npot_1_hp + 8);

    auto ta_xxxyy_x_1 = pbuffer.data(idx_npot_1_hp + 9);

    auto ta_xxxyy_y_1 = pbuffer.data(idx_npot_1_hp + 10);

    auto ta_xxxzz_x_1 = pbuffer.data(idx_npot_1_hp + 15);

    auto ta_xxxzz_z_1 = pbuffer.data(idx_npot_1_hp + 17);

    auto ta_xxyyy_x_1 = pbuffer.data(idx_npot_1_hp + 18);

    auto ta_xxyyy_y_1 = pbuffer.data(idx_npot_1_hp + 19);

    auto ta_xxzzz_x_1 = pbuffer.data(idx_npot_1_hp + 27);

    auto ta_xxzzz_z_1 = pbuffer.data(idx_npot_1_hp + 29);

    auto ta_xyyyy_x_1 = pbuffer.data(idx_npot_1_hp + 30);

    auto ta_xyyyy_y_1 = pbuffer.data(idx_npot_1_hp + 31);

    auto ta_xzzzz_x_1 = pbuffer.data(idx_npot_1_hp + 42);

    auto ta_xzzzz_z_1 = pbuffer.data(idx_npot_1_hp + 44);

    auto ta_yyyyy_x_1 = pbuffer.data(idx_npot_1_hp + 45);

    auto ta_yyyyy_y_1 = pbuffer.data(idx_npot_1_hp + 46);

    auto ta_yyyyy_z_1 = pbuffer.data(idx_npot_1_hp + 47);

    auto ta_yyyyz_y_1 = pbuffer.data(idx_npot_1_hp + 49);

    auto ta_yyyyz_z_1 = pbuffer.data(idx_npot_1_hp + 50);

    auto ta_yyyzz_y_1 = pbuffer.data(idx_npot_1_hp + 52);

    auto ta_yyyzz_z_1 = pbuffer.data(idx_npot_1_hp + 53);

    auto ta_yyzzz_y_1 = pbuffer.data(idx_npot_1_hp + 55);

    auto ta_yyzzz_z_1 = pbuffer.data(idx_npot_1_hp + 56);

    auto ta_yzzzz_y_1 = pbuffer.data(idx_npot_1_hp + 58);

    auto ta_yzzzz_z_1 = pbuffer.data(idx_npot_1_hp + 59);

    auto ta_zzzzz_x_1 = pbuffer.data(idx_npot_1_hp + 60);

    auto ta_zzzzz_y_1 = pbuffer.data(idx_npot_1_hp + 61);

    auto ta_zzzzz_z_1 = pbuffer.data(idx_npot_1_hp + 62);

    // Set up components of auxiliary buffer : HP

    auto ta1_x_xxxxx_x_0 = pbuffer.data(idx_npot_geom_010_0_hp);

    auto ta1_x_xxxxx_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 1);

    auto ta1_x_xxxxx_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 2);

    auto ta1_x_xxxxy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 3);

    auto ta1_x_xxxxy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 4);

    auto ta1_x_xxxxy_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 5);

    auto ta1_x_xxxxz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 6);

    auto ta1_x_xxxxz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 7);

    auto ta1_x_xxxxz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 8);

    auto ta1_x_xxxyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 9);

    auto ta1_x_xxxyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 10);

    auto ta1_x_xxxyy_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 11);

    auto ta1_x_xxxyz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 14);

    auto ta1_x_xxxzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 15);

    auto ta1_x_xxxzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 16);

    auto ta1_x_xxxzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 17);

    auto ta1_x_xxyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 18);

    auto ta1_x_xxyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 19);

    auto ta1_x_xxyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 20);

    auto ta1_x_xxyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 22);

    auto ta1_x_xxyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 23);

    auto ta1_x_xxyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 24);

    auto ta1_x_xxyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 26);

    auto ta1_x_xxzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 27);

    auto ta1_x_xxzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 28);

    auto ta1_x_xxzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 29);

    auto ta1_x_xyyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 30);

    auto ta1_x_xyyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 31);

    auto ta1_x_xyyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 36);

    auto ta1_x_xyzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 39);

    auto ta1_x_xzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 42);

    auto ta1_x_xzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 44);

    auto ta1_x_yyyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 45);

    auto ta1_x_yyyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 46);

    auto ta1_x_yyyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 47);

    auto ta1_x_yyyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 49);

    auto ta1_x_yyyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 50);

    auto ta1_x_yyyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 51);

    auto ta1_x_yyyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 52);

    auto ta1_x_yyyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 53);

    auto ta1_x_yyzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 54);

    auto ta1_x_yyzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 55);

    auto ta1_x_yyzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 56);

    auto ta1_x_yzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 57);

    auto ta1_x_yzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 58);

    auto ta1_x_yzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 59);

    auto ta1_x_zzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 60);

    auto ta1_x_zzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 61);

    auto ta1_x_zzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 62);

    auto ta1_y_xxxxx_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 63);

    auto ta1_y_xxxxx_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 64);

    auto ta1_y_xxxxx_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 65);

    auto ta1_y_xxxxy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 66);

    auto ta1_y_xxxxy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 67);

    auto ta1_y_xxxxz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 69);

    auto ta1_y_xxxxz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 71);

    auto ta1_y_xxxyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 72);

    auto ta1_y_xxxyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 73);

    auto ta1_y_xxxyy_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 74);

    auto ta1_y_xxxzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 78);

    auto ta1_y_xxxzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 79);

    auto ta1_y_xxxzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 80);

    auto ta1_y_xxyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 81);

    auto ta1_y_xxyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 82);

    auto ta1_y_xxyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 83);

    auto ta1_y_xxyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 84);

    auto ta1_y_xxyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 86);

    auto ta1_y_xxyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 88);

    auto ta1_y_xxzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 90);

    auto ta1_y_xxzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 91);

    auto ta1_y_xxzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 92);

    auto ta1_y_xyyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 93);

    auto ta1_y_xyyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 94);

    auto ta1_y_xyyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 95);

    auto ta1_y_xyyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 98);

    auto ta1_y_xyyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 100);

    auto ta1_y_xyyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 101);

    auto ta1_y_xyzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 103);

    auto ta1_y_xzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 105);

    auto ta1_y_xzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 106);

    auto ta1_y_xzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 107);

    auto ta1_y_yyyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 108);

    auto ta1_y_yyyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 109);

    auto ta1_y_yyyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 110);

    auto ta1_y_yyyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 111);

    auto ta1_y_yyyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 112);

    auto ta1_y_yyyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 113);

    auto ta1_y_yyyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 114);

    auto ta1_y_yyyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 115);

    auto ta1_y_yyyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 116);

    auto ta1_y_yyzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 117);

    auto ta1_y_yyzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 118);

    auto ta1_y_yyzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 119);

    auto ta1_y_yzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 121);

    auto ta1_y_yzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 122);

    auto ta1_y_zzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 123);

    auto ta1_y_zzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 124);

    auto ta1_y_zzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 125);

    auto ta1_z_xxxxx_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 126);

    auto ta1_z_xxxxx_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 127);

    auto ta1_z_xxxxx_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 128);

    auto ta1_z_xxxxy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 129);

    auto ta1_z_xxxxy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 130);

    auto ta1_z_xxxxz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 132);

    auto ta1_z_xxxxz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 134);

    auto ta1_z_xxxyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 135);

    auto ta1_z_xxxyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 136);

    auto ta1_z_xxxyy_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 137);

    auto ta1_z_xxxzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 141);

    auto ta1_z_xxxzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 142);

    auto ta1_z_xxxzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 143);

    auto ta1_z_xxyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 144);

    auto ta1_z_xxyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 145);

    auto ta1_z_xxyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 146);

    auto ta1_z_xxyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 149);

    auto ta1_z_xxyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 150);

    auto ta1_z_xxyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 151);

    auto ta1_z_xxzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 153);

    auto ta1_z_xxzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 154);

    auto ta1_z_xxzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 155);

    auto ta1_z_xyyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 156);

    auto ta1_z_xyyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 157);

    auto ta1_z_xyyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 158);

    auto ta1_z_xyyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 161);

    auto ta1_z_xyyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 163);

    auto ta1_z_xyyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 164);

    auto ta1_z_xyzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 166);

    auto ta1_z_xzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 168);

    auto ta1_z_xzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 169);

    auto ta1_z_xzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 170);

    auto ta1_z_yyyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 171);

    auto ta1_z_yyyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 172);

    auto ta1_z_yyyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 173);

    auto ta1_z_yyyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 175);

    auto ta1_z_yyyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 176);

    auto ta1_z_yyyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 177);

    auto ta1_z_yyyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 178);

    auto ta1_z_yyyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 179);

    auto ta1_z_yyzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 180);

    auto ta1_z_yyzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 181);

    auto ta1_z_yyzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 182);

    auto ta1_z_yzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 183);

    auto ta1_z_yzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 184);

    auto ta1_z_yzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 185);

    auto ta1_z_zzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 186);

    auto ta1_z_zzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 187);

    auto ta1_z_zzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 188);

    // Set up components of auxiliary buffer : HP

    auto ta1_x_xxxxx_x_1 = pbuffer.data(idx_npot_geom_010_1_hp);

    auto ta1_x_xxxxx_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 1);

    auto ta1_x_xxxxx_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 2);

    auto ta1_x_xxxxy_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 3);

    auto ta1_x_xxxxy_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 4);

    auto ta1_x_xxxxy_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 5);

    auto ta1_x_xxxxz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 6);

    auto ta1_x_xxxxz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 7);

    auto ta1_x_xxxxz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 8);

    auto ta1_x_xxxyy_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 9);

    auto ta1_x_xxxyy_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 10);

    auto ta1_x_xxxyy_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 11);

    auto ta1_x_xxxyz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 14);

    auto ta1_x_xxxzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 15);

    auto ta1_x_xxxzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 16);

    auto ta1_x_xxxzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 17);

    auto ta1_x_xxyyy_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 18);

    auto ta1_x_xxyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 19);

    auto ta1_x_xxyyy_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 20);

    auto ta1_x_xxyyz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 22);

    auto ta1_x_xxyyz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 23);

    auto ta1_x_xxyzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 24);

    auto ta1_x_xxyzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 26);

    auto ta1_x_xxzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 27);

    auto ta1_x_xxzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 28);

    auto ta1_x_xxzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 29);

    auto ta1_x_xyyyy_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 30);

    auto ta1_x_xyyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 31);

    auto ta1_x_xyyzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 36);

    auto ta1_x_xyzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 39);

    auto ta1_x_xzzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 42);

    auto ta1_x_xzzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 44);

    auto ta1_x_yyyyy_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 45);

    auto ta1_x_yyyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 46);

    auto ta1_x_yyyyy_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 47);

    auto ta1_x_yyyyz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 49);

    auto ta1_x_yyyyz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 50);

    auto ta1_x_yyyzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 51);

    auto ta1_x_yyyzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 52);

    auto ta1_x_yyyzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 53);

    auto ta1_x_yyzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 54);

    auto ta1_x_yyzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 55);

    auto ta1_x_yyzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 56);

    auto ta1_x_yzzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 57);

    auto ta1_x_yzzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 58);

    auto ta1_x_yzzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 59);

    auto ta1_x_zzzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 60);

    auto ta1_x_zzzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 61);

    auto ta1_x_zzzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 62);

    auto ta1_y_xxxxx_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 63);

    auto ta1_y_xxxxx_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 64);

    auto ta1_y_xxxxx_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 65);

    auto ta1_y_xxxxy_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 66);

    auto ta1_y_xxxxy_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 67);

    auto ta1_y_xxxxz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 69);

    auto ta1_y_xxxxz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 71);

    auto ta1_y_xxxyy_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 72);

    auto ta1_y_xxxyy_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 73);

    auto ta1_y_xxxyy_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 74);

    auto ta1_y_xxxzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 78);

    auto ta1_y_xxxzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 79);

    auto ta1_y_xxxzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 80);

    auto ta1_y_xxyyy_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 81);

    auto ta1_y_xxyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 82);

    auto ta1_y_xxyyy_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 83);

    auto ta1_y_xxyyz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 84);

    auto ta1_y_xxyyz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 86);

    auto ta1_y_xxyzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 88);

    auto ta1_y_xxzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 90);

    auto ta1_y_xxzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 91);

    auto ta1_y_xxzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 92);

    auto ta1_y_xyyyy_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 93);

    auto ta1_y_xyyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 94);

    auto ta1_y_xyyyy_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 95);

    auto ta1_y_xyyyz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 98);

    auto ta1_y_xyyzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 100);

    auto ta1_y_xyyzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 101);

    auto ta1_y_xyzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 103);

    auto ta1_y_xzzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 105);

    auto ta1_y_xzzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 106);

    auto ta1_y_xzzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 107);

    auto ta1_y_yyyyy_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 108);

    auto ta1_y_yyyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 109);

    auto ta1_y_yyyyy_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 110);

    auto ta1_y_yyyyz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 111);

    auto ta1_y_yyyyz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 112);

    auto ta1_y_yyyyz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 113);

    auto ta1_y_yyyzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 114);

    auto ta1_y_yyyzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 115);

    auto ta1_y_yyyzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 116);

    auto ta1_y_yyzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 117);

    auto ta1_y_yyzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 118);

    auto ta1_y_yyzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 119);

    auto ta1_y_yzzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 121);

    auto ta1_y_yzzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 122);

    auto ta1_y_zzzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 123);

    auto ta1_y_zzzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 124);

    auto ta1_y_zzzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 125);

    auto ta1_z_xxxxx_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 126);

    auto ta1_z_xxxxx_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 127);

    auto ta1_z_xxxxx_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 128);

    auto ta1_z_xxxxy_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 129);

    auto ta1_z_xxxxy_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 130);

    auto ta1_z_xxxxz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 132);

    auto ta1_z_xxxxz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 134);

    auto ta1_z_xxxyy_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 135);

    auto ta1_z_xxxyy_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 136);

    auto ta1_z_xxxyy_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 137);

    auto ta1_z_xxxzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 141);

    auto ta1_z_xxxzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 142);

    auto ta1_z_xxxzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 143);

    auto ta1_z_xxyyy_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 144);

    auto ta1_z_xxyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 145);

    auto ta1_z_xxyyy_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 146);

    auto ta1_z_xxyyz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 149);

    auto ta1_z_xxyzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 150);

    auto ta1_z_xxyzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 151);

    auto ta1_z_xxzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 153);

    auto ta1_z_xxzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 154);

    auto ta1_z_xxzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 155);

    auto ta1_z_xyyyy_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 156);

    auto ta1_z_xyyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 157);

    auto ta1_z_xyyyy_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 158);

    auto ta1_z_xyyyz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 161);

    auto ta1_z_xyyzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 163);

    auto ta1_z_xyyzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 164);

    auto ta1_z_xyzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 166);

    auto ta1_z_xzzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 168);

    auto ta1_z_xzzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 169);

    auto ta1_z_xzzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 170);

    auto ta1_z_yyyyy_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 171);

    auto ta1_z_yyyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 172);

    auto ta1_z_yyyyy_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 173);

    auto ta1_z_yyyyz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 175);

    auto ta1_z_yyyyz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 176);

    auto ta1_z_yyyzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 177);

    auto ta1_z_yyyzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 178);

    auto ta1_z_yyyzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 179);

    auto ta1_z_yyzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 180);

    auto ta1_z_yyzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 181);

    auto ta1_z_yyzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 182);

    auto ta1_z_yzzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 183);

    auto ta1_z_yzzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 184);

    auto ta1_z_yzzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 185);

    auto ta1_z_zzzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 186);

    auto ta1_z_zzzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 187);

    auto ta1_z_zzzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 188);

    // Set up 0-3 components of targeted buffer : IP

    auto ta1_x_xxxxxx_x_0 = pbuffer.data(idx_npot_geom_010_0_ip);

    auto ta1_x_xxxxxx_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 1);

    auto ta1_x_xxxxxx_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 2);

    #pragma omp simd aligned(pa_x, pc_x, ta1_x_xxxx_x_0, ta1_x_xxxx_x_1, ta1_x_xxxx_y_0, ta1_x_xxxx_y_1, ta1_x_xxxx_z_0, ta1_x_xxxx_z_1, ta1_x_xxxxx_0_0, ta1_x_xxxxx_0_1, ta1_x_xxxxx_x_0, ta1_x_xxxxx_x_1, ta1_x_xxxxx_y_0, ta1_x_xxxxx_y_1, ta1_x_xxxxx_z_0, ta1_x_xxxxx_z_1, ta1_x_xxxxxx_x_0, ta1_x_xxxxxx_y_0, ta1_x_xxxxxx_z_0, ta_xxxxx_x_1, ta_xxxxx_y_1, ta_xxxxx_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxxx_x_0[i] = 5.0 * ta1_x_xxxx_x_0[i] * fe_0 - 5.0 * ta1_x_xxxx_x_1[i] * fe_0 + ta1_x_xxxxx_0_0[i] * fe_0 - ta1_x_xxxxx_0_1[i] * fe_0 + ta_xxxxx_x_1[i] + ta1_x_xxxxx_x_0[i] * pa_x[i] - ta1_x_xxxxx_x_1[i] * pc_x[i];

        ta1_x_xxxxxx_y_0[i] = 5.0 * ta1_x_xxxx_y_0[i] * fe_0 - 5.0 * ta1_x_xxxx_y_1[i] * fe_0 + ta_xxxxx_y_1[i] + ta1_x_xxxxx_y_0[i] * pa_x[i] - ta1_x_xxxxx_y_1[i] * pc_x[i];

        ta1_x_xxxxxx_z_0[i] = 5.0 * ta1_x_xxxx_z_0[i] * fe_0 - 5.0 * ta1_x_xxxx_z_1[i] * fe_0 + ta_xxxxx_z_1[i] + ta1_x_xxxxx_z_0[i] * pa_x[i] - ta1_x_xxxxx_z_1[i] * pc_x[i];
    }

    // Set up 3-6 components of targeted buffer : IP

    auto ta1_x_xxxxxy_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 3);

    auto ta1_x_xxxxxy_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 4);

    auto ta1_x_xxxxxy_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 5);

    #pragma omp simd aligned(pa_y, pc_y, ta1_x_xxxxx_0_0, ta1_x_xxxxx_0_1, ta1_x_xxxxx_x_0, ta1_x_xxxxx_x_1, ta1_x_xxxxx_y_0, ta1_x_xxxxx_y_1, ta1_x_xxxxx_z_0, ta1_x_xxxxx_z_1, ta1_x_xxxxxy_x_0, ta1_x_xxxxxy_y_0, ta1_x_xxxxxy_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxxy_x_0[i] = ta1_x_xxxxx_x_0[i] * pa_y[i] - ta1_x_xxxxx_x_1[i] * pc_y[i];

        ta1_x_xxxxxy_y_0[i] = ta1_x_xxxxx_0_0[i] * fe_0 - ta1_x_xxxxx_0_1[i] * fe_0 + ta1_x_xxxxx_y_0[i] * pa_y[i] - ta1_x_xxxxx_y_1[i] * pc_y[i];

        ta1_x_xxxxxy_z_0[i] = ta1_x_xxxxx_z_0[i] * pa_y[i] - ta1_x_xxxxx_z_1[i] * pc_y[i];
    }

    // Set up 6-9 components of targeted buffer : IP

    auto ta1_x_xxxxxz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 6);

    auto ta1_x_xxxxxz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 7);

    auto ta1_x_xxxxxz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 8);

    #pragma omp simd aligned(pa_z, pc_z, ta1_x_xxxxx_0_0, ta1_x_xxxxx_0_1, ta1_x_xxxxx_x_0, ta1_x_xxxxx_x_1, ta1_x_xxxxx_y_0, ta1_x_xxxxx_y_1, ta1_x_xxxxx_z_0, ta1_x_xxxxx_z_1, ta1_x_xxxxxz_x_0, ta1_x_xxxxxz_y_0, ta1_x_xxxxxz_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxxz_x_0[i] = ta1_x_xxxxx_x_0[i] * pa_z[i] - ta1_x_xxxxx_x_1[i] * pc_z[i];

        ta1_x_xxxxxz_y_0[i] = ta1_x_xxxxx_y_0[i] * pa_z[i] - ta1_x_xxxxx_y_1[i] * pc_z[i];

        ta1_x_xxxxxz_z_0[i] = ta1_x_xxxxx_0_0[i] * fe_0 - ta1_x_xxxxx_0_1[i] * fe_0 + ta1_x_xxxxx_z_0[i] * pa_z[i] - ta1_x_xxxxx_z_1[i] * pc_z[i];
    }

    // Set up 9-12 components of targeted buffer : IP

    auto ta1_x_xxxxyy_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 9);

    auto ta1_x_xxxxyy_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 10);

    auto ta1_x_xxxxyy_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 11);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_xxxx_x_0, ta1_x_xxxx_x_1, ta1_x_xxxx_z_0, ta1_x_xxxx_z_1, ta1_x_xxxxy_x_0, ta1_x_xxxxy_x_1, ta1_x_xxxxy_z_0, ta1_x_xxxxy_z_1, ta1_x_xxxxyy_x_0, ta1_x_xxxxyy_y_0, ta1_x_xxxxyy_z_0, ta1_x_xxxyy_y_0, ta1_x_xxxyy_y_1, ta1_x_xxyy_y_0, ta1_x_xxyy_y_1, ta_xxxyy_y_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxyy_x_0[i] = ta1_x_xxxx_x_0[i] * fe_0 - ta1_x_xxxx_x_1[i] * fe_0 + ta1_x_xxxxy_x_0[i] * pa_y[i] - ta1_x_xxxxy_x_1[i] * pc_y[i];

        ta1_x_xxxxyy_y_0[i] = 3.0 * ta1_x_xxyy_y_0[i] * fe_0 - 3.0 * ta1_x_xxyy_y_1[i] * fe_0 + ta_xxxyy_y_1[i] + ta1_x_xxxyy_y_0[i] * pa_x[i] - ta1_x_xxxyy_y_1[i] * pc_x[i];

        ta1_x_xxxxyy_z_0[i] = ta1_x_xxxx_z_0[i] * fe_0 - ta1_x_xxxx_z_1[i] * fe_0 + ta1_x_xxxxy_z_0[i] * pa_y[i] - ta1_x_xxxxy_z_1[i] * pc_y[i];
    }

    // Set up 12-15 components of targeted buffer : IP

    auto ta1_x_xxxxyz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 12);

    auto ta1_x_xxxxyz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 13);

    auto ta1_x_xxxxyz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 14);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_xxxxy_y_0, ta1_x_xxxxy_y_1, ta1_x_xxxxyz_x_0, ta1_x_xxxxyz_y_0, ta1_x_xxxxyz_z_0, ta1_x_xxxxz_x_0, ta1_x_xxxxz_x_1, ta1_x_xxxxz_z_0, ta1_x_xxxxz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_x_xxxxyz_x_0[i] = ta1_x_xxxxz_x_0[i] * pa_y[i] - ta1_x_xxxxz_x_1[i] * pc_y[i];

        ta1_x_xxxxyz_y_0[i] = ta1_x_xxxxy_y_0[i] * pa_z[i] - ta1_x_xxxxy_y_1[i] * pc_z[i];

        ta1_x_xxxxyz_z_0[i] = ta1_x_xxxxz_z_0[i] * pa_y[i] - ta1_x_xxxxz_z_1[i] * pc_y[i];
    }

    // Set up 15-18 components of targeted buffer : IP

    auto ta1_x_xxxxzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 15);

    auto ta1_x_xxxxzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 16);

    auto ta1_x_xxxxzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 17);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_x_xxxx_x_0, ta1_x_xxxx_x_1, ta1_x_xxxx_y_0, ta1_x_xxxx_y_1, ta1_x_xxxxz_x_0, ta1_x_xxxxz_x_1, ta1_x_xxxxz_y_0, ta1_x_xxxxz_y_1, ta1_x_xxxxzz_x_0, ta1_x_xxxxzz_y_0, ta1_x_xxxxzz_z_0, ta1_x_xxxzz_z_0, ta1_x_xxxzz_z_1, ta1_x_xxzz_z_0, ta1_x_xxzz_z_1, ta_xxxzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxzz_x_0[i] = ta1_x_xxxx_x_0[i] * fe_0 - ta1_x_xxxx_x_1[i] * fe_0 + ta1_x_xxxxz_x_0[i] * pa_z[i] - ta1_x_xxxxz_x_1[i] * pc_z[i];

        ta1_x_xxxxzz_y_0[i] = ta1_x_xxxx_y_0[i] * fe_0 - ta1_x_xxxx_y_1[i] * fe_0 + ta1_x_xxxxz_y_0[i] * pa_z[i] - ta1_x_xxxxz_y_1[i] * pc_z[i];

        ta1_x_xxxxzz_z_0[i] = 3.0 * ta1_x_xxzz_z_0[i] * fe_0 - 3.0 * ta1_x_xxzz_z_1[i] * fe_0 + ta_xxxzz_z_1[i] + ta1_x_xxxzz_z_0[i] * pa_x[i] - ta1_x_xxxzz_z_1[i] * pc_x[i];
    }

    // Set up 18-21 components of targeted buffer : IP

    auto ta1_x_xxxyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 18);

    auto ta1_x_xxxyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 19);

    auto ta1_x_xxxyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 20);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_xxxy_x_0, ta1_x_xxxy_x_1, ta1_x_xxxy_z_0, ta1_x_xxxy_z_1, ta1_x_xxxyy_x_0, ta1_x_xxxyy_x_1, ta1_x_xxxyy_z_0, ta1_x_xxxyy_z_1, ta1_x_xxxyyy_x_0, ta1_x_xxxyyy_y_0, ta1_x_xxxyyy_z_0, ta1_x_xxyyy_y_0, ta1_x_xxyyy_y_1, ta1_x_xyyy_y_0, ta1_x_xyyy_y_1, ta_xxyyy_y_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxyyy_x_0[i] = 2.0 * ta1_x_xxxy_x_0[i] * fe_0 - 2.0 * ta1_x_xxxy_x_1[i] * fe_0 + ta1_x_xxxyy_x_0[i] * pa_y[i] - ta1_x_xxxyy_x_1[i] * pc_y[i];

        ta1_x_xxxyyy_y_0[i] = 2.0 * ta1_x_xyyy_y_0[i] * fe_0 - 2.0 * ta1_x_xyyy_y_1[i] * fe_0 + ta_xxyyy_y_1[i] + ta1_x_xxyyy_y_0[i] * pa_x[i] - ta1_x_xxyyy_y_1[i] * pc_x[i];

        ta1_x_xxxyyy_z_0[i] = 2.0 * ta1_x_xxxy_z_0[i] * fe_0 - 2.0 * ta1_x_xxxy_z_1[i] * fe_0 + ta1_x_xxxyy_z_0[i] * pa_y[i] - ta1_x_xxxyy_z_1[i] * pc_y[i];
    }

    // Set up 21-24 components of targeted buffer : IP

    auto ta1_x_xxxyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 21);

    auto ta1_x_xxxyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 22);

    auto ta1_x_xxxyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 23);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_xxxyy_x_0, ta1_x_xxxyy_x_1, ta1_x_xxxyy_y_0, ta1_x_xxxyy_y_1, ta1_x_xxxyyz_x_0, ta1_x_xxxyyz_y_0, ta1_x_xxxyyz_z_0, ta1_x_xxxyz_z_0, ta1_x_xxxyz_z_1, ta1_x_xxxz_z_0, ta1_x_xxxz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxyyz_x_0[i] = ta1_x_xxxyy_x_0[i] * pa_z[i] - ta1_x_xxxyy_x_1[i] * pc_z[i];

        ta1_x_xxxyyz_y_0[i] = ta1_x_xxxyy_y_0[i] * pa_z[i] - ta1_x_xxxyy_y_1[i] * pc_z[i];

        ta1_x_xxxyyz_z_0[i] = ta1_x_xxxz_z_0[i] * fe_0 - ta1_x_xxxz_z_1[i] * fe_0 + ta1_x_xxxyz_z_0[i] * pa_y[i] - ta1_x_xxxyz_z_1[i] * pc_y[i];
    }

    // Set up 24-27 components of targeted buffer : IP

    auto ta1_x_xxxyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 24);

    auto ta1_x_xxxyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 25);

    auto ta1_x_xxxyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 26);

    #pragma omp simd aligned(pa_y, pc_y, ta1_x_xxxyzz_x_0, ta1_x_xxxyzz_y_0, ta1_x_xxxyzz_z_0, ta1_x_xxxzz_0_0, ta1_x_xxxzz_0_1, ta1_x_xxxzz_x_0, ta1_x_xxxzz_x_1, ta1_x_xxxzz_y_0, ta1_x_xxxzz_y_1, ta1_x_xxxzz_z_0, ta1_x_xxxzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxyzz_x_0[i] = ta1_x_xxxzz_x_0[i] * pa_y[i] - ta1_x_xxxzz_x_1[i] * pc_y[i];

        ta1_x_xxxyzz_y_0[i] = ta1_x_xxxzz_0_0[i] * fe_0 - ta1_x_xxxzz_0_1[i] * fe_0 + ta1_x_xxxzz_y_0[i] * pa_y[i] - ta1_x_xxxzz_y_1[i] * pc_y[i];

        ta1_x_xxxyzz_z_0[i] = ta1_x_xxxzz_z_0[i] * pa_y[i] - ta1_x_xxxzz_z_1[i] * pc_y[i];
    }

    // Set up 27-30 components of targeted buffer : IP

    auto ta1_x_xxxzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 27);

    auto ta1_x_xxxzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 28);

    auto ta1_x_xxxzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 29);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_x_xxxz_x_0, ta1_x_xxxz_x_1, ta1_x_xxxz_y_0, ta1_x_xxxz_y_1, ta1_x_xxxzz_x_0, ta1_x_xxxzz_x_1, ta1_x_xxxzz_y_0, ta1_x_xxxzz_y_1, ta1_x_xxxzzz_x_0, ta1_x_xxxzzz_y_0, ta1_x_xxxzzz_z_0, ta1_x_xxzzz_z_0, ta1_x_xxzzz_z_1, ta1_x_xzzz_z_0, ta1_x_xzzz_z_1, ta_xxzzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxzzz_x_0[i] = 2.0 * ta1_x_xxxz_x_0[i] * fe_0 - 2.0 * ta1_x_xxxz_x_1[i] * fe_0 + ta1_x_xxxzz_x_0[i] * pa_z[i] - ta1_x_xxxzz_x_1[i] * pc_z[i];

        ta1_x_xxxzzz_y_0[i] = 2.0 * ta1_x_xxxz_y_0[i] * fe_0 - 2.0 * ta1_x_xxxz_y_1[i] * fe_0 + ta1_x_xxxzz_y_0[i] * pa_z[i] - ta1_x_xxxzz_y_1[i] * pc_z[i];

        ta1_x_xxxzzz_z_0[i] = 2.0 * ta1_x_xzzz_z_0[i] * fe_0 - 2.0 * ta1_x_xzzz_z_1[i] * fe_0 + ta_xxzzz_z_1[i] + ta1_x_xxzzz_z_0[i] * pa_x[i] - ta1_x_xxzzz_z_1[i] * pc_x[i];
    }

    // Set up 30-33 components of targeted buffer : IP

    auto ta1_x_xxyyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 30);

    auto ta1_x_xxyyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 31);

    auto ta1_x_xxyyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 32);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_xxyy_x_0, ta1_x_xxyy_x_1, ta1_x_xxyy_z_0, ta1_x_xxyy_z_1, ta1_x_xxyyy_x_0, ta1_x_xxyyy_x_1, ta1_x_xxyyy_z_0, ta1_x_xxyyy_z_1, ta1_x_xxyyyy_x_0, ta1_x_xxyyyy_y_0, ta1_x_xxyyyy_z_0, ta1_x_xyyyy_y_0, ta1_x_xyyyy_y_1, ta1_x_yyyy_y_0, ta1_x_yyyy_y_1, ta_xyyyy_y_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyyyy_x_0[i] = 3.0 * ta1_x_xxyy_x_0[i] * fe_0 - 3.0 * ta1_x_xxyy_x_1[i] * fe_0 + ta1_x_xxyyy_x_0[i] * pa_y[i] - ta1_x_xxyyy_x_1[i] * pc_y[i];

        ta1_x_xxyyyy_y_0[i] = ta1_x_yyyy_y_0[i] * fe_0 - ta1_x_yyyy_y_1[i] * fe_0 + ta_xyyyy_y_1[i] + ta1_x_xyyyy_y_0[i] * pa_x[i] - ta1_x_xyyyy_y_1[i] * pc_x[i];

        ta1_x_xxyyyy_z_0[i] = 3.0 * ta1_x_xxyy_z_0[i] * fe_0 - 3.0 * ta1_x_xxyy_z_1[i] * fe_0 + ta1_x_xxyyy_z_0[i] * pa_y[i] - ta1_x_xxyyy_z_1[i] * pc_y[i];
    }

    // Set up 33-36 components of targeted buffer : IP

    auto ta1_x_xxyyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 33);

    auto ta1_x_xxyyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 34);

    auto ta1_x_xxyyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 35);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_xxyyy_x_0, ta1_x_xxyyy_x_1, ta1_x_xxyyy_y_0, ta1_x_xxyyy_y_1, ta1_x_xxyyyz_x_0, ta1_x_xxyyyz_y_0, ta1_x_xxyyyz_z_0, ta1_x_xxyyz_z_0, ta1_x_xxyyz_z_1, ta1_x_xxyz_z_0, ta1_x_xxyz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyyyz_x_0[i] = ta1_x_xxyyy_x_0[i] * pa_z[i] - ta1_x_xxyyy_x_1[i] * pc_z[i];

        ta1_x_xxyyyz_y_0[i] = ta1_x_xxyyy_y_0[i] * pa_z[i] - ta1_x_xxyyy_y_1[i] * pc_z[i];

        ta1_x_xxyyyz_z_0[i] = 2.0 * ta1_x_xxyz_z_0[i] * fe_0 - 2.0 * ta1_x_xxyz_z_1[i] * fe_0 + ta1_x_xxyyz_z_0[i] * pa_y[i] - ta1_x_xxyyz_z_1[i] * pc_y[i];
    }

    // Set up 36-39 components of targeted buffer : IP

    auto ta1_x_xxyyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 36);

    auto ta1_x_xxyyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 37);

    auto ta1_x_xxyyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 38);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_xxyy_y_0, ta1_x_xxyy_y_1, ta1_x_xxyyz_y_0, ta1_x_xxyyz_y_1, ta1_x_xxyyzz_x_0, ta1_x_xxyyzz_y_0, ta1_x_xxyyzz_z_0, ta1_x_xxyzz_x_0, ta1_x_xxyzz_x_1, ta1_x_xxyzz_z_0, ta1_x_xxyzz_z_1, ta1_x_xxzz_x_0, ta1_x_xxzz_x_1, ta1_x_xxzz_z_0, ta1_x_xxzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyyzz_x_0[i] = ta1_x_xxzz_x_0[i] * fe_0 - ta1_x_xxzz_x_1[i] * fe_0 + ta1_x_xxyzz_x_0[i] * pa_y[i] - ta1_x_xxyzz_x_1[i] * pc_y[i];

        ta1_x_xxyyzz_y_0[i] = ta1_x_xxyy_y_0[i] * fe_0 - ta1_x_xxyy_y_1[i] * fe_0 + ta1_x_xxyyz_y_0[i] * pa_z[i] - ta1_x_xxyyz_y_1[i] * pc_z[i];

        ta1_x_xxyyzz_z_0[i] = ta1_x_xxzz_z_0[i] * fe_0 - ta1_x_xxzz_z_1[i] * fe_0 + ta1_x_xxyzz_z_0[i] * pa_y[i] - ta1_x_xxyzz_z_1[i] * pc_y[i];
    }

    // Set up 39-42 components of targeted buffer : IP

    auto ta1_x_xxyzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 39);

    auto ta1_x_xxyzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 40);

    auto ta1_x_xxyzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 41);

    #pragma omp simd aligned(pa_y, pc_y, ta1_x_xxyzzz_x_0, ta1_x_xxyzzz_y_0, ta1_x_xxyzzz_z_0, ta1_x_xxzzz_0_0, ta1_x_xxzzz_0_1, ta1_x_xxzzz_x_0, ta1_x_xxzzz_x_1, ta1_x_xxzzz_y_0, ta1_x_xxzzz_y_1, ta1_x_xxzzz_z_0, ta1_x_xxzzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyzzz_x_0[i] = ta1_x_xxzzz_x_0[i] * pa_y[i] - ta1_x_xxzzz_x_1[i] * pc_y[i];

        ta1_x_xxyzzz_y_0[i] = ta1_x_xxzzz_0_0[i] * fe_0 - ta1_x_xxzzz_0_1[i] * fe_0 + ta1_x_xxzzz_y_0[i] * pa_y[i] - ta1_x_xxzzz_y_1[i] * pc_y[i];

        ta1_x_xxyzzz_z_0[i] = ta1_x_xxzzz_z_0[i] * pa_y[i] - ta1_x_xxzzz_z_1[i] * pc_y[i];
    }

    // Set up 42-45 components of targeted buffer : IP

    auto ta1_x_xxzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 42);

    auto ta1_x_xxzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 43);

    auto ta1_x_xxzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 44);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_x_xxzz_x_0, ta1_x_xxzz_x_1, ta1_x_xxzz_y_0, ta1_x_xxzz_y_1, ta1_x_xxzzz_x_0, ta1_x_xxzzz_x_1, ta1_x_xxzzz_y_0, ta1_x_xxzzz_y_1, ta1_x_xxzzzz_x_0, ta1_x_xxzzzz_y_0, ta1_x_xxzzzz_z_0, ta1_x_xzzzz_z_0, ta1_x_xzzzz_z_1, ta1_x_zzzz_z_0, ta1_x_zzzz_z_1, ta_xzzzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxzzzz_x_0[i] = 3.0 * ta1_x_xxzz_x_0[i] * fe_0 - 3.0 * ta1_x_xxzz_x_1[i] * fe_0 + ta1_x_xxzzz_x_0[i] * pa_z[i] - ta1_x_xxzzz_x_1[i] * pc_z[i];

        ta1_x_xxzzzz_y_0[i] = 3.0 * ta1_x_xxzz_y_0[i] * fe_0 - 3.0 * ta1_x_xxzz_y_1[i] * fe_0 + ta1_x_xxzzz_y_0[i] * pa_z[i] - ta1_x_xxzzz_y_1[i] * pc_z[i];

        ta1_x_xxzzzz_z_0[i] = ta1_x_zzzz_z_0[i] * fe_0 - ta1_x_zzzz_z_1[i] * fe_0 + ta_xzzzz_z_1[i] + ta1_x_xzzzz_z_0[i] * pa_x[i] - ta1_x_xzzzz_z_1[i] * pc_x[i];
    }

    // Set up 45-48 components of targeted buffer : IP

    auto ta1_x_xyyyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 45);

    auto ta1_x_xyyyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 46);

    auto ta1_x_xyyyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 47);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_xyyy_x_0, ta1_x_xyyy_x_1, ta1_x_xyyyy_x_0, ta1_x_xyyyy_x_1, ta1_x_xyyyyy_x_0, ta1_x_xyyyyy_y_0, ta1_x_xyyyyy_z_0, ta1_x_yyyyy_y_0, ta1_x_yyyyy_y_1, ta1_x_yyyyy_z_0, ta1_x_yyyyy_z_1, ta_yyyyy_y_1, ta_yyyyy_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyyyy_x_0[i] = 4.0 * ta1_x_xyyy_x_0[i] * fe_0 - 4.0 * ta1_x_xyyy_x_1[i] * fe_0 + ta1_x_xyyyy_x_0[i] * pa_y[i] - ta1_x_xyyyy_x_1[i] * pc_y[i];

        ta1_x_xyyyyy_y_0[i] = ta_yyyyy_y_1[i] + ta1_x_yyyyy_y_0[i] * pa_x[i] - ta1_x_yyyyy_y_1[i] * pc_x[i];

        ta1_x_xyyyyy_z_0[i] = ta_yyyyy_z_1[i] + ta1_x_yyyyy_z_0[i] * pa_x[i] - ta1_x_yyyyy_z_1[i] * pc_x[i];
    }

    // Set up 48-51 components of targeted buffer : IP

    auto ta1_x_xyyyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 48);

    auto ta1_x_xyyyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 49);

    auto ta1_x_xyyyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 50);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_x_xyyyy_x_0, ta1_x_xyyyy_x_1, ta1_x_xyyyy_y_0, ta1_x_xyyyy_y_1, ta1_x_xyyyyz_x_0, ta1_x_xyyyyz_y_0, ta1_x_xyyyyz_z_0, ta1_x_yyyyz_z_0, ta1_x_yyyyz_z_1, ta_yyyyz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_x_xyyyyz_x_0[i] = ta1_x_xyyyy_x_0[i] * pa_z[i] - ta1_x_xyyyy_x_1[i] * pc_z[i];

        ta1_x_xyyyyz_y_0[i] = ta1_x_xyyyy_y_0[i] * pa_z[i] - ta1_x_xyyyy_y_1[i] * pc_z[i];

        ta1_x_xyyyyz_z_0[i] = ta_yyyyz_z_1[i] + ta1_x_yyyyz_z_0[i] * pa_x[i] - ta1_x_yyyyz_z_1[i] * pc_x[i];
    }

    // Set up 51-54 components of targeted buffer : IP

    auto ta1_x_xyyyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 51);

    auto ta1_x_xyyyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 52);

    auto ta1_x_xyyyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 53);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_xyyyzz_x_0, ta1_x_xyyyzz_y_0, ta1_x_xyyyzz_z_0, ta1_x_xyyzz_x_0, ta1_x_xyyzz_x_1, ta1_x_xyzz_x_0, ta1_x_xyzz_x_1, ta1_x_yyyzz_y_0, ta1_x_yyyzz_y_1, ta1_x_yyyzz_z_0, ta1_x_yyyzz_z_1, ta_yyyzz_y_1, ta_yyyzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyyzz_x_0[i] = 2.0 * ta1_x_xyzz_x_0[i] * fe_0 - 2.0 * ta1_x_xyzz_x_1[i] * fe_0 + ta1_x_xyyzz_x_0[i] * pa_y[i] - ta1_x_xyyzz_x_1[i] * pc_y[i];

        ta1_x_xyyyzz_y_0[i] = ta_yyyzz_y_1[i] + ta1_x_yyyzz_y_0[i] * pa_x[i] - ta1_x_yyyzz_y_1[i] * pc_x[i];

        ta1_x_xyyyzz_z_0[i] = ta_yyyzz_z_1[i] + ta1_x_yyyzz_z_0[i] * pa_x[i] - ta1_x_yyyzz_z_1[i] * pc_x[i];
    }

    // Set up 54-57 components of targeted buffer : IP

    auto ta1_x_xyyzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 54);

    auto ta1_x_xyyzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 55);

    auto ta1_x_xyyzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 56);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_xyyzzz_x_0, ta1_x_xyyzzz_y_0, ta1_x_xyyzzz_z_0, ta1_x_xyzzz_x_0, ta1_x_xyzzz_x_1, ta1_x_xzzz_x_0, ta1_x_xzzz_x_1, ta1_x_yyzzz_y_0, ta1_x_yyzzz_y_1, ta1_x_yyzzz_z_0, ta1_x_yyzzz_z_1, ta_yyzzz_y_1, ta_yyzzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyzzz_x_0[i] = ta1_x_xzzz_x_0[i] * fe_0 - ta1_x_xzzz_x_1[i] * fe_0 + ta1_x_xyzzz_x_0[i] * pa_y[i] - ta1_x_xyzzz_x_1[i] * pc_y[i];

        ta1_x_xyyzzz_y_0[i] = ta_yyzzz_y_1[i] + ta1_x_yyzzz_y_0[i] * pa_x[i] - ta1_x_yyzzz_y_1[i] * pc_x[i];

        ta1_x_xyyzzz_z_0[i] = ta_yyzzz_z_1[i] + ta1_x_yyzzz_z_0[i] * pa_x[i] - ta1_x_yyzzz_z_1[i] * pc_x[i];
    }

    // Set up 57-60 components of targeted buffer : IP

    auto ta1_x_xyzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 57);

    auto ta1_x_xyzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 58);

    auto ta1_x_xyzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 59);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_xyzzzz_x_0, ta1_x_xyzzzz_y_0, ta1_x_xyzzzz_z_0, ta1_x_xzzzz_x_0, ta1_x_xzzzz_x_1, ta1_x_xzzzz_z_0, ta1_x_xzzzz_z_1, ta1_x_yzzzz_y_0, ta1_x_yzzzz_y_1, ta_yzzzz_y_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_x_xyzzzz_x_0[i] = ta1_x_xzzzz_x_0[i] * pa_y[i] - ta1_x_xzzzz_x_1[i] * pc_y[i];

        ta1_x_xyzzzz_y_0[i] = ta_yzzzz_y_1[i] + ta1_x_yzzzz_y_0[i] * pa_x[i] - ta1_x_yzzzz_y_1[i] * pc_x[i];

        ta1_x_xyzzzz_z_0[i] = ta1_x_xzzzz_z_0[i] * pa_y[i] - ta1_x_xzzzz_z_1[i] * pc_y[i];
    }

    // Set up 60-63 components of targeted buffer : IP

    auto ta1_x_xzzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 60);

    auto ta1_x_xzzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 61);

    auto ta1_x_xzzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 62);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_x_xzzz_x_0, ta1_x_xzzz_x_1, ta1_x_xzzzz_x_0, ta1_x_xzzzz_x_1, ta1_x_xzzzzz_x_0, ta1_x_xzzzzz_y_0, ta1_x_xzzzzz_z_0, ta1_x_zzzzz_y_0, ta1_x_zzzzz_y_1, ta1_x_zzzzz_z_0, ta1_x_zzzzz_z_1, ta_zzzzz_y_1, ta_zzzzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xzzzzz_x_0[i] = 4.0 * ta1_x_xzzz_x_0[i] * fe_0 - 4.0 * ta1_x_xzzz_x_1[i] * fe_0 + ta1_x_xzzzz_x_0[i] * pa_z[i] - ta1_x_xzzzz_x_1[i] * pc_z[i];

        ta1_x_xzzzzz_y_0[i] = ta_zzzzz_y_1[i] + ta1_x_zzzzz_y_0[i] * pa_x[i] - ta1_x_zzzzz_y_1[i] * pc_x[i];

        ta1_x_xzzzzz_z_0[i] = ta_zzzzz_z_1[i] + ta1_x_zzzzz_z_0[i] * pa_x[i] - ta1_x_zzzzz_z_1[i] * pc_x[i];
    }

    // Set up 63-66 components of targeted buffer : IP

    auto ta1_x_yyyyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 63);

    auto ta1_x_yyyyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 64);

    auto ta1_x_yyyyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 65);

    #pragma omp simd aligned(pa_y, pc_y, ta1_x_yyyy_x_0, ta1_x_yyyy_x_1, ta1_x_yyyy_y_0, ta1_x_yyyy_y_1, ta1_x_yyyy_z_0, ta1_x_yyyy_z_1, ta1_x_yyyyy_0_0, ta1_x_yyyyy_0_1, ta1_x_yyyyy_x_0, ta1_x_yyyyy_x_1, ta1_x_yyyyy_y_0, ta1_x_yyyyy_y_1, ta1_x_yyyyy_z_0, ta1_x_yyyyy_z_1, ta1_x_yyyyyy_x_0, ta1_x_yyyyyy_y_0, ta1_x_yyyyyy_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyyyy_x_0[i] = 5.0 * ta1_x_yyyy_x_0[i] * fe_0 - 5.0 * ta1_x_yyyy_x_1[i] * fe_0 + ta1_x_yyyyy_x_0[i] * pa_y[i] - ta1_x_yyyyy_x_1[i] * pc_y[i];

        ta1_x_yyyyyy_y_0[i] = 5.0 * ta1_x_yyyy_y_0[i] * fe_0 - 5.0 * ta1_x_yyyy_y_1[i] * fe_0 + ta1_x_yyyyy_0_0[i] * fe_0 - ta1_x_yyyyy_0_1[i] * fe_0 + ta1_x_yyyyy_y_0[i] * pa_y[i] - ta1_x_yyyyy_y_1[i] * pc_y[i];

        ta1_x_yyyyyy_z_0[i] = 5.0 * ta1_x_yyyy_z_0[i] * fe_0 - 5.0 * ta1_x_yyyy_z_1[i] * fe_0 + ta1_x_yyyyy_z_0[i] * pa_y[i] - ta1_x_yyyyy_z_1[i] * pc_y[i];
    }

    // Set up 66-69 components of targeted buffer : IP

    auto ta1_x_yyyyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 66);

    auto ta1_x_yyyyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 67);

    auto ta1_x_yyyyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 68);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_yyyyy_x_0, ta1_x_yyyyy_x_1, ta1_x_yyyyy_y_0, ta1_x_yyyyy_y_1, ta1_x_yyyyyz_x_0, ta1_x_yyyyyz_y_0, ta1_x_yyyyyz_z_0, ta1_x_yyyyz_z_0, ta1_x_yyyyz_z_1, ta1_x_yyyz_z_0, ta1_x_yyyz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyyyz_x_0[i] = ta1_x_yyyyy_x_0[i] * pa_z[i] - ta1_x_yyyyy_x_1[i] * pc_z[i];

        ta1_x_yyyyyz_y_0[i] = ta1_x_yyyyy_y_0[i] * pa_z[i] - ta1_x_yyyyy_y_1[i] * pc_z[i];

        ta1_x_yyyyyz_z_0[i] = 4.0 * ta1_x_yyyz_z_0[i] * fe_0 - 4.0 * ta1_x_yyyz_z_1[i] * fe_0 + ta1_x_yyyyz_z_0[i] * pa_y[i] - ta1_x_yyyyz_z_1[i] * pc_y[i];
    }

    // Set up 69-72 components of targeted buffer : IP

    auto ta1_x_yyyyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 69);

    auto ta1_x_yyyyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 70);

    auto ta1_x_yyyyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 71);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_yyyy_y_0, ta1_x_yyyy_y_1, ta1_x_yyyyz_y_0, ta1_x_yyyyz_y_1, ta1_x_yyyyzz_x_0, ta1_x_yyyyzz_y_0, ta1_x_yyyyzz_z_0, ta1_x_yyyzz_x_0, ta1_x_yyyzz_x_1, ta1_x_yyyzz_z_0, ta1_x_yyyzz_z_1, ta1_x_yyzz_x_0, ta1_x_yyzz_x_1, ta1_x_yyzz_z_0, ta1_x_yyzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyyzz_x_0[i] = 3.0 * ta1_x_yyzz_x_0[i] * fe_0 - 3.0 * ta1_x_yyzz_x_1[i] * fe_0 + ta1_x_yyyzz_x_0[i] * pa_y[i] - ta1_x_yyyzz_x_1[i] * pc_y[i];

        ta1_x_yyyyzz_y_0[i] = ta1_x_yyyy_y_0[i] * fe_0 - ta1_x_yyyy_y_1[i] * fe_0 + ta1_x_yyyyz_y_0[i] * pa_z[i] - ta1_x_yyyyz_y_1[i] * pc_z[i];

        ta1_x_yyyyzz_z_0[i] = 3.0 * ta1_x_yyzz_z_0[i] * fe_0 - 3.0 * ta1_x_yyzz_z_1[i] * fe_0 + ta1_x_yyyzz_z_0[i] * pa_y[i] - ta1_x_yyyzz_z_1[i] * pc_y[i];
    }

    // Set up 72-75 components of targeted buffer : IP

    auto ta1_x_yyyzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 72);

    auto ta1_x_yyyzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 73);

    auto ta1_x_yyyzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 74);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_yyyz_y_0, ta1_x_yyyz_y_1, ta1_x_yyyzz_y_0, ta1_x_yyyzz_y_1, ta1_x_yyyzzz_x_0, ta1_x_yyyzzz_y_0, ta1_x_yyyzzz_z_0, ta1_x_yyzzz_x_0, ta1_x_yyzzz_x_1, ta1_x_yyzzz_z_0, ta1_x_yyzzz_z_1, ta1_x_yzzz_x_0, ta1_x_yzzz_x_1, ta1_x_yzzz_z_0, ta1_x_yzzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyzzz_x_0[i] = 2.0 * ta1_x_yzzz_x_0[i] * fe_0 - 2.0 * ta1_x_yzzz_x_1[i] * fe_0 + ta1_x_yyzzz_x_0[i] * pa_y[i] - ta1_x_yyzzz_x_1[i] * pc_y[i];

        ta1_x_yyyzzz_y_0[i] = 2.0 * ta1_x_yyyz_y_0[i] * fe_0 - 2.0 * ta1_x_yyyz_y_1[i] * fe_0 + ta1_x_yyyzz_y_0[i] * pa_z[i] - ta1_x_yyyzz_y_1[i] * pc_z[i];

        ta1_x_yyyzzz_z_0[i] = 2.0 * ta1_x_yzzz_z_0[i] * fe_0 - 2.0 * ta1_x_yzzz_z_1[i] * fe_0 + ta1_x_yyzzz_z_0[i] * pa_y[i] - ta1_x_yyzzz_z_1[i] * pc_y[i];
    }

    // Set up 75-78 components of targeted buffer : IP

    auto ta1_x_yyzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 75);

    auto ta1_x_yyzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 76);

    auto ta1_x_yyzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 77);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_yyzz_y_0, ta1_x_yyzz_y_1, ta1_x_yyzzz_y_0, ta1_x_yyzzz_y_1, ta1_x_yyzzzz_x_0, ta1_x_yyzzzz_y_0, ta1_x_yyzzzz_z_0, ta1_x_yzzzz_x_0, ta1_x_yzzzz_x_1, ta1_x_yzzzz_z_0, ta1_x_yzzzz_z_1, ta1_x_zzzz_x_0, ta1_x_zzzz_x_1, ta1_x_zzzz_z_0, ta1_x_zzzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyzzzz_x_0[i] = ta1_x_zzzz_x_0[i] * fe_0 - ta1_x_zzzz_x_1[i] * fe_0 + ta1_x_yzzzz_x_0[i] * pa_y[i] - ta1_x_yzzzz_x_1[i] * pc_y[i];

        ta1_x_yyzzzz_y_0[i] = 3.0 * ta1_x_yyzz_y_0[i] * fe_0 - 3.0 * ta1_x_yyzz_y_1[i] * fe_0 + ta1_x_yyzzz_y_0[i] * pa_z[i] - ta1_x_yyzzz_y_1[i] * pc_z[i];

        ta1_x_yyzzzz_z_0[i] = ta1_x_zzzz_z_0[i] * fe_0 - ta1_x_zzzz_z_1[i] * fe_0 + ta1_x_yzzzz_z_0[i] * pa_y[i] - ta1_x_yzzzz_z_1[i] * pc_y[i];
    }

    // Set up 78-81 components of targeted buffer : IP

    auto ta1_x_yzzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 78);

    auto ta1_x_yzzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 79);

    auto ta1_x_yzzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 80);

    #pragma omp simd aligned(pa_y, pc_y, ta1_x_yzzzzz_x_0, ta1_x_yzzzzz_y_0, ta1_x_yzzzzz_z_0, ta1_x_zzzzz_0_0, ta1_x_zzzzz_0_1, ta1_x_zzzzz_x_0, ta1_x_zzzzz_x_1, ta1_x_zzzzz_y_0, ta1_x_zzzzz_y_1, ta1_x_zzzzz_z_0, ta1_x_zzzzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yzzzzz_x_0[i] = ta1_x_zzzzz_x_0[i] * pa_y[i] - ta1_x_zzzzz_x_1[i] * pc_y[i];

        ta1_x_yzzzzz_y_0[i] = ta1_x_zzzzz_0_0[i] * fe_0 - ta1_x_zzzzz_0_1[i] * fe_0 + ta1_x_zzzzz_y_0[i] * pa_y[i] - ta1_x_zzzzz_y_1[i] * pc_y[i];

        ta1_x_yzzzzz_z_0[i] = ta1_x_zzzzz_z_0[i] * pa_y[i] - ta1_x_zzzzz_z_1[i] * pc_y[i];
    }

    // Set up 81-84 components of targeted buffer : IP

    auto ta1_x_zzzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 81);

    auto ta1_x_zzzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 82);

    auto ta1_x_zzzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 83);

    #pragma omp simd aligned(pa_z, pc_z, ta1_x_zzzz_x_0, ta1_x_zzzz_x_1, ta1_x_zzzz_y_0, ta1_x_zzzz_y_1, ta1_x_zzzz_z_0, ta1_x_zzzz_z_1, ta1_x_zzzzz_0_0, ta1_x_zzzzz_0_1, ta1_x_zzzzz_x_0, ta1_x_zzzzz_x_1, ta1_x_zzzzz_y_0, ta1_x_zzzzz_y_1, ta1_x_zzzzz_z_0, ta1_x_zzzzz_z_1, ta1_x_zzzzzz_x_0, ta1_x_zzzzzz_y_0, ta1_x_zzzzzz_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_zzzzzz_x_0[i] = 5.0 * ta1_x_zzzz_x_0[i] * fe_0 - 5.0 * ta1_x_zzzz_x_1[i] * fe_0 + ta1_x_zzzzz_x_0[i] * pa_z[i] - ta1_x_zzzzz_x_1[i] * pc_z[i];

        ta1_x_zzzzzz_y_0[i] = 5.0 * ta1_x_zzzz_y_0[i] * fe_0 - 5.0 * ta1_x_zzzz_y_1[i] * fe_0 + ta1_x_zzzzz_y_0[i] * pa_z[i] - ta1_x_zzzzz_y_1[i] * pc_z[i];

        ta1_x_zzzzzz_z_0[i] = 5.0 * ta1_x_zzzz_z_0[i] * fe_0 - 5.0 * ta1_x_zzzz_z_1[i] * fe_0 + ta1_x_zzzzz_0_0[i] * fe_0 - ta1_x_zzzzz_0_1[i] * fe_0 + ta1_x_zzzzz_z_0[i] * pa_z[i] - ta1_x_zzzzz_z_1[i] * pc_z[i];
    }

    // Set up 84-87 components of targeted buffer : IP

    auto ta1_y_xxxxxx_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 84);

    auto ta1_y_xxxxxx_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 85);

    auto ta1_y_xxxxxx_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 86);

    #pragma omp simd aligned(pa_x, pc_x, ta1_y_xxxx_x_0, ta1_y_xxxx_x_1, ta1_y_xxxx_y_0, ta1_y_xxxx_y_1, ta1_y_xxxx_z_0, ta1_y_xxxx_z_1, ta1_y_xxxxx_0_0, ta1_y_xxxxx_0_1, ta1_y_xxxxx_x_0, ta1_y_xxxxx_x_1, ta1_y_xxxxx_y_0, ta1_y_xxxxx_y_1, ta1_y_xxxxx_z_0, ta1_y_xxxxx_z_1, ta1_y_xxxxxx_x_0, ta1_y_xxxxxx_y_0, ta1_y_xxxxxx_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxxx_x_0[i] = 5.0 * ta1_y_xxxx_x_0[i] * fe_0 - 5.0 * ta1_y_xxxx_x_1[i] * fe_0 + ta1_y_xxxxx_0_0[i] * fe_0 - ta1_y_xxxxx_0_1[i] * fe_0 + ta1_y_xxxxx_x_0[i] * pa_x[i] - ta1_y_xxxxx_x_1[i] * pc_x[i];

        ta1_y_xxxxxx_y_0[i] = 5.0 * ta1_y_xxxx_y_0[i] * fe_0 - 5.0 * ta1_y_xxxx_y_1[i] * fe_0 + ta1_y_xxxxx_y_0[i] * pa_x[i] - ta1_y_xxxxx_y_1[i] * pc_x[i];

        ta1_y_xxxxxx_z_0[i] = 5.0 * ta1_y_xxxx_z_0[i] * fe_0 - 5.0 * ta1_y_xxxx_z_1[i] * fe_0 + ta1_y_xxxxx_z_0[i] * pa_x[i] - ta1_y_xxxxx_z_1[i] * pc_x[i];
    }

    // Set up 87-90 components of targeted buffer : IP

    auto ta1_y_xxxxxy_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 87);

    auto ta1_y_xxxxxy_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 88);

    auto ta1_y_xxxxxy_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 89);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_y_xxxxx_x_0, ta1_y_xxxxx_x_1, ta1_y_xxxxx_z_0, ta1_y_xxxxx_z_1, ta1_y_xxxxxy_x_0, ta1_y_xxxxxy_y_0, ta1_y_xxxxxy_z_0, ta1_y_xxxxy_y_0, ta1_y_xxxxy_y_1, ta1_y_xxxy_y_0, ta1_y_xxxy_y_1, ta_xxxxx_x_1, ta_xxxxx_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxxy_x_0[i] = ta_xxxxx_x_1[i] + ta1_y_xxxxx_x_0[i] * pa_y[i] - ta1_y_xxxxx_x_1[i] * pc_y[i];

        ta1_y_xxxxxy_y_0[i] = 4.0 * ta1_y_xxxy_y_0[i] * fe_0 - 4.0 * ta1_y_xxxy_y_1[i] * fe_0 + ta1_y_xxxxy_y_0[i] * pa_x[i] - ta1_y_xxxxy_y_1[i] * pc_x[i];

        ta1_y_xxxxxy_z_0[i] = ta_xxxxx_z_1[i] + ta1_y_xxxxx_z_0[i] * pa_y[i] - ta1_y_xxxxx_z_1[i] * pc_y[i];
    }

    // Set up 90-93 components of targeted buffer : IP

    auto ta1_y_xxxxxz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 90);

    auto ta1_y_xxxxxz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 91);

    auto ta1_y_xxxxxz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 92);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_xxxxx_x_0, ta1_y_xxxxx_x_1, ta1_y_xxxxx_y_0, ta1_y_xxxxx_y_1, ta1_y_xxxxxz_x_0, ta1_y_xxxxxz_y_0, ta1_y_xxxxxz_z_0, ta1_y_xxxxz_z_0, ta1_y_xxxxz_z_1, ta1_y_xxxz_z_0, ta1_y_xxxz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxxz_x_0[i] = ta1_y_xxxxx_x_0[i] * pa_z[i] - ta1_y_xxxxx_x_1[i] * pc_z[i];

        ta1_y_xxxxxz_y_0[i] = ta1_y_xxxxx_y_0[i] * pa_z[i] - ta1_y_xxxxx_y_1[i] * pc_z[i];

        ta1_y_xxxxxz_z_0[i] = 4.0 * ta1_y_xxxz_z_0[i] * fe_0 - 4.0 * ta1_y_xxxz_z_1[i] * fe_0 + ta1_y_xxxxz_z_0[i] * pa_x[i] - ta1_y_xxxxz_z_1[i] * pc_x[i];
    }

    // Set up 93-96 components of targeted buffer : IP

    auto ta1_y_xxxxyy_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 93);

    auto ta1_y_xxxxyy_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 94);

    auto ta1_y_xxxxyy_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 95);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_y_xxxx_x_0, ta1_y_xxxx_x_1, ta1_y_xxxxy_x_0, ta1_y_xxxxy_x_1, ta1_y_xxxxyy_x_0, ta1_y_xxxxyy_y_0, ta1_y_xxxxyy_z_0, ta1_y_xxxyy_y_0, ta1_y_xxxyy_y_1, ta1_y_xxxyy_z_0, ta1_y_xxxyy_z_1, ta1_y_xxyy_y_0, ta1_y_xxyy_y_1, ta1_y_xxyy_z_0, ta1_y_xxyy_z_1, ta_xxxxy_x_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxyy_x_0[i] = ta1_y_xxxx_x_0[i] * fe_0 - ta1_y_xxxx_x_1[i] * fe_0 + ta_xxxxy_x_1[i] + ta1_y_xxxxy_x_0[i] * pa_y[i] - ta1_y_xxxxy_x_1[i] * pc_y[i];

        ta1_y_xxxxyy_y_0[i] = 3.0 * ta1_y_xxyy_y_0[i] * fe_0 - 3.0 * ta1_y_xxyy_y_1[i] * fe_0 + ta1_y_xxxyy_y_0[i] * pa_x[i] - ta1_y_xxxyy_y_1[i] * pc_x[i];

        ta1_y_xxxxyy_z_0[i] = 3.0 * ta1_y_xxyy_z_0[i] * fe_0 - 3.0 * ta1_y_xxyy_z_1[i] * fe_0 + ta1_y_xxxyy_z_0[i] * pa_x[i] - ta1_y_xxxyy_z_1[i] * pc_x[i];
    }

    // Set up 96-99 components of targeted buffer : IP

    auto ta1_y_xxxxyz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 96);

    auto ta1_y_xxxxyz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 97);

    auto ta1_y_xxxxyz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 98);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_y_xxxxy_x_0, ta1_y_xxxxy_x_1, ta1_y_xxxxy_y_0, ta1_y_xxxxy_y_1, ta1_y_xxxxyz_x_0, ta1_y_xxxxyz_y_0, ta1_y_xxxxyz_z_0, ta1_y_xxxxz_z_0, ta1_y_xxxxz_z_1, ta_xxxxz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_y_xxxxyz_x_0[i] = ta1_y_xxxxy_x_0[i] * pa_z[i] - ta1_y_xxxxy_x_1[i] * pc_z[i];

        ta1_y_xxxxyz_y_0[i] = ta1_y_xxxxy_y_0[i] * pa_z[i] - ta1_y_xxxxy_y_1[i] * pc_z[i];

        ta1_y_xxxxyz_z_0[i] = ta_xxxxz_z_1[i] + ta1_y_xxxxz_z_0[i] * pa_y[i] - ta1_y_xxxxz_z_1[i] * pc_y[i];
    }

    // Set up 99-102 components of targeted buffer : IP

    auto ta1_y_xxxxzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 99);

    auto ta1_y_xxxxzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 100);

    auto ta1_y_xxxxzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 101);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_xxxx_x_0, ta1_y_xxxx_x_1, ta1_y_xxxxz_x_0, ta1_y_xxxxz_x_1, ta1_y_xxxxzz_x_0, ta1_y_xxxxzz_y_0, ta1_y_xxxxzz_z_0, ta1_y_xxxzz_y_0, ta1_y_xxxzz_y_1, ta1_y_xxxzz_z_0, ta1_y_xxxzz_z_1, ta1_y_xxzz_y_0, ta1_y_xxzz_y_1, ta1_y_xxzz_z_0, ta1_y_xxzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxzz_x_0[i] = ta1_y_xxxx_x_0[i] * fe_0 - ta1_y_xxxx_x_1[i] * fe_0 + ta1_y_xxxxz_x_0[i] * pa_z[i] - ta1_y_xxxxz_x_1[i] * pc_z[i];

        ta1_y_xxxxzz_y_0[i] = 3.0 * ta1_y_xxzz_y_0[i] * fe_0 - 3.0 * ta1_y_xxzz_y_1[i] * fe_0 + ta1_y_xxxzz_y_0[i] * pa_x[i] - ta1_y_xxxzz_y_1[i] * pc_x[i];

        ta1_y_xxxxzz_z_0[i] = 3.0 * ta1_y_xxzz_z_0[i] * fe_0 - 3.0 * ta1_y_xxzz_z_1[i] * fe_0 + ta1_y_xxxzz_z_0[i] * pa_x[i] - ta1_y_xxxzz_z_1[i] * pc_x[i];
    }

    // Set up 102-105 components of targeted buffer : IP

    auto ta1_y_xxxyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 102);

    auto ta1_y_xxxyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 103);

    auto ta1_y_xxxyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 104);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_y_xxxy_x_0, ta1_y_xxxy_x_1, ta1_y_xxxyy_x_0, ta1_y_xxxyy_x_1, ta1_y_xxxyyy_x_0, ta1_y_xxxyyy_y_0, ta1_y_xxxyyy_z_0, ta1_y_xxyyy_y_0, ta1_y_xxyyy_y_1, ta1_y_xxyyy_z_0, ta1_y_xxyyy_z_1, ta1_y_xyyy_y_0, ta1_y_xyyy_y_1, ta1_y_xyyy_z_0, ta1_y_xyyy_z_1, ta_xxxyy_x_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxyyy_x_0[i] = 2.0 * ta1_y_xxxy_x_0[i] * fe_0 - 2.0 * ta1_y_xxxy_x_1[i] * fe_0 + ta_xxxyy_x_1[i] + ta1_y_xxxyy_x_0[i] * pa_y[i] - ta1_y_xxxyy_x_1[i] * pc_y[i];

        ta1_y_xxxyyy_y_0[i] = 2.0 * ta1_y_xyyy_y_0[i] * fe_0 - 2.0 * ta1_y_xyyy_y_1[i] * fe_0 + ta1_y_xxyyy_y_0[i] * pa_x[i] - ta1_y_xxyyy_y_1[i] * pc_x[i];

        ta1_y_xxxyyy_z_0[i] = 2.0 * ta1_y_xyyy_z_0[i] * fe_0 - 2.0 * ta1_y_xyyy_z_1[i] * fe_0 + ta1_y_xxyyy_z_0[i] * pa_x[i] - ta1_y_xxyyy_z_1[i] * pc_x[i];
    }

    // Set up 105-108 components of targeted buffer : IP

    auto ta1_y_xxxyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 105);

    auto ta1_y_xxxyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 106);

    auto ta1_y_xxxyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 107);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_xxxyy_x_0, ta1_y_xxxyy_x_1, ta1_y_xxxyy_y_0, ta1_y_xxxyy_y_1, ta1_y_xxxyyz_x_0, ta1_y_xxxyyz_y_0, ta1_y_xxxyyz_z_0, ta1_y_xxyyz_z_0, ta1_y_xxyyz_z_1, ta1_y_xyyz_z_0, ta1_y_xyyz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxyyz_x_0[i] = ta1_y_xxxyy_x_0[i] * pa_z[i] - ta1_y_xxxyy_x_1[i] * pc_z[i];

        ta1_y_xxxyyz_y_0[i] = ta1_y_xxxyy_y_0[i] * pa_z[i] - ta1_y_xxxyy_y_1[i] * pc_z[i];

        ta1_y_xxxyyz_z_0[i] = 2.0 * ta1_y_xyyz_z_0[i] * fe_0 - 2.0 * ta1_y_xyyz_z_1[i] * fe_0 + ta1_y_xxyyz_z_0[i] * pa_x[i] - ta1_y_xxyyz_z_1[i] * pc_x[i];
    }

    // Set up 108-111 components of targeted buffer : IP

    auto ta1_y_xxxyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 108);

    auto ta1_y_xxxyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 109);

    auto ta1_y_xxxyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 110);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_y_xxxyzz_x_0, ta1_y_xxxyzz_y_0, ta1_y_xxxyzz_z_0, ta1_y_xxxzz_x_0, ta1_y_xxxzz_x_1, ta1_y_xxxzz_z_0, ta1_y_xxxzz_z_1, ta1_y_xxyzz_y_0, ta1_y_xxyzz_y_1, ta1_y_xyzz_y_0, ta1_y_xyzz_y_1, ta_xxxzz_x_1, ta_xxxzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxyzz_x_0[i] = ta_xxxzz_x_1[i] + ta1_y_xxxzz_x_0[i] * pa_y[i] - ta1_y_xxxzz_x_1[i] * pc_y[i];

        ta1_y_xxxyzz_y_0[i] = 2.0 * ta1_y_xyzz_y_0[i] * fe_0 - 2.0 * ta1_y_xyzz_y_1[i] * fe_0 + ta1_y_xxyzz_y_0[i] * pa_x[i] - ta1_y_xxyzz_y_1[i] * pc_x[i];

        ta1_y_xxxyzz_z_0[i] = ta_xxxzz_z_1[i] + ta1_y_xxxzz_z_0[i] * pa_y[i] - ta1_y_xxxzz_z_1[i] * pc_y[i];
    }

    // Set up 111-114 components of targeted buffer : IP

    auto ta1_y_xxxzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 111);

    auto ta1_y_xxxzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 112);

    auto ta1_y_xxxzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 113);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_xxxz_x_0, ta1_y_xxxz_x_1, ta1_y_xxxzz_x_0, ta1_y_xxxzz_x_1, ta1_y_xxxzzz_x_0, ta1_y_xxxzzz_y_0, ta1_y_xxxzzz_z_0, ta1_y_xxzzz_y_0, ta1_y_xxzzz_y_1, ta1_y_xxzzz_z_0, ta1_y_xxzzz_z_1, ta1_y_xzzz_y_0, ta1_y_xzzz_y_1, ta1_y_xzzz_z_0, ta1_y_xzzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxzzz_x_0[i] = 2.0 * ta1_y_xxxz_x_0[i] * fe_0 - 2.0 * ta1_y_xxxz_x_1[i] * fe_0 + ta1_y_xxxzz_x_0[i] * pa_z[i] - ta1_y_xxxzz_x_1[i] * pc_z[i];

        ta1_y_xxxzzz_y_0[i] = 2.0 * ta1_y_xzzz_y_0[i] * fe_0 - 2.0 * ta1_y_xzzz_y_1[i] * fe_0 + ta1_y_xxzzz_y_0[i] * pa_x[i] - ta1_y_xxzzz_y_1[i] * pc_x[i];

        ta1_y_xxxzzz_z_0[i] = 2.0 * ta1_y_xzzz_z_0[i] * fe_0 - 2.0 * ta1_y_xzzz_z_1[i] * fe_0 + ta1_y_xxzzz_z_0[i] * pa_x[i] - ta1_y_xxzzz_z_1[i] * pc_x[i];
    }

    // Set up 114-117 components of targeted buffer : IP

    auto ta1_y_xxyyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 114);

    auto ta1_y_xxyyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 115);

    auto ta1_y_xxyyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 116);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_y_xxyy_x_0, ta1_y_xxyy_x_1, ta1_y_xxyyy_x_0, ta1_y_xxyyy_x_1, ta1_y_xxyyyy_x_0, ta1_y_xxyyyy_y_0, ta1_y_xxyyyy_z_0, ta1_y_xyyyy_y_0, ta1_y_xyyyy_y_1, ta1_y_xyyyy_z_0, ta1_y_xyyyy_z_1, ta1_y_yyyy_y_0, ta1_y_yyyy_y_1, ta1_y_yyyy_z_0, ta1_y_yyyy_z_1, ta_xxyyy_x_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyyyy_x_0[i] = 3.0 * ta1_y_xxyy_x_0[i] * fe_0 - 3.0 * ta1_y_xxyy_x_1[i] * fe_0 + ta_xxyyy_x_1[i] + ta1_y_xxyyy_x_0[i] * pa_y[i] - ta1_y_xxyyy_x_1[i] * pc_y[i];

        ta1_y_xxyyyy_y_0[i] = ta1_y_yyyy_y_0[i] * fe_0 - ta1_y_yyyy_y_1[i] * fe_0 + ta1_y_xyyyy_y_0[i] * pa_x[i] - ta1_y_xyyyy_y_1[i] * pc_x[i];

        ta1_y_xxyyyy_z_0[i] = ta1_y_yyyy_z_0[i] * fe_0 - ta1_y_yyyy_z_1[i] * fe_0 + ta1_y_xyyyy_z_0[i] * pa_x[i] - ta1_y_xyyyy_z_1[i] * pc_x[i];
    }

    // Set up 117-120 components of targeted buffer : IP

    auto ta1_y_xxyyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 117);

    auto ta1_y_xxyyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 118);

    auto ta1_y_xxyyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 119);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_xxyyy_x_0, ta1_y_xxyyy_x_1, ta1_y_xxyyy_y_0, ta1_y_xxyyy_y_1, ta1_y_xxyyyz_x_0, ta1_y_xxyyyz_y_0, ta1_y_xxyyyz_z_0, ta1_y_xyyyz_z_0, ta1_y_xyyyz_z_1, ta1_y_yyyz_z_0, ta1_y_yyyz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyyyz_x_0[i] = ta1_y_xxyyy_x_0[i] * pa_z[i] - ta1_y_xxyyy_x_1[i] * pc_z[i];

        ta1_y_xxyyyz_y_0[i] = ta1_y_xxyyy_y_0[i] * pa_z[i] - ta1_y_xxyyy_y_1[i] * pc_z[i];

        ta1_y_xxyyyz_z_0[i] = ta1_y_yyyz_z_0[i] * fe_0 - ta1_y_yyyz_z_1[i] * fe_0 + ta1_y_xyyyz_z_0[i] * pa_x[i] - ta1_y_xyyyz_z_1[i] * pc_x[i];
    }

    // Set up 120-123 components of targeted buffer : IP

    auto ta1_y_xxyyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 120);

    auto ta1_y_xxyyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 121);

    auto ta1_y_xxyyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 122);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_xxyy_x_0, ta1_y_xxyy_x_1, ta1_y_xxyyz_x_0, ta1_y_xxyyz_x_1, ta1_y_xxyyzz_x_0, ta1_y_xxyyzz_y_0, ta1_y_xxyyzz_z_0, ta1_y_xyyzz_y_0, ta1_y_xyyzz_y_1, ta1_y_xyyzz_z_0, ta1_y_xyyzz_z_1, ta1_y_yyzz_y_0, ta1_y_yyzz_y_1, ta1_y_yyzz_z_0, ta1_y_yyzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyyzz_x_0[i] = ta1_y_xxyy_x_0[i] * fe_0 - ta1_y_xxyy_x_1[i] * fe_0 + ta1_y_xxyyz_x_0[i] * pa_z[i] - ta1_y_xxyyz_x_1[i] * pc_z[i];

        ta1_y_xxyyzz_y_0[i] = ta1_y_yyzz_y_0[i] * fe_0 - ta1_y_yyzz_y_1[i] * fe_0 + ta1_y_xyyzz_y_0[i] * pa_x[i] - ta1_y_xyyzz_y_1[i] * pc_x[i];

        ta1_y_xxyyzz_z_0[i] = ta1_y_yyzz_z_0[i] * fe_0 - ta1_y_yyzz_z_1[i] * fe_0 + ta1_y_xyyzz_z_0[i] * pa_x[i] - ta1_y_xyyzz_z_1[i] * pc_x[i];
    }

    // Set up 123-126 components of targeted buffer : IP

    auto ta1_y_xxyzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 123);

    auto ta1_y_xxyzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 124);

    auto ta1_y_xxyzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 125);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_y_xxyzzz_x_0, ta1_y_xxyzzz_y_0, ta1_y_xxyzzz_z_0, ta1_y_xxzzz_x_0, ta1_y_xxzzz_x_1, ta1_y_xxzzz_z_0, ta1_y_xxzzz_z_1, ta1_y_xyzzz_y_0, ta1_y_xyzzz_y_1, ta1_y_yzzz_y_0, ta1_y_yzzz_y_1, ta_xxzzz_x_1, ta_xxzzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyzzz_x_0[i] = ta_xxzzz_x_1[i] + ta1_y_xxzzz_x_0[i] * pa_y[i] - ta1_y_xxzzz_x_1[i] * pc_y[i];

        ta1_y_xxyzzz_y_0[i] = ta1_y_yzzz_y_0[i] * fe_0 - ta1_y_yzzz_y_1[i] * fe_0 + ta1_y_xyzzz_y_0[i] * pa_x[i] - ta1_y_xyzzz_y_1[i] * pc_x[i];

        ta1_y_xxyzzz_z_0[i] = ta_xxzzz_z_1[i] + ta1_y_xxzzz_z_0[i] * pa_y[i] - ta1_y_xxzzz_z_1[i] * pc_y[i];
    }

    // Set up 126-129 components of targeted buffer : IP

    auto ta1_y_xxzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 126);

    auto ta1_y_xxzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 127);

    auto ta1_y_xxzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 128);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_xxzz_x_0, ta1_y_xxzz_x_1, ta1_y_xxzzz_x_0, ta1_y_xxzzz_x_1, ta1_y_xxzzzz_x_0, ta1_y_xxzzzz_y_0, ta1_y_xxzzzz_z_0, ta1_y_xzzzz_y_0, ta1_y_xzzzz_y_1, ta1_y_xzzzz_z_0, ta1_y_xzzzz_z_1, ta1_y_zzzz_y_0, ta1_y_zzzz_y_1, ta1_y_zzzz_z_0, ta1_y_zzzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxzzzz_x_0[i] = 3.0 * ta1_y_xxzz_x_0[i] * fe_0 - 3.0 * ta1_y_xxzz_x_1[i] * fe_0 + ta1_y_xxzzz_x_0[i] * pa_z[i] - ta1_y_xxzzz_x_1[i] * pc_z[i];

        ta1_y_xxzzzz_y_0[i] = ta1_y_zzzz_y_0[i] * fe_0 - ta1_y_zzzz_y_1[i] * fe_0 + ta1_y_xzzzz_y_0[i] * pa_x[i] - ta1_y_xzzzz_y_1[i] * pc_x[i];

        ta1_y_xxzzzz_z_0[i] = ta1_y_zzzz_z_0[i] * fe_0 - ta1_y_zzzz_z_1[i] * fe_0 + ta1_y_xzzzz_z_0[i] * pa_x[i] - ta1_y_xzzzz_z_1[i] * pc_x[i];
    }

    // Set up 129-132 components of targeted buffer : IP

    auto ta1_y_xyyyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 129);

    auto ta1_y_xyyyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 130);

    auto ta1_y_xyyyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 131);

    #pragma omp simd aligned(pa_x, pc_x, ta1_y_xyyyyy_x_0, ta1_y_xyyyyy_y_0, ta1_y_xyyyyy_z_0, ta1_y_yyyyy_0_0, ta1_y_yyyyy_0_1, ta1_y_yyyyy_x_0, ta1_y_yyyyy_x_1, ta1_y_yyyyy_y_0, ta1_y_yyyyy_y_1, ta1_y_yyyyy_z_0, ta1_y_yyyyy_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyyyy_x_0[i] = ta1_y_yyyyy_0_0[i] * fe_0 - ta1_y_yyyyy_0_1[i] * fe_0 + ta1_y_yyyyy_x_0[i] * pa_x[i] - ta1_y_yyyyy_x_1[i] * pc_x[i];

        ta1_y_xyyyyy_y_0[i] = ta1_y_yyyyy_y_0[i] * pa_x[i] - ta1_y_yyyyy_y_1[i] * pc_x[i];

        ta1_y_xyyyyy_z_0[i] = ta1_y_yyyyy_z_0[i] * pa_x[i] - ta1_y_yyyyy_z_1[i] * pc_x[i];
    }

    // Set up 132-135 components of targeted buffer : IP

    auto ta1_y_xyyyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 132);

    auto ta1_y_xyyyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 133);

    auto ta1_y_xyyyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 134);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_xyyyy_x_0, ta1_y_xyyyy_x_1, ta1_y_xyyyyz_x_0, ta1_y_xyyyyz_y_0, ta1_y_xyyyyz_z_0, ta1_y_yyyyz_y_0, ta1_y_yyyyz_y_1, ta1_y_yyyyz_z_0, ta1_y_yyyyz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_y_xyyyyz_x_0[i] = ta1_y_xyyyy_x_0[i] * pa_z[i] - ta1_y_xyyyy_x_1[i] * pc_z[i];

        ta1_y_xyyyyz_y_0[i] = ta1_y_yyyyz_y_0[i] * pa_x[i] - ta1_y_yyyyz_y_1[i] * pc_x[i];

        ta1_y_xyyyyz_z_0[i] = ta1_y_yyyyz_z_0[i] * pa_x[i] - ta1_y_yyyyz_z_1[i] * pc_x[i];
    }

    // Set up 135-138 components of targeted buffer : IP

    auto ta1_y_xyyyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 135);

    auto ta1_y_xyyyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 136);

    auto ta1_y_xyyyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 137);

    #pragma omp simd aligned(pa_x, pc_x, ta1_y_xyyyzz_x_0, ta1_y_xyyyzz_y_0, ta1_y_xyyyzz_z_0, ta1_y_yyyzz_0_0, ta1_y_yyyzz_0_1, ta1_y_yyyzz_x_0, ta1_y_yyyzz_x_1, ta1_y_yyyzz_y_0, ta1_y_yyyzz_y_1, ta1_y_yyyzz_z_0, ta1_y_yyyzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyyzz_x_0[i] = ta1_y_yyyzz_0_0[i] * fe_0 - ta1_y_yyyzz_0_1[i] * fe_0 + ta1_y_yyyzz_x_0[i] * pa_x[i] - ta1_y_yyyzz_x_1[i] * pc_x[i];

        ta1_y_xyyyzz_y_0[i] = ta1_y_yyyzz_y_0[i] * pa_x[i] - ta1_y_yyyzz_y_1[i] * pc_x[i];

        ta1_y_xyyyzz_z_0[i] = ta1_y_yyyzz_z_0[i] * pa_x[i] - ta1_y_yyyzz_z_1[i] * pc_x[i];
    }

    // Set up 138-141 components of targeted buffer : IP

    auto ta1_y_xyyzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 138);

    auto ta1_y_xyyzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 139);

    auto ta1_y_xyyzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 140);

    #pragma omp simd aligned(pa_x, pc_x, ta1_y_xyyzzz_x_0, ta1_y_xyyzzz_y_0, ta1_y_xyyzzz_z_0, ta1_y_yyzzz_0_0, ta1_y_yyzzz_0_1, ta1_y_yyzzz_x_0, ta1_y_yyzzz_x_1, ta1_y_yyzzz_y_0, ta1_y_yyzzz_y_1, ta1_y_yyzzz_z_0, ta1_y_yyzzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyzzz_x_0[i] = ta1_y_yyzzz_0_0[i] * fe_0 - ta1_y_yyzzz_0_1[i] * fe_0 + ta1_y_yyzzz_x_0[i] * pa_x[i] - ta1_y_yyzzz_x_1[i] * pc_x[i];

        ta1_y_xyyzzz_y_0[i] = ta1_y_yyzzz_y_0[i] * pa_x[i] - ta1_y_yyzzz_y_1[i] * pc_x[i];

        ta1_y_xyyzzz_z_0[i] = ta1_y_yyzzz_z_0[i] * pa_x[i] - ta1_y_yyzzz_z_1[i] * pc_x[i];
    }

    // Set up 141-144 components of targeted buffer : IP

    auto ta1_y_xyzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 141);

    auto ta1_y_xyzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 142);

    auto ta1_y_xyzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 143);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_y_xyzzzz_x_0, ta1_y_xyzzzz_y_0, ta1_y_xyzzzz_z_0, ta1_y_xzzzz_x_0, ta1_y_xzzzz_x_1, ta1_y_yzzzz_y_0, ta1_y_yzzzz_y_1, ta1_y_yzzzz_z_0, ta1_y_yzzzz_z_1, ta_xzzzz_x_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_y_xyzzzz_x_0[i] = ta_xzzzz_x_1[i] + ta1_y_xzzzz_x_0[i] * pa_y[i] - ta1_y_xzzzz_x_1[i] * pc_y[i];

        ta1_y_xyzzzz_y_0[i] = ta1_y_yzzzz_y_0[i] * pa_x[i] - ta1_y_yzzzz_y_1[i] * pc_x[i];

        ta1_y_xyzzzz_z_0[i] = ta1_y_yzzzz_z_0[i] * pa_x[i] - ta1_y_yzzzz_z_1[i] * pc_x[i];
    }

    // Set up 144-147 components of targeted buffer : IP

    auto ta1_y_xzzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 144);

    auto ta1_y_xzzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 145);

    auto ta1_y_xzzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 146);

    #pragma omp simd aligned(pa_x, pc_x, ta1_y_xzzzzz_x_0, ta1_y_xzzzzz_y_0, ta1_y_xzzzzz_z_0, ta1_y_zzzzz_0_0, ta1_y_zzzzz_0_1, ta1_y_zzzzz_x_0, ta1_y_zzzzz_x_1, ta1_y_zzzzz_y_0, ta1_y_zzzzz_y_1, ta1_y_zzzzz_z_0, ta1_y_zzzzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xzzzzz_x_0[i] = ta1_y_zzzzz_0_0[i] * fe_0 - ta1_y_zzzzz_0_1[i] * fe_0 + ta1_y_zzzzz_x_0[i] * pa_x[i] - ta1_y_zzzzz_x_1[i] * pc_x[i];

        ta1_y_xzzzzz_y_0[i] = ta1_y_zzzzz_y_0[i] * pa_x[i] - ta1_y_zzzzz_y_1[i] * pc_x[i];

        ta1_y_xzzzzz_z_0[i] = ta1_y_zzzzz_z_0[i] * pa_x[i] - ta1_y_zzzzz_z_1[i] * pc_x[i];
    }

    // Set up 147-150 components of targeted buffer : IP

    auto ta1_y_yyyyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 147);

    auto ta1_y_yyyyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 148);

    auto ta1_y_yyyyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 149);

    #pragma omp simd aligned(pa_y, pc_y, ta1_y_yyyy_x_0, ta1_y_yyyy_x_1, ta1_y_yyyy_y_0, ta1_y_yyyy_y_1, ta1_y_yyyy_z_0, ta1_y_yyyy_z_1, ta1_y_yyyyy_0_0, ta1_y_yyyyy_0_1, ta1_y_yyyyy_x_0, ta1_y_yyyyy_x_1, ta1_y_yyyyy_y_0, ta1_y_yyyyy_y_1, ta1_y_yyyyy_z_0, ta1_y_yyyyy_z_1, ta1_y_yyyyyy_x_0, ta1_y_yyyyyy_y_0, ta1_y_yyyyyy_z_0, ta_yyyyy_x_1, ta_yyyyy_y_1, ta_yyyyy_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyyyy_x_0[i] = 5.0 * ta1_y_yyyy_x_0[i] * fe_0 - 5.0 * ta1_y_yyyy_x_1[i] * fe_0 + ta_yyyyy_x_1[i] + ta1_y_yyyyy_x_0[i] * pa_y[i] - ta1_y_yyyyy_x_1[i] * pc_y[i];

        ta1_y_yyyyyy_y_0[i] = 5.0 * ta1_y_yyyy_y_0[i] * fe_0 - 5.0 * ta1_y_yyyy_y_1[i] * fe_0 + ta1_y_yyyyy_0_0[i] * fe_0 - ta1_y_yyyyy_0_1[i] * fe_0 + ta_yyyyy_y_1[i] + ta1_y_yyyyy_y_0[i] * pa_y[i] - ta1_y_yyyyy_y_1[i] * pc_y[i];

        ta1_y_yyyyyy_z_0[i] = 5.0 * ta1_y_yyyy_z_0[i] * fe_0 - 5.0 * ta1_y_yyyy_z_1[i] * fe_0 + ta_yyyyy_z_1[i] + ta1_y_yyyyy_z_0[i] * pa_y[i] - ta1_y_yyyyy_z_1[i] * pc_y[i];
    }

    // Set up 150-153 components of targeted buffer : IP

    auto ta1_y_yyyyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 150);

    auto ta1_y_yyyyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 151);

    auto ta1_y_yyyyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 152);

    #pragma omp simd aligned(pa_z, pc_z, ta1_y_yyyyy_0_0, ta1_y_yyyyy_0_1, ta1_y_yyyyy_x_0, ta1_y_yyyyy_x_1, ta1_y_yyyyy_y_0, ta1_y_yyyyy_y_1, ta1_y_yyyyy_z_0, ta1_y_yyyyy_z_1, ta1_y_yyyyyz_x_0, ta1_y_yyyyyz_y_0, ta1_y_yyyyyz_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyyyz_x_0[i] = ta1_y_yyyyy_x_0[i] * pa_z[i] - ta1_y_yyyyy_x_1[i] * pc_z[i];

        ta1_y_yyyyyz_y_0[i] = ta1_y_yyyyy_y_0[i] * pa_z[i] - ta1_y_yyyyy_y_1[i] * pc_z[i];

        ta1_y_yyyyyz_z_0[i] = ta1_y_yyyyy_0_0[i] * fe_0 - ta1_y_yyyyy_0_1[i] * fe_0 + ta1_y_yyyyy_z_0[i] * pa_z[i] - ta1_y_yyyyy_z_1[i] * pc_z[i];
    }

    // Set up 153-156 components of targeted buffer : IP

    auto ta1_y_yyyyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 153);

    auto ta1_y_yyyyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 154);

    auto ta1_y_yyyyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 155);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_y_yyyy_x_0, ta1_y_yyyy_x_1, ta1_y_yyyy_y_0, ta1_y_yyyy_y_1, ta1_y_yyyyz_x_0, ta1_y_yyyyz_x_1, ta1_y_yyyyz_y_0, ta1_y_yyyyz_y_1, ta1_y_yyyyzz_x_0, ta1_y_yyyyzz_y_0, ta1_y_yyyyzz_z_0, ta1_y_yyyzz_z_0, ta1_y_yyyzz_z_1, ta1_y_yyzz_z_0, ta1_y_yyzz_z_1, ta_yyyzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyyzz_x_0[i] = ta1_y_yyyy_x_0[i] * fe_0 - ta1_y_yyyy_x_1[i] * fe_0 + ta1_y_yyyyz_x_0[i] * pa_z[i] - ta1_y_yyyyz_x_1[i] * pc_z[i];

        ta1_y_yyyyzz_y_0[i] = ta1_y_yyyy_y_0[i] * fe_0 - ta1_y_yyyy_y_1[i] * fe_0 + ta1_y_yyyyz_y_0[i] * pa_z[i] - ta1_y_yyyyz_y_1[i] * pc_z[i];

        ta1_y_yyyyzz_z_0[i] = 3.0 * ta1_y_yyzz_z_0[i] * fe_0 - 3.0 * ta1_y_yyzz_z_1[i] * fe_0 + ta_yyyzz_z_1[i] + ta1_y_yyyzz_z_0[i] * pa_y[i] - ta1_y_yyyzz_z_1[i] * pc_y[i];
    }

    // Set up 156-159 components of targeted buffer : IP

    auto ta1_y_yyyzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 156);

    auto ta1_y_yyyzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 157);

    auto ta1_y_yyyzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 158);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_y_yyyz_x_0, ta1_y_yyyz_x_1, ta1_y_yyyz_y_0, ta1_y_yyyz_y_1, ta1_y_yyyzz_x_0, ta1_y_yyyzz_x_1, ta1_y_yyyzz_y_0, ta1_y_yyyzz_y_1, ta1_y_yyyzzz_x_0, ta1_y_yyyzzz_y_0, ta1_y_yyyzzz_z_0, ta1_y_yyzzz_z_0, ta1_y_yyzzz_z_1, ta1_y_yzzz_z_0, ta1_y_yzzz_z_1, ta_yyzzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyzzz_x_0[i] = 2.0 * ta1_y_yyyz_x_0[i] * fe_0 - 2.0 * ta1_y_yyyz_x_1[i] * fe_0 + ta1_y_yyyzz_x_0[i] * pa_z[i] - ta1_y_yyyzz_x_1[i] * pc_z[i];

        ta1_y_yyyzzz_y_0[i] = 2.0 * ta1_y_yyyz_y_0[i] * fe_0 - 2.0 * ta1_y_yyyz_y_1[i] * fe_0 + ta1_y_yyyzz_y_0[i] * pa_z[i] - ta1_y_yyyzz_y_1[i] * pc_z[i];

        ta1_y_yyyzzz_z_0[i] = 2.0 * ta1_y_yzzz_z_0[i] * fe_0 - 2.0 * ta1_y_yzzz_z_1[i] * fe_0 + ta_yyzzz_z_1[i] + ta1_y_yyzzz_z_0[i] * pa_y[i] - ta1_y_yyzzz_z_1[i] * pc_y[i];
    }

    // Set up 159-162 components of targeted buffer : IP

    auto ta1_y_yyzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 159);

    auto ta1_y_yyzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 160);

    auto ta1_y_yyzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 161);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_y_yyzz_x_0, ta1_y_yyzz_x_1, ta1_y_yyzz_y_0, ta1_y_yyzz_y_1, ta1_y_yyzzz_x_0, ta1_y_yyzzz_x_1, ta1_y_yyzzz_y_0, ta1_y_yyzzz_y_1, ta1_y_yyzzzz_x_0, ta1_y_yyzzzz_y_0, ta1_y_yyzzzz_z_0, ta1_y_yzzzz_z_0, ta1_y_yzzzz_z_1, ta1_y_zzzz_z_0, ta1_y_zzzz_z_1, ta_yzzzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyzzzz_x_0[i] = 3.0 * ta1_y_yyzz_x_0[i] * fe_0 - 3.0 * ta1_y_yyzz_x_1[i] * fe_0 + ta1_y_yyzzz_x_0[i] * pa_z[i] - ta1_y_yyzzz_x_1[i] * pc_z[i];

        ta1_y_yyzzzz_y_0[i] = 3.0 * ta1_y_yyzz_y_0[i] * fe_0 - 3.0 * ta1_y_yyzz_y_1[i] * fe_0 + ta1_y_yyzzz_y_0[i] * pa_z[i] - ta1_y_yyzzz_y_1[i] * pc_z[i];

        ta1_y_yyzzzz_z_0[i] = ta1_y_zzzz_z_0[i] * fe_0 - ta1_y_zzzz_z_1[i] * fe_0 + ta_yzzzz_z_1[i] + ta1_y_yzzzz_z_0[i] * pa_y[i] - ta1_y_yzzzz_z_1[i] * pc_y[i];
    }

    // Set up 162-165 components of targeted buffer : IP

    auto ta1_y_yzzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 162);

    auto ta1_y_yzzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 163);

    auto ta1_y_yzzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 164);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_y_yzzz_y_0, ta1_y_yzzz_y_1, ta1_y_yzzzz_y_0, ta1_y_yzzzz_y_1, ta1_y_yzzzzz_x_0, ta1_y_yzzzzz_y_0, ta1_y_yzzzzz_z_0, ta1_y_zzzzz_x_0, ta1_y_zzzzz_x_1, ta1_y_zzzzz_z_0, ta1_y_zzzzz_z_1, ta_zzzzz_x_1, ta_zzzzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yzzzzz_x_0[i] = ta_zzzzz_x_1[i] + ta1_y_zzzzz_x_0[i] * pa_y[i] - ta1_y_zzzzz_x_1[i] * pc_y[i];

        ta1_y_yzzzzz_y_0[i] = 4.0 * ta1_y_yzzz_y_0[i] * fe_0 - 4.0 * ta1_y_yzzz_y_1[i] * fe_0 + ta1_y_yzzzz_y_0[i] * pa_z[i] - ta1_y_yzzzz_y_1[i] * pc_z[i];

        ta1_y_yzzzzz_z_0[i] = ta_zzzzz_z_1[i] + ta1_y_zzzzz_z_0[i] * pa_y[i] - ta1_y_zzzzz_z_1[i] * pc_y[i];
    }

    // Set up 165-168 components of targeted buffer : IP

    auto ta1_y_zzzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 165);

    auto ta1_y_zzzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 166);

    auto ta1_y_zzzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 167);

    #pragma omp simd aligned(pa_z, pc_z, ta1_y_zzzz_x_0, ta1_y_zzzz_x_1, ta1_y_zzzz_y_0, ta1_y_zzzz_y_1, ta1_y_zzzz_z_0, ta1_y_zzzz_z_1, ta1_y_zzzzz_0_0, ta1_y_zzzzz_0_1, ta1_y_zzzzz_x_0, ta1_y_zzzzz_x_1, ta1_y_zzzzz_y_0, ta1_y_zzzzz_y_1, ta1_y_zzzzz_z_0, ta1_y_zzzzz_z_1, ta1_y_zzzzzz_x_0, ta1_y_zzzzzz_y_0, ta1_y_zzzzzz_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_zzzzzz_x_0[i] = 5.0 * ta1_y_zzzz_x_0[i] * fe_0 - 5.0 * ta1_y_zzzz_x_1[i] * fe_0 + ta1_y_zzzzz_x_0[i] * pa_z[i] - ta1_y_zzzzz_x_1[i] * pc_z[i];

        ta1_y_zzzzzz_y_0[i] = 5.0 * ta1_y_zzzz_y_0[i] * fe_0 - 5.0 * ta1_y_zzzz_y_1[i] * fe_0 + ta1_y_zzzzz_y_0[i] * pa_z[i] - ta1_y_zzzzz_y_1[i] * pc_z[i];

        ta1_y_zzzzzz_z_0[i] = 5.0 * ta1_y_zzzz_z_0[i] * fe_0 - 5.0 * ta1_y_zzzz_z_1[i] * fe_0 + ta1_y_zzzzz_0_0[i] * fe_0 - ta1_y_zzzzz_0_1[i] * fe_0 + ta1_y_zzzzz_z_0[i] * pa_z[i] - ta1_y_zzzzz_z_1[i] * pc_z[i];
    }

    // Set up 168-171 components of targeted buffer : IP

    auto ta1_z_xxxxxx_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 168);

    auto ta1_z_xxxxxx_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 169);

    auto ta1_z_xxxxxx_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 170);

    #pragma omp simd aligned(pa_x, pc_x, ta1_z_xxxx_x_0, ta1_z_xxxx_x_1, ta1_z_xxxx_y_0, ta1_z_xxxx_y_1, ta1_z_xxxx_z_0, ta1_z_xxxx_z_1, ta1_z_xxxxx_0_0, ta1_z_xxxxx_0_1, ta1_z_xxxxx_x_0, ta1_z_xxxxx_x_1, ta1_z_xxxxx_y_0, ta1_z_xxxxx_y_1, ta1_z_xxxxx_z_0, ta1_z_xxxxx_z_1, ta1_z_xxxxxx_x_0, ta1_z_xxxxxx_y_0, ta1_z_xxxxxx_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxxx_x_0[i] = 5.0 * ta1_z_xxxx_x_0[i] * fe_0 - 5.0 * ta1_z_xxxx_x_1[i] * fe_0 + ta1_z_xxxxx_0_0[i] * fe_0 - ta1_z_xxxxx_0_1[i] * fe_0 + ta1_z_xxxxx_x_0[i] * pa_x[i] - ta1_z_xxxxx_x_1[i] * pc_x[i];

        ta1_z_xxxxxx_y_0[i] = 5.0 * ta1_z_xxxx_y_0[i] * fe_0 - 5.0 * ta1_z_xxxx_y_1[i] * fe_0 + ta1_z_xxxxx_y_0[i] * pa_x[i] - ta1_z_xxxxx_y_1[i] * pc_x[i];

        ta1_z_xxxxxx_z_0[i] = 5.0 * ta1_z_xxxx_z_0[i] * fe_0 - 5.0 * ta1_z_xxxx_z_1[i] * fe_0 + ta1_z_xxxxx_z_0[i] * pa_x[i] - ta1_z_xxxxx_z_1[i] * pc_x[i];
    }

    // Set up 171-174 components of targeted buffer : IP

    auto ta1_z_xxxxxy_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 171);

    auto ta1_z_xxxxxy_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 172);

    auto ta1_z_xxxxxy_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 173);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_xxxxx_x_0, ta1_z_xxxxx_x_1, ta1_z_xxxxx_z_0, ta1_z_xxxxx_z_1, ta1_z_xxxxxy_x_0, ta1_z_xxxxxy_y_0, ta1_z_xxxxxy_z_0, ta1_z_xxxxy_y_0, ta1_z_xxxxy_y_1, ta1_z_xxxy_y_0, ta1_z_xxxy_y_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxxy_x_0[i] = ta1_z_xxxxx_x_0[i] * pa_y[i] - ta1_z_xxxxx_x_1[i] * pc_y[i];

        ta1_z_xxxxxy_y_0[i] = 4.0 * ta1_z_xxxy_y_0[i] * fe_0 - 4.0 * ta1_z_xxxy_y_1[i] * fe_0 + ta1_z_xxxxy_y_0[i] * pa_x[i] - ta1_z_xxxxy_y_1[i] * pc_x[i];

        ta1_z_xxxxxy_z_0[i] = ta1_z_xxxxx_z_0[i] * pa_y[i] - ta1_z_xxxxx_z_1[i] * pc_y[i];
    }

    // Set up 174-177 components of targeted buffer : IP

    auto ta1_z_xxxxxz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 174);

    auto ta1_z_xxxxxz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 175);

    auto ta1_z_xxxxxz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 176);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_z_xxxxx_x_0, ta1_z_xxxxx_x_1, ta1_z_xxxxx_y_0, ta1_z_xxxxx_y_1, ta1_z_xxxxxz_x_0, ta1_z_xxxxxz_y_0, ta1_z_xxxxxz_z_0, ta1_z_xxxxz_z_0, ta1_z_xxxxz_z_1, ta1_z_xxxz_z_0, ta1_z_xxxz_z_1, ta_xxxxx_x_1, ta_xxxxx_y_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxxz_x_0[i] = ta_xxxxx_x_1[i] + ta1_z_xxxxx_x_0[i] * pa_z[i] - ta1_z_xxxxx_x_1[i] * pc_z[i];

        ta1_z_xxxxxz_y_0[i] = ta_xxxxx_y_1[i] + ta1_z_xxxxx_y_0[i] * pa_z[i] - ta1_z_xxxxx_y_1[i] * pc_z[i];

        ta1_z_xxxxxz_z_0[i] = 4.0 * ta1_z_xxxz_z_0[i] * fe_0 - 4.0 * ta1_z_xxxz_z_1[i] * fe_0 + ta1_z_xxxxz_z_0[i] * pa_x[i] - ta1_z_xxxxz_z_1[i] * pc_x[i];
    }

    // Set up 177-180 components of targeted buffer : IP

    auto ta1_z_xxxxyy_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 177);

    auto ta1_z_xxxxyy_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 178);

    auto ta1_z_xxxxyy_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 179);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_xxxx_x_0, ta1_z_xxxx_x_1, ta1_z_xxxxy_x_0, ta1_z_xxxxy_x_1, ta1_z_xxxxyy_x_0, ta1_z_xxxxyy_y_0, ta1_z_xxxxyy_z_0, ta1_z_xxxyy_y_0, ta1_z_xxxyy_y_1, ta1_z_xxxyy_z_0, ta1_z_xxxyy_z_1, ta1_z_xxyy_y_0, ta1_z_xxyy_y_1, ta1_z_xxyy_z_0, ta1_z_xxyy_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxyy_x_0[i] = ta1_z_xxxx_x_0[i] * fe_0 - ta1_z_xxxx_x_1[i] * fe_0 + ta1_z_xxxxy_x_0[i] * pa_y[i] - ta1_z_xxxxy_x_1[i] * pc_y[i];

        ta1_z_xxxxyy_y_0[i] = 3.0 * ta1_z_xxyy_y_0[i] * fe_0 - 3.0 * ta1_z_xxyy_y_1[i] * fe_0 + ta1_z_xxxyy_y_0[i] * pa_x[i] - ta1_z_xxxyy_y_1[i] * pc_x[i];

        ta1_z_xxxxyy_z_0[i] = 3.0 * ta1_z_xxyy_z_0[i] * fe_0 - 3.0 * ta1_z_xxyy_z_1[i] * fe_0 + ta1_z_xxxyy_z_0[i] * pa_x[i] - ta1_z_xxxyy_z_1[i] * pc_x[i];
    }

    // Set up 180-183 components of targeted buffer : IP

    auto ta1_z_xxxxyz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 180);

    auto ta1_z_xxxxyz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 181);

    auto ta1_z_xxxxyz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 182);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_z_xxxxy_y_0, ta1_z_xxxxy_y_1, ta1_z_xxxxyz_x_0, ta1_z_xxxxyz_y_0, ta1_z_xxxxyz_z_0, ta1_z_xxxxz_x_0, ta1_z_xxxxz_x_1, ta1_z_xxxxz_z_0, ta1_z_xxxxz_z_1, ta_xxxxy_y_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_z_xxxxyz_x_0[i] = ta1_z_xxxxz_x_0[i] * pa_y[i] - ta1_z_xxxxz_x_1[i] * pc_y[i];

        ta1_z_xxxxyz_y_0[i] = ta_xxxxy_y_1[i] + ta1_z_xxxxy_y_0[i] * pa_z[i] - ta1_z_xxxxy_y_1[i] * pc_z[i];

        ta1_z_xxxxyz_z_0[i] = ta1_z_xxxxz_z_0[i] * pa_y[i] - ta1_z_xxxxz_z_1[i] * pc_y[i];
    }

    // Set up 183-186 components of targeted buffer : IP

    auto ta1_z_xxxxzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 183);

    auto ta1_z_xxxxzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 184);

    auto ta1_z_xxxxzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 185);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_z_xxxx_x_0, ta1_z_xxxx_x_1, ta1_z_xxxxz_x_0, ta1_z_xxxxz_x_1, ta1_z_xxxxzz_x_0, ta1_z_xxxxzz_y_0, ta1_z_xxxxzz_z_0, ta1_z_xxxzz_y_0, ta1_z_xxxzz_y_1, ta1_z_xxxzz_z_0, ta1_z_xxxzz_z_1, ta1_z_xxzz_y_0, ta1_z_xxzz_y_1, ta1_z_xxzz_z_0, ta1_z_xxzz_z_1, ta_xxxxz_x_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxzz_x_0[i] = ta1_z_xxxx_x_0[i] * fe_0 - ta1_z_xxxx_x_1[i] * fe_0 + ta_xxxxz_x_1[i] + ta1_z_xxxxz_x_0[i] * pa_z[i] - ta1_z_xxxxz_x_1[i] * pc_z[i];

        ta1_z_xxxxzz_y_0[i] = 3.0 * ta1_z_xxzz_y_0[i] * fe_0 - 3.0 * ta1_z_xxzz_y_1[i] * fe_0 + ta1_z_xxxzz_y_0[i] * pa_x[i] - ta1_z_xxxzz_y_1[i] * pc_x[i];

        ta1_z_xxxxzz_z_0[i] = 3.0 * ta1_z_xxzz_z_0[i] * fe_0 - 3.0 * ta1_z_xxzz_z_1[i] * fe_0 + ta1_z_xxxzz_z_0[i] * pa_x[i] - ta1_z_xxxzz_z_1[i] * pc_x[i];
    }

    // Set up 186-189 components of targeted buffer : IP

    auto ta1_z_xxxyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 186);

    auto ta1_z_xxxyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 187);

    auto ta1_z_xxxyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 188);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_xxxy_x_0, ta1_z_xxxy_x_1, ta1_z_xxxyy_x_0, ta1_z_xxxyy_x_1, ta1_z_xxxyyy_x_0, ta1_z_xxxyyy_y_0, ta1_z_xxxyyy_z_0, ta1_z_xxyyy_y_0, ta1_z_xxyyy_y_1, ta1_z_xxyyy_z_0, ta1_z_xxyyy_z_1, ta1_z_xyyy_y_0, ta1_z_xyyy_y_1, ta1_z_xyyy_z_0, ta1_z_xyyy_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxyyy_x_0[i] = 2.0 * ta1_z_xxxy_x_0[i] * fe_0 - 2.0 * ta1_z_xxxy_x_1[i] * fe_0 + ta1_z_xxxyy_x_0[i] * pa_y[i] - ta1_z_xxxyy_x_1[i] * pc_y[i];

        ta1_z_xxxyyy_y_0[i] = 2.0 * ta1_z_xyyy_y_0[i] * fe_0 - 2.0 * ta1_z_xyyy_y_1[i] * fe_0 + ta1_z_xxyyy_y_0[i] * pa_x[i] - ta1_z_xxyyy_y_1[i] * pc_x[i];

        ta1_z_xxxyyy_z_0[i] = 2.0 * ta1_z_xyyy_z_0[i] * fe_0 - 2.0 * ta1_z_xyyy_z_1[i] * fe_0 + ta1_z_xxyyy_z_0[i] * pa_x[i] - ta1_z_xxyyy_z_1[i] * pc_x[i];
    }

    // Set up 189-192 components of targeted buffer : IP

    auto ta1_z_xxxyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 189);

    auto ta1_z_xxxyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 190);

    auto ta1_z_xxxyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 191);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_z_xxxyy_x_0, ta1_z_xxxyy_x_1, ta1_z_xxxyy_y_0, ta1_z_xxxyy_y_1, ta1_z_xxxyyz_x_0, ta1_z_xxxyyz_y_0, ta1_z_xxxyyz_z_0, ta1_z_xxyyz_z_0, ta1_z_xxyyz_z_1, ta1_z_xyyz_z_0, ta1_z_xyyz_z_1, ta_xxxyy_x_1, ta_xxxyy_y_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxyyz_x_0[i] = ta_xxxyy_x_1[i] + ta1_z_xxxyy_x_0[i] * pa_z[i] - ta1_z_xxxyy_x_1[i] * pc_z[i];

        ta1_z_xxxyyz_y_0[i] = ta_xxxyy_y_1[i] + ta1_z_xxxyy_y_0[i] * pa_z[i] - ta1_z_xxxyy_y_1[i] * pc_z[i];

        ta1_z_xxxyyz_z_0[i] = 2.0 * ta1_z_xyyz_z_0[i] * fe_0 - 2.0 * ta1_z_xyyz_z_1[i] * fe_0 + ta1_z_xxyyz_z_0[i] * pa_x[i] - ta1_z_xxyyz_z_1[i] * pc_x[i];
    }

    // Set up 192-195 components of targeted buffer : IP

    auto ta1_z_xxxyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 192);

    auto ta1_z_xxxyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 193);

    auto ta1_z_xxxyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 194);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_xxxyzz_x_0, ta1_z_xxxyzz_y_0, ta1_z_xxxyzz_z_0, ta1_z_xxxzz_x_0, ta1_z_xxxzz_x_1, ta1_z_xxxzz_z_0, ta1_z_xxxzz_z_1, ta1_z_xxyzz_y_0, ta1_z_xxyzz_y_1, ta1_z_xyzz_y_0, ta1_z_xyzz_y_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxyzz_x_0[i] = ta1_z_xxxzz_x_0[i] * pa_y[i] - ta1_z_xxxzz_x_1[i] * pc_y[i];

        ta1_z_xxxyzz_y_0[i] = 2.0 * ta1_z_xyzz_y_0[i] * fe_0 - 2.0 * ta1_z_xyzz_y_1[i] * fe_0 + ta1_z_xxyzz_y_0[i] * pa_x[i] - ta1_z_xxyzz_y_1[i] * pc_x[i];

        ta1_z_xxxyzz_z_0[i] = ta1_z_xxxzz_z_0[i] * pa_y[i] - ta1_z_xxxzz_z_1[i] * pc_y[i];
    }

    // Set up 195-198 components of targeted buffer : IP

    auto ta1_z_xxxzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 195);

    auto ta1_z_xxxzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 196);

    auto ta1_z_xxxzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 197);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_z_xxxz_x_0, ta1_z_xxxz_x_1, ta1_z_xxxzz_x_0, ta1_z_xxxzz_x_1, ta1_z_xxxzzz_x_0, ta1_z_xxxzzz_y_0, ta1_z_xxxzzz_z_0, ta1_z_xxzzz_y_0, ta1_z_xxzzz_y_1, ta1_z_xxzzz_z_0, ta1_z_xxzzz_z_1, ta1_z_xzzz_y_0, ta1_z_xzzz_y_1, ta1_z_xzzz_z_0, ta1_z_xzzz_z_1, ta_xxxzz_x_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxzzz_x_0[i] = 2.0 * ta1_z_xxxz_x_0[i] * fe_0 - 2.0 * ta1_z_xxxz_x_1[i] * fe_0 + ta_xxxzz_x_1[i] + ta1_z_xxxzz_x_0[i] * pa_z[i] - ta1_z_xxxzz_x_1[i] * pc_z[i];

        ta1_z_xxxzzz_y_0[i] = 2.0 * ta1_z_xzzz_y_0[i] * fe_0 - 2.0 * ta1_z_xzzz_y_1[i] * fe_0 + ta1_z_xxzzz_y_0[i] * pa_x[i] - ta1_z_xxzzz_y_1[i] * pc_x[i];

        ta1_z_xxxzzz_z_0[i] = 2.0 * ta1_z_xzzz_z_0[i] * fe_0 - 2.0 * ta1_z_xzzz_z_1[i] * fe_0 + ta1_z_xxzzz_z_0[i] * pa_x[i] - ta1_z_xxzzz_z_1[i] * pc_x[i];
    }

    // Set up 198-201 components of targeted buffer : IP

    auto ta1_z_xxyyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 198);

    auto ta1_z_xxyyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 199);

    auto ta1_z_xxyyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 200);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_xxyy_x_0, ta1_z_xxyy_x_1, ta1_z_xxyyy_x_0, ta1_z_xxyyy_x_1, ta1_z_xxyyyy_x_0, ta1_z_xxyyyy_y_0, ta1_z_xxyyyy_z_0, ta1_z_xyyyy_y_0, ta1_z_xyyyy_y_1, ta1_z_xyyyy_z_0, ta1_z_xyyyy_z_1, ta1_z_yyyy_y_0, ta1_z_yyyy_y_1, ta1_z_yyyy_z_0, ta1_z_yyyy_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyyyy_x_0[i] = 3.0 * ta1_z_xxyy_x_0[i] * fe_0 - 3.0 * ta1_z_xxyy_x_1[i] * fe_0 + ta1_z_xxyyy_x_0[i] * pa_y[i] - ta1_z_xxyyy_x_1[i] * pc_y[i];

        ta1_z_xxyyyy_y_0[i] = ta1_z_yyyy_y_0[i] * fe_0 - ta1_z_yyyy_y_1[i] * fe_0 + ta1_z_xyyyy_y_0[i] * pa_x[i] - ta1_z_xyyyy_y_1[i] * pc_x[i];

        ta1_z_xxyyyy_z_0[i] = ta1_z_yyyy_z_0[i] * fe_0 - ta1_z_yyyy_z_1[i] * fe_0 + ta1_z_xyyyy_z_0[i] * pa_x[i] - ta1_z_xyyyy_z_1[i] * pc_x[i];
    }

    // Set up 201-204 components of targeted buffer : IP

    auto ta1_z_xxyyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 201);

    auto ta1_z_xxyyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 202);

    auto ta1_z_xxyyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 203);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_z_xxyyy_x_0, ta1_z_xxyyy_x_1, ta1_z_xxyyy_y_0, ta1_z_xxyyy_y_1, ta1_z_xxyyyz_x_0, ta1_z_xxyyyz_y_0, ta1_z_xxyyyz_z_0, ta1_z_xyyyz_z_0, ta1_z_xyyyz_z_1, ta1_z_yyyz_z_0, ta1_z_yyyz_z_1, ta_xxyyy_x_1, ta_xxyyy_y_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyyyz_x_0[i] = ta_xxyyy_x_1[i] + ta1_z_xxyyy_x_0[i] * pa_z[i] - ta1_z_xxyyy_x_1[i] * pc_z[i];

        ta1_z_xxyyyz_y_0[i] = ta_xxyyy_y_1[i] + ta1_z_xxyyy_y_0[i] * pa_z[i] - ta1_z_xxyyy_y_1[i] * pc_z[i];

        ta1_z_xxyyyz_z_0[i] = ta1_z_yyyz_z_0[i] * fe_0 - ta1_z_yyyz_z_1[i] * fe_0 + ta1_z_xyyyz_z_0[i] * pa_x[i] - ta1_z_xyyyz_z_1[i] * pc_x[i];
    }

    // Set up 204-207 components of targeted buffer : IP

    auto ta1_z_xxyyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 204);

    auto ta1_z_xxyyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 205);

    auto ta1_z_xxyyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 206);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_xxyyzz_x_0, ta1_z_xxyyzz_y_0, ta1_z_xxyyzz_z_0, ta1_z_xxyzz_x_0, ta1_z_xxyzz_x_1, ta1_z_xxzz_x_0, ta1_z_xxzz_x_1, ta1_z_xyyzz_y_0, ta1_z_xyyzz_y_1, ta1_z_xyyzz_z_0, ta1_z_xyyzz_z_1, ta1_z_yyzz_y_0, ta1_z_yyzz_y_1, ta1_z_yyzz_z_0, ta1_z_yyzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyyzz_x_0[i] = ta1_z_xxzz_x_0[i] * fe_0 - ta1_z_xxzz_x_1[i] * fe_0 + ta1_z_xxyzz_x_0[i] * pa_y[i] - ta1_z_xxyzz_x_1[i] * pc_y[i];

        ta1_z_xxyyzz_y_0[i] = ta1_z_yyzz_y_0[i] * fe_0 - ta1_z_yyzz_y_1[i] * fe_0 + ta1_z_xyyzz_y_0[i] * pa_x[i] - ta1_z_xyyzz_y_1[i] * pc_x[i];

        ta1_z_xxyyzz_z_0[i] = ta1_z_yyzz_z_0[i] * fe_0 - ta1_z_yyzz_z_1[i] * fe_0 + ta1_z_xyyzz_z_0[i] * pa_x[i] - ta1_z_xyyzz_z_1[i] * pc_x[i];
    }

    // Set up 207-210 components of targeted buffer : IP

    auto ta1_z_xxyzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 207);

    auto ta1_z_xxyzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 208);

    auto ta1_z_xxyzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 209);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_xxyzzz_x_0, ta1_z_xxyzzz_y_0, ta1_z_xxyzzz_z_0, ta1_z_xxzzz_x_0, ta1_z_xxzzz_x_1, ta1_z_xxzzz_z_0, ta1_z_xxzzz_z_1, ta1_z_xyzzz_y_0, ta1_z_xyzzz_y_1, ta1_z_yzzz_y_0, ta1_z_yzzz_y_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyzzz_x_0[i] = ta1_z_xxzzz_x_0[i] * pa_y[i] - ta1_z_xxzzz_x_1[i] * pc_y[i];

        ta1_z_xxyzzz_y_0[i] = ta1_z_yzzz_y_0[i] * fe_0 - ta1_z_yzzz_y_1[i] * fe_0 + ta1_z_xyzzz_y_0[i] * pa_x[i] - ta1_z_xyzzz_y_1[i] * pc_x[i];

        ta1_z_xxyzzz_z_0[i] = ta1_z_xxzzz_z_0[i] * pa_y[i] - ta1_z_xxzzz_z_1[i] * pc_y[i];
    }

    // Set up 210-213 components of targeted buffer : IP

    auto ta1_z_xxzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 210);

    auto ta1_z_xxzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 211);

    auto ta1_z_xxzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 212);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_z_xxzz_x_0, ta1_z_xxzz_x_1, ta1_z_xxzzz_x_0, ta1_z_xxzzz_x_1, ta1_z_xxzzzz_x_0, ta1_z_xxzzzz_y_0, ta1_z_xxzzzz_z_0, ta1_z_xzzzz_y_0, ta1_z_xzzzz_y_1, ta1_z_xzzzz_z_0, ta1_z_xzzzz_z_1, ta1_z_zzzz_y_0, ta1_z_zzzz_y_1, ta1_z_zzzz_z_0, ta1_z_zzzz_z_1, ta_xxzzz_x_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxzzzz_x_0[i] = 3.0 * ta1_z_xxzz_x_0[i] * fe_0 - 3.0 * ta1_z_xxzz_x_1[i] * fe_0 + ta_xxzzz_x_1[i] + ta1_z_xxzzz_x_0[i] * pa_z[i] - ta1_z_xxzzz_x_1[i] * pc_z[i];

        ta1_z_xxzzzz_y_0[i] = ta1_z_zzzz_y_0[i] * fe_0 - ta1_z_zzzz_y_1[i] * fe_0 + ta1_z_xzzzz_y_0[i] * pa_x[i] - ta1_z_xzzzz_y_1[i] * pc_x[i];

        ta1_z_xxzzzz_z_0[i] = ta1_z_zzzz_z_0[i] * fe_0 - ta1_z_zzzz_z_1[i] * fe_0 + ta1_z_xzzzz_z_0[i] * pa_x[i] - ta1_z_xzzzz_z_1[i] * pc_x[i];
    }

    // Set up 213-216 components of targeted buffer : IP

    auto ta1_z_xyyyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 213);

    auto ta1_z_xyyyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 214);

    auto ta1_z_xyyyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 215);

    #pragma omp simd aligned(pa_x, pc_x, ta1_z_xyyyyy_x_0, ta1_z_xyyyyy_y_0, ta1_z_xyyyyy_z_0, ta1_z_yyyyy_0_0, ta1_z_yyyyy_0_1, ta1_z_yyyyy_x_0, ta1_z_yyyyy_x_1, ta1_z_yyyyy_y_0, ta1_z_yyyyy_y_1, ta1_z_yyyyy_z_0, ta1_z_yyyyy_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyyyy_x_0[i] = ta1_z_yyyyy_0_0[i] * fe_0 - ta1_z_yyyyy_0_1[i] * fe_0 + ta1_z_yyyyy_x_0[i] * pa_x[i] - ta1_z_yyyyy_x_1[i] * pc_x[i];

        ta1_z_xyyyyy_y_0[i] = ta1_z_yyyyy_y_0[i] * pa_x[i] - ta1_z_yyyyy_y_1[i] * pc_x[i];

        ta1_z_xyyyyy_z_0[i] = ta1_z_yyyyy_z_0[i] * pa_x[i] - ta1_z_yyyyy_z_1[i] * pc_x[i];
    }

    // Set up 216-219 components of targeted buffer : IP

    auto ta1_z_xyyyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 216);

    auto ta1_z_xyyyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 217);

    auto ta1_z_xyyyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 218);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_z_xyyyy_x_0, ta1_z_xyyyy_x_1, ta1_z_xyyyyz_x_0, ta1_z_xyyyyz_y_0, ta1_z_xyyyyz_z_0, ta1_z_yyyyz_y_0, ta1_z_yyyyz_y_1, ta1_z_yyyyz_z_0, ta1_z_yyyyz_z_1, ta_xyyyy_x_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_z_xyyyyz_x_0[i] = ta_xyyyy_x_1[i] + ta1_z_xyyyy_x_0[i] * pa_z[i] - ta1_z_xyyyy_x_1[i] * pc_z[i];

        ta1_z_xyyyyz_y_0[i] = ta1_z_yyyyz_y_0[i] * pa_x[i] - ta1_z_yyyyz_y_1[i] * pc_x[i];

        ta1_z_xyyyyz_z_0[i] = ta1_z_yyyyz_z_0[i] * pa_x[i] - ta1_z_yyyyz_z_1[i] * pc_x[i];
    }

    // Set up 219-222 components of targeted buffer : IP

    auto ta1_z_xyyyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 219);

    auto ta1_z_xyyyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 220);

    auto ta1_z_xyyyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 221);

    #pragma omp simd aligned(pa_x, pc_x, ta1_z_xyyyzz_x_0, ta1_z_xyyyzz_y_0, ta1_z_xyyyzz_z_0, ta1_z_yyyzz_0_0, ta1_z_yyyzz_0_1, ta1_z_yyyzz_x_0, ta1_z_yyyzz_x_1, ta1_z_yyyzz_y_0, ta1_z_yyyzz_y_1, ta1_z_yyyzz_z_0, ta1_z_yyyzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyyzz_x_0[i] = ta1_z_yyyzz_0_0[i] * fe_0 - ta1_z_yyyzz_0_1[i] * fe_0 + ta1_z_yyyzz_x_0[i] * pa_x[i] - ta1_z_yyyzz_x_1[i] * pc_x[i];

        ta1_z_xyyyzz_y_0[i] = ta1_z_yyyzz_y_0[i] * pa_x[i] - ta1_z_yyyzz_y_1[i] * pc_x[i];

        ta1_z_xyyyzz_z_0[i] = ta1_z_yyyzz_z_0[i] * pa_x[i] - ta1_z_yyyzz_z_1[i] * pc_x[i];
    }

    // Set up 222-225 components of targeted buffer : IP

    auto ta1_z_xyyzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 222);

    auto ta1_z_xyyzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 223);

    auto ta1_z_xyyzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 224);

    #pragma omp simd aligned(pa_x, pc_x, ta1_z_xyyzzz_x_0, ta1_z_xyyzzz_y_0, ta1_z_xyyzzz_z_0, ta1_z_yyzzz_0_0, ta1_z_yyzzz_0_1, ta1_z_yyzzz_x_0, ta1_z_yyzzz_x_1, ta1_z_yyzzz_y_0, ta1_z_yyzzz_y_1, ta1_z_yyzzz_z_0, ta1_z_yyzzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyzzz_x_0[i] = ta1_z_yyzzz_0_0[i] * fe_0 - ta1_z_yyzzz_0_1[i] * fe_0 + ta1_z_yyzzz_x_0[i] * pa_x[i] - ta1_z_yyzzz_x_1[i] * pc_x[i];

        ta1_z_xyyzzz_y_0[i] = ta1_z_yyzzz_y_0[i] * pa_x[i] - ta1_z_yyzzz_y_1[i] * pc_x[i];

        ta1_z_xyyzzz_z_0[i] = ta1_z_yyzzz_z_0[i] * pa_x[i] - ta1_z_yyzzz_z_1[i] * pc_x[i];
    }

    // Set up 225-228 components of targeted buffer : IP

    auto ta1_z_xyzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 225);

    auto ta1_z_xyzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 226);

    auto ta1_z_xyzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 227);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_xyzzzz_x_0, ta1_z_xyzzzz_y_0, ta1_z_xyzzzz_z_0, ta1_z_xzzzz_x_0, ta1_z_xzzzz_x_1, ta1_z_yzzzz_y_0, ta1_z_yzzzz_y_1, ta1_z_yzzzz_z_0, ta1_z_yzzzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_z_xyzzzz_x_0[i] = ta1_z_xzzzz_x_0[i] * pa_y[i] - ta1_z_xzzzz_x_1[i] * pc_y[i];

        ta1_z_xyzzzz_y_0[i] = ta1_z_yzzzz_y_0[i] * pa_x[i] - ta1_z_yzzzz_y_1[i] * pc_x[i];

        ta1_z_xyzzzz_z_0[i] = ta1_z_yzzzz_z_0[i] * pa_x[i] - ta1_z_yzzzz_z_1[i] * pc_x[i];
    }

    // Set up 228-231 components of targeted buffer : IP

    auto ta1_z_xzzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 228);

    auto ta1_z_xzzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 229);

    auto ta1_z_xzzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 230);

    #pragma omp simd aligned(pa_x, pc_x, ta1_z_xzzzzz_x_0, ta1_z_xzzzzz_y_0, ta1_z_xzzzzz_z_0, ta1_z_zzzzz_0_0, ta1_z_zzzzz_0_1, ta1_z_zzzzz_x_0, ta1_z_zzzzz_x_1, ta1_z_zzzzz_y_0, ta1_z_zzzzz_y_1, ta1_z_zzzzz_z_0, ta1_z_zzzzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xzzzzz_x_0[i] = ta1_z_zzzzz_0_0[i] * fe_0 - ta1_z_zzzzz_0_1[i] * fe_0 + ta1_z_zzzzz_x_0[i] * pa_x[i] - ta1_z_zzzzz_x_1[i] * pc_x[i];

        ta1_z_xzzzzz_y_0[i] = ta1_z_zzzzz_y_0[i] * pa_x[i] - ta1_z_zzzzz_y_1[i] * pc_x[i];

        ta1_z_xzzzzz_z_0[i] = ta1_z_zzzzz_z_0[i] * pa_x[i] - ta1_z_zzzzz_z_1[i] * pc_x[i];
    }

    // Set up 231-234 components of targeted buffer : IP

    auto ta1_z_yyyyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 231);

    auto ta1_z_yyyyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 232);

    auto ta1_z_yyyyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 233);

    #pragma omp simd aligned(pa_y, pc_y, ta1_z_yyyy_x_0, ta1_z_yyyy_x_1, ta1_z_yyyy_y_0, ta1_z_yyyy_y_1, ta1_z_yyyy_z_0, ta1_z_yyyy_z_1, ta1_z_yyyyy_0_0, ta1_z_yyyyy_0_1, ta1_z_yyyyy_x_0, ta1_z_yyyyy_x_1, ta1_z_yyyyy_y_0, ta1_z_yyyyy_y_1, ta1_z_yyyyy_z_0, ta1_z_yyyyy_z_1, ta1_z_yyyyyy_x_0, ta1_z_yyyyyy_y_0, ta1_z_yyyyyy_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyyyy_x_0[i] = 5.0 * ta1_z_yyyy_x_0[i] * fe_0 - 5.0 * ta1_z_yyyy_x_1[i] * fe_0 + ta1_z_yyyyy_x_0[i] * pa_y[i] - ta1_z_yyyyy_x_1[i] * pc_y[i];

        ta1_z_yyyyyy_y_0[i] = 5.0 * ta1_z_yyyy_y_0[i] * fe_0 - 5.0 * ta1_z_yyyy_y_1[i] * fe_0 + ta1_z_yyyyy_0_0[i] * fe_0 - ta1_z_yyyyy_0_1[i] * fe_0 + ta1_z_yyyyy_y_0[i] * pa_y[i] - ta1_z_yyyyy_y_1[i] * pc_y[i];

        ta1_z_yyyyyy_z_0[i] = 5.0 * ta1_z_yyyy_z_0[i] * fe_0 - 5.0 * ta1_z_yyyy_z_1[i] * fe_0 + ta1_z_yyyyy_z_0[i] * pa_y[i] - ta1_z_yyyyy_z_1[i] * pc_y[i];
    }

    // Set up 234-237 components of targeted buffer : IP

    auto ta1_z_yyyyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 234);

    auto ta1_z_yyyyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 235);

    auto ta1_z_yyyyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 236);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_z_yyyyy_x_0, ta1_z_yyyyy_x_1, ta1_z_yyyyy_y_0, ta1_z_yyyyy_y_1, ta1_z_yyyyyz_x_0, ta1_z_yyyyyz_y_0, ta1_z_yyyyyz_z_0, ta1_z_yyyyz_z_0, ta1_z_yyyyz_z_1, ta1_z_yyyz_z_0, ta1_z_yyyz_z_1, ta_yyyyy_x_1, ta_yyyyy_y_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyyyz_x_0[i] = ta_yyyyy_x_1[i] + ta1_z_yyyyy_x_0[i] * pa_z[i] - ta1_z_yyyyy_x_1[i] * pc_z[i];

        ta1_z_yyyyyz_y_0[i] = ta_yyyyy_y_1[i] + ta1_z_yyyyy_y_0[i] * pa_z[i] - ta1_z_yyyyy_y_1[i] * pc_z[i];

        ta1_z_yyyyyz_z_0[i] = 4.0 * ta1_z_yyyz_z_0[i] * fe_0 - 4.0 * ta1_z_yyyz_z_1[i] * fe_0 + ta1_z_yyyyz_z_0[i] * pa_y[i] - ta1_z_yyyyz_z_1[i] * pc_y[i];
    }

    // Set up 237-240 components of targeted buffer : IP

    auto ta1_z_yyyyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 237);

    auto ta1_z_yyyyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 238);

    auto ta1_z_yyyyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 239);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_z_yyyy_y_0, ta1_z_yyyy_y_1, ta1_z_yyyyz_y_0, ta1_z_yyyyz_y_1, ta1_z_yyyyzz_x_0, ta1_z_yyyyzz_y_0, ta1_z_yyyyzz_z_0, ta1_z_yyyzz_x_0, ta1_z_yyyzz_x_1, ta1_z_yyyzz_z_0, ta1_z_yyyzz_z_1, ta1_z_yyzz_x_0, ta1_z_yyzz_x_1, ta1_z_yyzz_z_0, ta1_z_yyzz_z_1, ta_yyyyz_y_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyyzz_x_0[i] = 3.0 * ta1_z_yyzz_x_0[i] * fe_0 - 3.0 * ta1_z_yyzz_x_1[i] * fe_0 + ta1_z_yyyzz_x_0[i] * pa_y[i] - ta1_z_yyyzz_x_1[i] * pc_y[i];

        ta1_z_yyyyzz_y_0[i] = ta1_z_yyyy_y_0[i] * fe_0 - ta1_z_yyyy_y_1[i] * fe_0 + ta_yyyyz_y_1[i] + ta1_z_yyyyz_y_0[i] * pa_z[i] - ta1_z_yyyyz_y_1[i] * pc_z[i];

        ta1_z_yyyyzz_z_0[i] = 3.0 * ta1_z_yyzz_z_0[i] * fe_0 - 3.0 * ta1_z_yyzz_z_1[i] * fe_0 + ta1_z_yyyzz_z_0[i] * pa_y[i] - ta1_z_yyyzz_z_1[i] * pc_y[i];
    }

    // Set up 240-243 components of targeted buffer : IP

    auto ta1_z_yyyzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 240);

    auto ta1_z_yyyzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 241);

    auto ta1_z_yyyzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 242);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_z_yyyz_y_0, ta1_z_yyyz_y_1, ta1_z_yyyzz_y_0, ta1_z_yyyzz_y_1, ta1_z_yyyzzz_x_0, ta1_z_yyyzzz_y_0, ta1_z_yyyzzz_z_0, ta1_z_yyzzz_x_0, ta1_z_yyzzz_x_1, ta1_z_yyzzz_z_0, ta1_z_yyzzz_z_1, ta1_z_yzzz_x_0, ta1_z_yzzz_x_1, ta1_z_yzzz_z_0, ta1_z_yzzz_z_1, ta_yyyzz_y_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyzzz_x_0[i] = 2.0 * ta1_z_yzzz_x_0[i] * fe_0 - 2.0 * ta1_z_yzzz_x_1[i] * fe_0 + ta1_z_yyzzz_x_0[i] * pa_y[i] - ta1_z_yyzzz_x_1[i] * pc_y[i];

        ta1_z_yyyzzz_y_0[i] = 2.0 * ta1_z_yyyz_y_0[i] * fe_0 - 2.0 * ta1_z_yyyz_y_1[i] * fe_0 + ta_yyyzz_y_1[i] + ta1_z_yyyzz_y_0[i] * pa_z[i] - ta1_z_yyyzz_y_1[i] * pc_z[i];

        ta1_z_yyyzzz_z_0[i] = 2.0 * ta1_z_yzzz_z_0[i] * fe_0 - 2.0 * ta1_z_yzzz_z_1[i] * fe_0 + ta1_z_yyzzz_z_0[i] * pa_y[i] - ta1_z_yyzzz_z_1[i] * pc_y[i];
    }

    // Set up 243-246 components of targeted buffer : IP

    auto ta1_z_yyzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 243);

    auto ta1_z_yyzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 244);

    auto ta1_z_yyzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 245);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_z_yyzz_y_0, ta1_z_yyzz_y_1, ta1_z_yyzzz_y_0, ta1_z_yyzzz_y_1, ta1_z_yyzzzz_x_0, ta1_z_yyzzzz_y_0, ta1_z_yyzzzz_z_0, ta1_z_yzzzz_x_0, ta1_z_yzzzz_x_1, ta1_z_yzzzz_z_0, ta1_z_yzzzz_z_1, ta1_z_zzzz_x_0, ta1_z_zzzz_x_1, ta1_z_zzzz_z_0, ta1_z_zzzz_z_1, ta_yyzzz_y_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyzzzz_x_0[i] = ta1_z_zzzz_x_0[i] * fe_0 - ta1_z_zzzz_x_1[i] * fe_0 + ta1_z_yzzzz_x_0[i] * pa_y[i] - ta1_z_yzzzz_x_1[i] * pc_y[i];

        ta1_z_yyzzzz_y_0[i] = 3.0 * ta1_z_yyzz_y_0[i] * fe_0 - 3.0 * ta1_z_yyzz_y_1[i] * fe_0 + ta_yyzzz_y_1[i] + ta1_z_yyzzz_y_0[i] * pa_z[i] - ta1_z_yyzzz_y_1[i] * pc_z[i];

        ta1_z_yyzzzz_z_0[i] = ta1_z_zzzz_z_0[i] * fe_0 - ta1_z_zzzz_z_1[i] * fe_0 + ta1_z_yzzzz_z_0[i] * pa_y[i] - ta1_z_yzzzz_z_1[i] * pc_y[i];
    }

    // Set up 246-249 components of targeted buffer : IP

    auto ta1_z_yzzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 246);

    auto ta1_z_yzzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 247);

    auto ta1_z_yzzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 248);

    #pragma omp simd aligned(pa_y, pc_y, ta1_z_yzzzzz_x_0, ta1_z_yzzzzz_y_0, ta1_z_yzzzzz_z_0, ta1_z_zzzzz_0_0, ta1_z_zzzzz_0_1, ta1_z_zzzzz_x_0, ta1_z_zzzzz_x_1, ta1_z_zzzzz_y_0, ta1_z_zzzzz_y_1, ta1_z_zzzzz_z_0, ta1_z_zzzzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yzzzzz_x_0[i] = ta1_z_zzzzz_x_0[i] * pa_y[i] - ta1_z_zzzzz_x_1[i] * pc_y[i];

        ta1_z_yzzzzz_y_0[i] = ta1_z_zzzzz_0_0[i] * fe_0 - ta1_z_zzzzz_0_1[i] * fe_0 + ta1_z_zzzzz_y_0[i] * pa_y[i] - ta1_z_zzzzz_y_1[i] * pc_y[i];

        ta1_z_yzzzzz_z_0[i] = ta1_z_zzzzz_z_0[i] * pa_y[i] - ta1_z_zzzzz_z_1[i] * pc_y[i];
    }

    // Set up 249-252 components of targeted buffer : IP

    auto ta1_z_zzzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_ip + 249);

    auto ta1_z_zzzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_ip + 250);

    auto ta1_z_zzzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_ip + 251);

    #pragma omp simd aligned(pa_z, pc_z, ta1_z_zzzz_x_0, ta1_z_zzzz_x_1, ta1_z_zzzz_y_0, ta1_z_zzzz_y_1, ta1_z_zzzz_z_0, ta1_z_zzzz_z_1, ta1_z_zzzzz_0_0, ta1_z_zzzzz_0_1, ta1_z_zzzzz_x_0, ta1_z_zzzzz_x_1, ta1_z_zzzzz_y_0, ta1_z_zzzzz_y_1, ta1_z_zzzzz_z_0, ta1_z_zzzzz_z_1, ta1_z_zzzzzz_x_0, ta1_z_zzzzzz_y_0, ta1_z_zzzzzz_z_0, ta_zzzzz_x_1, ta_zzzzz_y_1, ta_zzzzz_z_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_zzzzzz_x_0[i] = 5.0 * ta1_z_zzzz_x_0[i] * fe_0 - 5.0 * ta1_z_zzzz_x_1[i] * fe_0 + ta_zzzzz_x_1[i] + ta1_z_zzzzz_x_0[i] * pa_z[i] - ta1_z_zzzzz_x_1[i] * pc_z[i];

        ta1_z_zzzzzz_y_0[i] = 5.0 * ta1_z_zzzz_y_0[i] * fe_0 - 5.0 * ta1_z_zzzz_y_1[i] * fe_0 + ta_zzzzz_y_1[i] + ta1_z_zzzzz_y_0[i] * pa_z[i] - ta1_z_zzzzz_y_1[i] * pc_z[i];

        ta1_z_zzzzzz_z_0[i] = 5.0 * ta1_z_zzzz_z_0[i] * fe_0 - 5.0 * ta1_z_zzzz_z_1[i] * fe_0 + ta1_z_zzzzz_0_0[i] * fe_0 - ta1_z_zzzzz_0_1[i] * fe_0 + ta_zzzzz_z_1[i] + ta1_z_zzzzz_z_0[i] * pa_z[i] - ta1_z_zzzzz_z_1[i] * pc_z[i];
    }

}

} // npotrec namespace

