#include "ProjectedCorePotentialPrimRecDGForG.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_dg_g(CSimdArray<double>& pbuffer, 
                                        const size_t idx_dg_g_0_0_0,
                                        const size_t idx_sg_g_0_0_0,
                                        const size_t idx_pg_g_0_0_0,
                                        const size_t idx_pf_f_0_0_1,
                                        const size_t idx_pg_f_0_0_1,
                                        const size_t idx_sg_g_1_0_0,
                                        const size_t idx_pg_g_1_0_0,
                                        const size_t idx_sg_d_1_0_1,
                                        const size_t idx_pg_d_1_0_1,
                                        const size_t idx_pf_p_1_1_1,
                                        const size_t idx_pg_p_1_1_1,
                                        const size_t idx_sg_s_2_1_1,
                                        const size_t idx_pg_s_2_1_1,
                                        const int p,
                                        const size_t idx_sg_g_0_0_1,
                                        const size_t idx_pg_g_0_0_1,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_b,
                                        const TPoint<double>& r_a,
                                        const double a_exp,
                                        const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents on ket side

    auto b_exps = factors.data(0);

    // Set up B center coordinates

    auto rb_x = factors.data(idx_b);

    auto rb_y = factors.data(idx_b + 1);

    auto rb_z = factors.data(idx_b + 2);

    // Set up A center coordinates

    const auto xyz = r_a.coordinates();

    const auto a_x = xyz[0];

    const auto a_y = xyz[1];

    const auto a_z = xyz[2];

    // Set up components of auxiliary buffer : SG

    auto tg_0_xxxx_g_0_0_0 = pbuffer.data(idx_sg_g_0_0_0);

    auto tg_0_xxxy_g_0_0_0 = pbuffer.data(idx_sg_g_0_0_0 + 1);

    auto tg_0_xxxz_g_0_0_0 = pbuffer.data(idx_sg_g_0_0_0 + 2);

    auto tg_0_xxyy_g_0_0_0 = pbuffer.data(idx_sg_g_0_0_0 + 3);

    auto tg_0_xxyz_g_0_0_0 = pbuffer.data(idx_sg_g_0_0_0 + 4);

    auto tg_0_xxzz_g_0_0_0 = pbuffer.data(idx_sg_g_0_0_0 + 5);

    auto tg_0_xyyy_g_0_0_0 = pbuffer.data(idx_sg_g_0_0_0 + 6);

    auto tg_0_xyyz_g_0_0_0 = pbuffer.data(idx_sg_g_0_0_0 + 7);

    auto tg_0_xyzz_g_0_0_0 = pbuffer.data(idx_sg_g_0_0_0 + 8);

    auto tg_0_xzzz_g_0_0_0 = pbuffer.data(idx_sg_g_0_0_0 + 9);

    auto tg_0_yyyy_g_0_0_0 = pbuffer.data(idx_sg_g_0_0_0 + 10);

    auto tg_0_yyyz_g_0_0_0 = pbuffer.data(idx_sg_g_0_0_0 + 11);

    auto tg_0_yyzz_g_0_0_0 = pbuffer.data(idx_sg_g_0_0_0 + 12);

    auto tg_0_yzzz_g_0_0_0 = pbuffer.data(idx_sg_g_0_0_0 + 13);

    auto tg_0_zzzz_g_0_0_0 = pbuffer.data(idx_sg_g_0_0_0 + 14);

    // Set up components of auxiliary buffer : PG

    auto tg_x_xxxx_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0);

    auto tg_x_xxxy_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 1);

    auto tg_x_xxxz_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 2);

    auto tg_x_xxyy_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 3);

    auto tg_x_xxyz_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 4);

    auto tg_x_xxzz_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 5);

    auto tg_x_xyyy_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 6);

    auto tg_x_xyyz_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 7);

    auto tg_x_xyzz_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 8);

    auto tg_x_xzzz_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 9);

    auto tg_x_yyyy_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 10);

    auto tg_x_yyyz_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 11);

    auto tg_x_yyzz_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 12);

    auto tg_x_yzzz_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 13);

    auto tg_x_zzzz_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 14);

    auto tg_y_xxxx_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 15);

    auto tg_y_xxxy_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 16);

    auto tg_y_xxxz_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 17);

    auto tg_y_xxyy_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 18);

    auto tg_y_xxyz_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 19);

    auto tg_y_xxzz_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 20);

    auto tg_y_xyyy_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 21);

    auto tg_y_xyyz_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 22);

    auto tg_y_xyzz_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 23);

    auto tg_y_xzzz_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 24);

    auto tg_y_yyyy_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 25);

    auto tg_y_yyyz_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 26);

    auto tg_y_yyzz_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 27);

    auto tg_y_yzzz_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 28);

    auto tg_y_zzzz_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 29);

    auto tg_z_xxxx_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 30);

    auto tg_z_xxxy_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 31);

    auto tg_z_xxxz_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 32);

    auto tg_z_xxyy_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 33);

    auto tg_z_xxyz_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 34);

    auto tg_z_xxzz_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 35);

    auto tg_z_xyyy_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 36);

    auto tg_z_xyyz_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 37);

    auto tg_z_xyzz_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 38);

    auto tg_z_xzzz_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 39);

    auto tg_z_yyyy_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 40);

    auto tg_z_yyyz_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 41);

    auto tg_z_yyzz_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 42);

    auto tg_z_yzzz_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 43);

    auto tg_z_zzzz_g_0_0_0 = pbuffer.data(idx_pg_g_0_0_0 + 44);

    // Set up components of auxiliary buffer : PF

    auto tg_x_xxx_f_0_0_1 = pbuffer.data(idx_pf_f_0_0_1);

    auto tg_x_xxy_f_0_0_1 = pbuffer.data(idx_pf_f_0_0_1 + 1);

    auto tg_x_xxz_f_0_0_1 = pbuffer.data(idx_pf_f_0_0_1 + 2);

    auto tg_x_xyy_f_0_0_1 = pbuffer.data(idx_pf_f_0_0_1 + 3);

    auto tg_x_xyz_f_0_0_1 = pbuffer.data(idx_pf_f_0_0_1 + 4);

    auto tg_x_xzz_f_0_0_1 = pbuffer.data(idx_pf_f_0_0_1 + 5);

    auto tg_x_yyy_f_0_0_1 = pbuffer.data(idx_pf_f_0_0_1 + 6);

    auto tg_x_yyz_f_0_0_1 = pbuffer.data(idx_pf_f_0_0_1 + 7);

    auto tg_x_yzz_f_0_0_1 = pbuffer.data(idx_pf_f_0_0_1 + 8);

    auto tg_x_zzz_f_0_0_1 = pbuffer.data(idx_pf_f_0_0_1 + 9);

    auto tg_y_xxx_f_0_0_1 = pbuffer.data(idx_pf_f_0_0_1 + 10);

    auto tg_y_xxy_f_0_0_1 = pbuffer.data(idx_pf_f_0_0_1 + 11);

    auto tg_y_xxz_f_0_0_1 = pbuffer.data(idx_pf_f_0_0_1 + 12);

    auto tg_y_xyy_f_0_0_1 = pbuffer.data(idx_pf_f_0_0_1 + 13);

    auto tg_y_xyz_f_0_0_1 = pbuffer.data(idx_pf_f_0_0_1 + 14);

    auto tg_y_xzz_f_0_0_1 = pbuffer.data(idx_pf_f_0_0_1 + 15);

    auto tg_y_yyy_f_0_0_1 = pbuffer.data(idx_pf_f_0_0_1 + 16);

    auto tg_y_yyz_f_0_0_1 = pbuffer.data(idx_pf_f_0_0_1 + 17);

    auto tg_y_yzz_f_0_0_1 = pbuffer.data(idx_pf_f_0_0_1 + 18);

    auto tg_y_zzz_f_0_0_1 = pbuffer.data(idx_pf_f_0_0_1 + 19);

    auto tg_z_xxx_f_0_0_1 = pbuffer.data(idx_pf_f_0_0_1 + 20);

    auto tg_z_xxy_f_0_0_1 = pbuffer.data(idx_pf_f_0_0_1 + 21);

    auto tg_z_xxz_f_0_0_1 = pbuffer.data(idx_pf_f_0_0_1 + 22);

    auto tg_z_xyy_f_0_0_1 = pbuffer.data(idx_pf_f_0_0_1 + 23);

    auto tg_z_xyz_f_0_0_1 = pbuffer.data(idx_pf_f_0_0_1 + 24);

    auto tg_z_xzz_f_0_0_1 = pbuffer.data(idx_pf_f_0_0_1 + 25);

    auto tg_z_yyy_f_0_0_1 = pbuffer.data(idx_pf_f_0_0_1 + 26);

    auto tg_z_yyz_f_0_0_1 = pbuffer.data(idx_pf_f_0_0_1 + 27);

    auto tg_z_yzz_f_0_0_1 = pbuffer.data(idx_pf_f_0_0_1 + 28);

    auto tg_z_zzz_f_0_0_1 = pbuffer.data(idx_pf_f_0_0_1 + 29);

    // Set up components of auxiliary buffer : PG

    auto tg_x_xxxx_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1);

    auto tg_x_xxxy_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 1);

    auto tg_x_xxxz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 2);

    auto tg_x_xxyy_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 3);

    auto tg_x_xxyz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 4);

    auto tg_x_xxzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 5);

    auto tg_x_xyyy_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 6);

    auto tg_x_xyyz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 7);

    auto tg_x_xyzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 8);

    auto tg_x_xzzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 9);

    auto tg_x_yyyy_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 10);

    auto tg_x_yyyz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 11);

    auto tg_x_yyzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 12);

    auto tg_x_yzzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 13);

    auto tg_x_zzzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 14);

    auto tg_y_xxxx_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 15);

    auto tg_y_xxxy_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 16);

    auto tg_y_xxxz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 17);

    auto tg_y_xxyy_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 18);

    auto tg_y_xxyz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 19);

    auto tg_y_xxzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 20);

    auto tg_y_xyyy_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 21);

    auto tg_y_xyyz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 22);

    auto tg_y_xyzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 23);

    auto tg_y_xzzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 24);

    auto tg_y_yyyy_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 25);

    auto tg_y_yyyz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 26);

    auto tg_y_yyzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 27);

    auto tg_y_yzzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 28);

    auto tg_y_zzzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 29);

    auto tg_z_xxxx_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 30);

    auto tg_z_xxxy_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 31);

    auto tg_z_xxxz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 32);

    auto tg_z_xxyy_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 33);

    auto tg_z_xxyz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 34);

    auto tg_z_xxzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 35);

    auto tg_z_xyyy_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 36);

    auto tg_z_xyyz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 37);

    auto tg_z_xyzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 38);

    auto tg_z_xzzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 39);

    auto tg_z_yyyy_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 40);

    auto tg_z_yyyz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 41);

    auto tg_z_yyzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 42);

    auto tg_z_yzzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 43);

    auto tg_z_zzzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 44);

    // Set up components of auxiliary buffer : SG

    auto tg_0_xxxx_g_1_0_0 = pbuffer.data(idx_sg_g_1_0_0);

    auto tg_0_xxxy_g_1_0_0 = pbuffer.data(idx_sg_g_1_0_0 + 1);

    auto tg_0_xxxz_g_1_0_0 = pbuffer.data(idx_sg_g_1_0_0 + 2);

    auto tg_0_xxyy_g_1_0_0 = pbuffer.data(idx_sg_g_1_0_0 + 3);

    auto tg_0_xxyz_g_1_0_0 = pbuffer.data(idx_sg_g_1_0_0 + 4);

    auto tg_0_xxzz_g_1_0_0 = pbuffer.data(idx_sg_g_1_0_0 + 5);

    auto tg_0_xyyy_g_1_0_0 = pbuffer.data(idx_sg_g_1_0_0 + 6);

    auto tg_0_xyyz_g_1_0_0 = pbuffer.data(idx_sg_g_1_0_0 + 7);

    auto tg_0_xyzz_g_1_0_0 = pbuffer.data(idx_sg_g_1_0_0 + 8);

    auto tg_0_xzzz_g_1_0_0 = pbuffer.data(idx_sg_g_1_0_0 + 9);

    auto tg_0_yyyy_g_1_0_0 = pbuffer.data(idx_sg_g_1_0_0 + 10);

    auto tg_0_yyyz_g_1_0_0 = pbuffer.data(idx_sg_g_1_0_0 + 11);

    auto tg_0_yyzz_g_1_0_0 = pbuffer.data(idx_sg_g_1_0_0 + 12);

    auto tg_0_yzzz_g_1_0_0 = pbuffer.data(idx_sg_g_1_0_0 + 13);

    auto tg_0_zzzz_g_1_0_0 = pbuffer.data(idx_sg_g_1_0_0 + 14);

    // Set up components of auxiliary buffer : PG

    auto tg_x_xxxx_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0);

    auto tg_x_xxxy_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 1);

    auto tg_x_xxxz_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 2);

    auto tg_x_xxyy_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 3);

    auto tg_x_xxyz_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 4);

    auto tg_x_xxzz_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 5);

    auto tg_x_xyyy_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 6);

    auto tg_x_xyyz_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 7);

    auto tg_x_xyzz_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 8);

    auto tg_x_xzzz_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 9);

    auto tg_x_yyyy_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 10);

    auto tg_x_yyyz_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 11);

    auto tg_x_yyzz_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 12);

    auto tg_x_yzzz_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 13);

    auto tg_x_zzzz_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 14);

    auto tg_y_xxxx_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 15);

    auto tg_y_xxxy_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 16);

    auto tg_y_xxxz_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 17);

    auto tg_y_xxyy_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 18);

    auto tg_y_xxyz_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 19);

    auto tg_y_xxzz_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 20);

    auto tg_y_xyyy_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 21);

    auto tg_y_xyyz_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 22);

    auto tg_y_xyzz_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 23);

    auto tg_y_xzzz_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 24);

    auto tg_y_yyyy_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 25);

    auto tg_y_yyyz_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 26);

    auto tg_y_yyzz_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 27);

    auto tg_y_yzzz_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 28);

    auto tg_y_zzzz_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 29);

    auto tg_z_xxxx_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 30);

    auto tg_z_xxxy_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 31);

    auto tg_z_xxxz_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 32);

    auto tg_z_xxyy_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 33);

    auto tg_z_xxyz_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 34);

    auto tg_z_xxzz_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 35);

    auto tg_z_xyyy_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 36);

    auto tg_z_xyyz_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 37);

    auto tg_z_xyzz_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 38);

    auto tg_z_xzzz_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 39);

    auto tg_z_yyyy_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 40);

    auto tg_z_yyyz_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 41);

    auto tg_z_yyzz_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 42);

    auto tg_z_yzzz_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 43);

    auto tg_z_zzzz_g_1_0_0 = pbuffer.data(idx_pg_g_1_0_0 + 44);

    // Set up components of auxiliary buffer : SG

    auto tg_0_xxxx_d_1_0_1 = pbuffer.data(idx_sg_d_1_0_1);

    auto tg_0_xxxy_d_1_0_1 = pbuffer.data(idx_sg_d_1_0_1 + 1);

    auto tg_0_xxxz_d_1_0_1 = pbuffer.data(idx_sg_d_1_0_1 + 2);

    auto tg_0_xxyy_d_1_0_1 = pbuffer.data(idx_sg_d_1_0_1 + 3);

    auto tg_0_xxyz_d_1_0_1 = pbuffer.data(idx_sg_d_1_0_1 + 4);

    auto tg_0_xxzz_d_1_0_1 = pbuffer.data(idx_sg_d_1_0_1 + 5);

    auto tg_0_xyyy_d_1_0_1 = pbuffer.data(idx_sg_d_1_0_1 + 6);

    auto tg_0_xyyz_d_1_0_1 = pbuffer.data(idx_sg_d_1_0_1 + 7);

    auto tg_0_xyzz_d_1_0_1 = pbuffer.data(idx_sg_d_1_0_1 + 8);

    auto tg_0_xzzz_d_1_0_1 = pbuffer.data(idx_sg_d_1_0_1 + 9);

    auto tg_0_yyyy_d_1_0_1 = pbuffer.data(idx_sg_d_1_0_1 + 10);

    auto tg_0_yyyz_d_1_0_1 = pbuffer.data(idx_sg_d_1_0_1 + 11);

    auto tg_0_yyzz_d_1_0_1 = pbuffer.data(idx_sg_d_1_0_1 + 12);

    auto tg_0_yzzz_d_1_0_1 = pbuffer.data(idx_sg_d_1_0_1 + 13);

    auto tg_0_zzzz_d_1_0_1 = pbuffer.data(idx_sg_d_1_0_1 + 14);

    // Set up components of auxiliary buffer : PG

    auto tg_x_xxxx_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1);

    auto tg_x_xxxy_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 1);

    auto tg_x_xxxz_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 2);

    auto tg_x_xxyy_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 3);

    auto tg_x_xxyz_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 4);

    auto tg_x_xxzz_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 5);

    auto tg_x_xyyy_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 6);

    auto tg_x_xyyz_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 7);

    auto tg_x_xyzz_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 8);

    auto tg_x_xzzz_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 9);

    auto tg_x_yyyy_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 10);

    auto tg_x_yyyz_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 11);

    auto tg_x_yyzz_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 12);

    auto tg_x_yzzz_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 13);

    auto tg_x_zzzz_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 14);

    auto tg_y_xxxx_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 15);

    auto tg_y_xxxy_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 16);

    auto tg_y_xxxz_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 17);

    auto tg_y_xxyy_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 18);

    auto tg_y_xxyz_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 19);

    auto tg_y_xxzz_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 20);

    auto tg_y_xyyy_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 21);

    auto tg_y_xyyz_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 22);

    auto tg_y_xyzz_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 23);

    auto tg_y_xzzz_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 24);

    auto tg_y_yyyy_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 25);

    auto tg_y_yyyz_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 26);

    auto tg_y_yyzz_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 27);

    auto tg_y_yzzz_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 28);

    auto tg_y_zzzz_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 29);

    auto tg_z_xxxx_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 30);

    auto tg_z_xxxy_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 31);

    auto tg_z_xxxz_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 32);

    auto tg_z_xxyy_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 33);

    auto tg_z_xxyz_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 34);

    auto tg_z_xxzz_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 35);

    auto tg_z_xyyy_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 36);

    auto tg_z_xyyz_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 37);

    auto tg_z_xyzz_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 38);

    auto tg_z_xzzz_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 39);

    auto tg_z_yyyy_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 40);

    auto tg_z_yyyz_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 41);

    auto tg_z_yyzz_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 42);

    auto tg_z_yzzz_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 43);

    auto tg_z_zzzz_d_1_0_1 = pbuffer.data(idx_pg_d_1_0_1 + 44);

    // Set up components of auxiliary buffer : PF

    auto tg_x_xxx_p_1_1_1 = pbuffer.data(idx_pf_p_1_1_1);

    auto tg_x_xxy_p_1_1_1 = pbuffer.data(idx_pf_p_1_1_1 + 1);

    auto tg_x_xxz_p_1_1_1 = pbuffer.data(idx_pf_p_1_1_1 + 2);

    auto tg_x_xyy_p_1_1_1 = pbuffer.data(idx_pf_p_1_1_1 + 3);

    auto tg_x_xyz_p_1_1_1 = pbuffer.data(idx_pf_p_1_1_1 + 4);

    auto tg_x_xzz_p_1_1_1 = pbuffer.data(idx_pf_p_1_1_1 + 5);

    auto tg_x_yyy_p_1_1_1 = pbuffer.data(idx_pf_p_1_1_1 + 6);

    auto tg_x_yyz_p_1_1_1 = pbuffer.data(idx_pf_p_1_1_1 + 7);

    auto tg_x_yzz_p_1_1_1 = pbuffer.data(idx_pf_p_1_1_1 + 8);

    auto tg_x_zzz_p_1_1_1 = pbuffer.data(idx_pf_p_1_1_1 + 9);

    auto tg_y_xxx_p_1_1_1 = pbuffer.data(idx_pf_p_1_1_1 + 10);

    auto tg_y_xxy_p_1_1_1 = pbuffer.data(idx_pf_p_1_1_1 + 11);

    auto tg_y_xxz_p_1_1_1 = pbuffer.data(idx_pf_p_1_1_1 + 12);

    auto tg_y_xyy_p_1_1_1 = pbuffer.data(idx_pf_p_1_1_1 + 13);

    auto tg_y_xyz_p_1_1_1 = pbuffer.data(idx_pf_p_1_1_1 + 14);

    auto tg_y_xzz_p_1_1_1 = pbuffer.data(idx_pf_p_1_1_1 + 15);

    auto tg_y_yyy_p_1_1_1 = pbuffer.data(idx_pf_p_1_1_1 + 16);

    auto tg_y_yyz_p_1_1_1 = pbuffer.data(idx_pf_p_1_1_1 + 17);

    auto tg_y_yzz_p_1_1_1 = pbuffer.data(idx_pf_p_1_1_1 + 18);

    auto tg_y_zzz_p_1_1_1 = pbuffer.data(idx_pf_p_1_1_1 + 19);

    auto tg_z_xxx_p_1_1_1 = pbuffer.data(idx_pf_p_1_1_1 + 20);

    auto tg_z_xxy_p_1_1_1 = pbuffer.data(idx_pf_p_1_1_1 + 21);

    auto tg_z_xxz_p_1_1_1 = pbuffer.data(idx_pf_p_1_1_1 + 22);

    auto tg_z_xyy_p_1_1_1 = pbuffer.data(idx_pf_p_1_1_1 + 23);

    auto tg_z_xyz_p_1_1_1 = pbuffer.data(idx_pf_p_1_1_1 + 24);

    auto tg_z_xzz_p_1_1_1 = pbuffer.data(idx_pf_p_1_1_1 + 25);

    auto tg_z_yyy_p_1_1_1 = pbuffer.data(idx_pf_p_1_1_1 + 26);

    auto tg_z_yyz_p_1_1_1 = pbuffer.data(idx_pf_p_1_1_1 + 27);

    auto tg_z_yzz_p_1_1_1 = pbuffer.data(idx_pf_p_1_1_1 + 28);

    auto tg_z_zzz_p_1_1_1 = pbuffer.data(idx_pf_p_1_1_1 + 29);

    // Set up components of auxiliary buffer : PG

    auto tg_x_xxxx_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1);

    auto tg_x_xxxy_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 1);

    auto tg_x_xxxz_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 2);

    auto tg_x_xxyy_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 3);

    auto tg_x_xxyz_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 4);

    auto tg_x_xxzz_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 5);

    auto tg_x_xyyy_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 6);

    auto tg_x_xyyz_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 7);

    auto tg_x_xyzz_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 8);

    auto tg_x_xzzz_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 9);

    auto tg_x_yyyy_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 10);

    auto tg_x_yyyz_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 11);

    auto tg_x_yyzz_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 12);

    auto tg_x_yzzz_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 13);

    auto tg_x_zzzz_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 14);

    auto tg_y_xxxx_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 15);

    auto tg_y_xxxy_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 16);

    auto tg_y_xxxz_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 17);

    auto tg_y_xxyy_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 18);

    auto tg_y_xxyz_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 19);

    auto tg_y_xxzz_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 20);

    auto tg_y_xyyy_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 21);

    auto tg_y_xyyz_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 22);

    auto tg_y_xyzz_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 23);

    auto tg_y_xzzz_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 24);

    auto tg_y_yyyy_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 25);

    auto tg_y_yyyz_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 26);

    auto tg_y_yyzz_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 27);

    auto tg_y_yzzz_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 28);

    auto tg_y_zzzz_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 29);

    auto tg_z_xxxx_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 30);

    auto tg_z_xxxy_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 31);

    auto tg_z_xxxz_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 32);

    auto tg_z_xxyy_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 33);

    auto tg_z_xxyz_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 34);

    auto tg_z_xxzz_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 35);

    auto tg_z_xyyy_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 36);

    auto tg_z_xyyz_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 37);

    auto tg_z_xyzz_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 38);

    auto tg_z_xzzz_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 39);

    auto tg_z_yyyy_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 40);

    auto tg_z_yyyz_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 41);

    auto tg_z_yyzz_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 42);

    auto tg_z_yzzz_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 43);

    auto tg_z_zzzz_p_1_1_1 = pbuffer.data(idx_pg_p_1_1_1 + 44);

    // Set up components of auxiliary buffer : SG

    auto tg_0_xxxx_s_2_1_1 = pbuffer.data(idx_sg_s_2_1_1);

    auto tg_0_xxxy_s_2_1_1 = pbuffer.data(idx_sg_s_2_1_1 + 1);

    auto tg_0_xxxz_s_2_1_1 = pbuffer.data(idx_sg_s_2_1_1 + 2);

    auto tg_0_xxyy_s_2_1_1 = pbuffer.data(idx_sg_s_2_1_1 + 3);

    auto tg_0_xxyz_s_2_1_1 = pbuffer.data(idx_sg_s_2_1_1 + 4);

    auto tg_0_xxzz_s_2_1_1 = pbuffer.data(idx_sg_s_2_1_1 + 5);

    auto tg_0_xyyy_s_2_1_1 = pbuffer.data(idx_sg_s_2_1_1 + 6);

    auto tg_0_xyyz_s_2_1_1 = pbuffer.data(idx_sg_s_2_1_1 + 7);

    auto tg_0_xyzz_s_2_1_1 = pbuffer.data(idx_sg_s_2_1_1 + 8);

    auto tg_0_xzzz_s_2_1_1 = pbuffer.data(idx_sg_s_2_1_1 + 9);

    auto tg_0_yyyy_s_2_1_1 = pbuffer.data(idx_sg_s_2_1_1 + 10);

    auto tg_0_yyyz_s_2_1_1 = pbuffer.data(idx_sg_s_2_1_1 + 11);

    auto tg_0_yyzz_s_2_1_1 = pbuffer.data(idx_sg_s_2_1_1 + 12);

    auto tg_0_yzzz_s_2_1_1 = pbuffer.data(idx_sg_s_2_1_1 + 13);

    auto tg_0_zzzz_s_2_1_1 = pbuffer.data(idx_sg_s_2_1_1 + 14);

    // Set up components of auxiliary buffer : PG

    auto tg_x_xxxx_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1);

    auto tg_x_xxxy_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 1);

    auto tg_x_xxxz_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 2);

    auto tg_x_xxyy_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 3);

    auto tg_x_xxyz_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 4);

    auto tg_x_xxzz_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 5);

    auto tg_x_xyyy_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 6);

    auto tg_x_xyyz_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 7);

    auto tg_x_xyzz_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 8);

    auto tg_x_xzzz_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 9);

    auto tg_x_yyyy_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 10);

    auto tg_x_yyyz_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 11);

    auto tg_x_yyzz_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 12);

    auto tg_x_yzzz_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 13);

    auto tg_x_zzzz_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 14);

    auto tg_y_xxxx_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 15);

    auto tg_y_xxxy_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 16);

    auto tg_y_xxxz_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 17);

    auto tg_y_xxyy_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 18);

    auto tg_y_xxyz_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 19);

    auto tg_y_xxzz_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 20);

    auto tg_y_xyyy_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 21);

    auto tg_y_xyyz_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 22);

    auto tg_y_xyzz_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 23);

    auto tg_y_xzzz_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 24);

    auto tg_y_yyyy_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 25);

    auto tg_y_yyyz_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 26);

    auto tg_y_yyzz_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 27);

    auto tg_y_yzzz_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 28);

    auto tg_y_zzzz_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 29);

    auto tg_z_xxxx_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 30);

    auto tg_z_xxxy_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 31);

    auto tg_z_xxxz_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 32);

    auto tg_z_xxyy_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 33);

    auto tg_z_xxyz_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 34);

    auto tg_z_xxzz_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 35);

    auto tg_z_xyyy_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 36);

    auto tg_z_xyyz_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 37);

    auto tg_z_xyzz_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 38);

    auto tg_z_xzzz_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 39);

    auto tg_z_yyyy_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 40);

    auto tg_z_yyyz_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 41);

    auto tg_z_yyzz_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 42);

    auto tg_z_yzzz_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 43);

    auto tg_z_zzzz_s_2_1_1 = pbuffer.data(idx_pg_s_2_1_1 + 44);

    // Set up components of targeted buffer : DG

    auto tg_xx_xxxx_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0);

    auto tg_xx_xxxy_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 1);

    auto tg_xx_xxxz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 2);

    auto tg_xx_xxyy_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 3);

    auto tg_xx_xxyz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 4);

    auto tg_xx_xxzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 5);

    auto tg_xx_xyyy_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 6);

    auto tg_xx_xyyz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 7);

    auto tg_xx_xyzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 8);

    auto tg_xx_xzzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 9);

    auto tg_xx_yyyy_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 10);

    auto tg_xx_yyyz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 11);

    auto tg_xx_yyzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 12);

    auto tg_xx_yzzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 13);

    auto tg_xx_zzzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 14);

    auto tg_xy_xxxx_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 15);

    auto tg_xy_xxxy_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 16);

    auto tg_xy_xxxz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 17);

    auto tg_xy_xxyy_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 18);

    auto tg_xy_xxyz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 19);

    auto tg_xy_xxzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 20);

    auto tg_xy_xyyy_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 21);

    auto tg_xy_xyyz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 22);

    auto tg_xy_xyzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 23);

    auto tg_xy_xzzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 24);

    auto tg_xy_yyyy_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 25);

    auto tg_xy_yyyz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 26);

    auto tg_xy_yyzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 27);

    auto tg_xy_yzzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 28);

    auto tg_xy_zzzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 29);

    auto tg_xz_xxxx_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 30);

    auto tg_xz_xxxy_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 31);

    auto tg_xz_xxxz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 32);

    auto tg_xz_xxyy_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 33);

    auto tg_xz_xxyz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 34);

    auto tg_xz_xxzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 35);

    auto tg_xz_xyyy_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 36);

    auto tg_xz_xyyz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 37);

    auto tg_xz_xyzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 38);

    auto tg_xz_xzzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 39);

    auto tg_xz_yyyy_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 40);

    auto tg_xz_yyyz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 41);

    auto tg_xz_yyzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 42);

    auto tg_xz_yzzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 43);

    auto tg_xz_zzzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 44);

    auto tg_yy_xxxx_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 45);

    auto tg_yy_xxxy_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 46);

    auto tg_yy_xxxz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 47);

    auto tg_yy_xxyy_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 48);

    auto tg_yy_xxyz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 49);

    auto tg_yy_xxzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 50);

    auto tg_yy_xyyy_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 51);

    auto tg_yy_xyyz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 52);

    auto tg_yy_xyzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 53);

    auto tg_yy_xzzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 54);

    auto tg_yy_yyyy_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 55);

    auto tg_yy_yyyz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 56);

    auto tg_yy_yyzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 57);

    auto tg_yy_yzzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 58);

    auto tg_yy_zzzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 59);

    auto tg_yz_xxxx_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 60);

    auto tg_yz_xxxy_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 61);

    auto tg_yz_xxxz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 62);

    auto tg_yz_xxyy_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 63);

    auto tg_yz_xxyz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 64);

    auto tg_yz_xxzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 65);

    auto tg_yz_xyyy_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 66);

    auto tg_yz_xyyz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 67);

    auto tg_yz_xyzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 68);

    auto tg_yz_xzzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 69);

    auto tg_yz_yyyy_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 70);

    auto tg_yz_yyyz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 71);

    auto tg_yz_yyzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 72);

    auto tg_yz_yzzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 73);

    auto tg_yz_zzzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 74);

    auto tg_zz_xxxx_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 75);

    auto tg_zz_xxxy_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 76);

    auto tg_zz_xxxz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 77);

    auto tg_zz_xxyy_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 78);

    auto tg_zz_xxyz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 79);

    auto tg_zz_xxzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 80);

    auto tg_zz_xyyy_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 81);

    auto tg_zz_xyyz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 82);

    auto tg_zz_xyzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 83);

    auto tg_zz_xzzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 84);

    auto tg_zz_yyyy_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 85);

    auto tg_zz_yyyz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 86);

    auto tg_zz_yyzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 87);

    auto tg_zz_yzzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 88);

    auto tg_zz_zzzz_g_0_0_0 = pbuffer.data(idx_dg_g_0_0_0 + 89);

    #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_0_xxxx_d_1_0_1, tg_0_xxxx_g_0_0_0, tg_0_xxxx_g_1_0_0, tg_0_xxxx_s_2_1_1, tg_0_xxxy_d_1_0_1, tg_0_xxxy_g_0_0_0, tg_0_xxxy_g_1_0_0, tg_0_xxxy_s_2_1_1, tg_0_xxxz_d_1_0_1, tg_0_xxxz_g_0_0_0, tg_0_xxxz_g_1_0_0, tg_0_xxxz_s_2_1_1, tg_0_xxyy_d_1_0_1, tg_0_xxyy_g_0_0_0, tg_0_xxyy_g_1_0_0, tg_0_xxyy_s_2_1_1, tg_0_xxyz_d_1_0_1, tg_0_xxyz_g_0_0_0, tg_0_xxyz_g_1_0_0, tg_0_xxyz_s_2_1_1, tg_0_xxzz_d_1_0_1, tg_0_xxzz_g_0_0_0, tg_0_xxzz_g_1_0_0, tg_0_xxzz_s_2_1_1, tg_0_xyyy_d_1_0_1, tg_0_xyyy_g_0_0_0, tg_0_xyyy_g_1_0_0, tg_0_xyyy_s_2_1_1, tg_0_xyyz_d_1_0_1, tg_0_xyyz_g_0_0_0, tg_0_xyyz_g_1_0_0, tg_0_xyyz_s_2_1_1, tg_0_xyzz_d_1_0_1, tg_0_xyzz_g_0_0_0, tg_0_xyzz_g_1_0_0, tg_0_xyzz_s_2_1_1, tg_0_xzzz_d_1_0_1, tg_0_xzzz_g_0_0_0, tg_0_xzzz_g_1_0_0, tg_0_xzzz_s_2_1_1, tg_0_yyyy_d_1_0_1, tg_0_yyyy_g_0_0_0, tg_0_yyyy_g_1_0_0, tg_0_yyyy_s_2_1_1, tg_0_yyyz_d_1_0_1, tg_0_yyyz_g_0_0_0, tg_0_yyyz_g_1_0_0, tg_0_yyyz_s_2_1_1, tg_0_yyzz_d_1_0_1, tg_0_yyzz_g_0_0_0, tg_0_yyzz_g_1_0_0, tg_0_yyzz_s_2_1_1, tg_0_yzzz_d_1_0_1, tg_0_yzzz_g_0_0_0, tg_0_yzzz_g_1_0_0, tg_0_yzzz_s_2_1_1, tg_0_zzzz_d_1_0_1, tg_0_zzzz_g_0_0_0, tg_0_zzzz_g_1_0_0, tg_0_zzzz_s_2_1_1, tg_x_xxx_f_0_0_1, tg_x_xxx_p_1_1_1, tg_x_xxxx_d_1_0_1, tg_x_xxxx_f_0_0_1, tg_x_xxxx_g_0_0_0, tg_x_xxxx_g_1_0_0, tg_x_xxxx_p_1_1_1, tg_x_xxxx_s_2_1_1, tg_x_xxxy_d_1_0_1, tg_x_xxxy_f_0_0_1, tg_x_xxxy_g_0_0_0, tg_x_xxxy_g_1_0_0, tg_x_xxxy_p_1_1_1, tg_x_xxxy_s_2_1_1, tg_x_xxxz_d_1_0_1, tg_x_xxxz_f_0_0_1, tg_x_xxxz_g_0_0_0, tg_x_xxxz_g_1_0_0, tg_x_xxxz_p_1_1_1, tg_x_xxxz_s_2_1_1, tg_x_xxy_f_0_0_1, tg_x_xxy_p_1_1_1, tg_x_xxyy_d_1_0_1, tg_x_xxyy_f_0_0_1, tg_x_xxyy_g_0_0_0, tg_x_xxyy_g_1_0_0, tg_x_xxyy_p_1_1_1, tg_x_xxyy_s_2_1_1, tg_x_xxyz_d_1_0_1, tg_x_xxyz_f_0_0_1, tg_x_xxyz_g_0_0_0, tg_x_xxyz_g_1_0_0, tg_x_xxyz_p_1_1_1, tg_x_xxyz_s_2_1_1, tg_x_xxz_f_0_0_1, tg_x_xxz_p_1_1_1, tg_x_xxzz_d_1_0_1, tg_x_xxzz_f_0_0_1, tg_x_xxzz_g_0_0_0, tg_x_xxzz_g_1_0_0, tg_x_xxzz_p_1_1_1, tg_x_xxzz_s_2_1_1, tg_x_xyy_f_0_0_1, tg_x_xyy_p_1_1_1, tg_x_xyyy_d_1_0_1, tg_x_xyyy_f_0_0_1, tg_x_xyyy_g_0_0_0, tg_x_xyyy_g_1_0_0, tg_x_xyyy_p_1_1_1, tg_x_xyyy_s_2_1_1, tg_x_xyyz_d_1_0_1, tg_x_xyyz_f_0_0_1, tg_x_xyyz_g_0_0_0, tg_x_xyyz_g_1_0_0, tg_x_xyyz_p_1_1_1, tg_x_xyyz_s_2_1_1, tg_x_xyz_f_0_0_1, tg_x_xyz_p_1_1_1, tg_x_xyzz_d_1_0_1, tg_x_xyzz_f_0_0_1, tg_x_xyzz_g_0_0_0, tg_x_xyzz_g_1_0_0, tg_x_xyzz_p_1_1_1, tg_x_xyzz_s_2_1_1, tg_x_xzz_f_0_0_1, tg_x_xzz_p_1_1_1, tg_x_xzzz_d_1_0_1, tg_x_xzzz_f_0_0_1, tg_x_xzzz_g_0_0_0, tg_x_xzzz_g_1_0_0, tg_x_xzzz_p_1_1_1, tg_x_xzzz_s_2_1_1, tg_x_yyy_f_0_0_1, tg_x_yyy_p_1_1_1, tg_x_yyyy_d_1_0_1, tg_x_yyyy_f_0_0_1, tg_x_yyyy_g_0_0_0, tg_x_yyyy_g_1_0_0, tg_x_yyyy_p_1_1_1, tg_x_yyyy_s_2_1_1, tg_x_yyyz_d_1_0_1, tg_x_yyyz_f_0_0_1, tg_x_yyyz_g_0_0_0, tg_x_yyyz_g_1_0_0, tg_x_yyyz_p_1_1_1, tg_x_yyyz_s_2_1_1, tg_x_yyz_f_0_0_1, tg_x_yyz_p_1_1_1, tg_x_yyzz_d_1_0_1, tg_x_yyzz_f_0_0_1, tg_x_yyzz_g_0_0_0, tg_x_yyzz_g_1_0_0, tg_x_yyzz_p_1_1_1, tg_x_yyzz_s_2_1_1, tg_x_yzz_f_0_0_1, tg_x_yzz_p_1_1_1, tg_x_yzzz_d_1_0_1, tg_x_yzzz_f_0_0_1, tg_x_yzzz_g_0_0_0, tg_x_yzzz_g_1_0_0, tg_x_yzzz_p_1_1_1, tg_x_yzzz_s_2_1_1, tg_x_zzz_f_0_0_1, tg_x_zzz_p_1_1_1, tg_x_zzzz_d_1_0_1, tg_x_zzzz_f_0_0_1, tg_x_zzzz_g_0_0_0, tg_x_zzzz_g_1_0_0, tg_x_zzzz_p_1_1_1, tg_x_zzzz_s_2_1_1, tg_xx_xxxx_g_0_0_0, tg_xx_xxxy_g_0_0_0, tg_xx_xxxz_g_0_0_0, tg_xx_xxyy_g_0_0_0, tg_xx_xxyz_g_0_0_0, tg_xx_xxzz_g_0_0_0, tg_xx_xyyy_g_0_0_0, tg_xx_xyyz_g_0_0_0, tg_xx_xyzz_g_0_0_0, tg_xx_xzzz_g_0_0_0, tg_xx_yyyy_g_0_0_0, tg_xx_yyyz_g_0_0_0, tg_xx_yyzz_g_0_0_0, tg_xx_yzzz_g_0_0_0, tg_xx_zzzz_g_0_0_0, tg_xy_xxxx_g_0_0_0, tg_xy_xxxy_g_0_0_0, tg_xy_xxxz_g_0_0_0, tg_xy_xxyy_g_0_0_0, tg_xy_xxyz_g_0_0_0, tg_xy_xxzz_g_0_0_0, tg_xy_xyyy_g_0_0_0, tg_xy_xyyz_g_0_0_0, tg_xy_xyzz_g_0_0_0, tg_xy_xzzz_g_0_0_0, tg_xy_yyyy_g_0_0_0, tg_xy_yyyz_g_0_0_0, tg_xy_yyzz_g_0_0_0, tg_xy_yzzz_g_0_0_0, tg_xy_zzzz_g_0_0_0, tg_xz_xxxx_g_0_0_0, tg_xz_xxxy_g_0_0_0, tg_xz_xxxz_g_0_0_0, tg_xz_xxyy_g_0_0_0, tg_xz_xxyz_g_0_0_0, tg_xz_xxzz_g_0_0_0, tg_xz_xyyy_g_0_0_0, tg_xz_xyyz_g_0_0_0, tg_xz_xyzz_g_0_0_0, tg_xz_xzzz_g_0_0_0, tg_xz_yyyy_g_0_0_0, tg_xz_yyyz_g_0_0_0, tg_xz_yyzz_g_0_0_0, tg_xz_yzzz_g_0_0_0, tg_xz_zzzz_g_0_0_0, tg_y_xxx_f_0_0_1, tg_y_xxx_p_1_1_1, tg_y_xxxx_d_1_0_1, tg_y_xxxx_f_0_0_1, tg_y_xxxx_g_0_0_0, tg_y_xxxx_g_1_0_0, tg_y_xxxx_p_1_1_1, tg_y_xxxx_s_2_1_1, tg_y_xxxy_d_1_0_1, tg_y_xxxy_f_0_0_1, tg_y_xxxy_g_0_0_0, tg_y_xxxy_g_1_0_0, tg_y_xxxy_p_1_1_1, tg_y_xxxy_s_2_1_1, tg_y_xxxz_d_1_0_1, tg_y_xxxz_f_0_0_1, tg_y_xxxz_g_0_0_0, tg_y_xxxz_g_1_0_0, tg_y_xxxz_p_1_1_1, tg_y_xxxz_s_2_1_1, tg_y_xxy_f_0_0_1, tg_y_xxy_p_1_1_1, tg_y_xxyy_d_1_0_1, tg_y_xxyy_f_0_0_1, tg_y_xxyy_g_0_0_0, tg_y_xxyy_g_1_0_0, tg_y_xxyy_p_1_1_1, tg_y_xxyy_s_2_1_1, tg_y_xxyz_d_1_0_1, tg_y_xxyz_f_0_0_1, tg_y_xxyz_g_0_0_0, tg_y_xxyz_g_1_0_0, tg_y_xxyz_p_1_1_1, tg_y_xxyz_s_2_1_1, tg_y_xxz_f_0_0_1, tg_y_xxz_p_1_1_1, tg_y_xxzz_d_1_0_1, tg_y_xxzz_f_0_0_1, tg_y_xxzz_g_0_0_0, tg_y_xxzz_g_1_0_0, tg_y_xxzz_p_1_1_1, tg_y_xxzz_s_2_1_1, tg_y_xyy_f_0_0_1, tg_y_xyy_p_1_1_1, tg_y_xyyy_d_1_0_1, tg_y_xyyy_f_0_0_1, tg_y_xyyy_g_0_0_0, tg_y_xyyy_g_1_0_0, tg_y_xyyy_p_1_1_1, tg_y_xyyy_s_2_1_1, tg_y_xyyz_d_1_0_1, tg_y_xyyz_f_0_0_1, tg_y_xyyz_g_0_0_0, tg_y_xyyz_g_1_0_0, tg_y_xyyz_p_1_1_1, tg_y_xyyz_s_2_1_1, tg_y_xyz_f_0_0_1, tg_y_xyz_p_1_1_1, tg_y_xyzz_d_1_0_1, tg_y_xyzz_f_0_0_1, tg_y_xyzz_g_0_0_0, tg_y_xyzz_g_1_0_0, tg_y_xyzz_p_1_1_1, tg_y_xyzz_s_2_1_1, tg_y_xzz_f_0_0_1, tg_y_xzz_p_1_1_1, tg_y_xzzz_d_1_0_1, tg_y_xzzz_f_0_0_1, tg_y_xzzz_g_0_0_0, tg_y_xzzz_g_1_0_0, tg_y_xzzz_p_1_1_1, tg_y_xzzz_s_2_1_1, tg_y_yyy_f_0_0_1, tg_y_yyy_p_1_1_1, tg_y_yyyy_d_1_0_1, tg_y_yyyy_f_0_0_1, tg_y_yyyy_g_0_0_0, tg_y_yyyy_g_1_0_0, tg_y_yyyy_p_1_1_1, tg_y_yyyy_s_2_1_1, tg_y_yyyz_d_1_0_1, tg_y_yyyz_f_0_0_1, tg_y_yyyz_g_0_0_0, tg_y_yyyz_g_1_0_0, tg_y_yyyz_p_1_1_1, tg_y_yyyz_s_2_1_1, tg_y_yyz_f_0_0_1, tg_y_yyz_p_1_1_1, tg_y_yyzz_d_1_0_1, tg_y_yyzz_f_0_0_1, tg_y_yyzz_g_0_0_0, tg_y_yyzz_g_1_0_0, tg_y_yyzz_p_1_1_1, tg_y_yyzz_s_2_1_1, tg_y_yzz_f_0_0_1, tg_y_yzz_p_1_1_1, tg_y_yzzz_d_1_0_1, tg_y_yzzz_f_0_0_1, tg_y_yzzz_g_0_0_0, tg_y_yzzz_g_1_0_0, tg_y_yzzz_p_1_1_1, tg_y_yzzz_s_2_1_1, tg_y_zzz_f_0_0_1, tg_y_zzz_p_1_1_1, tg_y_zzzz_d_1_0_1, tg_y_zzzz_f_0_0_1, tg_y_zzzz_g_0_0_0, tg_y_zzzz_g_1_0_0, tg_y_zzzz_p_1_1_1, tg_y_zzzz_s_2_1_1, tg_yy_xxxx_g_0_0_0, tg_yy_xxxy_g_0_0_0, tg_yy_xxxz_g_0_0_0, tg_yy_xxyy_g_0_0_0, tg_yy_xxyz_g_0_0_0, tg_yy_xxzz_g_0_0_0, tg_yy_xyyy_g_0_0_0, tg_yy_xyyz_g_0_0_0, tg_yy_xyzz_g_0_0_0, tg_yy_xzzz_g_0_0_0, tg_yy_yyyy_g_0_0_0, tg_yy_yyyz_g_0_0_0, tg_yy_yyzz_g_0_0_0, tg_yy_yzzz_g_0_0_0, tg_yy_zzzz_g_0_0_0, tg_yz_xxxx_g_0_0_0, tg_yz_xxxy_g_0_0_0, tg_yz_xxxz_g_0_0_0, tg_yz_xxyy_g_0_0_0, tg_yz_xxyz_g_0_0_0, tg_yz_xxzz_g_0_0_0, tg_yz_xyyy_g_0_0_0, tg_yz_xyyz_g_0_0_0, tg_yz_xyzz_g_0_0_0, tg_yz_xzzz_g_0_0_0, tg_yz_yyyy_g_0_0_0, tg_yz_yyyz_g_0_0_0, tg_yz_yyzz_g_0_0_0, tg_yz_yzzz_g_0_0_0, tg_yz_zzzz_g_0_0_0, tg_z_xxx_f_0_0_1, tg_z_xxx_p_1_1_1, tg_z_xxxx_d_1_0_1, tg_z_xxxx_f_0_0_1, tg_z_xxxx_g_0_0_0, tg_z_xxxx_g_1_0_0, tg_z_xxxx_p_1_1_1, tg_z_xxxx_s_2_1_1, tg_z_xxxy_d_1_0_1, tg_z_xxxy_f_0_0_1, tg_z_xxxy_g_0_0_0, tg_z_xxxy_g_1_0_0, tg_z_xxxy_p_1_1_1, tg_z_xxxy_s_2_1_1, tg_z_xxxz_d_1_0_1, tg_z_xxxz_f_0_0_1, tg_z_xxxz_g_0_0_0, tg_z_xxxz_g_1_0_0, tg_z_xxxz_p_1_1_1, tg_z_xxxz_s_2_1_1, tg_z_xxy_f_0_0_1, tg_z_xxy_p_1_1_1, tg_z_xxyy_d_1_0_1, tg_z_xxyy_f_0_0_1, tg_z_xxyy_g_0_0_0, tg_z_xxyy_g_1_0_0, tg_z_xxyy_p_1_1_1, tg_z_xxyy_s_2_1_1, tg_z_xxyz_d_1_0_1, tg_z_xxyz_f_0_0_1, tg_z_xxyz_g_0_0_0, tg_z_xxyz_g_1_0_0, tg_z_xxyz_p_1_1_1, tg_z_xxyz_s_2_1_1, tg_z_xxz_f_0_0_1, tg_z_xxz_p_1_1_1, tg_z_xxzz_d_1_0_1, tg_z_xxzz_f_0_0_1, tg_z_xxzz_g_0_0_0, tg_z_xxzz_g_1_0_0, tg_z_xxzz_p_1_1_1, tg_z_xxzz_s_2_1_1, tg_z_xyy_f_0_0_1, tg_z_xyy_p_1_1_1, tg_z_xyyy_d_1_0_1, tg_z_xyyy_f_0_0_1, tg_z_xyyy_g_0_0_0, tg_z_xyyy_g_1_0_0, tg_z_xyyy_p_1_1_1, tg_z_xyyy_s_2_1_1, tg_z_xyyz_d_1_0_1, tg_z_xyyz_f_0_0_1, tg_z_xyyz_g_0_0_0, tg_z_xyyz_g_1_0_0, tg_z_xyyz_p_1_1_1, tg_z_xyyz_s_2_1_1, tg_z_xyz_f_0_0_1, tg_z_xyz_p_1_1_1, tg_z_xyzz_d_1_0_1, tg_z_xyzz_f_0_0_1, tg_z_xyzz_g_0_0_0, tg_z_xyzz_g_1_0_0, tg_z_xyzz_p_1_1_1, tg_z_xyzz_s_2_1_1, tg_z_xzz_f_0_0_1, tg_z_xzz_p_1_1_1, tg_z_xzzz_d_1_0_1, tg_z_xzzz_f_0_0_1, tg_z_xzzz_g_0_0_0, tg_z_xzzz_g_1_0_0, tg_z_xzzz_p_1_1_1, tg_z_xzzz_s_2_1_1, tg_z_yyy_f_0_0_1, tg_z_yyy_p_1_1_1, tg_z_yyyy_d_1_0_1, tg_z_yyyy_f_0_0_1, tg_z_yyyy_g_0_0_0, tg_z_yyyy_g_1_0_0, tg_z_yyyy_p_1_1_1, tg_z_yyyy_s_2_1_1, tg_z_yyyz_d_1_0_1, tg_z_yyyz_f_0_0_1, tg_z_yyyz_g_0_0_0, tg_z_yyyz_g_1_0_0, tg_z_yyyz_p_1_1_1, tg_z_yyyz_s_2_1_1, tg_z_yyz_f_0_0_1, tg_z_yyz_p_1_1_1, tg_z_yyzz_d_1_0_1, tg_z_yyzz_f_0_0_1, tg_z_yyzz_g_0_0_0, tg_z_yyzz_g_1_0_0, tg_z_yyzz_p_1_1_1, tg_z_yyzz_s_2_1_1, tg_z_yzz_f_0_0_1, tg_z_yzz_p_1_1_1, tg_z_yzzz_d_1_0_1, tg_z_yzzz_f_0_0_1, tg_z_yzzz_g_0_0_0, tg_z_yzzz_g_1_0_0, tg_z_yzzz_p_1_1_1, tg_z_yzzz_s_2_1_1, tg_z_zzz_f_0_0_1, tg_z_zzz_p_1_1_1, tg_z_zzzz_d_1_0_1, tg_z_zzzz_f_0_0_1, tg_z_zzzz_g_0_0_0, tg_z_zzzz_g_1_0_0, tg_z_zzzz_p_1_1_1, tg_z_zzzz_s_2_1_1, tg_zz_xxxx_g_0_0_0, tg_zz_xxxy_g_0_0_0, tg_zz_xxxz_g_0_0_0, tg_zz_xxyy_g_0_0_0, tg_zz_xxyz_g_0_0_0, tg_zz_xxzz_g_0_0_0, tg_zz_xyyy_g_0_0_0, tg_zz_xyyz_g_0_0_0, tg_zz_xyzz_g_0_0_0, tg_zz_xzzz_g_0_0_0, tg_zz_yyyy_g_0_0_0, tg_zz_yyyz_g_0_0_0, tg_zz_yyzz_g_0_0_0, tg_zz_yzzz_g_0_0_0, tg_zz_zzzz_g_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double faz_0 = -(b_exps[i] + c_exp) * fzi_0;

        const double fbzi_0 = b_exps[i] * fzi_0;

        const double f2abz_0 = 2.0 * a_exp * b_exps[i] * fzi_0;

            const double fbi_0 = 1.0 / b_exps[i];

            const double fai_0 = 1.0 / a_exp;

        tg_xx_xxxx_g_0_0_0[i] = -9.0 / 2.0 * tg_0_xxxx_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_xxxx_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_xxxx_g_0_0_0[i] * fzi_0 + tg_0_xxxx_g_1_0_0[i] * fbzi_0 * fbzi_0 + 18.0 * tg_x_xxx_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 18.0 * tg_x_xxx_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_x_xxxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_x_xxxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_x_xxxx_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_x_xxxx_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_xxxx_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxxx_g_0_0_0[i] * a_x * faz_0;

        tg_xx_xxxy_g_0_0_0[i] = -9.0 / 2.0 * tg_0_xxxy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_xxxy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_xxxy_g_0_0_0[i] * fzi_0 + tg_0_xxxy_g_1_0_0[i] * fbzi_0 * fbzi_0 + 27.0 / 2.0 * tg_x_xxy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_x_xxy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_x_xxxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_x_xxxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_x_xxxy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_x_xxxy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_xxxy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxxy_g_0_0_0[i] * a_x * faz_0;

        tg_xx_xxxz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_xxxz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_xxxz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_xxxz_g_0_0_0[i] * fzi_0 + tg_0_xxxz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 27.0 / 2.0 * tg_x_xxz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_x_xxz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_x_xxxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_x_xxxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_x_xxxz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_x_xxxz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_xxxz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxxz_g_0_0_0[i] * a_x * faz_0;

        tg_xx_xxyy_g_0_0_0[i] = -9.0 / 2.0 * tg_0_xxyy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_xxyy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_xxyy_g_0_0_0[i] * fzi_0 + tg_0_xxyy_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 * tg_x_xyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_x_xyy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_x_xxyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_x_xxyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_x_xxyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_x_xxyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_xxyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxyy_g_0_0_0[i] * a_x * faz_0;

        tg_xx_xxyz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_xxyz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_xxyz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_xxyz_g_0_0_0[i] * fzi_0 + tg_0_xxyz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 * tg_x_xyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_x_xyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_x_xxyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_x_xxyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_x_xxyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_x_xxyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_xxyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxyz_g_0_0_0[i] * a_x * faz_0;

        tg_xx_xxzz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_xxzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_xxzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_xxzz_g_0_0_0[i] * fzi_0 + tg_0_xxzz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 * tg_x_xzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_x_xzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_x_xxzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_x_xxzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_x_xxzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_x_xxzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_xxzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxzz_g_0_0_0[i] * a_x * faz_0;

        tg_xx_xyyy_g_0_0_0[i] = -9.0 / 2.0 * tg_0_xyyy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_xyyy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_xyyy_g_0_0_0[i] * fzi_0 + tg_0_xyyy_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 / 2.0 * tg_x_yyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_x_yyy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_x_xyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_x_xyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_x_xyyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_x_xyyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_xyyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xyyy_g_0_0_0[i] * a_x * faz_0;

        tg_xx_xyyz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_xyyz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_xyyz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_xyyz_g_0_0_0[i] * fzi_0 + tg_0_xyyz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 / 2.0 * tg_x_yyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_x_yyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_x_xyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_x_xyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_x_xyyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_x_xyyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_xyyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xyyz_g_0_0_0[i] * a_x * faz_0;

        tg_xx_xyzz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_xyzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_xyzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_xyzz_g_0_0_0[i] * fzi_0 + tg_0_xyzz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 / 2.0 * tg_x_yzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_x_yzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_x_xyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_x_xyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_x_xyzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_x_xyzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_xyzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xyzz_g_0_0_0[i] * a_x * faz_0;

        tg_xx_xzzz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_xzzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_xzzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_xzzz_g_0_0_0[i] * fzi_0 + tg_0_xzzz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 / 2.0 * tg_x_zzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_x_zzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_x_xzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_x_xzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_x_xzzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_x_xzzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_xzzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xzzz_g_0_0_0[i] * a_x * faz_0;

        tg_xx_yyyy_g_0_0_0[i] = -9.0 / 2.0 * tg_0_yyyy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_yyyy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_yyyy_g_0_0_0[i] * fzi_0 + tg_0_yyyy_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_x_yyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_x_yyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_x_yyyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_x_yyyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_yyyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_yyyy_g_0_0_0[i] * a_x * faz_0;

        tg_xx_yyyz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_yyyz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_yyyz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_yyyz_g_0_0_0[i] * fzi_0 + tg_0_yyyz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_x_yyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_x_yyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_x_yyyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_x_yyyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_yyyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_yyyz_g_0_0_0[i] * a_x * faz_0;

        tg_xx_yyzz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_yyzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_yyzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_yyzz_g_0_0_0[i] * fzi_0 + tg_0_yyzz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_x_yyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_x_yyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_x_yyzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_x_yyzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_yyzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_yyzz_g_0_0_0[i] * a_x * faz_0;

        tg_xx_yzzz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_yzzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_yzzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_yzzz_g_0_0_0[i] * fzi_0 + tg_0_yzzz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_x_yzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_x_yzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_x_yzzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_x_yzzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_yzzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_yzzz_g_0_0_0[i] * a_x * faz_0;

        tg_xx_zzzz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_zzzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_zzzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_zzzz_g_0_0_0[i] * fzi_0 + tg_0_zzzz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_x_zzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_x_zzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_x_zzzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_x_zzzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_zzzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_zzzz_g_0_0_0[i] * a_x * faz_0;

        tg_xy_xxxx_g_0_0_0[i] = -9.0 * tg_x_xxxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_x_xxxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_x_xxxx_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_x_xxxx_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_x_xxxx_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_x_xxxx_g_0_0_0[i] * a_y * faz_0;

        tg_xy_xxxy_g_0_0_0[i] = 27.0 / 2.0 * tg_y_xxy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_y_xxy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_y_xxxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_y_xxxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_y_xxxy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_y_xxxy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_y_xxxy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxy_g_0_0_0[i] * a_x * faz_0;

        tg_xy_xxxz_g_0_0_0[i] = -9.0 * tg_x_xxxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_x_xxxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_x_xxxz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_x_xxxz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_x_xxxz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_x_xxxz_g_0_0_0[i] * a_y * faz_0;

        tg_xy_xxyy_g_0_0_0[i] = 9.0 * tg_y_xyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_y_xyy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_y_xxyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_y_xxyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_y_xxyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_y_xxyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_y_xxyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xxyy_g_0_0_0[i] * a_x * faz_0;

        tg_xy_xxyz_g_0_0_0[i] = 9.0 * tg_y_xyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_y_xyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_y_xxyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_y_xxyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_y_xxyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_y_xxyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_y_xxyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xxyz_g_0_0_0[i] * a_x * faz_0;

        tg_xy_xxzz_g_0_0_0[i] = -9.0 * tg_x_xxzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_x_xxzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_x_xxzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_x_xxzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_x_xxzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_x_xxzz_g_0_0_0[i] * a_y * faz_0;

        tg_xy_xyyy_g_0_0_0[i] = 9.0 / 2.0 * tg_y_yyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_y_yyy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_y_xyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_y_xyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_y_xyyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_y_xyyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_y_xyyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xyyy_g_0_0_0[i] * a_x * faz_0;

        tg_xy_xyyz_g_0_0_0[i] = 9.0 / 2.0 * tg_y_yyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_y_yyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_y_xyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_y_xyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_y_xyyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_y_xyyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_y_xyyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xyyz_g_0_0_0[i] * a_x * faz_0;

        tg_xy_xyzz_g_0_0_0[i] = 9.0 / 2.0 * tg_y_yzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_y_yzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_y_xyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_y_xyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_y_xyzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_y_xyzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_y_xyzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xyzz_g_0_0_0[i] * a_x * faz_0;

        tg_xy_xzzz_g_0_0_0[i] = -9.0 * tg_x_xzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_x_xzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_x_xzzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_x_xzzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_x_xzzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_x_xzzz_g_0_0_0[i] * a_y * faz_0;

        tg_xy_yyyy_g_0_0_0[i] = -9.0 * tg_y_yyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_y_yyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_y_yyyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_y_yyyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_y_yyyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_yyyy_g_0_0_0[i] * a_x * faz_0;

        tg_xy_yyyz_g_0_0_0[i] = -9.0 * tg_y_yyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_y_yyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_y_yyyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_y_yyyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_y_yyyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_yyyz_g_0_0_0[i] * a_x * faz_0;

        tg_xy_yyzz_g_0_0_0[i] = -9.0 * tg_y_yyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_y_yyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_y_yyzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_y_yyzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_y_yyzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_yyzz_g_0_0_0[i] * a_x * faz_0;

        tg_xy_yzzz_g_0_0_0[i] = -9.0 * tg_y_yzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_y_yzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_y_yzzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_y_yzzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_y_yzzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_yzzz_g_0_0_0[i] * a_x * faz_0;

        tg_xy_zzzz_g_0_0_0[i] = -9.0 * tg_y_zzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_y_zzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_y_zzzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_y_zzzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_y_zzzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_zzzz_g_0_0_0[i] * a_x * faz_0;

        tg_xz_xxxx_g_0_0_0[i] = -9.0 * tg_x_xxxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_x_xxxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_x_xxxx_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_x_xxxx_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_x_xxxx_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_x_xxxx_g_0_0_0[i] * a_z * faz_0;

        tg_xz_xxxy_g_0_0_0[i] = -9.0 * tg_x_xxxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_x_xxxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_x_xxxy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_x_xxxy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_x_xxxy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_x_xxxy_g_0_0_0[i] * a_z * faz_0;

        tg_xz_xxxz_g_0_0_0[i] = 27.0 / 2.0 * tg_z_xxz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_z_xxz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_z_xxxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_z_xxxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_z_xxxz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_z_xxxz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_z_xxxz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxz_g_0_0_0[i] * a_x * faz_0;

        tg_xz_xxyy_g_0_0_0[i] = -9.0 * tg_x_xxyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_x_xxyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_x_xxyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_x_xxyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_x_xxyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_x_xxyy_g_0_0_0[i] * a_z * faz_0;

        tg_xz_xxyz_g_0_0_0[i] = 9.0 * tg_z_xyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_z_xyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_z_xxyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_z_xxyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_z_xxyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_z_xxyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_z_xxyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xxyz_g_0_0_0[i] * a_x * faz_0;

        tg_xz_xxzz_g_0_0_0[i] = 9.0 * tg_z_xzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_z_xzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_z_xxzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_z_xxzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_z_xxzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_z_xxzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_z_xxzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xxzz_g_0_0_0[i] * a_x * faz_0;

        tg_xz_xyyy_g_0_0_0[i] = -9.0 * tg_x_xyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_x_xyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_x_xyyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_x_xyyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_x_xyyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_x_xyyy_g_0_0_0[i] * a_z * faz_0;

        tg_xz_xyyz_g_0_0_0[i] = 9.0 / 2.0 * tg_z_yyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_z_yyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_z_xyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_z_xyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_z_xyyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_z_xyyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_z_xyyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xyyz_g_0_0_0[i] * a_x * faz_0;

        tg_xz_xyzz_g_0_0_0[i] = 9.0 / 2.0 * tg_z_yzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_z_yzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_z_xyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_z_xyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_z_xyzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_z_xyzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_z_xyzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xyzz_g_0_0_0[i] * a_x * faz_0;

        tg_xz_xzzz_g_0_0_0[i] = 9.0 / 2.0 * tg_z_zzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_z_zzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_z_xzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_z_xzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_z_xzzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_z_xzzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_z_xzzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xzzz_g_0_0_0[i] * a_x * faz_0;

        tg_xz_yyyy_g_0_0_0[i] = -9.0 * tg_z_yyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_z_yyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_z_yyyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_z_yyyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_z_yyyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_yyyy_g_0_0_0[i] * a_x * faz_0;

        tg_xz_yyyz_g_0_0_0[i] = -9.0 * tg_z_yyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_z_yyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_z_yyyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_z_yyyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_z_yyyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_yyyz_g_0_0_0[i] * a_x * faz_0;

        tg_xz_yyzz_g_0_0_0[i] = -9.0 * tg_z_yyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_z_yyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_z_yyzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_z_yyzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_z_yyzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_yyzz_g_0_0_0[i] * a_x * faz_0;

        tg_xz_yzzz_g_0_0_0[i] = -9.0 * tg_z_yzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_z_yzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_z_yzzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_z_yzzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_z_yzzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_yzzz_g_0_0_0[i] * a_x * faz_0;

        tg_xz_zzzz_g_0_0_0[i] = -9.0 * tg_z_zzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_z_zzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_z_zzzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_z_zzzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_z_zzzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_zzzz_g_0_0_0[i] * a_x * faz_0;

        tg_yy_xxxx_g_0_0_0[i] = -9.0 / 2.0 * tg_0_xxxx_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_xxxx_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_xxxx_g_0_0_0[i] * fzi_0 + tg_0_xxxx_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_y_xxxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_y_xxxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_y_xxxx_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_y_xxxx_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_xxxx_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxx_g_0_0_0[i] * a_y * faz_0;

        tg_yy_xxxy_g_0_0_0[i] = -9.0 / 2.0 * tg_0_xxxy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_xxxy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_xxxy_g_0_0_0[i] * fzi_0 + tg_0_xxxy_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 / 2.0 * tg_y_xxx_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_y_xxx_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_y_xxxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_y_xxxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_y_xxxy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_y_xxxy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_xxxy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxy_g_0_0_0[i] * a_y * faz_0;

        tg_yy_xxxz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_xxxz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_xxxz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_xxxz_g_0_0_0[i] * fzi_0 + tg_0_xxxz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_y_xxxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_y_xxxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_y_xxxz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_y_xxxz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_xxxz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxz_g_0_0_0[i] * a_y * faz_0;

        tg_yy_xxyy_g_0_0_0[i] = -9.0 / 2.0 * tg_0_xxyy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_xxyy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_xxyy_g_0_0_0[i] * fzi_0 + tg_0_xxyy_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 * tg_y_xxy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_y_xxy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_y_xxyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_y_xxyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_y_xxyy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_y_xxyy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_xxyy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxyy_g_0_0_0[i] * a_y * faz_0;

        tg_yy_xxyz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_xxyz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_xxyz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_xxyz_g_0_0_0[i] * fzi_0 + tg_0_xxyz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 / 2.0 * tg_y_xxz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_y_xxz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_y_xxyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_y_xxyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_y_xxyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_y_xxyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_xxyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxyz_g_0_0_0[i] * a_y * faz_0;

        tg_yy_xxzz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_xxzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_xxzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_xxzz_g_0_0_0[i] * fzi_0 + tg_0_xxzz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_y_xxzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_y_xxzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_y_xxzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_y_xxzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_xxzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxzz_g_0_0_0[i] * a_y * faz_0;

        tg_yy_xyyy_g_0_0_0[i] = -9.0 / 2.0 * tg_0_xyyy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_xyyy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_xyyy_g_0_0_0[i] * fzi_0 + tg_0_xyyy_g_1_0_0[i] * fbzi_0 * fbzi_0 + 27.0 / 2.0 * tg_y_xyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_y_xyy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_y_xyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_y_xyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_y_xyyy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_y_xyyy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_xyyy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xyyy_g_0_0_0[i] * a_y * faz_0;

        tg_yy_xyyz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_xyyz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_xyyz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_xyyz_g_0_0_0[i] * fzi_0 + tg_0_xyyz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 * tg_y_xyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_y_xyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_y_xyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_y_xyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_y_xyyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_y_xyyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_xyyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xyyz_g_0_0_0[i] * a_y * faz_0;

        tg_yy_xyzz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_xyzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_xyzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_xyzz_g_0_0_0[i] * fzi_0 + tg_0_xyzz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 / 2.0 * tg_y_xzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_y_xzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_y_xyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_y_xyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_y_xyzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_y_xyzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_xyzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xyzz_g_0_0_0[i] * a_y * faz_0;

        tg_yy_xzzz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_xzzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_xzzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_xzzz_g_0_0_0[i] * fzi_0 + tg_0_xzzz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_y_xzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_y_xzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_y_xzzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_y_xzzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_xzzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xzzz_g_0_0_0[i] * a_y * faz_0;

        tg_yy_yyyy_g_0_0_0[i] = -9.0 / 2.0 * tg_0_yyyy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_yyyy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_yyyy_g_0_0_0[i] * fzi_0 + tg_0_yyyy_g_1_0_0[i] * fbzi_0 * fbzi_0 + 18.0 * tg_y_yyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 18.0 * tg_y_yyy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_y_yyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_y_yyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_y_yyyy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_y_yyyy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_yyyy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_yyyy_g_0_0_0[i] * a_y * faz_0;

        tg_yy_yyyz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_yyyz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_yyyz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_yyyz_g_0_0_0[i] * fzi_0 + tg_0_yyyz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 27.0 / 2.0 * tg_y_yyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_y_yyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_y_yyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_y_yyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_y_yyyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_y_yyyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_yyyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_yyyz_g_0_0_0[i] * a_y * faz_0;

        tg_yy_yyzz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_yyzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_yyzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_yyzz_g_0_0_0[i] * fzi_0 + tg_0_yyzz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 * tg_y_yzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_y_yzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_y_yyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_y_yyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_y_yyzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_y_yyzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_yyzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_yyzz_g_0_0_0[i] * a_y * faz_0;

        tg_yy_yzzz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_yzzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_yzzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_yzzz_g_0_0_0[i] * fzi_0 + tg_0_yzzz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 / 2.0 * tg_y_zzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_y_zzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_y_yzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_y_yzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_y_yzzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_y_yzzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_yzzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_yzzz_g_0_0_0[i] * a_y * faz_0;

        tg_yy_zzzz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_zzzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_zzzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_zzzz_g_0_0_0[i] * fzi_0 + tg_0_zzzz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_y_zzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_y_zzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_y_zzzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_y_zzzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_zzzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_zzzz_g_0_0_0[i] * a_y * faz_0;

        tg_yz_xxxx_g_0_0_0[i] = -9.0 * tg_z_xxxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_z_xxxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_z_xxxx_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_z_xxxx_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_z_xxxx_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxx_g_0_0_0[i] * a_y * faz_0;

        tg_yz_xxxy_g_0_0_0[i] = -9.0 * tg_y_xxxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_y_xxxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_y_xxxy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_y_xxxy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_y_xxxy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_y_xxxy_g_0_0_0[i] * a_z * faz_0;

        tg_yz_xxxz_g_0_0_0[i] = -9.0 * tg_z_xxxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_z_xxxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_z_xxxz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_z_xxxz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_z_xxxz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxz_g_0_0_0[i] * a_y * faz_0;

        tg_yz_xxyy_g_0_0_0[i] = -9.0 * tg_y_xxyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_y_xxyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_y_xxyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_y_xxyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_y_xxyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_y_xxyy_g_0_0_0[i] * a_z * faz_0;

        tg_yz_xxyz_g_0_0_0[i] = 9.0 / 2.0 * tg_z_xxz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_z_xxz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_z_xxyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_z_xxyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_z_xxyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_z_xxyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_z_xxyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xxyz_g_0_0_0[i] * a_y * faz_0;

        tg_yz_xxzz_g_0_0_0[i] = -9.0 * tg_z_xxzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_z_xxzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_z_xxzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_z_xxzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_z_xxzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xxzz_g_0_0_0[i] * a_y * faz_0;

        tg_yz_xyyy_g_0_0_0[i] = -9.0 * tg_y_xyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_y_xyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_y_xyyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_y_xyyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_y_xyyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_y_xyyy_g_0_0_0[i] * a_z * faz_0;

        tg_yz_xyyz_g_0_0_0[i] = 9.0 * tg_z_xyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_z_xyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_z_xyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_z_xyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_z_xyyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_z_xyyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_z_xyyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xyyz_g_0_0_0[i] * a_y * faz_0;

        tg_yz_xyzz_g_0_0_0[i] = 9.0 / 2.0 * tg_z_xzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_z_xzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_z_xyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_z_xyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_z_xyzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_z_xyzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_z_xyzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xyzz_g_0_0_0[i] * a_y * faz_0;

        tg_yz_xzzz_g_0_0_0[i] = -9.0 * tg_z_xzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_z_xzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_z_xzzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_z_xzzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_z_xzzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xzzz_g_0_0_0[i] * a_y * faz_0;

        tg_yz_yyyy_g_0_0_0[i] = -9.0 * tg_y_yyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_y_yyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_y_yyyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_y_yyyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_y_yyyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_y_yyyy_g_0_0_0[i] * a_z * faz_0;

        tg_yz_yyyz_g_0_0_0[i] = 27.0 / 2.0 * tg_z_yyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_z_yyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_z_yyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_z_yyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_z_yyyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_z_yyyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_z_yyyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_yyyz_g_0_0_0[i] * a_y * faz_0;

        tg_yz_yyzz_g_0_0_0[i] = 9.0 * tg_z_yzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_z_yzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_z_yyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_z_yyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_z_yyzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_z_yyzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_z_yyzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_yyzz_g_0_0_0[i] * a_y * faz_0;

        tg_yz_yzzz_g_0_0_0[i] = 9.0 / 2.0 * tg_z_zzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_z_zzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_z_yzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_z_yzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_z_yzzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_z_yzzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_z_yzzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_yzzz_g_0_0_0[i] * a_y * faz_0;

        tg_yz_zzzz_g_0_0_0[i] = -9.0 * tg_z_zzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_z_zzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_z_zzzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_z_zzzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_z_zzzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_zzzz_g_0_0_0[i] * a_y * faz_0;

        tg_zz_xxxx_g_0_0_0[i] = -9.0 / 2.0 * tg_0_xxxx_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_xxxx_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_xxxx_g_0_0_0[i] * fzi_0 + tg_0_xxxx_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_z_xxxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_z_xxxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_z_xxxx_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_z_xxxx_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_xxxx_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxx_g_0_0_0[i] * a_z * faz_0;

        tg_zz_xxxy_g_0_0_0[i] = -9.0 / 2.0 * tg_0_xxxy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_xxxy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_xxxy_g_0_0_0[i] * fzi_0 + tg_0_xxxy_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_z_xxxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_z_xxxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_z_xxxy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_z_xxxy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_xxxy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxy_g_0_0_0[i] * a_z * faz_0;

        tg_zz_xxxz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_xxxz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_xxxz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_xxxz_g_0_0_0[i] * fzi_0 + tg_0_xxxz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 / 2.0 * tg_z_xxx_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_z_xxx_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_z_xxxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_z_xxxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_z_xxxz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_z_xxxz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_xxxz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxxz_g_0_0_0[i] * a_z * faz_0;

        tg_zz_xxyy_g_0_0_0[i] = -9.0 / 2.0 * tg_0_xxyy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_xxyy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_xxyy_g_0_0_0[i] * fzi_0 + tg_0_xxyy_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_z_xxyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_z_xxyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_z_xxyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_z_xxyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_xxyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxyy_g_0_0_0[i] * a_z * faz_0;

        tg_zz_xxyz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_xxyz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_xxyz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_xxyz_g_0_0_0[i] * fzi_0 + tg_0_xxyz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 / 2.0 * tg_z_xxy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_z_xxy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_z_xxyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_z_xxyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_z_xxyz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_z_xxyz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_xxyz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxyz_g_0_0_0[i] * a_z * faz_0;

        tg_zz_xxzz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_xxzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_xxzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_xxzz_g_0_0_0[i] * fzi_0 + tg_0_xxzz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 * tg_z_xxz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_z_xxz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_z_xxzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_z_xxzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_z_xxzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_z_xxzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_xxzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxzz_g_0_0_0[i] * a_z * faz_0;

        tg_zz_xyyy_g_0_0_0[i] = -9.0 / 2.0 * tg_0_xyyy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_xyyy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_xyyy_g_0_0_0[i] * fzi_0 + tg_0_xyyy_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_z_xyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_z_xyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_z_xyyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_z_xyyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_xyyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xyyy_g_0_0_0[i] * a_z * faz_0;

        tg_zz_xyyz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_xyyz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_xyyz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_xyyz_g_0_0_0[i] * fzi_0 + tg_0_xyyz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 / 2.0 * tg_z_xyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_z_xyy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_z_xyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_z_xyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_z_xyyz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_z_xyyz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_xyyz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xyyz_g_0_0_0[i] * a_z * faz_0;

        tg_zz_xyzz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_xyzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_xyzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_xyzz_g_0_0_0[i] * fzi_0 + tg_0_xyzz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 * tg_z_xyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_z_xyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_z_xyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_z_xyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_z_xyzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_z_xyzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_xyzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xyzz_g_0_0_0[i] * a_z * faz_0;

        tg_zz_xzzz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_xzzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_xzzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_xzzz_g_0_0_0[i] * fzi_0 + tg_0_xzzz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 27.0 / 2.0 * tg_z_xzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_z_xzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_z_xzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_z_xzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_z_xzzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_z_xzzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_xzzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xzzz_g_0_0_0[i] * a_z * faz_0;

        tg_zz_yyyy_g_0_0_0[i] = -9.0 / 2.0 * tg_0_yyyy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_yyyy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_yyyy_g_0_0_0[i] * fzi_0 + tg_0_yyyy_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_z_yyyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_z_yyyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_z_yyyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_z_yyyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_yyyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_yyyy_g_0_0_0[i] * a_z * faz_0;

        tg_zz_yyyz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_yyyz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_yyyz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_yyyz_g_0_0_0[i] * fzi_0 + tg_0_yyyz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 / 2.0 * tg_z_yyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_z_yyy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_z_yyyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_z_yyyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_z_yyyz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_z_yyyz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_yyyz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_yyyz_g_0_0_0[i] * a_z * faz_0;

        tg_zz_yyzz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_yyzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_yyzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_yyzz_g_0_0_0[i] * fzi_0 + tg_0_yyzz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 * tg_z_yyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_z_yyz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_z_yyzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_z_yyzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_z_yyzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_z_yyzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_yyzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_yyzz_g_0_0_0[i] * a_z * faz_0;

        tg_zz_yzzz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_yzzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_yzzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_yzzz_g_0_0_0[i] * fzi_0 + tg_0_yzzz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 27.0 / 2.0 * tg_z_yzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_z_yzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_z_yzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_z_yzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_z_yzzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_z_yzzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_yzzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_yzzz_g_0_0_0[i] * a_z * faz_0;

        tg_zz_zzzz_g_0_0_0[i] = -9.0 / 2.0 * tg_0_zzzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_0_zzzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_0_zzzz_g_0_0_0[i] * fzi_0 + tg_0_zzzz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 18.0 * tg_z_zzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 18.0 * tg_z_zzz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_z_zzzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_z_zzzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_z_zzzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_z_zzzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_zzzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_zzzz_g_0_0_0[i] * a_z * faz_0;
    }

    if (p > 0)
    {
        const double fp_0 = (double)p;

        // Set up components of auxiliary buffer : SG

        auto tg_0_xxxx_g_0_0_1 = pbuffer.data(idx_sg_g_0_0_1);

        auto tg_0_xxxy_g_0_0_1 = pbuffer.data(idx_sg_g_0_0_1 + 1);

        auto tg_0_xxxz_g_0_0_1 = pbuffer.data(idx_sg_g_0_0_1 + 2);

        auto tg_0_xxyy_g_0_0_1 = pbuffer.data(idx_sg_g_0_0_1 + 3);

        auto tg_0_xxyz_g_0_0_1 = pbuffer.data(idx_sg_g_0_0_1 + 4);

        auto tg_0_xxzz_g_0_0_1 = pbuffer.data(idx_sg_g_0_0_1 + 5);

        auto tg_0_xyyy_g_0_0_1 = pbuffer.data(idx_sg_g_0_0_1 + 6);

        auto tg_0_xyyz_g_0_0_1 = pbuffer.data(idx_sg_g_0_0_1 + 7);

        auto tg_0_xyzz_g_0_0_1 = pbuffer.data(idx_sg_g_0_0_1 + 8);

        auto tg_0_xzzz_g_0_0_1 = pbuffer.data(idx_sg_g_0_0_1 + 9);

        auto tg_0_yyyy_g_0_0_1 = pbuffer.data(idx_sg_g_0_0_1 + 10);

        auto tg_0_yyyz_g_0_0_1 = pbuffer.data(idx_sg_g_0_0_1 + 11);

        auto tg_0_yyzz_g_0_0_1 = pbuffer.data(idx_sg_g_0_0_1 + 12);

        auto tg_0_yzzz_g_0_0_1 = pbuffer.data(idx_sg_g_0_0_1 + 13);

        auto tg_0_zzzz_g_0_0_1 = pbuffer.data(idx_sg_g_0_0_1 + 14);

        // Set up components of auxiliary buffer : PG

        auto tg_x_xxxx_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1);

        auto tg_x_xxxy_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 1);

        auto tg_x_xxxz_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 2);

        auto tg_x_xxyy_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 3);

        auto tg_x_xxyz_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 4);

        auto tg_x_xxzz_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 5);

        auto tg_x_xyyy_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 6);

        auto tg_x_xyyz_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 7);

        auto tg_x_xyzz_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 8);

        auto tg_x_xzzz_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 9);

        auto tg_x_yyyy_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 10);

        auto tg_x_yyyz_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 11);

        auto tg_x_yyzz_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 12);

        auto tg_x_yzzz_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 13);

        auto tg_x_zzzz_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 14);

        auto tg_y_xxxx_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 15);

        auto tg_y_xxxy_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 16);

        auto tg_y_xxxz_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 17);

        auto tg_y_xxyy_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 18);

        auto tg_y_xxyz_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 19);

        auto tg_y_xxzz_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 20);

        auto tg_y_xyyy_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 21);

        auto tg_y_xyyz_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 22);

        auto tg_y_xyzz_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 23);

        auto tg_y_xzzz_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 24);

        auto tg_y_yyyy_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 25);

        auto tg_y_yyyz_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 26);

        auto tg_y_yyzz_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 27);

        auto tg_y_yzzz_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 28);

        auto tg_y_zzzz_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 29);

        auto tg_z_xxxx_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 30);

        auto tg_z_xxxy_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 31);

        auto tg_z_xxxz_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 32);

        auto tg_z_xxyy_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 33);

        auto tg_z_xxyz_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 34);

        auto tg_z_xxzz_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 35);

        auto tg_z_xyyy_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 36);

        auto tg_z_xyyz_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 37);

        auto tg_z_xyzz_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 38);

        auto tg_z_xzzz_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 39);

        auto tg_z_yyyy_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 40);

        auto tg_z_yyyz_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 41);

        auto tg_z_yyzz_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 42);

        auto tg_z_yzzz_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 43);

        auto tg_z_zzzz_g_0_0_1 = pbuffer.data(idx_pg_g_0_0_1 + 44);

        #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_0_xxxx_g_0_0_1, tg_0_xxxy_g_0_0_1, tg_0_xxxz_g_0_0_1, tg_0_xxyy_g_0_0_1, tg_0_xxyz_g_0_0_1, tg_0_xxzz_g_0_0_1, tg_0_xyyy_g_0_0_1, tg_0_xyyz_g_0_0_1, tg_0_xyzz_g_0_0_1, tg_0_xzzz_g_0_0_1, tg_0_yyyy_g_0_0_1, tg_0_yyyz_g_0_0_1, tg_0_yyzz_g_0_0_1, tg_0_yzzz_g_0_0_1, tg_0_zzzz_g_0_0_1, tg_x_xxxx_g_0_0_1, tg_x_xxxy_g_0_0_1, tg_x_xxxz_g_0_0_1, tg_x_xxyy_g_0_0_1, tg_x_xxyz_g_0_0_1, tg_x_xxzz_g_0_0_1, tg_x_xyyy_g_0_0_1, tg_x_xyyz_g_0_0_1, tg_x_xyzz_g_0_0_1, tg_x_xzzz_g_0_0_1, tg_x_yyyy_g_0_0_1, tg_x_yyyz_g_0_0_1, tg_x_yyzz_g_0_0_1, tg_x_yzzz_g_0_0_1, tg_x_zzzz_g_0_0_1, tg_xx_xxxx_g_0_0_0, tg_xx_xxxy_g_0_0_0, tg_xx_xxxz_g_0_0_0, tg_xx_xxyy_g_0_0_0, tg_xx_xxyz_g_0_0_0, tg_xx_xxzz_g_0_0_0, tg_xx_xyyy_g_0_0_0, tg_xx_xyyz_g_0_0_0, tg_xx_xyzz_g_0_0_0, tg_xx_xzzz_g_0_0_0, tg_xx_yyyy_g_0_0_0, tg_xx_yyyz_g_0_0_0, tg_xx_yyzz_g_0_0_0, tg_xx_yzzz_g_0_0_0, tg_xx_zzzz_g_0_0_0, tg_xy_xxxx_g_0_0_0, tg_xy_xxxy_g_0_0_0, tg_xy_xxxz_g_0_0_0, tg_xy_xxyy_g_0_0_0, tg_xy_xxyz_g_0_0_0, tg_xy_xxzz_g_0_0_0, tg_xy_xyyy_g_0_0_0, tg_xy_xyyz_g_0_0_0, tg_xy_xyzz_g_0_0_0, tg_xy_xzzz_g_0_0_0, tg_xy_yyyy_g_0_0_0, tg_xy_yyyz_g_0_0_0, tg_xy_yyzz_g_0_0_0, tg_xy_yzzz_g_0_0_0, tg_xy_zzzz_g_0_0_0, tg_xz_xxxx_g_0_0_0, tg_xz_xxxy_g_0_0_0, tg_xz_xxxz_g_0_0_0, tg_xz_xxyy_g_0_0_0, tg_xz_xxyz_g_0_0_0, tg_xz_xxzz_g_0_0_0, tg_xz_xyyy_g_0_0_0, tg_xz_xyyz_g_0_0_0, tg_xz_xyzz_g_0_0_0, tg_xz_xzzz_g_0_0_0, tg_xz_yyyy_g_0_0_0, tg_xz_yyyz_g_0_0_0, tg_xz_yyzz_g_0_0_0, tg_xz_yzzz_g_0_0_0, tg_xz_zzzz_g_0_0_0, tg_y_xxxx_g_0_0_1, tg_y_xxxy_g_0_0_1, tg_y_xxxz_g_0_0_1, tg_y_xxyy_g_0_0_1, tg_y_xxyz_g_0_0_1, tg_y_xxzz_g_0_0_1, tg_y_xyyy_g_0_0_1, tg_y_xyyz_g_0_0_1, tg_y_xyzz_g_0_0_1, tg_y_xzzz_g_0_0_1, tg_y_yyyy_g_0_0_1, tg_y_yyyz_g_0_0_1, tg_y_yyzz_g_0_0_1, tg_y_yzzz_g_0_0_1, tg_y_zzzz_g_0_0_1, tg_yy_xxxx_g_0_0_0, tg_yy_xxxy_g_0_0_0, tg_yy_xxxz_g_0_0_0, tg_yy_xxyy_g_0_0_0, tg_yy_xxyz_g_0_0_0, tg_yy_xxzz_g_0_0_0, tg_yy_xyyy_g_0_0_0, tg_yy_xyyz_g_0_0_0, tg_yy_xyzz_g_0_0_0, tg_yy_xzzz_g_0_0_0, tg_yy_yyyy_g_0_0_0, tg_yy_yyyz_g_0_0_0, tg_yy_yyzz_g_0_0_0, tg_yy_yzzz_g_0_0_0, tg_yy_zzzz_g_0_0_0, tg_yz_xxxx_g_0_0_0, tg_yz_xxxy_g_0_0_0, tg_yz_xxxz_g_0_0_0, tg_yz_xxyy_g_0_0_0, tg_yz_xxyz_g_0_0_0, tg_yz_xxzz_g_0_0_0, tg_yz_xyyy_g_0_0_0, tg_yz_xyyz_g_0_0_0, tg_yz_xyzz_g_0_0_0, tg_yz_xzzz_g_0_0_0, tg_yz_yyyy_g_0_0_0, tg_yz_yyyz_g_0_0_0, tg_yz_yyzz_g_0_0_0, tg_yz_yzzz_g_0_0_0, tg_yz_zzzz_g_0_0_0, tg_z_xxxx_g_0_0_1, tg_z_xxxy_g_0_0_1, tg_z_xxxz_g_0_0_1, tg_z_xxyy_g_0_0_1, tg_z_xxyz_g_0_0_1, tg_z_xxzz_g_0_0_1, tg_z_xyyy_g_0_0_1, tg_z_xyyz_g_0_0_1, tg_z_xyzz_g_0_0_1, tg_z_xzzz_g_0_0_1, tg_z_yyyy_g_0_0_1, tg_z_yyyz_g_0_0_1, tg_z_yyzz_g_0_0_1, tg_z_yzzz_g_0_0_1, tg_z_zzzz_g_0_0_1, tg_zz_xxxx_g_0_0_0, tg_zz_xxxy_g_0_0_0, tg_zz_xxxz_g_0_0_0, tg_zz_xxyy_g_0_0_0, tg_zz_xxyz_g_0_0_0, tg_zz_xxzz_g_0_0_0, tg_zz_xyyy_g_0_0_0, tg_zz_xyyz_g_0_0_0, tg_zz_xyzz_g_0_0_0, tg_zz_xzzz_g_0_0_0, tg_zz_yyyy_g_0_0_0, tg_zz_yyyz_g_0_0_0, tg_zz_yyzz_g_0_0_0, tg_zz_yzzz_g_0_0_0, tg_zz_zzzz_g_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fai_0 = 1.0 / a_exp;

            tg_xx_xxxx_g_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxx_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxxx_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xxxy_g_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxxy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xxxz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxxz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xxyy_g_0_0_0[i] += 1.0 / 2.0 * tg_0_xxyy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xxyz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_xxyz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xxzz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_xxzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xyyy_g_0_0_0[i] += 1.0 / 2.0 * tg_0_xyyy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xyyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xyyz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_xyyz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xyyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xyzz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_xyzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xyzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xzzz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_xzzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xzzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_yyyy_g_0_0_0[i] += 1.0 / 2.0 * tg_0_yyyy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_yyyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_yyyz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_yyyz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_yyyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_yyzz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_yyzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_yyzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_yzzz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_yzzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_yzzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_zzzz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_zzzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_zzzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxxx_g_0_0_0[i] += tg_y_xxxx_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxxy_g_0_0_0[i] += tg_y_xxxy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxxz_g_0_0_0[i] += tg_y_xxxz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxyy_g_0_0_0[i] += tg_y_xxyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxyz_g_0_0_0[i] += tg_y_xxyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxzz_g_0_0_0[i] += tg_y_xxzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xyyy_g_0_0_0[i] += tg_y_xyyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xyyz_g_0_0_0[i] += tg_y_xyyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xyzz_g_0_0_0[i] += tg_y_xyzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xzzz_g_0_0_0[i] += tg_y_xzzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_yyyy_g_0_0_0[i] += tg_y_yyyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_yyyz_g_0_0_0[i] += tg_y_yyyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_yyzz_g_0_0_0[i] += tg_y_yyzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_yzzz_g_0_0_0[i] += tg_y_yzzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_zzzz_g_0_0_0[i] += tg_y_zzzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxxx_g_0_0_0[i] += tg_z_xxxx_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxxy_g_0_0_0[i] += tg_z_xxxy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxxz_g_0_0_0[i] += tg_z_xxxz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxyy_g_0_0_0[i] += tg_z_xxyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxyz_g_0_0_0[i] += tg_z_xxyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxzz_g_0_0_0[i] += tg_z_xxzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xyyy_g_0_0_0[i] += tg_z_xyyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xyyz_g_0_0_0[i] += tg_z_xyyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xyzz_g_0_0_0[i] += tg_z_xyzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xzzz_g_0_0_0[i] += tg_z_xzzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_yyyy_g_0_0_0[i] += tg_z_yyyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_yyyz_g_0_0_0[i] += tg_z_yyyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_yyzz_g_0_0_0[i] += tg_z_yyzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_yzzz_g_0_0_0[i] += tg_z_yzzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_zzzz_g_0_0_0[i] += tg_z_zzzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_yy_xxxx_g_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxx_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxxx_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xxxy_g_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxxy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xxxz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxxz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xxyy_g_0_0_0[i] += 1.0 / 2.0 * tg_0_xxyy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xxyz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_xxyz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xxzz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_xxzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xyyy_g_0_0_0[i] += 1.0 / 2.0 * tg_0_xyyy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xyyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xyyz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_xyyz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xyyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xyzz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_xyzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xyzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xzzz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_xzzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xzzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_yyyy_g_0_0_0[i] += 1.0 / 2.0 * tg_0_yyyy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_yyyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_yyyz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_yyyz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_yyyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_yyzz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_yyzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_yyzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_yzzz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_yzzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_yzzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_zzzz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_zzzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_zzzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxxx_g_0_0_0[i] += tg_z_xxxx_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxxy_g_0_0_0[i] += tg_z_xxxy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxxz_g_0_0_0[i] += tg_z_xxxz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxyy_g_0_0_0[i] += tg_z_xxyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxyz_g_0_0_0[i] += tg_z_xxyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxzz_g_0_0_0[i] += tg_z_xxzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xyyy_g_0_0_0[i] += tg_z_xyyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xyyz_g_0_0_0[i] += tg_z_xyyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xyzz_g_0_0_0[i] += tg_z_xyzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xzzz_g_0_0_0[i] += tg_z_xzzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_yyyy_g_0_0_0[i] += tg_z_yyyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_yyyz_g_0_0_0[i] += tg_z_yyyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_yyzz_g_0_0_0[i] += tg_z_yyzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_yzzz_g_0_0_0[i] += tg_z_yzzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_zzzz_g_0_0_0[i] += tg_z_zzzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_zz_xxxx_g_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxx_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxxx_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xxxy_g_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxxy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xxxz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_xxxz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxxz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xxyy_g_0_0_0[i] += 1.0 / 2.0 * tg_0_xxyy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxyy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xxyz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_xxyz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxyz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xxzz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_xxzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xyyy_g_0_0_0[i] += 1.0 / 2.0 * tg_0_xyyy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xyyy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xyyz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_xyyz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xyyz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xyzz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_xyzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xyzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xzzz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_xzzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xzzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_yyyy_g_0_0_0[i] += 1.0 / 2.0 * tg_0_yyyy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_yyyy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_yyyz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_yyyz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_yyyz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_yyzz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_yyzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_yyzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_yzzz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_yzzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_yzzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_zzzz_g_0_0_0[i] += 1.0 / 2.0 * tg_0_zzzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_zzzz_g_0_0_1[i] * fai_0 * a_z * fp_0;
        }
    }
}

} // t2pecp namespace

