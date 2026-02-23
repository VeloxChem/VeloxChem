#include "ProjectedCorePotentialPrimRecPGForF.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_pg_f(CSimdArray<double>& pbuffer, 
                                        const size_t idx_pg_f_0_0_0,
                                        const size_t idx_sg_f_0_0_0,
                                        const size_t idx_sf_d_0_0_1,
                                        const size_t idx_sg_d_0_0_1,
                                        const size_t idx_sg_f_1_0_0,
                                        const size_t idx_sg_p_1_0_1,
                                        const size_t idx_sf_s_1_1_1,
                                        const size_t idx_sg_s_1_1_1,
                                        const int p,
                                        const size_t idx_sg_f_0_0_1,
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

    auto tg_0_xxxx_f_0_0_0 = pbuffer.data(idx_sg_f_0_0_0);

    auto tg_0_xxxy_f_0_0_0 = pbuffer.data(idx_sg_f_0_0_0 + 1);

    auto tg_0_xxxz_f_0_0_0 = pbuffer.data(idx_sg_f_0_0_0 + 2);

    auto tg_0_xxyy_f_0_0_0 = pbuffer.data(idx_sg_f_0_0_0 + 3);

    auto tg_0_xxyz_f_0_0_0 = pbuffer.data(idx_sg_f_0_0_0 + 4);

    auto tg_0_xxzz_f_0_0_0 = pbuffer.data(idx_sg_f_0_0_0 + 5);

    auto tg_0_xyyy_f_0_0_0 = pbuffer.data(idx_sg_f_0_0_0 + 6);

    auto tg_0_xyyz_f_0_0_0 = pbuffer.data(idx_sg_f_0_0_0 + 7);

    auto tg_0_xyzz_f_0_0_0 = pbuffer.data(idx_sg_f_0_0_0 + 8);

    auto tg_0_xzzz_f_0_0_0 = pbuffer.data(idx_sg_f_0_0_0 + 9);

    auto tg_0_yyyy_f_0_0_0 = pbuffer.data(idx_sg_f_0_0_0 + 10);

    auto tg_0_yyyz_f_0_0_0 = pbuffer.data(idx_sg_f_0_0_0 + 11);

    auto tg_0_yyzz_f_0_0_0 = pbuffer.data(idx_sg_f_0_0_0 + 12);

    auto tg_0_yzzz_f_0_0_0 = pbuffer.data(idx_sg_f_0_0_0 + 13);

    auto tg_0_zzzz_f_0_0_0 = pbuffer.data(idx_sg_f_0_0_0 + 14);

    // Set up components of auxiliary buffer : SF

    auto tg_0_xxx_d_0_0_1 = pbuffer.data(idx_sf_d_0_0_1);

    auto tg_0_xxy_d_0_0_1 = pbuffer.data(idx_sf_d_0_0_1 + 1);

    auto tg_0_xxz_d_0_0_1 = pbuffer.data(idx_sf_d_0_0_1 + 2);

    auto tg_0_xyy_d_0_0_1 = pbuffer.data(idx_sf_d_0_0_1 + 3);

    auto tg_0_xyz_d_0_0_1 = pbuffer.data(idx_sf_d_0_0_1 + 4);

    auto tg_0_xzz_d_0_0_1 = pbuffer.data(idx_sf_d_0_0_1 + 5);

    auto tg_0_yyy_d_0_0_1 = pbuffer.data(idx_sf_d_0_0_1 + 6);

    auto tg_0_yyz_d_0_0_1 = pbuffer.data(idx_sf_d_0_0_1 + 7);

    auto tg_0_yzz_d_0_0_1 = pbuffer.data(idx_sf_d_0_0_1 + 8);

    auto tg_0_zzz_d_0_0_1 = pbuffer.data(idx_sf_d_0_0_1 + 9);

    // Set up components of auxiliary buffer : SG

    auto tg_0_xxxx_d_0_0_1 = pbuffer.data(idx_sg_d_0_0_1);

    auto tg_0_xxxy_d_0_0_1 = pbuffer.data(idx_sg_d_0_0_1 + 1);

    auto tg_0_xxxz_d_0_0_1 = pbuffer.data(idx_sg_d_0_0_1 + 2);

    auto tg_0_xxyy_d_0_0_1 = pbuffer.data(idx_sg_d_0_0_1 + 3);

    auto tg_0_xxyz_d_0_0_1 = pbuffer.data(idx_sg_d_0_0_1 + 4);

    auto tg_0_xxzz_d_0_0_1 = pbuffer.data(idx_sg_d_0_0_1 + 5);

    auto tg_0_xyyy_d_0_0_1 = pbuffer.data(idx_sg_d_0_0_1 + 6);

    auto tg_0_xyyz_d_0_0_1 = pbuffer.data(idx_sg_d_0_0_1 + 7);

    auto tg_0_xyzz_d_0_0_1 = pbuffer.data(idx_sg_d_0_0_1 + 8);

    auto tg_0_xzzz_d_0_0_1 = pbuffer.data(idx_sg_d_0_0_1 + 9);

    auto tg_0_yyyy_d_0_0_1 = pbuffer.data(idx_sg_d_0_0_1 + 10);

    auto tg_0_yyyz_d_0_0_1 = pbuffer.data(idx_sg_d_0_0_1 + 11);

    auto tg_0_yyzz_d_0_0_1 = pbuffer.data(idx_sg_d_0_0_1 + 12);

    auto tg_0_yzzz_d_0_0_1 = pbuffer.data(idx_sg_d_0_0_1 + 13);

    auto tg_0_zzzz_d_0_0_1 = pbuffer.data(idx_sg_d_0_0_1 + 14);

    // Set up components of auxiliary buffer : SG

    auto tg_0_xxxx_f_1_0_0 = pbuffer.data(idx_sg_f_1_0_0);

    auto tg_0_xxxy_f_1_0_0 = pbuffer.data(idx_sg_f_1_0_0 + 1);

    auto tg_0_xxxz_f_1_0_0 = pbuffer.data(idx_sg_f_1_0_0 + 2);

    auto tg_0_xxyy_f_1_0_0 = pbuffer.data(idx_sg_f_1_0_0 + 3);

    auto tg_0_xxyz_f_1_0_0 = pbuffer.data(idx_sg_f_1_0_0 + 4);

    auto tg_0_xxzz_f_1_0_0 = pbuffer.data(idx_sg_f_1_0_0 + 5);

    auto tg_0_xyyy_f_1_0_0 = pbuffer.data(idx_sg_f_1_0_0 + 6);

    auto tg_0_xyyz_f_1_0_0 = pbuffer.data(idx_sg_f_1_0_0 + 7);

    auto tg_0_xyzz_f_1_0_0 = pbuffer.data(idx_sg_f_1_0_0 + 8);

    auto tg_0_xzzz_f_1_0_0 = pbuffer.data(idx_sg_f_1_0_0 + 9);

    auto tg_0_yyyy_f_1_0_0 = pbuffer.data(idx_sg_f_1_0_0 + 10);

    auto tg_0_yyyz_f_1_0_0 = pbuffer.data(idx_sg_f_1_0_0 + 11);

    auto tg_0_yyzz_f_1_0_0 = pbuffer.data(idx_sg_f_1_0_0 + 12);

    auto tg_0_yzzz_f_1_0_0 = pbuffer.data(idx_sg_f_1_0_0 + 13);

    auto tg_0_zzzz_f_1_0_0 = pbuffer.data(idx_sg_f_1_0_0 + 14);

    // Set up components of auxiliary buffer : SG

    auto tg_0_xxxx_p_1_0_1 = pbuffer.data(idx_sg_p_1_0_1);

    auto tg_0_xxxy_p_1_0_1 = pbuffer.data(idx_sg_p_1_0_1 + 1);

    auto tg_0_xxxz_p_1_0_1 = pbuffer.data(idx_sg_p_1_0_1 + 2);

    auto tg_0_xxyy_p_1_0_1 = pbuffer.data(idx_sg_p_1_0_1 + 3);

    auto tg_0_xxyz_p_1_0_1 = pbuffer.data(idx_sg_p_1_0_1 + 4);

    auto tg_0_xxzz_p_1_0_1 = pbuffer.data(idx_sg_p_1_0_1 + 5);

    auto tg_0_xyyy_p_1_0_1 = pbuffer.data(idx_sg_p_1_0_1 + 6);

    auto tg_0_xyyz_p_1_0_1 = pbuffer.data(idx_sg_p_1_0_1 + 7);

    auto tg_0_xyzz_p_1_0_1 = pbuffer.data(idx_sg_p_1_0_1 + 8);

    auto tg_0_xzzz_p_1_0_1 = pbuffer.data(idx_sg_p_1_0_1 + 9);

    auto tg_0_yyyy_p_1_0_1 = pbuffer.data(idx_sg_p_1_0_1 + 10);

    auto tg_0_yyyz_p_1_0_1 = pbuffer.data(idx_sg_p_1_0_1 + 11);

    auto tg_0_yyzz_p_1_0_1 = pbuffer.data(idx_sg_p_1_0_1 + 12);

    auto tg_0_yzzz_p_1_0_1 = pbuffer.data(idx_sg_p_1_0_1 + 13);

    auto tg_0_zzzz_p_1_0_1 = pbuffer.data(idx_sg_p_1_0_1 + 14);

    // Set up components of auxiliary buffer : SF

    auto tg_0_xxx_s_1_1_1 = pbuffer.data(idx_sf_s_1_1_1);

    auto tg_0_xxy_s_1_1_1 = pbuffer.data(idx_sf_s_1_1_1 + 1);

    auto tg_0_xxz_s_1_1_1 = pbuffer.data(idx_sf_s_1_1_1 + 2);

    auto tg_0_xyy_s_1_1_1 = pbuffer.data(idx_sf_s_1_1_1 + 3);

    auto tg_0_xyz_s_1_1_1 = pbuffer.data(idx_sf_s_1_1_1 + 4);

    auto tg_0_xzz_s_1_1_1 = pbuffer.data(idx_sf_s_1_1_1 + 5);

    auto tg_0_yyy_s_1_1_1 = pbuffer.data(idx_sf_s_1_1_1 + 6);

    auto tg_0_yyz_s_1_1_1 = pbuffer.data(idx_sf_s_1_1_1 + 7);

    auto tg_0_yzz_s_1_1_1 = pbuffer.data(idx_sf_s_1_1_1 + 8);

    auto tg_0_zzz_s_1_1_1 = pbuffer.data(idx_sf_s_1_1_1 + 9);

    // Set up components of auxiliary buffer : SG

    auto tg_0_xxxx_s_1_1_1 = pbuffer.data(idx_sg_s_1_1_1);

    auto tg_0_xxxy_s_1_1_1 = pbuffer.data(idx_sg_s_1_1_1 + 1);

    auto tg_0_xxxz_s_1_1_1 = pbuffer.data(idx_sg_s_1_1_1 + 2);

    auto tg_0_xxyy_s_1_1_1 = pbuffer.data(idx_sg_s_1_1_1 + 3);

    auto tg_0_xxyz_s_1_1_1 = pbuffer.data(idx_sg_s_1_1_1 + 4);

    auto tg_0_xxzz_s_1_1_1 = pbuffer.data(idx_sg_s_1_1_1 + 5);

    auto tg_0_xyyy_s_1_1_1 = pbuffer.data(idx_sg_s_1_1_1 + 6);

    auto tg_0_xyyz_s_1_1_1 = pbuffer.data(idx_sg_s_1_1_1 + 7);

    auto tg_0_xyzz_s_1_1_1 = pbuffer.data(idx_sg_s_1_1_1 + 8);

    auto tg_0_xzzz_s_1_1_1 = pbuffer.data(idx_sg_s_1_1_1 + 9);

    auto tg_0_yyyy_s_1_1_1 = pbuffer.data(idx_sg_s_1_1_1 + 10);

    auto tg_0_yyyz_s_1_1_1 = pbuffer.data(idx_sg_s_1_1_1 + 11);

    auto tg_0_yyzz_s_1_1_1 = pbuffer.data(idx_sg_s_1_1_1 + 12);

    auto tg_0_yzzz_s_1_1_1 = pbuffer.data(idx_sg_s_1_1_1 + 13);

    auto tg_0_zzzz_s_1_1_1 = pbuffer.data(idx_sg_s_1_1_1 + 14);

    // Set up components of targeted buffer : PG

    auto tg_x_xxxx_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0);

    auto tg_x_xxxy_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 1);

    auto tg_x_xxxz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 2);

    auto tg_x_xxyy_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 3);

    auto tg_x_xxyz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 4);

    auto tg_x_xxzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 5);

    auto tg_x_xyyy_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 6);

    auto tg_x_xyyz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 7);

    auto tg_x_xyzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 8);

    auto tg_x_xzzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 9);

    auto tg_x_yyyy_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 10);

    auto tg_x_yyyz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 11);

    auto tg_x_yyzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 12);

    auto tg_x_yzzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 13);

    auto tg_x_zzzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 14);

    auto tg_y_xxxx_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 15);

    auto tg_y_xxxy_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 16);

    auto tg_y_xxxz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 17);

    auto tg_y_xxyy_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 18);

    auto tg_y_xxyz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 19);

    auto tg_y_xxzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 20);

    auto tg_y_xyyy_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 21);

    auto tg_y_xyyz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 22);

    auto tg_y_xyzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 23);

    auto tg_y_xzzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 24);

    auto tg_y_yyyy_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 25);

    auto tg_y_yyyz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 26);

    auto tg_y_yyzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 27);

    auto tg_y_yzzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 28);

    auto tg_y_zzzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 29);

    auto tg_z_xxxx_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 30);

    auto tg_z_xxxy_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 31);

    auto tg_z_xxxz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 32);

    auto tg_z_xxyy_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 33);

    auto tg_z_xxyz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 34);

    auto tg_z_xxzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 35);

    auto tg_z_xyyy_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 36);

    auto tg_z_xyyz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 37);

    auto tg_z_xyzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 38);

    auto tg_z_xzzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 39);

    auto tg_z_yyyy_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 40);

    auto tg_z_yyyz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 41);

    auto tg_z_yyzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 42);

    auto tg_z_yzzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 43);

    auto tg_z_zzzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 44);

    #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_0_xxxx_d_0_0_1, tg_0_xxxx_f_0_0_0, tg_0_xxxx_f_1_0_0, tg_0_xxxx_p_1_0_1, tg_0_xxxx_s_1_1_1, tg_0_xxxy_d_0_0_1, tg_0_xxxy_f_0_0_0, tg_0_xxxy_f_1_0_0, tg_0_xxxy_p_1_0_1, tg_0_xxxy_s_1_1_1, tg_0_xxxz_d_0_0_1, tg_0_xxxz_f_0_0_0, tg_0_xxxz_f_1_0_0, tg_0_xxxz_p_1_0_1, tg_0_xxxz_s_1_1_1, tg_0_xxyy_d_0_0_1, tg_0_xxyy_f_0_0_0, tg_0_xxyy_f_1_0_0, tg_0_xxyy_p_1_0_1, tg_0_xxyy_s_1_1_1, tg_0_xxyz_d_0_0_1, tg_0_xxyz_f_0_0_0, tg_0_xxyz_f_1_0_0, tg_0_xxyz_p_1_0_1, tg_0_xxyz_s_1_1_1, tg_0_xxzz_d_0_0_1, tg_0_xxzz_f_0_0_0, tg_0_xxzz_f_1_0_0, tg_0_xxzz_p_1_0_1, tg_0_xxzz_s_1_1_1, tg_0_xyyy_d_0_0_1, tg_0_xyyy_f_0_0_0, tg_0_xyyy_f_1_0_0, tg_0_xyyy_p_1_0_1, tg_0_xyyy_s_1_1_1, tg_0_xyyz_d_0_0_1, tg_0_xyyz_f_0_0_0, tg_0_xyyz_f_1_0_0, tg_0_xyyz_p_1_0_1, tg_0_xyyz_s_1_1_1, tg_0_xyzz_d_0_0_1, tg_0_xyzz_f_0_0_0, tg_0_xyzz_f_1_0_0, tg_0_xyzz_p_1_0_1, tg_0_xyzz_s_1_1_1, tg_0_xzzz_d_0_0_1, tg_0_xzzz_f_0_0_0, tg_0_xzzz_f_1_0_0, tg_0_xzzz_p_1_0_1, tg_0_xzzz_s_1_1_1, tg_0_yyyy_d_0_0_1, tg_0_yyyy_f_0_0_0, tg_0_yyyy_f_1_0_0, tg_0_yyyy_p_1_0_1, tg_0_yyyy_s_1_1_1, tg_0_yyyz_d_0_0_1, tg_0_yyyz_f_0_0_0, tg_0_yyyz_f_1_0_0, tg_0_yyyz_p_1_0_1, tg_0_yyyz_s_1_1_1, tg_0_yyzz_d_0_0_1, tg_0_yyzz_f_0_0_0, tg_0_yyzz_f_1_0_0, tg_0_yyzz_p_1_0_1, tg_0_yyzz_s_1_1_1, tg_0_yzzz_d_0_0_1, tg_0_yzzz_f_0_0_0, tg_0_yzzz_f_1_0_0, tg_0_yzzz_p_1_0_1, tg_0_yzzz_s_1_1_1, tg_0_zzzz_d_0_0_1, tg_0_zzzz_f_0_0_0, tg_0_zzzz_f_1_0_0, tg_0_zzzz_p_1_0_1, tg_0_zzzz_s_1_1_1, tg_x_xxxx_f_0_0_0, tg_x_xxxy_f_0_0_0, tg_x_xxxz_f_0_0_0, tg_x_xxyy_f_0_0_0, tg_x_xxyz_f_0_0_0, tg_x_xxzz_f_0_0_0, tg_x_xyyy_f_0_0_0, tg_x_xyyz_f_0_0_0, tg_x_xyzz_f_0_0_0, tg_x_xzzz_f_0_0_0, tg_x_yyyy_f_0_0_0, tg_x_yyyz_f_0_0_0, tg_x_yyzz_f_0_0_0, tg_x_yzzz_f_0_0_0, tg_x_zzzz_f_0_0_0, tg_y_xxxx_f_0_0_0, tg_y_xxxy_f_0_0_0, tg_y_xxxz_f_0_0_0, tg_y_xxyy_f_0_0_0, tg_y_xxyz_f_0_0_0, tg_y_xxzz_f_0_0_0, tg_y_xyyy_f_0_0_0, tg_y_xyyz_f_0_0_0, tg_y_xyzz_f_0_0_0, tg_y_xzzz_f_0_0_0, tg_y_yyyy_f_0_0_0, tg_y_yyyz_f_0_0_0, tg_y_yyzz_f_0_0_0, tg_y_yzzz_f_0_0_0, tg_y_zzzz_f_0_0_0, tg_z_xxxx_f_0_0_0, tg_z_xxxy_f_0_0_0, tg_z_xxxz_f_0_0_0, tg_z_xxyy_f_0_0_0, tg_z_xxyz_f_0_0_0, tg_z_xxzz_f_0_0_0, tg_z_xyyy_f_0_0_0, tg_z_xyyz_f_0_0_0, tg_z_xyzz_f_0_0_0, tg_z_xzzz_f_0_0_0, tg_z_yyyy_f_0_0_0, tg_z_yyyz_f_0_0_0, tg_z_yyzz_f_0_0_0, tg_z_yzzz_f_0_0_0, tg_z_zzzz_f_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double faz_0 = -(b_exps[i] + c_exp) * fzi_0;

        const double fbzi_0 = b_exps[i] * fzi_0;

        const double f2abz_0 = 2.0 * a_exp * b_exps[i] * fzi_0;

            const double fbi_0 = 1.0 / b_exps[i];

        tg_x_xxxx_f_0_0_0[i] = 14.0 * tg_0_xxxx_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_0_xxxx_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_0_xxxx_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 14.0 * tg_0_xxxx_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_xxxx_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxxx_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxx_f_0_0_0[i] * a_x * faz_0;

        tg_x_xxxy_f_0_0_0[i] = 21.0 / 2.0 * tg_0_xxxy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_0_xxxy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_0_xxxy_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 21.0 / 2.0 * tg_0_xxxy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_xxxy_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxxy_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxy_f_0_0_0[i] * a_x * faz_0;

        tg_x_xxxz_f_0_0_0[i] = 21.0 / 2.0 * tg_0_xxxz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_0_xxxz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_0_xxxz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 21.0 / 2.0 * tg_0_xxxz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_xxxz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxxz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxz_f_0_0_0[i] * a_x * faz_0;

        tg_x_xxyy_f_0_0_0[i] = 7.0 * tg_0_xxyy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_0_xxyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_0_xxyy_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_0_xxyy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_xxyy_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxyy_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyy_f_0_0_0[i] * a_x * faz_0;

        tg_x_xxyz_f_0_0_0[i] = 7.0 * tg_0_xxyz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_0_xxyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_0_xxyz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_0_xxyz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_xxyz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxyz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyz_f_0_0_0[i] * a_x * faz_0;

        tg_x_xxzz_f_0_0_0[i] = 7.0 * tg_0_xxzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_0_xxzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_0_xxzz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_0_xxzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_xxzz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxzz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxzz_f_0_0_0[i] * a_x * faz_0;

        tg_x_xyyy_f_0_0_0[i] = 7.0 / 2.0 * tg_0_xyyy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_0_xyyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_0_xyyy_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 / 2.0 * tg_0_xyyy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_xyyy_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xyyy_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyy_f_0_0_0[i] * a_x * faz_0;

        tg_x_xyyz_f_0_0_0[i] = 7.0 / 2.0 * tg_0_xyyz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_0_xyyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_0_xyyz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 / 2.0 * tg_0_xyyz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_xyyz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xyyz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyz_f_0_0_0[i] * a_x * faz_0;

        tg_x_xyzz_f_0_0_0[i] = 7.0 / 2.0 * tg_0_xyzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_0_xyzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_0_xyzz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 / 2.0 * tg_0_xyzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_xyzz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xyzz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xyzz_f_0_0_0[i] * a_x * faz_0;

        tg_x_xzzz_f_0_0_0[i] = 7.0 / 2.0 * tg_0_xzzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_0_xzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_0_xzzz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 / 2.0 * tg_0_xzzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_xzzz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xzzz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xzzz_f_0_0_0[i] * a_x * faz_0;

        tg_x_yyyy_f_0_0_0[i] = 7.0 * tg_0_yyyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_0_yyyy_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_0_yyyy_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_yyyy_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyy_f_0_0_0[i] * a_x * faz_0;

        tg_x_yyyz_f_0_0_0[i] = 7.0 * tg_0_yyyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_0_yyyz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_0_yyyz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_yyyz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyz_f_0_0_0[i] * a_x * faz_0;

        tg_x_yyzz_f_0_0_0[i] = 7.0 * tg_0_yyzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_0_yyzz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_0_yyzz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_yyzz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_yyzz_f_0_0_0[i] * a_x * faz_0;

        tg_x_yzzz_f_0_0_0[i] = 7.0 * tg_0_yzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_0_yzzz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_0_yzzz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_yzzz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_yzzz_f_0_0_0[i] * a_x * faz_0;

        tg_x_zzzz_f_0_0_0[i] = 7.0 * tg_0_zzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_0_zzzz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_0_zzzz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_zzzz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_zzzz_f_0_0_0[i] * a_x * faz_0;

        tg_y_xxxx_f_0_0_0[i] = 7.0 * tg_0_xxxx_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_0_xxxx_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_0_xxxx_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxxx_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxx_f_0_0_0[i] * a_y * faz_0;

        tg_y_xxxy_f_0_0_0[i] = 7.0 / 2.0 * tg_0_xxxy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_0_xxxy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_0_xxxy_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 / 2.0 * tg_0_xxxy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_xxxy_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxxy_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxy_f_0_0_0[i] * a_y * faz_0;

        tg_y_xxxz_f_0_0_0[i] = 7.0 * tg_0_xxxz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_0_xxxz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_0_xxxz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxxz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxz_f_0_0_0[i] * a_y * faz_0;

        tg_y_xxyy_f_0_0_0[i] = 7.0 * tg_0_xxyy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_0_xxyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_0_xxyy_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_0_xxyy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_xxyy_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxyy_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyy_f_0_0_0[i] * a_y * faz_0;

        tg_y_xxyz_f_0_0_0[i] = 7.0 / 2.0 * tg_0_xxyz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_0_xxyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_0_xxyz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 / 2.0 * tg_0_xxyz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_xxyz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxyz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyz_f_0_0_0[i] * a_y * faz_0;

        tg_y_xxzz_f_0_0_0[i] = 7.0 * tg_0_xxzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_0_xxzz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_0_xxzz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxzz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxzz_f_0_0_0[i] * a_y * faz_0;

        tg_y_xyyy_f_0_0_0[i] = 21.0 / 2.0 * tg_0_xyyy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_0_xyyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_0_xyyy_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 21.0 / 2.0 * tg_0_xyyy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_xyyy_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xyyy_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyy_f_0_0_0[i] * a_y * faz_0;

        tg_y_xyyz_f_0_0_0[i] = 7.0 * tg_0_xyyz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_0_xyyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_0_xyyz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_0_xyyz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_xyyz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xyyz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyz_f_0_0_0[i] * a_y * faz_0;

        tg_y_xyzz_f_0_0_0[i] = 7.0 / 2.0 * tg_0_xyzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_0_xyzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_0_xyzz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 / 2.0 * tg_0_xyzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_xyzz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xyzz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xyzz_f_0_0_0[i] * a_y * faz_0;

        tg_y_xzzz_f_0_0_0[i] = 7.0 * tg_0_xzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_0_xzzz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_0_xzzz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xzzz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xzzz_f_0_0_0[i] * a_y * faz_0;

        tg_y_yyyy_f_0_0_0[i] = 14.0 * tg_0_yyyy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_0_yyyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_0_yyyy_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 14.0 * tg_0_yyyy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_yyyy_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_yyyy_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyy_f_0_0_0[i] * a_y * faz_0;

        tg_y_yyyz_f_0_0_0[i] = 21.0 / 2.0 * tg_0_yyyz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_0_yyyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_0_yyyz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 21.0 / 2.0 * tg_0_yyyz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_yyyz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_yyyz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyz_f_0_0_0[i] * a_y * faz_0;

        tg_y_yyzz_f_0_0_0[i] = 7.0 * tg_0_yyzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_0_yyzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_0_yyzz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_0_yyzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_yyzz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_yyzz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_yyzz_f_0_0_0[i] * a_y * faz_0;

        tg_y_yzzz_f_0_0_0[i] = 7.0 / 2.0 * tg_0_yzzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_0_yzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_0_yzzz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 / 2.0 * tg_0_yzzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_yzzz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_yzzz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_yzzz_f_0_0_0[i] * a_y * faz_0;

        tg_y_zzzz_f_0_0_0[i] = 7.0 * tg_0_zzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_0_zzzz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_0_zzzz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_zzzz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_zzzz_f_0_0_0[i] * a_y * faz_0;

        tg_z_xxxx_f_0_0_0[i] = 7.0 * tg_0_xxxx_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_0_xxxx_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_0_xxxx_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxxx_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxx_f_0_0_0[i] * a_z * faz_0;

        tg_z_xxxy_f_0_0_0[i] = 7.0 * tg_0_xxxy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_0_xxxy_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_0_xxxy_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxxy_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxy_f_0_0_0[i] * a_z * faz_0;

        tg_z_xxxz_f_0_0_0[i] = 7.0 / 2.0 * tg_0_xxxz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_0_xxxz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_0_xxxz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 / 2.0 * tg_0_xxxz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_xxxz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxxz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxz_f_0_0_0[i] * a_z * faz_0;

        tg_z_xxyy_f_0_0_0[i] = 7.0 * tg_0_xxyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_0_xxyy_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_0_xxyy_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxyy_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyy_f_0_0_0[i] * a_z * faz_0;

        tg_z_xxyz_f_0_0_0[i] = 7.0 / 2.0 * tg_0_xxyz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_0_xxyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_0_xxyz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 / 2.0 * tg_0_xxyz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_xxyz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxyz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyz_f_0_0_0[i] * a_z * faz_0;

        tg_z_xxzz_f_0_0_0[i] = 7.0 * tg_0_xxzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_0_xxzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_0_xxzz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_0_xxzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_xxzz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxzz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxzz_f_0_0_0[i] * a_z * faz_0;

        tg_z_xyyy_f_0_0_0[i] = 7.0 * tg_0_xyyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_0_xyyy_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_0_xyyy_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xyyy_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyy_f_0_0_0[i] * a_z * faz_0;

        tg_z_xyyz_f_0_0_0[i] = 7.0 / 2.0 * tg_0_xyyz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_0_xyyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_0_xyyz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 / 2.0 * tg_0_xyyz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_xyyz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xyyz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyz_f_0_0_0[i] * a_z * faz_0;

        tg_z_xyzz_f_0_0_0[i] = 7.0 * tg_0_xyzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_0_xyzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_0_xyzz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_0_xyzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_xyzz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xyzz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xyzz_f_0_0_0[i] * a_z * faz_0;

        tg_z_xzzz_f_0_0_0[i] = 21.0 / 2.0 * tg_0_xzzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_0_xzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_0_xzzz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 21.0 / 2.0 * tg_0_xzzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_xzzz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xzzz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xzzz_f_0_0_0[i] * a_z * faz_0;

        tg_z_yyyy_f_0_0_0[i] = 7.0 * tg_0_yyyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_0_yyyy_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_0_yyyy_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_yyyy_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyy_f_0_0_0[i] * a_z * faz_0;

        tg_z_yyyz_f_0_0_0[i] = 7.0 / 2.0 * tg_0_yyyz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_0_yyyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_0_yyyz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 / 2.0 * tg_0_yyyz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_yyyz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_yyyz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyz_f_0_0_0[i] * a_z * faz_0;

        tg_z_yyzz_f_0_0_0[i] = 7.0 * tg_0_yyzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_0_yyzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_0_yyzz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_0_yyzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_yyzz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_yyzz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_yyzz_f_0_0_0[i] * a_z * faz_0;

        tg_z_yzzz_f_0_0_0[i] = 21.0 / 2.0 * tg_0_yzzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_0_yzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_0_yzzz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 21.0 / 2.0 * tg_0_yzzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_yzzz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_yzzz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_yzzz_f_0_0_0[i] * a_z * faz_0;

        tg_z_zzzz_f_0_0_0[i] = 14.0 * tg_0_zzzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_0_zzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_0_zzzz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 14.0 * tg_0_zzzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_zzzz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_zzzz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_zzzz_f_0_0_0[i] * a_z * faz_0;
    }

    if (p > 0)
    {
        const double fp_0 = (double)p;

        // Set up components of auxiliary buffer : SG

        auto tg_0_xxxx_f_0_0_1 = pbuffer.data(idx_sg_f_0_0_1);

        auto tg_0_xxxy_f_0_0_1 = pbuffer.data(idx_sg_f_0_0_1 + 1);

        auto tg_0_xxxz_f_0_0_1 = pbuffer.data(idx_sg_f_0_0_1 + 2);

        auto tg_0_xxyy_f_0_0_1 = pbuffer.data(idx_sg_f_0_0_1 + 3);

        auto tg_0_xxyz_f_0_0_1 = pbuffer.data(idx_sg_f_0_0_1 + 4);

        auto tg_0_xxzz_f_0_0_1 = pbuffer.data(idx_sg_f_0_0_1 + 5);

        auto tg_0_xyyy_f_0_0_1 = pbuffer.data(idx_sg_f_0_0_1 + 6);

        auto tg_0_xyyz_f_0_0_1 = pbuffer.data(idx_sg_f_0_0_1 + 7);

        auto tg_0_xyzz_f_0_0_1 = pbuffer.data(idx_sg_f_0_0_1 + 8);

        auto tg_0_xzzz_f_0_0_1 = pbuffer.data(idx_sg_f_0_0_1 + 9);

        auto tg_0_yyyy_f_0_0_1 = pbuffer.data(idx_sg_f_0_0_1 + 10);

        auto tg_0_yyyz_f_0_0_1 = pbuffer.data(idx_sg_f_0_0_1 + 11);

        auto tg_0_yyzz_f_0_0_1 = pbuffer.data(idx_sg_f_0_0_1 + 12);

        auto tg_0_yzzz_f_0_0_1 = pbuffer.data(idx_sg_f_0_0_1 + 13);

        auto tg_0_zzzz_f_0_0_1 = pbuffer.data(idx_sg_f_0_0_1 + 14);

        #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_0_xxxx_f_0_0_1, tg_0_xxxy_f_0_0_1, tg_0_xxxz_f_0_0_1, tg_0_xxyy_f_0_0_1, tg_0_xxyz_f_0_0_1, tg_0_xxzz_f_0_0_1, tg_0_xyyy_f_0_0_1, tg_0_xyyz_f_0_0_1, tg_0_xyzz_f_0_0_1, tg_0_xzzz_f_0_0_1, tg_0_yyyy_f_0_0_1, tg_0_yyyz_f_0_0_1, tg_0_yyzz_f_0_0_1, tg_0_yzzz_f_0_0_1, tg_0_zzzz_f_0_0_1, tg_x_xxxx_f_0_0_0, tg_x_xxxy_f_0_0_0, tg_x_xxxz_f_0_0_0, tg_x_xxyy_f_0_0_0, tg_x_xxyz_f_0_0_0, tg_x_xxzz_f_0_0_0, tg_x_xyyy_f_0_0_0, tg_x_xyyz_f_0_0_0, tg_x_xyzz_f_0_0_0, tg_x_xzzz_f_0_0_0, tg_x_yyyy_f_0_0_0, tg_x_yyyz_f_0_0_0, tg_x_yyzz_f_0_0_0, tg_x_yzzz_f_0_0_0, tg_x_zzzz_f_0_0_0, tg_y_xxxx_f_0_0_0, tg_y_xxxy_f_0_0_0, tg_y_xxxz_f_0_0_0, tg_y_xxyy_f_0_0_0, tg_y_xxyz_f_0_0_0, tg_y_xxzz_f_0_0_0, tg_y_xyyy_f_0_0_0, tg_y_xyyz_f_0_0_0, tg_y_xyzz_f_0_0_0, tg_y_xzzz_f_0_0_0, tg_y_yyyy_f_0_0_0, tg_y_yyyz_f_0_0_0, tg_y_yyzz_f_0_0_0, tg_y_yzzz_f_0_0_0, tg_y_zzzz_f_0_0_0, tg_z_xxxx_f_0_0_0, tg_z_xxxy_f_0_0_0, tg_z_xxxz_f_0_0_0, tg_z_xxyy_f_0_0_0, tg_z_xxyz_f_0_0_0, tg_z_xxzz_f_0_0_0, tg_z_xyyy_f_0_0_0, tg_z_xyyz_f_0_0_0, tg_z_xyzz_f_0_0_0, tg_z_xzzz_f_0_0_0, tg_z_yyyy_f_0_0_0, tg_z_yyyz_f_0_0_0, tg_z_yyzz_f_0_0_0, tg_z_yzzz_f_0_0_0, tg_z_zzzz_f_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fai_0 = 1.0 / a_exp;

            tg_x_xxxx_f_0_0_0[i] = tg_0_xxxx_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxxy_f_0_0_0[i] = tg_0_xxxy_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxxz_f_0_0_0[i] = tg_0_xxxz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxyy_f_0_0_0[i] = tg_0_xxyy_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxyz_f_0_0_0[i] = tg_0_xxyz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxzz_f_0_0_0[i] = tg_0_xxzz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xyyy_f_0_0_0[i] = tg_0_xyyy_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xyyz_f_0_0_0[i] = tg_0_xyyz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xyzz_f_0_0_0[i] = tg_0_xyzz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xzzz_f_0_0_0[i] = tg_0_xzzz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_yyyy_f_0_0_0[i] = tg_0_yyyy_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_yyyz_f_0_0_0[i] = tg_0_yyyz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_yyzz_f_0_0_0[i] = tg_0_yyzz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_yzzz_f_0_0_0[i] = tg_0_yzzz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_zzzz_f_0_0_0[i] = tg_0_zzzz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_y_xxxx_f_0_0_0[i] = tg_0_xxxx_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxxy_f_0_0_0[i] = tg_0_xxxy_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxxz_f_0_0_0[i] = tg_0_xxxz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxyy_f_0_0_0[i] = tg_0_xxyy_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxyz_f_0_0_0[i] = tg_0_xxyz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxzz_f_0_0_0[i] = tg_0_xxzz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xyyy_f_0_0_0[i] = tg_0_xyyy_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xyyz_f_0_0_0[i] = tg_0_xyyz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xyzz_f_0_0_0[i] = tg_0_xyzz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xzzz_f_0_0_0[i] = tg_0_xzzz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_yyyy_f_0_0_0[i] = tg_0_yyyy_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_yyyz_f_0_0_0[i] = tg_0_yyyz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_yyzz_f_0_0_0[i] = tg_0_yyzz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_yzzz_f_0_0_0[i] = tg_0_yzzz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_zzzz_f_0_0_0[i] = tg_0_zzzz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_z_xxxx_f_0_0_0[i] = tg_0_xxxx_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxxy_f_0_0_0[i] = tg_0_xxxy_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxxz_f_0_0_0[i] = tg_0_xxxz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxyy_f_0_0_0[i] = tg_0_xxyy_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxyz_f_0_0_0[i] = tg_0_xxyz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxzz_f_0_0_0[i] = tg_0_xxzz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xyyy_f_0_0_0[i] = tg_0_xyyy_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xyyz_f_0_0_0[i] = tg_0_xyyz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xyzz_f_0_0_0[i] = tg_0_xyzz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xzzz_f_0_0_0[i] = tg_0_xzzz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_yyyy_f_0_0_0[i] = tg_0_yyyy_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_yyyz_f_0_0_0[i] = tg_0_yyyz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_yyzz_f_0_0_0[i] = tg_0_yyzz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_yzzz_f_0_0_0[i] = tg_0_yzzz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_zzzz_f_0_0_0[i] = tg_0_zzzz_f_0_0_1[i] * fai_0 * a_z * fp_0;
        }
    }
}

} // t2pecp namespace

