#include "ProjectedCorePotentialPrimRecPGForD.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_pg_d(CSimdArray<double>& pbuffer, 
                                        const size_t idx_pg_d_0_0_0,
                                        const size_t idx_sg_d_0_0_0,
                                        const size_t idx_sf_p_0_0_1,
                                        const size_t idx_sg_p_0_0_1,
                                        const size_t idx_sg_d_1_0_0,
                                        const size_t idx_sg_s_1_0_1,
                                        const int p,
                                        const size_t idx_sg_d_0_0_1,
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

    auto tg_0_xxxx_d_0_0_0 = pbuffer.data(idx_sg_d_0_0_0);

    auto tg_0_xxxy_d_0_0_0 = pbuffer.data(idx_sg_d_0_0_0 + 1);

    auto tg_0_xxxz_d_0_0_0 = pbuffer.data(idx_sg_d_0_0_0 + 2);

    auto tg_0_xxyy_d_0_0_0 = pbuffer.data(idx_sg_d_0_0_0 + 3);

    auto tg_0_xxyz_d_0_0_0 = pbuffer.data(idx_sg_d_0_0_0 + 4);

    auto tg_0_xxzz_d_0_0_0 = pbuffer.data(idx_sg_d_0_0_0 + 5);

    auto tg_0_xyyy_d_0_0_0 = pbuffer.data(idx_sg_d_0_0_0 + 6);

    auto tg_0_xyyz_d_0_0_0 = pbuffer.data(idx_sg_d_0_0_0 + 7);

    auto tg_0_xyzz_d_0_0_0 = pbuffer.data(idx_sg_d_0_0_0 + 8);

    auto tg_0_xzzz_d_0_0_0 = pbuffer.data(idx_sg_d_0_0_0 + 9);

    auto tg_0_yyyy_d_0_0_0 = pbuffer.data(idx_sg_d_0_0_0 + 10);

    auto tg_0_yyyz_d_0_0_0 = pbuffer.data(idx_sg_d_0_0_0 + 11);

    auto tg_0_yyzz_d_0_0_0 = pbuffer.data(idx_sg_d_0_0_0 + 12);

    auto tg_0_yzzz_d_0_0_0 = pbuffer.data(idx_sg_d_0_0_0 + 13);

    auto tg_0_zzzz_d_0_0_0 = pbuffer.data(idx_sg_d_0_0_0 + 14);

    // Set up components of auxiliary buffer : SF

    auto tg_0_xxx_p_0_0_1 = pbuffer.data(idx_sf_p_0_0_1);

    auto tg_0_xxy_p_0_0_1 = pbuffer.data(idx_sf_p_0_0_1 + 1);

    auto tg_0_xxz_p_0_0_1 = pbuffer.data(idx_sf_p_0_0_1 + 2);

    auto tg_0_xyy_p_0_0_1 = pbuffer.data(idx_sf_p_0_0_1 + 3);

    auto tg_0_xyz_p_0_0_1 = pbuffer.data(idx_sf_p_0_0_1 + 4);

    auto tg_0_xzz_p_0_0_1 = pbuffer.data(idx_sf_p_0_0_1 + 5);

    auto tg_0_yyy_p_0_0_1 = pbuffer.data(idx_sf_p_0_0_1 + 6);

    auto tg_0_yyz_p_0_0_1 = pbuffer.data(idx_sf_p_0_0_1 + 7);

    auto tg_0_yzz_p_0_0_1 = pbuffer.data(idx_sf_p_0_0_1 + 8);

    auto tg_0_zzz_p_0_0_1 = pbuffer.data(idx_sf_p_0_0_1 + 9);

    // Set up components of auxiliary buffer : SG

    auto tg_0_xxxx_p_0_0_1 = pbuffer.data(idx_sg_p_0_0_1);

    auto tg_0_xxxy_p_0_0_1 = pbuffer.data(idx_sg_p_0_0_1 + 1);

    auto tg_0_xxxz_p_0_0_1 = pbuffer.data(idx_sg_p_0_0_1 + 2);

    auto tg_0_xxyy_p_0_0_1 = pbuffer.data(idx_sg_p_0_0_1 + 3);

    auto tg_0_xxyz_p_0_0_1 = pbuffer.data(idx_sg_p_0_0_1 + 4);

    auto tg_0_xxzz_p_0_0_1 = pbuffer.data(idx_sg_p_0_0_1 + 5);

    auto tg_0_xyyy_p_0_0_1 = pbuffer.data(idx_sg_p_0_0_1 + 6);

    auto tg_0_xyyz_p_0_0_1 = pbuffer.data(idx_sg_p_0_0_1 + 7);

    auto tg_0_xyzz_p_0_0_1 = pbuffer.data(idx_sg_p_0_0_1 + 8);

    auto tg_0_xzzz_p_0_0_1 = pbuffer.data(idx_sg_p_0_0_1 + 9);

    auto tg_0_yyyy_p_0_0_1 = pbuffer.data(idx_sg_p_0_0_1 + 10);

    auto tg_0_yyyz_p_0_0_1 = pbuffer.data(idx_sg_p_0_0_1 + 11);

    auto tg_0_yyzz_p_0_0_1 = pbuffer.data(idx_sg_p_0_0_1 + 12);

    auto tg_0_yzzz_p_0_0_1 = pbuffer.data(idx_sg_p_0_0_1 + 13);

    auto tg_0_zzzz_p_0_0_1 = pbuffer.data(idx_sg_p_0_0_1 + 14);

    // Set up components of auxiliary buffer : SG

    auto tg_0_xxxx_d_1_0_0 = pbuffer.data(idx_sg_d_1_0_0);

    auto tg_0_xxxy_d_1_0_0 = pbuffer.data(idx_sg_d_1_0_0 + 1);

    auto tg_0_xxxz_d_1_0_0 = pbuffer.data(idx_sg_d_1_0_0 + 2);

    auto tg_0_xxyy_d_1_0_0 = pbuffer.data(idx_sg_d_1_0_0 + 3);

    auto tg_0_xxyz_d_1_0_0 = pbuffer.data(idx_sg_d_1_0_0 + 4);

    auto tg_0_xxzz_d_1_0_0 = pbuffer.data(idx_sg_d_1_0_0 + 5);

    auto tg_0_xyyy_d_1_0_0 = pbuffer.data(idx_sg_d_1_0_0 + 6);

    auto tg_0_xyyz_d_1_0_0 = pbuffer.data(idx_sg_d_1_0_0 + 7);

    auto tg_0_xyzz_d_1_0_0 = pbuffer.data(idx_sg_d_1_0_0 + 8);

    auto tg_0_xzzz_d_1_0_0 = pbuffer.data(idx_sg_d_1_0_0 + 9);

    auto tg_0_yyyy_d_1_0_0 = pbuffer.data(idx_sg_d_1_0_0 + 10);

    auto tg_0_yyyz_d_1_0_0 = pbuffer.data(idx_sg_d_1_0_0 + 11);

    auto tg_0_yyzz_d_1_0_0 = pbuffer.data(idx_sg_d_1_0_0 + 12);

    auto tg_0_yzzz_d_1_0_0 = pbuffer.data(idx_sg_d_1_0_0 + 13);

    auto tg_0_zzzz_d_1_0_0 = pbuffer.data(idx_sg_d_1_0_0 + 14);

    // Set up components of auxiliary buffer : SG

    auto tg_0_xxxx_s_1_0_1 = pbuffer.data(idx_sg_s_1_0_1);

    auto tg_0_xxxy_s_1_0_1 = pbuffer.data(idx_sg_s_1_0_1 + 1);

    auto tg_0_xxxz_s_1_0_1 = pbuffer.data(idx_sg_s_1_0_1 + 2);

    auto tg_0_xxyy_s_1_0_1 = pbuffer.data(idx_sg_s_1_0_1 + 3);

    auto tg_0_xxyz_s_1_0_1 = pbuffer.data(idx_sg_s_1_0_1 + 4);

    auto tg_0_xxzz_s_1_0_1 = pbuffer.data(idx_sg_s_1_0_1 + 5);

    auto tg_0_xyyy_s_1_0_1 = pbuffer.data(idx_sg_s_1_0_1 + 6);

    auto tg_0_xyyz_s_1_0_1 = pbuffer.data(idx_sg_s_1_0_1 + 7);

    auto tg_0_xyzz_s_1_0_1 = pbuffer.data(idx_sg_s_1_0_1 + 8);

    auto tg_0_xzzz_s_1_0_1 = pbuffer.data(idx_sg_s_1_0_1 + 9);

    auto tg_0_yyyy_s_1_0_1 = pbuffer.data(idx_sg_s_1_0_1 + 10);

    auto tg_0_yyyz_s_1_0_1 = pbuffer.data(idx_sg_s_1_0_1 + 11);

    auto tg_0_yyzz_s_1_0_1 = pbuffer.data(idx_sg_s_1_0_1 + 12);

    auto tg_0_yzzz_s_1_0_1 = pbuffer.data(idx_sg_s_1_0_1 + 13);

    auto tg_0_zzzz_s_1_0_1 = pbuffer.data(idx_sg_s_1_0_1 + 14);

    // Set up components of targeted buffer : PG

    auto tg_x_xxxx_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0);

    auto tg_x_xxxy_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 1);

    auto tg_x_xxxz_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 2);

    auto tg_x_xxyy_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 3);

    auto tg_x_xxyz_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 4);

    auto tg_x_xxzz_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 5);

    auto tg_x_xyyy_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 6);

    auto tg_x_xyyz_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 7);

    auto tg_x_xyzz_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 8);

    auto tg_x_xzzz_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 9);

    auto tg_x_yyyy_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 10);

    auto tg_x_yyyz_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 11);

    auto tg_x_yyzz_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 12);

    auto tg_x_yzzz_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 13);

    auto tg_x_zzzz_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 14);

    auto tg_y_xxxx_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 15);

    auto tg_y_xxxy_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 16);

    auto tg_y_xxxz_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 17);

    auto tg_y_xxyy_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 18);

    auto tg_y_xxyz_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 19);

    auto tg_y_xxzz_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 20);

    auto tg_y_xyyy_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 21);

    auto tg_y_xyyz_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 22);

    auto tg_y_xyzz_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 23);

    auto tg_y_xzzz_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 24);

    auto tg_y_yyyy_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 25);

    auto tg_y_yyyz_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 26);

    auto tg_y_yyzz_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 27);

    auto tg_y_yzzz_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 28);

    auto tg_y_zzzz_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 29);

    auto tg_z_xxxx_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 30);

    auto tg_z_xxxy_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 31);

    auto tg_z_xxxz_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 32);

    auto tg_z_xxyy_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 33);

    auto tg_z_xxyz_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 34);

    auto tg_z_xxzz_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 35);

    auto tg_z_xyyy_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 36);

    auto tg_z_xyyz_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 37);

    auto tg_z_xyzz_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 38);

    auto tg_z_xzzz_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 39);

    auto tg_z_yyyy_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 40);

    auto tg_z_yyyz_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 41);

    auto tg_z_yyzz_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 42);

    auto tg_z_yzzz_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 43);

    auto tg_z_zzzz_d_0_0_0 = pbuffer.data(idx_pg_d_0_0_0 + 44);

    #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_0_xxxx_d_0_0_0, tg_0_xxxx_d_1_0_0, tg_0_xxxx_p_0_0_1, tg_0_xxxx_s_1_0_1, tg_0_xxxy_d_0_0_0, tg_0_xxxy_d_1_0_0, tg_0_xxxy_p_0_0_1, tg_0_xxxy_s_1_0_1, tg_0_xxxz_d_0_0_0, tg_0_xxxz_d_1_0_0, tg_0_xxxz_p_0_0_1, tg_0_xxxz_s_1_0_1, tg_0_xxyy_d_0_0_0, tg_0_xxyy_d_1_0_0, tg_0_xxyy_p_0_0_1, tg_0_xxyy_s_1_0_1, tg_0_xxyz_d_0_0_0, tg_0_xxyz_d_1_0_0, tg_0_xxyz_p_0_0_1, tg_0_xxyz_s_1_0_1, tg_0_xxzz_d_0_0_0, tg_0_xxzz_d_1_0_0, tg_0_xxzz_p_0_0_1, tg_0_xxzz_s_1_0_1, tg_0_xyyy_d_0_0_0, tg_0_xyyy_d_1_0_0, tg_0_xyyy_p_0_0_1, tg_0_xyyy_s_1_0_1, tg_0_xyyz_d_0_0_0, tg_0_xyyz_d_1_0_0, tg_0_xyyz_p_0_0_1, tg_0_xyyz_s_1_0_1, tg_0_xyzz_d_0_0_0, tg_0_xyzz_d_1_0_0, tg_0_xyzz_p_0_0_1, tg_0_xyzz_s_1_0_1, tg_0_xzzz_d_0_0_0, tg_0_xzzz_d_1_0_0, tg_0_xzzz_p_0_0_1, tg_0_xzzz_s_1_0_1, tg_0_yyyy_d_0_0_0, tg_0_yyyy_d_1_0_0, tg_0_yyyy_p_0_0_1, tg_0_yyyy_s_1_0_1, tg_0_yyyz_d_0_0_0, tg_0_yyyz_d_1_0_0, tg_0_yyyz_p_0_0_1, tg_0_yyyz_s_1_0_1, tg_0_yyzz_d_0_0_0, tg_0_yyzz_d_1_0_0, tg_0_yyzz_p_0_0_1, tg_0_yyzz_s_1_0_1, tg_0_yzzz_d_0_0_0, tg_0_yzzz_d_1_0_0, tg_0_yzzz_p_0_0_1, tg_0_yzzz_s_1_0_1, tg_0_zzzz_d_0_0_0, tg_0_zzzz_d_1_0_0, tg_0_zzzz_p_0_0_1, tg_0_zzzz_s_1_0_1, tg_x_xxxx_d_0_0_0, tg_x_xxxy_d_0_0_0, tg_x_xxxz_d_0_0_0, tg_x_xxyy_d_0_0_0, tg_x_xxyz_d_0_0_0, tg_x_xxzz_d_0_0_0, tg_x_xyyy_d_0_0_0, tg_x_xyyz_d_0_0_0, tg_x_xyzz_d_0_0_0, tg_x_xzzz_d_0_0_0, tg_x_yyyy_d_0_0_0, tg_x_yyyz_d_0_0_0, tg_x_yyzz_d_0_0_0, tg_x_yzzz_d_0_0_0, tg_x_zzzz_d_0_0_0, tg_y_xxxx_d_0_0_0, tg_y_xxxy_d_0_0_0, tg_y_xxxz_d_0_0_0, tg_y_xxyy_d_0_0_0, tg_y_xxyz_d_0_0_0, tg_y_xxzz_d_0_0_0, tg_y_xyyy_d_0_0_0, tg_y_xyyz_d_0_0_0, tg_y_xyzz_d_0_0_0, tg_y_xzzz_d_0_0_0, tg_y_yyyy_d_0_0_0, tg_y_yyyz_d_0_0_0, tg_y_yyzz_d_0_0_0, tg_y_yzzz_d_0_0_0, tg_y_zzzz_d_0_0_0, tg_z_xxxx_d_0_0_0, tg_z_xxxy_d_0_0_0, tg_z_xxxz_d_0_0_0, tg_z_xxyy_d_0_0_0, tg_z_xxyz_d_0_0_0, tg_z_xxzz_d_0_0_0, tg_z_xyyy_d_0_0_0, tg_z_xyyz_d_0_0_0, tg_z_xyzz_d_0_0_0, tg_z_xzzz_d_0_0_0, tg_z_yyyy_d_0_0_0, tg_z_yyyz_d_0_0_0, tg_z_yyzz_d_0_0_0, tg_z_yzzz_d_0_0_0, tg_z_zzzz_d_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double faz_0 = -(b_exps[i] + c_exp) * fzi_0;

        const double fbzi_0 = b_exps[i] * fzi_0;

        const double f2abz_0 = 2.0 * a_exp * b_exps[i] * fzi_0;

            const double fbi_0 = 1.0 / b_exps[i];

        tg_x_xxxx_d_0_0_0[i] = -5.0 * tg_0_xxxx_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 10.0 * tg_0_xxxx_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_xxxx_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxxx_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxx_d_0_0_0[i] * a_x * faz_0;

        tg_x_xxxy_d_0_0_0[i] = -5.0 * tg_0_xxxy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 15.0 / 2.0 * tg_0_xxxy_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_xxxy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxxy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxy_d_0_0_0[i] * a_x * faz_0;

        tg_x_xxxz_d_0_0_0[i] = -5.0 * tg_0_xxxz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 15.0 / 2.0 * tg_0_xxxz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_xxxz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxxz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxz_d_0_0_0[i] * a_x * faz_0;

        tg_x_xxyy_d_0_0_0[i] = -5.0 * tg_0_xxyy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_0_xxyy_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_xxyy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxyy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyy_d_0_0_0[i] * a_x * faz_0;

        tg_x_xxyz_d_0_0_0[i] = -5.0 * tg_0_xxyz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_0_xxyz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_xxyz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxyz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyz_d_0_0_0[i] * a_x * faz_0;

        tg_x_xxzz_d_0_0_0[i] = -5.0 * tg_0_xxzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_0_xxzz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_xxzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xxzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xxzz_d_0_0_0[i] * a_x * faz_0;

        tg_x_xyyy_d_0_0_0[i] = -5.0 * tg_0_xyyy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 / 2.0 * tg_0_xyyy_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_xyyy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xyyy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyy_d_0_0_0[i] * a_x * faz_0;

        tg_x_xyyz_d_0_0_0[i] = -5.0 * tg_0_xyyz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 / 2.0 * tg_0_xyyz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_xyyz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xyyz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyz_d_0_0_0[i] * a_x * faz_0;

        tg_x_xyzz_d_0_0_0[i] = -5.0 * tg_0_xyzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 / 2.0 * tg_0_xyzz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_xyzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xyzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xyzz_d_0_0_0[i] * a_x * faz_0;

        tg_x_xzzz_d_0_0_0[i] = -5.0 * tg_0_xzzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 / 2.0 * tg_0_xzzz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_xzzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xzzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xzzz_d_0_0_0[i] * a_x * faz_0;

        tg_x_yyyy_d_0_0_0[i] = -5.0 * tg_0_yyyy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_0_yyyy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_yyyy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyy_d_0_0_0[i] * a_x * faz_0;

        tg_x_yyyz_d_0_0_0[i] = -5.0 * tg_0_yyyz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_0_yyyz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_yyyz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyz_d_0_0_0[i] * a_x * faz_0;

        tg_x_yyzz_d_0_0_0[i] = -5.0 * tg_0_yyzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_0_yyzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_yyzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_yyzz_d_0_0_0[i] * a_x * faz_0;

        tg_x_yzzz_d_0_0_0[i] = -5.0 * tg_0_yzzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_0_yzzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_yzzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_yzzz_d_0_0_0[i] * a_x * faz_0;

        tg_x_zzzz_d_0_0_0[i] = -5.0 * tg_0_zzzz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_0_zzzz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_zzzz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_zzzz_d_0_0_0[i] * a_x * faz_0;

        tg_y_xxxx_d_0_0_0[i] = -5.0 * tg_0_xxxx_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_0_xxxx_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxxx_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxx_d_0_0_0[i] * a_y * faz_0;

        tg_y_xxxy_d_0_0_0[i] = -5.0 * tg_0_xxxy_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 / 2.0 * tg_0_xxxy_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_xxxy_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxxy_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxy_d_0_0_0[i] * a_y * faz_0;

        tg_y_xxxz_d_0_0_0[i] = -5.0 * tg_0_xxxz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_0_xxxz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxxz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxz_d_0_0_0[i] * a_y * faz_0;

        tg_y_xxyy_d_0_0_0[i] = -5.0 * tg_0_xxyy_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_0_xxyy_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_xxyy_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxyy_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyy_d_0_0_0[i] * a_y * faz_0;

        tg_y_xxyz_d_0_0_0[i] = -5.0 * tg_0_xxyz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 / 2.0 * tg_0_xxyz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_xxyz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxyz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyz_d_0_0_0[i] * a_y * faz_0;

        tg_y_xxzz_d_0_0_0[i] = -5.0 * tg_0_xxzz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_0_xxzz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xxzz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xxzz_d_0_0_0[i] * a_y * faz_0;

        tg_y_xyyy_d_0_0_0[i] = -5.0 * tg_0_xyyy_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 15.0 / 2.0 * tg_0_xyyy_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_xyyy_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xyyy_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyy_d_0_0_0[i] * a_y * faz_0;

        tg_y_xyyz_d_0_0_0[i] = -5.0 * tg_0_xyyz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_0_xyyz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_xyyz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xyyz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyz_d_0_0_0[i] * a_y * faz_0;

        tg_y_xyzz_d_0_0_0[i] = -5.0 * tg_0_xyzz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 / 2.0 * tg_0_xyzz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_xyzz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xyzz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xyzz_d_0_0_0[i] * a_y * faz_0;

        tg_y_xzzz_d_0_0_0[i] = -5.0 * tg_0_xzzz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_0_xzzz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xzzz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xzzz_d_0_0_0[i] * a_y * faz_0;

        tg_y_yyyy_d_0_0_0[i] = -5.0 * tg_0_yyyy_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 10.0 * tg_0_yyyy_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_yyyy_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_yyyy_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyy_d_0_0_0[i] * a_y * faz_0;

        tg_y_yyyz_d_0_0_0[i] = -5.0 * tg_0_yyyz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 15.0 / 2.0 * tg_0_yyyz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_yyyz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_yyyz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyz_d_0_0_0[i] * a_y * faz_0;

        tg_y_yyzz_d_0_0_0[i] = -5.0 * tg_0_yyzz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_0_yyzz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_yyzz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_yyzz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_yyzz_d_0_0_0[i] * a_y * faz_0;

        tg_y_yzzz_d_0_0_0[i] = -5.0 * tg_0_yzzz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 / 2.0 * tg_0_yzzz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_yzzz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_yzzz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_yzzz_d_0_0_0[i] * a_y * faz_0;

        tg_y_zzzz_d_0_0_0[i] = -5.0 * tg_0_zzzz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_0_zzzz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_zzzz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_zzzz_d_0_0_0[i] * a_y * faz_0;

        tg_z_xxxx_d_0_0_0[i] = -5.0 * tg_0_xxxx_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_0_xxxx_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxxx_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxx_d_0_0_0[i] * a_z * faz_0;

        tg_z_xxxy_d_0_0_0[i] = -5.0 * tg_0_xxxy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_0_xxxy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxxy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxy_d_0_0_0[i] * a_z * faz_0;

        tg_z_xxxz_d_0_0_0[i] = -5.0 * tg_0_xxxz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 / 2.0 * tg_0_xxxz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_xxxz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxxz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxxz_d_0_0_0[i] * a_z * faz_0;

        tg_z_xxyy_d_0_0_0[i] = -5.0 * tg_0_xxyy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_0_xxyy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxyy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyy_d_0_0_0[i] * a_z * faz_0;

        tg_z_xxyz_d_0_0_0[i] = -5.0 * tg_0_xxyz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 / 2.0 * tg_0_xxyz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_xxyz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxyz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxyz_d_0_0_0[i] * a_z * faz_0;

        tg_z_xxzz_d_0_0_0[i] = -5.0 * tg_0_xxzz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_0_xxzz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_xxzz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xxzz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xxzz_d_0_0_0[i] * a_z * faz_0;

        tg_z_xyyy_d_0_0_0[i] = -5.0 * tg_0_xyyy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_0_xyyy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xyyy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyy_d_0_0_0[i] * a_z * faz_0;

        tg_z_xyyz_d_0_0_0[i] = -5.0 * tg_0_xyyz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 / 2.0 * tg_0_xyyz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_xyyz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xyyz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xyyz_d_0_0_0[i] * a_z * faz_0;

        tg_z_xyzz_d_0_0_0[i] = -5.0 * tg_0_xyzz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_0_xyzz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_xyzz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xyzz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xyzz_d_0_0_0[i] * a_z * faz_0;

        tg_z_xzzz_d_0_0_0[i] = -5.0 * tg_0_xzzz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 15.0 / 2.0 * tg_0_xzzz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_xzzz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xzzz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xzzz_d_0_0_0[i] * a_z * faz_0;

        tg_z_yyyy_d_0_0_0[i] = -5.0 * tg_0_yyyy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_0_yyyy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_yyyy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyy_d_0_0_0[i] * a_z * faz_0;

        tg_z_yyyz_d_0_0_0[i] = -5.0 * tg_0_yyyz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 / 2.0 * tg_0_yyyz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_yyyz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_yyyz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_yyyz_d_0_0_0[i] * a_z * faz_0;

        tg_z_yyzz_d_0_0_0[i] = -5.0 * tg_0_yyzz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_0_yyzz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_yyzz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_yyzz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_yyzz_d_0_0_0[i] * a_z * faz_0;

        tg_z_yzzz_d_0_0_0[i] = -5.0 * tg_0_yzzz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 15.0 / 2.0 * tg_0_yzzz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_yzzz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_yzzz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_yzzz_d_0_0_0[i] * a_z * faz_0;

        tg_z_zzzz_d_0_0_0[i] = -5.0 * tg_0_zzzz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 10.0 * tg_0_zzzz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_zzzz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_zzzz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_zzzz_d_0_0_0[i] * a_z * faz_0;
    }

    if (p > 0)
    {
        const double fp_0 = (double)p;

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

        #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_0_xxxx_d_0_0_1, tg_0_xxxy_d_0_0_1, tg_0_xxxz_d_0_0_1, tg_0_xxyy_d_0_0_1, tg_0_xxyz_d_0_0_1, tg_0_xxzz_d_0_0_1, tg_0_xyyy_d_0_0_1, tg_0_xyyz_d_0_0_1, tg_0_xyzz_d_0_0_1, tg_0_xzzz_d_0_0_1, tg_0_yyyy_d_0_0_1, tg_0_yyyz_d_0_0_1, tg_0_yyzz_d_0_0_1, tg_0_yzzz_d_0_0_1, tg_0_zzzz_d_0_0_1, tg_x_xxxx_d_0_0_0, tg_x_xxxy_d_0_0_0, tg_x_xxxz_d_0_0_0, tg_x_xxyy_d_0_0_0, tg_x_xxyz_d_0_0_0, tg_x_xxzz_d_0_0_0, tg_x_xyyy_d_0_0_0, tg_x_xyyz_d_0_0_0, tg_x_xyzz_d_0_0_0, tg_x_xzzz_d_0_0_0, tg_x_yyyy_d_0_0_0, tg_x_yyyz_d_0_0_0, tg_x_yyzz_d_0_0_0, tg_x_yzzz_d_0_0_0, tg_x_zzzz_d_0_0_0, tg_y_xxxx_d_0_0_0, tg_y_xxxy_d_0_0_0, tg_y_xxxz_d_0_0_0, tg_y_xxyy_d_0_0_0, tg_y_xxyz_d_0_0_0, tg_y_xxzz_d_0_0_0, tg_y_xyyy_d_0_0_0, tg_y_xyyz_d_0_0_0, tg_y_xyzz_d_0_0_0, tg_y_xzzz_d_0_0_0, tg_y_yyyy_d_0_0_0, tg_y_yyyz_d_0_0_0, tg_y_yyzz_d_0_0_0, tg_y_yzzz_d_0_0_0, tg_y_zzzz_d_0_0_0, tg_z_xxxx_d_0_0_0, tg_z_xxxy_d_0_0_0, tg_z_xxxz_d_0_0_0, tg_z_xxyy_d_0_0_0, tg_z_xxyz_d_0_0_0, tg_z_xxzz_d_0_0_0, tg_z_xyyy_d_0_0_0, tg_z_xyyz_d_0_0_0, tg_z_xyzz_d_0_0_0, tg_z_xzzz_d_0_0_0, tg_z_yyyy_d_0_0_0, tg_z_yyyz_d_0_0_0, tg_z_yyzz_d_0_0_0, tg_z_yzzz_d_0_0_0, tg_z_zzzz_d_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fai_0 = 1.0 / a_exp;

            tg_x_xxxx_d_0_0_0[i] = tg_0_xxxx_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxxy_d_0_0_0[i] = tg_0_xxxy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxxz_d_0_0_0[i] = tg_0_xxxz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxyy_d_0_0_0[i] = tg_0_xxyy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxyz_d_0_0_0[i] = tg_0_xxyz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xxzz_d_0_0_0[i] = tg_0_xxzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xyyy_d_0_0_0[i] = tg_0_xyyy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xyyz_d_0_0_0[i] = tg_0_xyyz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xyzz_d_0_0_0[i] = tg_0_xyzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xzzz_d_0_0_0[i] = tg_0_xzzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_yyyy_d_0_0_0[i] = tg_0_yyyy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_yyyz_d_0_0_0[i] = tg_0_yyyz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_yyzz_d_0_0_0[i] = tg_0_yyzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_yzzz_d_0_0_0[i] = tg_0_yzzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_zzzz_d_0_0_0[i] = tg_0_zzzz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_y_xxxx_d_0_0_0[i] = tg_0_xxxx_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxxy_d_0_0_0[i] = tg_0_xxxy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxxz_d_0_0_0[i] = tg_0_xxxz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxyy_d_0_0_0[i] = tg_0_xxyy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxyz_d_0_0_0[i] = tg_0_xxyz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xxzz_d_0_0_0[i] = tg_0_xxzz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xyyy_d_0_0_0[i] = tg_0_xyyy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xyyz_d_0_0_0[i] = tg_0_xyyz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xyzz_d_0_0_0[i] = tg_0_xyzz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xzzz_d_0_0_0[i] = tg_0_xzzz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_yyyy_d_0_0_0[i] = tg_0_yyyy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_yyyz_d_0_0_0[i] = tg_0_yyyz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_yyzz_d_0_0_0[i] = tg_0_yyzz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_yzzz_d_0_0_0[i] = tg_0_yzzz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_zzzz_d_0_0_0[i] = tg_0_zzzz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_z_xxxx_d_0_0_0[i] = tg_0_xxxx_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxxy_d_0_0_0[i] = tg_0_xxxy_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxxz_d_0_0_0[i] = tg_0_xxxz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxyy_d_0_0_0[i] = tg_0_xxyy_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxyz_d_0_0_0[i] = tg_0_xxyz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xxzz_d_0_0_0[i] = tg_0_xxzz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xyyy_d_0_0_0[i] = tg_0_xyyy_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xyyz_d_0_0_0[i] = tg_0_xyyz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xyzz_d_0_0_0[i] = tg_0_xyzz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xzzz_d_0_0_0[i] = tg_0_xzzz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_yyyy_d_0_0_0[i] = tg_0_yyyy_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_yyyz_d_0_0_0[i] = tg_0_yyyz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_yyzz_d_0_0_0[i] = tg_0_yyzz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_yzzz_d_0_0_0[i] = tg_0_yzzz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_zzzz_d_0_0_0[i] = tg_0_zzzz_d_0_0_1[i] * fai_0 * a_z * fp_0;
        }
    }
}

} // t2pecp namespace

