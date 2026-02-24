#include "ProjectedCorePotentialPrimRecDPForF.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_dp_f(CSimdArray<double>& pbuffer, 
                                        const size_t idx_dp_f_0_0_0,
                                        const size_t idx_sp_f_0_0_0,
                                        const size_t idx_pp_f_0_0_0,
                                        const size_t idx_ps_d_0_0_1,
                                        const size_t idx_pp_d_0_0_1,
                                        const size_t idx_sp_f_1_0_0,
                                        const size_t idx_pp_f_1_0_0,
                                        const size_t idx_sp_p_1_0_1,
                                        const size_t idx_pp_p_1_0_1,
                                        const size_t idx_ps_s_1_1_1,
                                        const size_t idx_pp_s_1_1_1,
                                        const int p,
                                        const size_t idx_sp_f_0_0_1,
                                        const size_t idx_pp_f_0_0_1,
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

    // Set up components of auxiliary buffer : SP

    auto tg_0_x_f_0_0_0 = pbuffer.data(idx_sp_f_0_0_0);

    auto tg_0_y_f_0_0_0 = pbuffer.data(idx_sp_f_0_0_0 + 1);

    auto tg_0_z_f_0_0_0 = pbuffer.data(idx_sp_f_0_0_0 + 2);

    // Set up components of auxiliary buffer : PP

    auto tg_x_x_f_0_0_0 = pbuffer.data(idx_pp_f_0_0_0);

    auto tg_x_y_f_0_0_0 = pbuffer.data(idx_pp_f_0_0_0 + 1);

    auto tg_x_z_f_0_0_0 = pbuffer.data(idx_pp_f_0_0_0 + 2);

    auto tg_y_x_f_0_0_0 = pbuffer.data(idx_pp_f_0_0_0 + 3);

    auto tg_y_y_f_0_0_0 = pbuffer.data(idx_pp_f_0_0_0 + 4);

    auto tg_y_z_f_0_0_0 = pbuffer.data(idx_pp_f_0_0_0 + 5);

    auto tg_z_x_f_0_0_0 = pbuffer.data(idx_pp_f_0_0_0 + 6);

    auto tg_z_y_f_0_0_0 = pbuffer.data(idx_pp_f_0_0_0 + 7);

    auto tg_z_z_f_0_0_0 = pbuffer.data(idx_pp_f_0_0_0 + 8);

    // Set up components of auxiliary buffer : PS

    auto tg_x_0_d_0_0_1 = pbuffer.data(idx_ps_d_0_0_1);

    auto tg_y_0_d_0_0_1 = pbuffer.data(idx_ps_d_0_0_1 + 1);

    auto tg_z_0_d_0_0_1 = pbuffer.data(idx_ps_d_0_0_1 + 2);

    // Set up components of auxiliary buffer : PP

    auto tg_x_x_d_0_0_1 = pbuffer.data(idx_pp_d_0_0_1);

    auto tg_x_y_d_0_0_1 = pbuffer.data(idx_pp_d_0_0_1 + 1);

    auto tg_x_z_d_0_0_1 = pbuffer.data(idx_pp_d_0_0_1 + 2);

    auto tg_y_x_d_0_0_1 = pbuffer.data(idx_pp_d_0_0_1 + 3);

    auto tg_y_y_d_0_0_1 = pbuffer.data(idx_pp_d_0_0_1 + 4);

    auto tg_y_z_d_0_0_1 = pbuffer.data(idx_pp_d_0_0_1 + 5);

    auto tg_z_x_d_0_0_1 = pbuffer.data(idx_pp_d_0_0_1 + 6);

    auto tg_z_y_d_0_0_1 = pbuffer.data(idx_pp_d_0_0_1 + 7);

    auto tg_z_z_d_0_0_1 = pbuffer.data(idx_pp_d_0_0_1 + 8);

    // Set up components of auxiliary buffer : SP

    auto tg_0_x_f_1_0_0 = pbuffer.data(idx_sp_f_1_0_0);

    auto tg_0_y_f_1_0_0 = pbuffer.data(idx_sp_f_1_0_0 + 1);

    auto tg_0_z_f_1_0_0 = pbuffer.data(idx_sp_f_1_0_0 + 2);

    // Set up components of auxiliary buffer : PP

    auto tg_x_x_f_1_0_0 = pbuffer.data(idx_pp_f_1_0_0);

    auto tg_x_y_f_1_0_0 = pbuffer.data(idx_pp_f_1_0_0 + 1);

    auto tg_x_z_f_1_0_0 = pbuffer.data(idx_pp_f_1_0_0 + 2);

    auto tg_y_x_f_1_0_0 = pbuffer.data(idx_pp_f_1_0_0 + 3);

    auto tg_y_y_f_1_0_0 = pbuffer.data(idx_pp_f_1_0_0 + 4);

    auto tg_y_z_f_1_0_0 = pbuffer.data(idx_pp_f_1_0_0 + 5);

    auto tg_z_x_f_1_0_0 = pbuffer.data(idx_pp_f_1_0_0 + 6);

    auto tg_z_y_f_1_0_0 = pbuffer.data(idx_pp_f_1_0_0 + 7);

    auto tg_z_z_f_1_0_0 = pbuffer.data(idx_pp_f_1_0_0 + 8);

    // Set up components of auxiliary buffer : SP

    auto tg_0_x_p_1_0_1 = pbuffer.data(idx_sp_p_1_0_1);

    auto tg_0_y_p_1_0_1 = pbuffer.data(idx_sp_p_1_0_1 + 1);

    auto tg_0_z_p_1_0_1 = pbuffer.data(idx_sp_p_1_0_1 + 2);

    // Set up components of auxiliary buffer : PP

    auto tg_x_x_p_1_0_1 = pbuffer.data(idx_pp_p_1_0_1);

    auto tg_x_y_p_1_0_1 = pbuffer.data(idx_pp_p_1_0_1 + 1);

    auto tg_x_z_p_1_0_1 = pbuffer.data(idx_pp_p_1_0_1 + 2);

    auto tg_y_x_p_1_0_1 = pbuffer.data(idx_pp_p_1_0_1 + 3);

    auto tg_y_y_p_1_0_1 = pbuffer.data(idx_pp_p_1_0_1 + 4);

    auto tg_y_z_p_1_0_1 = pbuffer.data(idx_pp_p_1_0_1 + 5);

    auto tg_z_x_p_1_0_1 = pbuffer.data(idx_pp_p_1_0_1 + 6);

    auto tg_z_y_p_1_0_1 = pbuffer.data(idx_pp_p_1_0_1 + 7);

    auto tg_z_z_p_1_0_1 = pbuffer.data(idx_pp_p_1_0_1 + 8);

    // Set up components of auxiliary buffer : PS

    auto tg_x_0_s_1_1_1 = pbuffer.data(idx_ps_s_1_1_1);

    auto tg_y_0_s_1_1_1 = pbuffer.data(idx_ps_s_1_1_1 + 1);

    auto tg_z_0_s_1_1_1 = pbuffer.data(idx_ps_s_1_1_1 + 2);

    // Set up components of auxiliary buffer : PP

    auto tg_x_x_s_1_1_1 = pbuffer.data(idx_pp_s_1_1_1);

    auto tg_x_y_s_1_1_1 = pbuffer.data(idx_pp_s_1_1_1 + 1);

    auto tg_x_z_s_1_1_1 = pbuffer.data(idx_pp_s_1_1_1 + 2);

    auto tg_y_x_s_1_1_1 = pbuffer.data(idx_pp_s_1_1_1 + 3);

    auto tg_y_y_s_1_1_1 = pbuffer.data(idx_pp_s_1_1_1 + 4);

    auto tg_y_z_s_1_1_1 = pbuffer.data(idx_pp_s_1_1_1 + 5);

    auto tg_z_x_s_1_1_1 = pbuffer.data(idx_pp_s_1_1_1 + 6);

    auto tg_z_y_s_1_1_1 = pbuffer.data(idx_pp_s_1_1_1 + 7);

    auto tg_z_z_s_1_1_1 = pbuffer.data(idx_pp_s_1_1_1 + 8);

    // Set up components of targeted buffer : DP

    auto tg_xx_x_f_0_0_0 = pbuffer.data(idx_dp_f_0_0_0);

    auto tg_xx_y_f_0_0_0 = pbuffer.data(idx_dp_f_0_0_0 + 1);

    auto tg_xx_z_f_0_0_0 = pbuffer.data(idx_dp_f_0_0_0 + 2);

    auto tg_xy_x_f_0_0_0 = pbuffer.data(idx_dp_f_0_0_0 + 3);

    auto tg_xy_y_f_0_0_0 = pbuffer.data(idx_dp_f_0_0_0 + 4);

    auto tg_xy_z_f_0_0_0 = pbuffer.data(idx_dp_f_0_0_0 + 5);

    auto tg_xz_x_f_0_0_0 = pbuffer.data(idx_dp_f_0_0_0 + 6);

    auto tg_xz_y_f_0_0_0 = pbuffer.data(idx_dp_f_0_0_0 + 7);

    auto tg_xz_z_f_0_0_0 = pbuffer.data(idx_dp_f_0_0_0 + 8);

    auto tg_yy_x_f_0_0_0 = pbuffer.data(idx_dp_f_0_0_0 + 9);

    auto tg_yy_y_f_0_0_0 = pbuffer.data(idx_dp_f_0_0_0 + 10);

    auto tg_yy_z_f_0_0_0 = pbuffer.data(idx_dp_f_0_0_0 + 11);

    auto tg_yz_x_f_0_0_0 = pbuffer.data(idx_dp_f_0_0_0 + 12);

    auto tg_yz_y_f_0_0_0 = pbuffer.data(idx_dp_f_0_0_0 + 13);

    auto tg_yz_z_f_0_0_0 = pbuffer.data(idx_dp_f_0_0_0 + 14);

    auto tg_zz_x_f_0_0_0 = pbuffer.data(idx_dp_f_0_0_0 + 15);

    auto tg_zz_y_f_0_0_0 = pbuffer.data(idx_dp_f_0_0_0 + 16);

    auto tg_zz_z_f_0_0_0 = pbuffer.data(idx_dp_f_0_0_0 + 17);

    #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_0_x_f_0_0_0, tg_0_x_f_1_0_0, tg_0_x_p_1_0_1, tg_0_y_f_0_0_0, tg_0_y_f_1_0_0, tg_0_y_p_1_0_1, tg_0_z_f_0_0_0, tg_0_z_f_1_0_0, tg_0_z_p_1_0_1, tg_x_0_d_0_0_1, tg_x_0_s_1_1_1, tg_x_x_d_0_0_1, tg_x_x_f_0_0_0, tg_x_x_f_1_0_0, tg_x_x_p_1_0_1, tg_x_x_s_1_1_1, tg_x_y_d_0_0_1, tg_x_y_f_0_0_0, tg_x_y_f_1_0_0, tg_x_y_p_1_0_1, tg_x_y_s_1_1_1, tg_x_z_d_0_0_1, tg_x_z_f_0_0_0, tg_x_z_f_1_0_0, tg_x_z_p_1_0_1, tg_x_z_s_1_1_1, tg_xx_x_f_0_0_0, tg_xx_y_f_0_0_0, tg_xx_z_f_0_0_0, tg_xy_x_f_0_0_0, tg_xy_y_f_0_0_0, tg_xy_z_f_0_0_0, tg_xz_x_f_0_0_0, tg_xz_y_f_0_0_0, tg_xz_z_f_0_0_0, tg_y_0_d_0_0_1, tg_y_0_s_1_1_1, tg_y_x_d_0_0_1, tg_y_x_f_0_0_0, tg_y_x_f_1_0_0, tg_y_x_p_1_0_1, tg_y_x_s_1_1_1, tg_y_y_d_0_0_1, tg_y_y_f_0_0_0, tg_y_y_f_1_0_0, tg_y_y_p_1_0_1, tg_y_y_s_1_1_1, tg_y_z_d_0_0_1, tg_y_z_f_0_0_0, tg_y_z_f_1_0_0, tg_y_z_p_1_0_1, tg_y_z_s_1_1_1, tg_yy_x_f_0_0_0, tg_yy_y_f_0_0_0, tg_yy_z_f_0_0_0, tg_yz_x_f_0_0_0, tg_yz_y_f_0_0_0, tg_yz_z_f_0_0_0, tg_z_0_d_0_0_1, tg_z_0_s_1_1_1, tg_z_x_d_0_0_1, tg_z_x_f_0_0_0, tg_z_x_f_1_0_0, tg_z_x_p_1_0_1, tg_z_x_s_1_1_1, tg_z_y_d_0_0_1, tg_z_y_f_0_0_0, tg_z_y_f_1_0_0, tg_z_y_p_1_0_1, tg_z_y_s_1_1_1, tg_z_z_d_0_0_1, tg_z_z_f_0_0_0, tg_z_z_f_1_0_0, tg_z_z_p_1_0_1, tg_z_z_s_1_1_1, tg_zz_x_f_0_0_0, tg_zz_y_f_0_0_0, tg_zz_z_f_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double faz_0 = -(b_exps[i] + c_exp) * fzi_0;

        const double fbzi_0 = b_exps[i] * fzi_0;

        const double f2abz_0 = 2.0 * a_exp * b_exps[i] * fzi_0;

            const double fbi_0 = 1.0 / b_exps[i];

            const double fai_0 = 1.0 / a_exp;

        tg_xx_x_f_0_0_0[i] = -7.0 / 2.0 * tg_0_x_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_x_f_0_0_0[i] * fzi_0 + tg_0_x_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 / 2.0 * tg_x_0_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_x_0_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_x_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_x_x_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_x_x_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_x_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_x_f_0_0_0[i] * a_x * faz_0;

        tg_xx_y_f_0_0_0[i] = -7.0 / 2.0 * tg_0_y_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_y_f_0_0_0[i] * fzi_0 + tg_0_y_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_x_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_x_y_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_x_y_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_y_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_y_f_0_0_0[i] * a_x * faz_0;

        tg_xx_z_f_0_0_0[i] = -7.0 / 2.0 * tg_0_z_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_z_f_0_0_0[i] * fzi_0 + tg_0_z_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_x_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_x_z_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_x_z_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_x_z_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_z_f_0_0_0[i] * a_x * faz_0;

        tg_xy_x_f_0_0_0[i] = 7.0 * tg_x_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_x_x_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_x_x_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_x_x_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_x_x_f_0_0_0[i] * a_y * faz_0;

        tg_xy_y_f_0_0_0[i] = 7.0 * tg_y_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_y_y_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_y_y_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_y_y_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_y_f_0_0_0[i] * a_x * faz_0;

        tg_xy_z_f_0_0_0[i] = 7.0 * tg_y_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_y_z_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_y_z_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_y_z_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_z_f_0_0_0[i] * a_x * faz_0;

        tg_xz_x_f_0_0_0[i] = 7.0 * tg_x_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_x_x_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_x_x_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_x_x_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_x_x_f_0_0_0[i] * a_z * faz_0;

        tg_xz_y_f_0_0_0[i] = 7.0 * tg_z_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_z_y_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_z_y_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_z_y_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_y_f_0_0_0[i] * a_x * faz_0;

        tg_xz_z_f_0_0_0[i] = 7.0 * tg_z_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_z_z_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_z_z_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_z_z_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_z_f_0_0_0[i] * a_x * faz_0;

        tg_yy_x_f_0_0_0[i] = -7.0 / 2.0 * tg_0_x_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_x_f_0_0_0[i] * fzi_0 + tg_0_x_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_y_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_y_x_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_y_x_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_x_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_x_f_0_0_0[i] * a_y * faz_0;

        tg_yy_y_f_0_0_0[i] = -7.0 / 2.0 * tg_0_y_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_y_f_0_0_0[i] * fzi_0 + tg_0_y_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 / 2.0 * tg_y_0_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_y_0_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_y_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_y_y_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_y_y_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_y_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_y_f_0_0_0[i] * a_y * faz_0;

        tg_yy_z_f_0_0_0[i] = -7.0 / 2.0 * tg_0_z_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_z_f_0_0_0[i] * fzi_0 + tg_0_z_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_y_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_y_z_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_y_z_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_y_z_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_z_f_0_0_0[i] * a_y * faz_0;

        tg_yz_x_f_0_0_0[i] = 7.0 * tg_z_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_z_x_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_z_x_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_z_x_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_x_f_0_0_0[i] * a_y * faz_0;

        tg_yz_y_f_0_0_0[i] = 7.0 * tg_y_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_y_y_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_y_y_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_y_y_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_y_y_f_0_0_0[i] * a_z * faz_0;

        tg_yz_z_f_0_0_0[i] = 7.0 * tg_z_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_z_z_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_z_z_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_z_z_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_z_f_0_0_0[i] * a_y * faz_0;

        tg_zz_x_f_0_0_0[i] = -7.0 / 2.0 * tg_0_x_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_x_f_0_0_0[i] * fzi_0 + tg_0_x_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_z_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_z_x_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_z_x_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_x_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_x_f_0_0_0[i] * a_z * faz_0;

        tg_zz_y_f_0_0_0[i] = -7.0 / 2.0 * tg_0_y_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_y_f_0_0_0[i] * fzi_0 + tg_0_y_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_z_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_z_y_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_z_y_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_y_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_y_f_0_0_0[i] * a_z * faz_0;

        tg_zz_z_f_0_0_0[i] = -7.0 / 2.0 * tg_0_z_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_0_z_f_0_0_0[i] * fzi_0 + tg_0_z_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 / 2.0 * tg_z_0_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_z_0_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_z_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_z_z_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_z_z_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_z_z_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_z_f_0_0_0[i] * a_z * faz_0;
    }

    if (p > 0)
    {
        const double fp_0 = (double)p;

        // Set up components of auxiliary buffer : SP

        auto tg_0_x_f_0_0_1 = pbuffer.data(idx_sp_f_0_0_1);

        auto tg_0_y_f_0_0_1 = pbuffer.data(idx_sp_f_0_0_1 + 1);

        auto tg_0_z_f_0_0_1 = pbuffer.data(idx_sp_f_0_0_1 + 2);

        // Set up components of auxiliary buffer : PP

        auto tg_x_x_f_0_0_1 = pbuffer.data(idx_pp_f_0_0_1);

        auto tg_x_y_f_0_0_1 = pbuffer.data(idx_pp_f_0_0_1 + 1);

        auto tg_x_z_f_0_0_1 = pbuffer.data(idx_pp_f_0_0_1 + 2);

        auto tg_y_x_f_0_0_1 = pbuffer.data(idx_pp_f_0_0_1 + 3);

        auto tg_y_y_f_0_0_1 = pbuffer.data(idx_pp_f_0_0_1 + 4);

        auto tg_y_z_f_0_0_1 = pbuffer.data(idx_pp_f_0_0_1 + 5);

        auto tg_z_x_f_0_0_1 = pbuffer.data(idx_pp_f_0_0_1 + 6);

        auto tg_z_y_f_0_0_1 = pbuffer.data(idx_pp_f_0_0_1 + 7);

        auto tg_z_z_f_0_0_1 = pbuffer.data(idx_pp_f_0_0_1 + 8);

        #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_0_x_f_0_0_1, tg_0_y_f_0_0_1, tg_0_z_f_0_0_1, tg_x_x_f_0_0_1, tg_x_y_f_0_0_1, tg_x_z_f_0_0_1, tg_xx_x_f_0_0_0, tg_xx_y_f_0_0_0, tg_xx_z_f_0_0_0, tg_xy_x_f_0_0_0, tg_xy_y_f_0_0_0, tg_xy_z_f_0_0_0, tg_xz_x_f_0_0_0, tg_xz_y_f_0_0_0, tg_xz_z_f_0_0_0, tg_y_x_f_0_0_1, tg_y_y_f_0_0_1, tg_y_z_f_0_0_1, tg_yy_x_f_0_0_0, tg_yy_y_f_0_0_0, tg_yy_z_f_0_0_0, tg_yz_x_f_0_0_0, tg_yz_y_f_0_0_0, tg_yz_z_f_0_0_0, tg_z_x_f_0_0_1, tg_z_y_f_0_0_1, tg_z_z_f_0_0_1, tg_zz_x_f_0_0_0, tg_zz_y_f_0_0_0, tg_zz_z_f_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fai_0 = 1.0 / a_exp;

            tg_xx_x_f_0_0_0[i] += 1.0 / 2.0 * tg_0_x_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_x_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_y_f_0_0_0[i] += 1.0 / 2.0 * tg_0_y_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_y_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_z_f_0_0_0[i] += 1.0 / 2.0 * tg_0_z_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_z_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_x_f_0_0_0[i] += tg_y_x_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_y_f_0_0_0[i] += tg_y_y_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_z_f_0_0_0[i] += tg_y_z_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_x_f_0_0_0[i] += tg_z_x_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_y_f_0_0_0[i] += tg_z_y_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_z_f_0_0_0[i] += tg_z_z_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_yy_x_f_0_0_0[i] += 1.0 / 2.0 * tg_0_x_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_x_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_y_f_0_0_0[i] += 1.0 / 2.0 * tg_0_y_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_y_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_z_f_0_0_0[i] += 1.0 / 2.0 * tg_0_z_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_z_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_x_f_0_0_0[i] += tg_z_x_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_y_f_0_0_0[i] += tg_z_y_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_z_f_0_0_0[i] += tg_z_z_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_zz_x_f_0_0_0[i] += 1.0 / 2.0 * tg_0_x_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_x_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_y_f_0_0_0[i] += 1.0 / 2.0 * tg_0_y_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_y_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_z_f_0_0_0[i] += 1.0 / 2.0 * tg_0_z_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_z_f_0_0_1[i] * fai_0 * a_z * fp_0;
        }
    }
}

} // t2pecp namespace

