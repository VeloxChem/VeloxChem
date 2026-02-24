#include "ProjectedCorePotentialPrimRecPDForD.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_pd_d(CSimdArray<double>& pbuffer, 
                                        const size_t idx_pd_d_0_0_0,
                                        const size_t idx_sd_d_0_0_0,
                                        const size_t idx_sp_p_0_0_1,
                                        const size_t idx_sd_p_0_0_1,
                                        const size_t idx_sd_d_1_0_0,
                                        const size_t idx_sd_s_1_0_1,
                                        const int p,
                                        const size_t idx_sd_d_0_0_1,
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

    // Set up components of auxiliary buffer : SD

    auto tg_0_xx_d_0_0_0 = pbuffer.data(idx_sd_d_0_0_0);

    auto tg_0_xy_d_0_0_0 = pbuffer.data(idx_sd_d_0_0_0 + 1);

    auto tg_0_xz_d_0_0_0 = pbuffer.data(idx_sd_d_0_0_0 + 2);

    auto tg_0_yy_d_0_0_0 = pbuffer.data(idx_sd_d_0_0_0 + 3);

    auto tg_0_yz_d_0_0_0 = pbuffer.data(idx_sd_d_0_0_0 + 4);

    auto tg_0_zz_d_0_0_0 = pbuffer.data(idx_sd_d_0_0_0 + 5);

    // Set up components of auxiliary buffer : SP

    auto tg_0_x_p_0_0_1 = pbuffer.data(idx_sp_p_0_0_1);

    auto tg_0_y_p_0_0_1 = pbuffer.data(idx_sp_p_0_0_1 + 1);

    auto tg_0_z_p_0_0_1 = pbuffer.data(idx_sp_p_0_0_1 + 2);

    // Set up components of auxiliary buffer : SD

    auto tg_0_xx_p_0_0_1 = pbuffer.data(idx_sd_p_0_0_1);

    auto tg_0_xy_p_0_0_1 = pbuffer.data(idx_sd_p_0_0_1 + 1);

    auto tg_0_xz_p_0_0_1 = pbuffer.data(idx_sd_p_0_0_1 + 2);

    auto tg_0_yy_p_0_0_1 = pbuffer.data(idx_sd_p_0_0_1 + 3);

    auto tg_0_yz_p_0_0_1 = pbuffer.data(idx_sd_p_0_0_1 + 4);

    auto tg_0_zz_p_0_0_1 = pbuffer.data(idx_sd_p_0_0_1 + 5);

    // Set up components of auxiliary buffer : SD

    auto tg_0_xx_d_1_0_0 = pbuffer.data(idx_sd_d_1_0_0);

    auto tg_0_xy_d_1_0_0 = pbuffer.data(idx_sd_d_1_0_0 + 1);

    auto tg_0_xz_d_1_0_0 = pbuffer.data(idx_sd_d_1_0_0 + 2);

    auto tg_0_yy_d_1_0_0 = pbuffer.data(idx_sd_d_1_0_0 + 3);

    auto tg_0_yz_d_1_0_0 = pbuffer.data(idx_sd_d_1_0_0 + 4);

    auto tg_0_zz_d_1_0_0 = pbuffer.data(idx_sd_d_1_0_0 + 5);

    // Set up components of auxiliary buffer : SD

    auto tg_0_xx_s_1_0_1 = pbuffer.data(idx_sd_s_1_0_1);

    auto tg_0_xy_s_1_0_1 = pbuffer.data(idx_sd_s_1_0_1 + 1);

    auto tg_0_xz_s_1_0_1 = pbuffer.data(idx_sd_s_1_0_1 + 2);

    auto tg_0_yy_s_1_0_1 = pbuffer.data(idx_sd_s_1_0_1 + 3);

    auto tg_0_yz_s_1_0_1 = pbuffer.data(idx_sd_s_1_0_1 + 4);

    auto tg_0_zz_s_1_0_1 = pbuffer.data(idx_sd_s_1_0_1 + 5);

    // Set up components of targeted buffer : PD

    auto tg_x_xx_d_0_0_0 = pbuffer.data(idx_pd_d_0_0_0);

    auto tg_x_xy_d_0_0_0 = pbuffer.data(idx_pd_d_0_0_0 + 1);

    auto tg_x_xz_d_0_0_0 = pbuffer.data(idx_pd_d_0_0_0 + 2);

    auto tg_x_yy_d_0_0_0 = pbuffer.data(idx_pd_d_0_0_0 + 3);

    auto tg_x_yz_d_0_0_0 = pbuffer.data(idx_pd_d_0_0_0 + 4);

    auto tg_x_zz_d_0_0_0 = pbuffer.data(idx_pd_d_0_0_0 + 5);

    auto tg_y_xx_d_0_0_0 = pbuffer.data(idx_pd_d_0_0_0 + 6);

    auto tg_y_xy_d_0_0_0 = pbuffer.data(idx_pd_d_0_0_0 + 7);

    auto tg_y_xz_d_0_0_0 = pbuffer.data(idx_pd_d_0_0_0 + 8);

    auto tg_y_yy_d_0_0_0 = pbuffer.data(idx_pd_d_0_0_0 + 9);

    auto tg_y_yz_d_0_0_0 = pbuffer.data(idx_pd_d_0_0_0 + 10);

    auto tg_y_zz_d_0_0_0 = pbuffer.data(idx_pd_d_0_0_0 + 11);

    auto tg_z_xx_d_0_0_0 = pbuffer.data(idx_pd_d_0_0_0 + 12);

    auto tg_z_xy_d_0_0_0 = pbuffer.data(idx_pd_d_0_0_0 + 13);

    auto tg_z_xz_d_0_0_0 = pbuffer.data(idx_pd_d_0_0_0 + 14);

    auto tg_z_yy_d_0_0_0 = pbuffer.data(idx_pd_d_0_0_0 + 15);

    auto tg_z_yz_d_0_0_0 = pbuffer.data(idx_pd_d_0_0_0 + 16);

    auto tg_z_zz_d_0_0_0 = pbuffer.data(idx_pd_d_0_0_0 + 17);

    #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_0_xx_d_0_0_0, tg_0_xx_d_1_0_0, tg_0_xx_p_0_0_1, tg_0_xx_s_1_0_1, tg_0_xy_d_0_0_0, tg_0_xy_d_1_0_0, tg_0_xy_p_0_0_1, tg_0_xy_s_1_0_1, tg_0_xz_d_0_0_0, tg_0_xz_d_1_0_0, tg_0_xz_p_0_0_1, tg_0_xz_s_1_0_1, tg_0_yy_d_0_0_0, tg_0_yy_d_1_0_0, tg_0_yy_p_0_0_1, tg_0_yy_s_1_0_1, tg_0_yz_d_0_0_0, tg_0_yz_d_1_0_0, tg_0_yz_p_0_0_1, tg_0_yz_s_1_0_1, tg_0_zz_d_0_0_0, tg_0_zz_d_1_0_0, tg_0_zz_p_0_0_1, tg_0_zz_s_1_0_1, tg_x_xx_d_0_0_0, tg_x_xy_d_0_0_0, tg_x_xz_d_0_0_0, tg_x_yy_d_0_0_0, tg_x_yz_d_0_0_0, tg_x_zz_d_0_0_0, tg_y_xx_d_0_0_0, tg_y_xy_d_0_0_0, tg_y_xz_d_0_0_0, tg_y_yy_d_0_0_0, tg_y_yz_d_0_0_0, tg_y_zz_d_0_0_0, tg_z_xx_d_0_0_0, tg_z_xy_d_0_0_0, tg_z_xz_d_0_0_0, tg_z_yy_d_0_0_0, tg_z_yz_d_0_0_0, tg_z_zz_d_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double faz_0 = -(b_exps[i] + c_exp) * fzi_0;

        const double fbzi_0 = b_exps[i] * fzi_0;

        const double f2abz_0 = 2.0 * a_exp * b_exps[i] * fzi_0;

            const double fbi_0 = 1.0 / b_exps[i];

        tg_x_xx_d_0_0_0[i] = -5.0 * tg_0_xx_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_0_xx_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_xx_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xx_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xx_d_0_0_0[i] * a_x * faz_0;

        tg_x_xy_d_0_0_0[i] = -5.0 * tg_0_xy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 / 2.0 * tg_0_xy_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_xy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xy_d_0_0_0[i] * a_x * faz_0;

        tg_x_xz_d_0_0_0[i] = -5.0 * tg_0_xz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 / 2.0 * tg_0_xz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_xz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_xz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_xz_d_0_0_0[i] * a_x * faz_0;

        tg_x_yy_d_0_0_0[i] = -5.0 * tg_0_yy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_0_yy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_yy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_yy_d_0_0_0[i] * a_x * faz_0;

        tg_x_yz_d_0_0_0[i] = -5.0 * tg_0_yz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_0_yz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_yz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_yz_d_0_0_0[i] * a_x * faz_0;

        tg_x_zz_d_0_0_0[i] = -5.0 * tg_0_zz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_0_zz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_zz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_zz_d_0_0_0[i] * a_x * faz_0;

        tg_y_xx_d_0_0_0[i] = -5.0 * tg_0_xx_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_0_xx_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xx_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xx_d_0_0_0[i] * a_y * faz_0;

        tg_y_xy_d_0_0_0[i] = -5.0 * tg_0_xy_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 / 2.0 * tg_0_xy_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_xy_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xy_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xy_d_0_0_0[i] * a_y * faz_0;

        tg_y_xz_d_0_0_0[i] = -5.0 * tg_0_xz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_0_xz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_xz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_xz_d_0_0_0[i] * a_y * faz_0;

        tg_y_yy_d_0_0_0[i] = -5.0 * tg_0_yy_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_0_yy_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_yy_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_yy_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_yy_d_0_0_0[i] * a_y * faz_0;

        tg_y_yz_d_0_0_0[i] = -5.0 * tg_0_yz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 / 2.0 * tg_0_yz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_yz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_yz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_yz_d_0_0_0[i] * a_y * faz_0;

        tg_y_zz_d_0_0_0[i] = -5.0 * tg_0_zz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_0_zz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_zz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_zz_d_0_0_0[i] * a_y * faz_0;

        tg_z_xx_d_0_0_0[i] = -5.0 * tg_0_xx_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_0_xx_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xx_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xx_d_0_0_0[i] * a_z * faz_0;

        tg_z_xy_d_0_0_0[i] = -5.0 * tg_0_xy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_0_xy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xy_d_0_0_0[i] * a_z * faz_0;

        tg_z_xz_d_0_0_0[i] = -5.0 * tg_0_xz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 / 2.0 * tg_0_xz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_xz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_xz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_xz_d_0_0_0[i] * a_z * faz_0;

        tg_z_yy_d_0_0_0[i] = -5.0 * tg_0_yy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_0_yy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_yy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_yy_d_0_0_0[i] * a_z * faz_0;

        tg_z_yz_d_0_0_0[i] = -5.0 * tg_0_yz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 / 2.0 * tg_0_yz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_yz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_yz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_yz_d_0_0_0[i] * a_z * faz_0;

        tg_z_zz_d_0_0_0[i] = -5.0 * tg_0_zz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_0_zz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_0_zz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_zz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_zz_d_0_0_0[i] * a_z * faz_0;
    }

    if (p > 0)
    {
        const double fp_0 = (double)p;

        // Set up components of auxiliary buffer : SD

        auto tg_0_xx_d_0_0_1 = pbuffer.data(idx_sd_d_0_0_1);

        auto tg_0_xy_d_0_0_1 = pbuffer.data(idx_sd_d_0_0_1 + 1);

        auto tg_0_xz_d_0_0_1 = pbuffer.data(idx_sd_d_0_0_1 + 2);

        auto tg_0_yy_d_0_0_1 = pbuffer.data(idx_sd_d_0_0_1 + 3);

        auto tg_0_yz_d_0_0_1 = pbuffer.data(idx_sd_d_0_0_1 + 4);

        auto tg_0_zz_d_0_0_1 = pbuffer.data(idx_sd_d_0_0_1 + 5);

        #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_0_xx_d_0_0_1, tg_0_xy_d_0_0_1, tg_0_xz_d_0_0_1, tg_0_yy_d_0_0_1, tg_0_yz_d_0_0_1, tg_0_zz_d_0_0_1, tg_x_xx_d_0_0_0, tg_x_xy_d_0_0_0, tg_x_xz_d_0_0_0, tg_x_yy_d_0_0_0, tg_x_yz_d_0_0_0, tg_x_zz_d_0_0_0, tg_y_xx_d_0_0_0, tg_y_xy_d_0_0_0, tg_y_xz_d_0_0_0, tg_y_yy_d_0_0_0, tg_y_yz_d_0_0_0, tg_y_zz_d_0_0_0, tg_z_xx_d_0_0_0, tg_z_xy_d_0_0_0, tg_z_xz_d_0_0_0, tg_z_yy_d_0_0_0, tg_z_yz_d_0_0_0, tg_z_zz_d_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fai_0 = 1.0 / a_exp;

            tg_x_xx_d_0_0_0[i] += tg_0_xx_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xy_d_0_0_0[i] += tg_0_xy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_xz_d_0_0_0[i] += tg_0_xz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_yy_d_0_0_0[i] += tg_0_yy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_yz_d_0_0_0[i] += tg_0_yz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_zz_d_0_0_0[i] += tg_0_zz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_y_xx_d_0_0_0[i] += tg_0_xx_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xy_d_0_0_0[i] += tg_0_xy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_xz_d_0_0_0[i] += tg_0_xz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_yy_d_0_0_0[i] += tg_0_yy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_yz_d_0_0_0[i] += tg_0_yz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_zz_d_0_0_0[i] += tg_0_zz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_z_xx_d_0_0_0[i] += tg_0_xx_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xy_d_0_0_0[i] += tg_0_xy_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_xz_d_0_0_0[i] += tg_0_xz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_yy_d_0_0_0[i] += tg_0_yy_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_yz_d_0_0_0[i] += tg_0_yz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_zz_d_0_0_0[i] += tg_0_zz_d_0_0_1[i] * fai_0 * a_z * fp_0;
        }
    }
}

} // t2pecp namespace

