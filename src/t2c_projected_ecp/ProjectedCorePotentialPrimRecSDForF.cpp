#include "ProjectedCorePotentialPrimRecSDForF.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_sd_f(CSimdArray<double>& pbuffer, 
                                        const size_t idx_sd_f_0_0_0,
                                        const size_t idx_ss_f_0_0_0,
                                        const size_t idx_sp_f_0_0_0,
                                        const size_t idx_sp_d_0_0_1,
                                        const size_t idx_ss_f_0_1_0,
                                        const size_t idx_sp_f_0_1_0,
                                        const size_t idx_ss_p_0_1_1,
                                        const size_t idx_sp_p_0_1_1,
                                        const size_t idx_sp_s_1_1_1,
                                        const int m,
                                        const size_t idx_ss_f_0_0_1,
                                        const size_t idx_sp_f_0_0_1,
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

    // Set up components of auxiliary buffer : SS

    auto tg_0_0_f_0_0_0 = pbuffer.data(idx_ss_f_0_0_0);

    // Set up components of auxiliary buffer : SP

    auto tg_0_x_f_0_0_0 = pbuffer.data(idx_sp_f_0_0_0);

    auto tg_0_y_f_0_0_0 = pbuffer.data(idx_sp_f_0_0_0 + 1);

    auto tg_0_z_f_0_0_0 = pbuffer.data(idx_sp_f_0_0_0 + 2);

    // Set up components of auxiliary buffer : SP

    auto tg_0_x_d_0_0_1 = pbuffer.data(idx_sp_d_0_0_1);

    auto tg_0_y_d_0_0_1 = pbuffer.data(idx_sp_d_0_0_1 + 1);

    auto tg_0_z_d_0_0_1 = pbuffer.data(idx_sp_d_0_0_1 + 2);

    // Set up components of auxiliary buffer : SS

    auto tg_0_0_f_0_1_0 = pbuffer.data(idx_ss_f_0_1_0);

    // Set up components of auxiliary buffer : SP

    auto tg_0_x_f_0_1_0 = pbuffer.data(idx_sp_f_0_1_0);

    auto tg_0_y_f_0_1_0 = pbuffer.data(idx_sp_f_0_1_0 + 1);

    auto tg_0_z_f_0_1_0 = pbuffer.data(idx_sp_f_0_1_0 + 2);

    // Set up components of auxiliary buffer : SS

    auto tg_0_0_p_0_1_1 = pbuffer.data(idx_ss_p_0_1_1);

    // Set up components of auxiliary buffer : SP

    auto tg_0_x_p_0_1_1 = pbuffer.data(idx_sp_p_0_1_1);

    auto tg_0_y_p_0_1_1 = pbuffer.data(idx_sp_p_0_1_1 + 1);

    auto tg_0_z_p_0_1_1 = pbuffer.data(idx_sp_p_0_1_1 + 2);

    // Set up components of auxiliary buffer : SP

    auto tg_0_x_s_1_1_1 = pbuffer.data(idx_sp_s_1_1_1);

    auto tg_0_y_s_1_1_1 = pbuffer.data(idx_sp_s_1_1_1 + 1);

    auto tg_0_z_s_1_1_1 = pbuffer.data(idx_sp_s_1_1_1 + 2);

    // Set up components of targeted buffer : SD

    auto tg_0_xx_f_0_0_0 = pbuffer.data(idx_sd_f_0_0_0);

    auto tg_0_xy_f_0_0_0 = pbuffer.data(idx_sd_f_0_0_0 + 1);

    auto tg_0_xz_f_0_0_0 = pbuffer.data(idx_sd_f_0_0_0 + 2);

    auto tg_0_yy_f_0_0_0 = pbuffer.data(idx_sd_f_0_0_0 + 3);

    auto tg_0_yz_f_0_0_0 = pbuffer.data(idx_sd_f_0_0_0 + 4);

    auto tg_0_zz_f_0_0_0 = pbuffer.data(idx_sd_f_0_0_0 + 5);

    #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_0_0_f_0_0_0, tg_0_0_f_0_1_0, tg_0_x_d_0_0_1, tg_0_x_f_0_0_0, tg_0_x_f_0_1_0, tg_0_x_p_0_1_1, tg_0_x_s_1_1_1, tg_0_xx_f_0_0_0, tg_0_xy_f_0_0_0, tg_0_xz_f_0_0_0, tg_0_y_d_0_0_1, tg_0_y_f_0_0_0, tg_0_y_f_0_1_0, tg_0_y_p_0_1_1, tg_0_y_s_1_1_1, tg_0_yy_f_0_0_0, tg_0_yz_f_0_0_0, tg_0_z_d_0_0_1, tg_0_z_f_0_0_0, tg_0_z_f_0_1_0, tg_0_z_p_0_1_1, tg_0_z_s_1_1_1, tg_0_zz_f_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double fbz_0 = -(a_exp + c_exp) * fzi_0;

        const double fazi_0 = a_exp * fzi_0;

        const double fb_0 = b_exps[i];

        const double f2abz_0 = 2.0 * a_exp * b_exps[i] * fzi_0;

        tg_0_xx_f_0_0_0[i] = tg_0_0_f_0_0_0[i] * fzi_0 + tg_0_0_f_0_1_0[i] * fazi_0 * fazi_0 + 7.0 * tg_0_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * a_x * fazi_0 - 7.0 * tg_0_x_p_0_1_1[i] * f2abz_0 * rb_x[i] * fazi_0 + 7.0 * tg_0_x_d_0_0_1[i] * a_x * fazi_0 + 2.0 * tg_0_x_f_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_x_f_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xy_f_0_0_0[i] = 7.0 * tg_0_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * a_x * fazi_0 - 7.0 * tg_0_y_p_0_1_1[i] * f2abz_0 * rb_x[i] * fazi_0 + 7.0 * tg_0_y_d_0_0_1[i] * a_x * fazi_0 + 2.0 * tg_0_y_f_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_y_f_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_xz_f_0_0_0[i] = 7.0 * tg_0_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * a_x * fazi_0 - 7.0 * tg_0_z_p_0_1_1[i] * f2abz_0 * rb_x[i] * fazi_0 + 7.0 * tg_0_z_d_0_0_1[i] * a_x * fazi_0 + 2.0 * tg_0_z_f_0_1_0[i] * rb_x[i] * fazi_0 * fazi_0 * fb_0 + tg_0_z_f_0_0_0[i] * rb_x[i] * fbz_0;

        tg_0_yy_f_0_0_0[i] = tg_0_0_f_0_0_0[i] * fzi_0 + tg_0_0_f_0_1_0[i] * fazi_0 * fazi_0 + 7.0 * tg_0_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * a_y * fazi_0 - 7.0 * tg_0_y_p_0_1_1[i] * f2abz_0 * rb_y[i] * fazi_0 + 7.0 * tg_0_y_d_0_0_1[i] * a_y * fazi_0 + 2.0 * tg_0_y_f_0_1_0[i] * rb_y[i] * fazi_0 * fazi_0 * fb_0 + tg_0_y_f_0_0_0[i] * rb_y[i] * fbz_0;

        tg_0_yz_f_0_0_0[i] = 7.0 * tg_0_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * a_y * fazi_0 - 7.0 * tg_0_z_p_0_1_1[i] * f2abz_0 * rb_y[i] * fazi_0 + 7.0 * tg_0_z_d_0_0_1[i] * a_y * fazi_0 + 2.0 * tg_0_z_f_0_1_0[i] * rb_y[i] * fazi_0 * fazi_0 * fb_0 + tg_0_z_f_0_0_0[i] * rb_y[i] * fbz_0;

        tg_0_zz_f_0_0_0[i] = tg_0_0_f_0_0_0[i] * fzi_0 + tg_0_0_f_0_1_0[i] * fazi_0 * fazi_0 + 7.0 * tg_0_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * a_z * fazi_0 - 7.0 * tg_0_z_p_0_1_1[i] * f2abz_0 * rb_z[i] * fazi_0 + 7.0 * tg_0_z_d_0_0_1[i] * a_z * fazi_0 + 2.0 * tg_0_z_f_0_1_0[i] * rb_z[i] * fazi_0 * fazi_0 * fb_0 + tg_0_z_f_0_0_0[i] * rb_z[i] * fbz_0;
    }

    if (m > 0)
    {
        const double fm_0 = (double)m;

        // Set up components of auxiliary buffer : SS

        auto tg_0_0_f_0_0_1 = pbuffer.data(idx_ss_f_0_0_1);

        // Set up components of auxiliary buffer : SP

        auto tg_0_x_f_0_0_1 = pbuffer.data(idx_sp_f_0_0_1);

        auto tg_0_y_f_0_0_1 = pbuffer.data(idx_sp_f_0_0_1 + 1);

        auto tg_0_z_f_0_0_1 = pbuffer.data(idx_sp_f_0_0_1 + 2);

        #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_0_0_f_0_0_1, tg_0_x_f_0_0_1, tg_0_xx_f_0_0_0, tg_0_xy_f_0_0_0, tg_0_xz_f_0_0_0, tg_0_y_f_0_0_1, tg_0_yy_f_0_0_0, tg_0_yz_f_0_0_0, tg_0_z_f_0_0_1, tg_0_zz_f_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fbi_0 = 1.0 / b_exps[i];

            tg_0_xx_f_0_0_0[i] += 1.0 / 2.0 * tg_0_0_f_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_x_f_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xy_f_0_0_0[i] += tg_0_y_f_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_xz_f_0_0_0[i] += tg_0_z_f_0_0_1[i] * fbi_0 * rb_x[i] * fm_0;

            tg_0_yy_f_0_0_0[i] += 1.0 / 2.0 * tg_0_0_f_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_y_f_0_0_1[i] * fbi_0 * rb_y[i] * fm_0;

            tg_0_yz_f_0_0_0[i] += tg_0_z_f_0_0_1[i] * fbi_0 * rb_y[i] * fm_0;

            tg_0_zz_f_0_0_0[i] += 1.0 / 2.0 * tg_0_0_f_0_0_1[i] * fbi_0 * fbi_0 * fm_0 + tg_0_z_f_0_0_1[i] * fbi_0 * rb_z[i] * fm_0;
        }
    }
}

} // t2pecp namespace

