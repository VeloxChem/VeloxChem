#include "ProjectedCorePotentialPrimRecPPForF.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_pp_f(CSimdArray<double>& pbuffer, 
                                        const size_t idx_pp_f_0_0_0,
                                        const size_t idx_sp_f_0_0_0,
                                        const size_t idx_ss_d_0_0_1,
                                        const size_t idx_sp_d_0_0_1,
                                        const size_t idx_sp_f_1_0_0,
                                        const size_t idx_sp_p_1_0_1,
                                        const size_t idx_ss_s_1_1_1,
                                        const size_t idx_sp_s_1_1_1,
                                        const int p,
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

    // Set up components of auxiliary buffer : SP

    auto tg_0_x_f_0_0_0 = pbuffer.data(idx_sp_f_0_0_0);

    auto tg_0_y_f_0_0_0 = pbuffer.data(idx_sp_f_0_0_0 + 1);

    auto tg_0_z_f_0_0_0 = pbuffer.data(idx_sp_f_0_0_0 + 2);

    // Set up components of auxiliary buffer : SS

    auto tg_0_0_d_0_0_1 = pbuffer.data(idx_ss_d_0_0_1);

    // Set up components of auxiliary buffer : SP

    auto tg_0_x_d_0_0_1 = pbuffer.data(idx_sp_d_0_0_1);

    auto tg_0_y_d_0_0_1 = pbuffer.data(idx_sp_d_0_0_1 + 1);

    auto tg_0_z_d_0_0_1 = pbuffer.data(idx_sp_d_0_0_1 + 2);

    // Set up components of auxiliary buffer : SP

    auto tg_0_x_f_1_0_0 = pbuffer.data(idx_sp_f_1_0_0);

    auto tg_0_y_f_1_0_0 = pbuffer.data(idx_sp_f_1_0_0 + 1);

    auto tg_0_z_f_1_0_0 = pbuffer.data(idx_sp_f_1_0_0 + 2);

    // Set up components of auxiliary buffer : SP

    auto tg_0_x_p_1_0_1 = pbuffer.data(idx_sp_p_1_0_1);

    auto tg_0_y_p_1_0_1 = pbuffer.data(idx_sp_p_1_0_1 + 1);

    auto tg_0_z_p_1_0_1 = pbuffer.data(idx_sp_p_1_0_1 + 2);

    // Set up components of auxiliary buffer : SS

    auto tg_0_0_s_1_1_1 = pbuffer.data(idx_ss_s_1_1_1);

    // Set up components of auxiliary buffer : SP

    auto tg_0_x_s_1_1_1 = pbuffer.data(idx_sp_s_1_1_1);

    auto tg_0_y_s_1_1_1 = pbuffer.data(idx_sp_s_1_1_1 + 1);

    auto tg_0_z_s_1_1_1 = pbuffer.data(idx_sp_s_1_1_1 + 2);

    // Set up components of targeted buffer : PP

    auto tg_x_x_f_0_0_0 = pbuffer.data(idx_pp_f_0_0_0);

    auto tg_x_y_f_0_0_0 = pbuffer.data(idx_pp_f_0_0_0 + 1);

    auto tg_x_z_f_0_0_0 = pbuffer.data(idx_pp_f_0_0_0 + 2);

    auto tg_y_x_f_0_0_0 = pbuffer.data(idx_pp_f_0_0_0 + 3);

    auto tg_y_y_f_0_0_0 = pbuffer.data(idx_pp_f_0_0_0 + 4);

    auto tg_y_z_f_0_0_0 = pbuffer.data(idx_pp_f_0_0_0 + 5);

    auto tg_z_x_f_0_0_0 = pbuffer.data(idx_pp_f_0_0_0 + 6);

    auto tg_z_y_f_0_0_0 = pbuffer.data(idx_pp_f_0_0_0 + 7);

    auto tg_z_z_f_0_0_0 = pbuffer.data(idx_pp_f_0_0_0 + 8);

    #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_0_0_d_0_0_1, tg_0_0_s_1_1_1, tg_0_x_d_0_0_1, tg_0_x_f_0_0_0, tg_0_x_f_1_0_0, tg_0_x_p_1_0_1, tg_0_x_s_1_1_1, tg_0_y_d_0_0_1, tg_0_y_f_0_0_0, tg_0_y_f_1_0_0, tg_0_y_p_1_0_1, tg_0_y_s_1_1_1, tg_0_z_d_0_0_1, tg_0_z_f_0_0_0, tg_0_z_f_1_0_0, tg_0_z_p_1_0_1, tg_0_z_s_1_1_1, tg_x_x_f_0_0_0, tg_x_y_f_0_0_0, tg_x_z_f_0_0_0, tg_y_x_f_0_0_0, tg_y_y_f_0_0_0, tg_y_z_f_0_0_0, tg_z_x_f_0_0_0, tg_z_y_f_0_0_0, tg_z_z_f_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double faz_0 = -(b_exps[i] + c_exp) * fzi_0;

        const double fbzi_0 = b_exps[i] * fzi_0;

        const double f2abz_0 = 2.0 * a_exp * b_exps[i] * fzi_0;

            const double fbi_0 = 1.0 / b_exps[i];

        tg_x_x_f_0_0_0[i] = 7.0 / 2.0 * tg_0_0_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_0_0_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_0_x_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_0_x_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_x_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_x_f_0_0_0[i] * a_x * faz_0;

        tg_x_y_f_0_0_0[i] = 7.0 * tg_0_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_0_y_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_0_y_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_y_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_y_f_0_0_0[i] * a_x * faz_0;

        tg_x_z_f_0_0_0[i] = 7.0 * tg_0_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_0_z_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_0_z_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_0_z_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_z_f_0_0_0[i] * a_x * faz_0;

        tg_y_x_f_0_0_0[i] = 7.0 * tg_0_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_0_x_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_0_x_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_x_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_x_f_0_0_0[i] * a_y * faz_0;

        tg_y_y_f_0_0_0[i] = 7.0 / 2.0 * tg_0_0_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_0_0_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_0_y_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_0_y_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_y_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_y_f_0_0_0[i] * a_y * faz_0;

        tg_y_z_f_0_0_0[i] = 7.0 * tg_0_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_0_z_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_0_z_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_0_z_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_z_f_0_0_0[i] * a_y * faz_0;

        tg_z_x_f_0_0_0[i] = 7.0 * tg_0_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_0_x_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_0_x_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_x_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_x_f_0_0_0[i] * a_z * faz_0;

        tg_z_y_f_0_0_0[i] = 7.0 * tg_0_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_0_y_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_0_y_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_y_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_y_f_0_0_0[i] * a_z * faz_0;

        tg_z_z_f_0_0_0[i] = 7.0 / 2.0 * tg_0_0_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_0_0_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_0_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_0_z_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_0_z_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_0_z_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_z_f_0_0_0[i] * a_z * faz_0;
    }

    if (p > 0)
    {
        const double fp_0 = (double)p;

        // Set up components of auxiliary buffer : SP

        auto tg_0_x_f_0_0_1 = pbuffer.data(idx_sp_f_0_0_1);

        auto tg_0_y_f_0_0_1 = pbuffer.data(idx_sp_f_0_0_1 + 1);

        auto tg_0_z_f_0_0_1 = pbuffer.data(idx_sp_f_0_0_1 + 2);

        #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_0_x_f_0_0_1, tg_0_y_f_0_0_1, tg_0_z_f_0_0_1, tg_x_x_f_0_0_0, tg_x_y_f_0_0_0, tg_x_z_f_0_0_0, tg_y_x_f_0_0_0, tg_y_y_f_0_0_0, tg_y_z_f_0_0_0, tg_z_x_f_0_0_0, tg_z_y_f_0_0_0, tg_z_z_f_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fai_0 = 1.0 / a_exp;

            tg_x_x_f_0_0_0[i] += tg_0_x_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_y_f_0_0_0[i] += tg_0_y_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_x_z_f_0_0_0[i] += tg_0_z_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_y_x_f_0_0_0[i] += tg_0_x_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_y_f_0_0_0[i] += tg_0_y_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_y_z_f_0_0_0[i] += tg_0_z_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_z_x_f_0_0_0[i] += tg_0_x_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_y_f_0_0_0[i] += tg_0_y_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_z_z_f_0_0_0[i] += tg_0_z_f_0_0_1[i] * fai_0 * a_z * fp_0;
        }
    }
}

} // t2pecp namespace

