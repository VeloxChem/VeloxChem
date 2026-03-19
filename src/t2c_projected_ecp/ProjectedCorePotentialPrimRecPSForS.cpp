#include "ProjectedCorePotentialPrimRecPSForS.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_ps_s(CSimdArray<double>& pbuffer, 
                                        const size_t idx_ps_s_0_0_0,
                                        const size_t idx_ss_s_0_0_0,
                                        const size_t idx_ss_s_1_0_0,
                                        const int p,
                                        const size_t idx_ss_s_0_0_1,
                                        const CSimdArray<double>& factors,
                                        const TPoint<double>& r_a,
                                        const double a_exp,
                                        const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents on ket side

    auto b_exps = factors.data(0);

    // Set up A center coordinates

    const auto xyz = r_a.coordinates();

    const auto a_x = xyz[0];

    const auto a_y = xyz[1];

    const auto a_z = xyz[2];

    // Set up components of auxiliary buffer : SS

    auto tg_0_0_s_0_0_0 = pbuffer.data(idx_ss_s_0_0_0);

    // Set up components of auxiliary buffer : SS

    auto tg_0_0_s_1_0_0 = pbuffer.data(idx_ss_s_1_0_0);

    // Set up components of targeted buffer : PS

    auto tg_x_0_s_0_0_0 = pbuffer.data(idx_ps_s_0_0_0);

    auto tg_y_0_s_0_0_0 = pbuffer.data(idx_ps_s_0_0_0 + 1);

    auto tg_z_0_s_0_0_0 = pbuffer.data(idx_ps_s_0_0_0 + 2);

    #pragma omp simd aligned(b_exps, tg_0_0_s_0_0_0, tg_0_0_s_1_0_0, tg_x_0_s_0_0_0, tg_y_0_s_0_0_0, tg_z_0_s_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double faz_0 = -(b_exps[i] + c_exp) * fzi_0;

        const double fbzi_0 = b_exps[i] * fzi_0;

        tg_x_0_s_0_0_0[i] = 2.0 * tg_0_0_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_0_0_s_0_0_0[i] * a_x * faz_0;

        tg_y_0_s_0_0_0[i] = 2.0 * tg_0_0_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_0_0_s_0_0_0[i] * a_y * faz_0;

        tg_z_0_s_0_0_0[i] = 2.0 * tg_0_0_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_0_0_s_0_0_0[i] * a_z * faz_0;
    }

    if (p > 0)
    {
        const double fp_0 = (double)p;

        // Set up components of auxiliary buffer : SS

        auto tg_0_0_s_0_0_1 = pbuffer.data(idx_ss_s_0_0_1);

        #pragma omp simd aligned(b_exps, tg_0_0_s_0_0_1, tg_x_0_s_0_0_0, tg_y_0_s_0_0_0, tg_z_0_s_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fai_0 = 1.0 / a_exp;

            tg_x_0_s_0_0_0[i] += tg_0_0_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_y_0_s_0_0_0[i] += tg_0_0_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_z_0_s_0_0_0[i] += tg_0_0_s_0_0_1[i] * fai_0 * a_z * fp_0;
        }
    }
}

} // t2pecp namespace

