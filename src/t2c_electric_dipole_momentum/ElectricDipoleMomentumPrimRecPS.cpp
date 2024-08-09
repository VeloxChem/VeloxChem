#include "ElectricDipoleMomentumPrimRecPS.hpp"

namespace diprec { // diprec namespace

auto
comp_prim_electric_dipole_momentum_ps(CSimdArray<double>& pbuffer, 
                                      const size_t idx_dip_ps,
                                      const size_t idx_ovl_ss,
                                      const size_t idx_dip_ss,
                                      const CSimdArray<double>& factors,
                                      const size_t idx_rpa,
                                      const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : SS

    auto ts_0_0 = pbuffer.data(idx_ovl_ss);

    // Set up components of auxiliary buffer : SS

    auto tr_x_0_0 = pbuffer.data(idx_dip_ss);

    auto tr_y_0_0 = pbuffer.data(idx_dip_ss + 1);

    auto tr_z_0_0 = pbuffer.data(idx_dip_ss + 2);

    // Set up components of targeted buffer : PS

    auto tr_x_x_0 = pbuffer.data(idx_dip_ps);

    auto tr_x_y_0 = pbuffer.data(idx_dip_ps + 1);

    auto tr_x_z_0 = pbuffer.data(idx_dip_ps + 2);

    auto tr_y_x_0 = pbuffer.data(idx_dip_ps + 3);

    auto tr_y_y_0 = pbuffer.data(idx_dip_ps + 4);

    auto tr_y_z_0 = pbuffer.data(idx_dip_ps + 5);

    auto tr_z_x_0 = pbuffer.data(idx_dip_ps + 6);

    auto tr_z_y_0 = pbuffer.data(idx_dip_ps + 7);

    auto tr_z_z_0 = pbuffer.data(idx_dip_ps + 8);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_x_0_0, tr_x_x_0, tr_x_y_0, tr_x_z_0, tr_y_0_0, tr_y_x_0, tr_y_y_0, tr_y_z_0, tr_z_0_0, tr_z_x_0, tr_z_y_0, tr_z_z_0, ts_0_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_x_0[i] = ts_0_0[i] * fe_0 + tr_x_0_0[i] * pa_x[i];

        tr_x_y_0[i] = tr_x_0_0[i] * pa_y[i];

        tr_x_z_0[i] = tr_x_0_0[i] * pa_z[i];

        tr_y_x_0[i] = tr_y_0_0[i] * pa_x[i];

        tr_y_y_0[i] = ts_0_0[i] * fe_0 + tr_y_0_0[i] * pa_y[i];

        tr_y_z_0[i] = tr_y_0_0[i] * pa_z[i];

        tr_z_x_0[i] = tr_z_0_0[i] * pa_x[i];

        tr_z_y_0[i] = tr_z_0_0[i] * pa_y[i];

        tr_z_z_0[i] = ts_0_0[i] * fe_0 + tr_z_0_0[i] * pa_z[i];
    }
}

} // diprec namespace

