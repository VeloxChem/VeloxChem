#include "ElectricDipoleMomentumPrimRecSP.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_sp(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_sp,
                                      const size_t              idx_ovl_ss,
                                      const size_t              idx_dip_ss,
                                      const CSimdArray<double>& factors,
                                      const size_t              idx_rpb,
                                      const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PB) distances

    auto pb_x = factors.data(idx_rpb);

    auto pb_y = factors.data(idx_rpb + 1);

    auto pb_z = factors.data(idx_rpb + 2);

    // Set up components of auxiliary buffer : SS

    auto ts_0_0 = pbuffer.data(idx_ovl_ss);

    // Set up components of auxiliary buffer : SS

    auto tr_x_0_0 = pbuffer.data(idx_dip_ss);

    auto tr_y_0_0 = pbuffer.data(idx_dip_ss + 1);

    auto tr_z_0_0 = pbuffer.data(idx_dip_ss + 2);

    // Set up components of targeted buffer : SP

    auto tr_x_0_x = pbuffer.data(idx_dip_sp);

    auto tr_x_0_y = pbuffer.data(idx_dip_sp + 1);

    auto tr_x_0_z = pbuffer.data(idx_dip_sp + 2);

    auto tr_y_0_x = pbuffer.data(idx_dip_sp + 3);

    auto tr_y_0_y = pbuffer.data(idx_dip_sp + 4);

    auto tr_y_0_z = pbuffer.data(idx_dip_sp + 5);

    auto tr_z_0_x = pbuffer.data(idx_dip_sp + 6);

    auto tr_z_0_y = pbuffer.data(idx_dip_sp + 7);

    auto tr_z_0_z = pbuffer.data(idx_dip_sp + 8);

#pragma omp simd aligned(pb_x,         \
                             pb_y,     \
                             pb_z,     \
                             tr_x_0_0, \
                             tr_x_0_x, \
                             tr_x_0_y, \
                             tr_x_0_z, \
                             tr_y_0_0, \
                             tr_y_0_x, \
                             tr_y_0_y, \
                             tr_y_0_z, \
                             tr_z_0_0, \
                             tr_z_0_x, \
                             tr_z_0_y, \
                             tr_z_0_z, \
                             ts_0_0,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_0_x[i] = ts_0_0[i] * fe_0 + tr_x_0_0[i] * pb_x[i];

        tr_x_0_y[i] = tr_x_0_0[i] * pb_y[i];

        tr_x_0_z[i] = tr_x_0_0[i] * pb_z[i];

        tr_y_0_x[i] = tr_y_0_0[i] * pb_x[i];

        tr_y_0_y[i] = ts_0_0[i] * fe_0 + tr_y_0_0[i] * pb_y[i];

        tr_y_0_z[i] = tr_y_0_0[i] * pb_z[i];

        tr_z_0_x[i] = tr_z_0_0[i] * pb_x[i];

        tr_z_0_y[i] = tr_z_0_0[i] * pb_y[i];

        tr_z_0_z[i] = ts_0_0[i] * fe_0 + tr_z_0_0[i] * pb_z[i];
    }
}

}  // namespace diprec
