#include "ThreeCenterOverlapPrimRecDS.hpp"

namespace t3ovlrec { // t3ovlrec namespace

auto
comp_prim_overlap_ds(CSimdArray<double>& pbuffer, 
                     const size_t idx_ds,
                     const size_t idx_ss,
                     const size_t idx_ps,
                     const CSimdArray<double>& factors,
                     const size_t idx_rga,
                     const double a_exp,
                     const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(GA) distances

    auto ga_x = factors.data(idx_rga);

    auto ga_y = factors.data(idx_rga + 1);

    auto ga_z = factors.data(idx_rga + 2);

    // Set up components of auxiliary buffer : SS

    auto ts_0_0 = pbuffer.data(idx_ss);

    // Set up components of auxiliary buffer : PS

    auto ts_x_0 = pbuffer.data(idx_ps);

    auto ts_y_0 = pbuffer.data(idx_ps + 1);

    auto ts_z_0 = pbuffer.data(idx_ps + 2);

    // Set up components of targeted buffer : DS

    auto ts_xx_0 = pbuffer.data(idx_ds);

    auto ts_xy_0 = pbuffer.data(idx_ds + 1);

    auto ts_xz_0 = pbuffer.data(idx_ds + 2);

    auto ts_yy_0 = pbuffer.data(idx_ds + 3);

    auto ts_yz_0 = pbuffer.data(idx_ds + 4);

    auto ts_zz_0 = pbuffer.data(idx_ds + 5);

    #pragma omp simd aligned(ga_x, ga_y, ga_z, ts_0_0, ts_x_0, ts_xx_0, ts_xy_0, ts_xz_0, ts_y_0, ts_yy_0, ts_yz_0, ts_z_0, ts_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xx_0[i] = ts_0_0[i] * gfe_0 + ts_x_0[i] * ga_x[i];

        ts_xy_0[i] = ts_y_0[i] * ga_x[i];

        ts_xz_0[i] = ts_z_0[i] * ga_x[i];

        ts_yy_0[i] = ts_0_0[i] * gfe_0 + ts_y_0[i] * ga_y[i];

        ts_yz_0[i] = ts_z_0[i] * ga_y[i];

        ts_zz_0[i] = ts_0_0[i] * gfe_0 + ts_z_0[i] * ga_z[i];
    }
}

} // t3ovlrec namespace

