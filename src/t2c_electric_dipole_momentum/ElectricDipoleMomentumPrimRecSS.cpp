#include "ElectricDipoleMomentumPrimRecSS.hpp"


namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_ss(CSimdArray<double>& pbuffer,
                                      const size_t idx_dip_ss,
                                      const size_t idx_ovl_ss,
                                      const CSimdArray<double>& factors,
                                      const size_t idx_rpc) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();
    
    // Set up R(PC) distances

    auto pc_x = factors.data(idx_rpc);

    auto pc_y = factors.data(idx_rpc + 1);

    auto pc_z = factors.data(idx_rpc + 2);

    /// Set up components of auxiliary buffer : SS

    auto ts_0_0 = pbuffer.data(idx_ovl_ss);

    /// Set up components of targeted buffer : SS

    auto tr_x_0_0 = pbuffer.data(idx_dip_ss);

    auto tr_y_0_0 = pbuffer.data(idx_dip_ss + 1);

    auto tr_z_0_0 = pbuffer.data(idx_dip_ss + 2);

#pragma omp simd aligned(tr_x_0_0, tr_y_0_0, tr_z_0_0, ts_0_0, pc_x, pc_y, pc_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fact = ts_0_0[i];

        tr_x_0_0[i] = fact * pc_x[i];

        tr_y_0_0[i] = fact * pc_y[i];

        tr_z_0_0[i] = fact * pc_z[i];
    }
}

}  // namespace diprec
