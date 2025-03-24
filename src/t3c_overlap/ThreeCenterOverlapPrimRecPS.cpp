#include "ThreeCenterOverlapPrimRecPS.hpp"

namespace t3ovlrec { // t3ovlrec namespace

auto
comp_prim_overlap_ps(CSimdArray<double>& pbuffer, 
                     const size_t idx_ps,
                     const size_t idx_ss,
                     const CSimdArray<double>& factors,
                     const size_t idx_rga) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(GA) distances

    auto ga_x = factors.data(idx_rga);

    auto ga_y = factors.data(idx_rga + 1);

    auto ga_z = factors.data(idx_rga + 2);

    // Set up components of auxiliary buffer : SS

    auto ts_0_0 = pbuffer.data(idx_ss);

    // Set up components of targeted buffer : PS

    auto ts_x_0 = pbuffer.data(idx_ps);

    auto ts_y_0 = pbuffer.data(idx_ps + 1);

    auto ts_z_0 = pbuffer.data(idx_ps + 2);

    #pragma omp simd aligned(ga_x, ga_y, ga_z, ts_0_0, ts_x_0, ts_y_0, ts_z_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ts_x_0[i] = ts_0_0[i] * ga_x[i];

        ts_y_0[i] = ts_0_0[i] * ga_y[i];

        ts_z_0[i] = ts_0_0[i] * ga_z[i];
    }
}

} // t3ovlrec namespace

