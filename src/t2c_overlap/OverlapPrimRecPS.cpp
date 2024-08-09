#include "OverlapPrimRecPS.hpp"

namespace ovlrec {  // ovlrec namespace

auto
comp_prim_overlap_ps(CSimdArray<double>&       pbuffer,
                     const size_t              idx_ovl_ps,
                     const size_t              idx_ovl_ss,
                     const CSimdArray<double>& factors,
                     const size_t              idx_rpa) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : SS

    auto ts_0_0 = pbuffer.data(idx_ovl_ss);

    // Set up components of targeted buffer : PS

    auto ts_x_0 = pbuffer.data(idx_ovl_ps);

    auto ts_y_0 = pbuffer.data(idx_ovl_ps + 1);

    auto ts_z_0 = pbuffer.data(idx_ovl_ps + 2);

#pragma omp simd aligned(pa_x, pa_y, pa_z, ts_0_0, ts_x_0, ts_y_0, ts_z_0 : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ts_x_0[i] = ts_0_0[i] * pa_x[i];

        ts_y_0[i] = ts_0_0[i] * pa_y[i];

        ts_z_0[i] = ts_0_0[i] * pa_z[i];
    }
}

}  // namespace ovlrec
