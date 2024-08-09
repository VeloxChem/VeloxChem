#include "OverlapPrimRecSP.hpp"

namespace ovlrec {  // ovlrec namespace

auto
comp_prim_overlap_sp(CSimdArray<double>&       pbuffer,
                     const size_t              idx_ovl_sp,
                     const size_t              idx_ovl_ss,
                     const CSimdArray<double>& factors,
                     const size_t              idx_rpb) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(PB) distances

    auto pb_x = factors.data(idx_rpb);

    auto pb_y = factors.data(idx_rpb + 1);

    auto pb_z = factors.data(idx_rpb + 2);

    // Set up components of auxiliary buffer : SS

    auto ts_0_0 = pbuffer.data(idx_ovl_ss);

    // Set up components of targeted buffer : SP

    auto ts_0_x = pbuffer.data(idx_ovl_sp);

    auto ts_0_y = pbuffer.data(idx_ovl_sp + 1);

    auto ts_0_z = pbuffer.data(idx_ovl_sp + 2);

#pragma omp simd aligned(pb_x, pb_y, pb_z, ts_0_0, ts_0_x, ts_0_y, ts_0_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ts_0_x[i] = ts_0_0[i] * pb_x[i];

        ts_0_y[i] = ts_0_0[i] * pb_y[i];

        ts_0_z[i] = ts_0_0[i] * pb_z[i];
    }
}

}  // namespace ovlrec
