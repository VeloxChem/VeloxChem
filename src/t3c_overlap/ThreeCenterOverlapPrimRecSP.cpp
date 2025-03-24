#include "ThreeCenterOverlapPrimRecSP.hpp"

namespace t3ovlrec { // t3ovlrec namespace

auto
comp_prim_overlap_sp(CSimdArray<double>& pbuffer, 
                     const size_t idx_sp,
                     const size_t idx_ss,
                     const CSimdArray<double>& factors,
                     const size_t idx_rgb) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(GB) distances

    auto gb_x = factors.data(idx_rgb);

    auto gb_y = factors.data(idx_rgb + 1);

    auto gb_z = factors.data(idx_rgb + 2);

    // Set up components of auxiliary buffer : SS

    auto ts_0_0 = pbuffer.data(idx_ss);

    // Set up components of targeted buffer : SP

    auto ts_0_x = pbuffer.data(idx_sp);

    auto ts_0_y = pbuffer.data(idx_sp + 1);

    auto ts_0_z = pbuffer.data(idx_sp + 2);

    #pragma omp simd aligned(gb_x, gb_y, gb_z, ts_0_0, ts_0_x, ts_0_y, ts_0_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ts_0_x[i] = ts_0_0[i] * gb_x[i];

        ts_0_y[i] = ts_0_0[i] * gb_y[i];

        ts_0_z[i] = ts_0_0[i] * gb_z[i];
    }
}

} // t3ovlrec namespace

