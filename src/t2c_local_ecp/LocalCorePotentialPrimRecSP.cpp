#include "LocalCorePotentialPrimRecSP.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_sp(CSimdArray<double>& pbuffer, 
                                  const size_t idx_sp,
                                  const size_t idx_ss,
                                  const CSimdArray<double>& factors) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(RB) distances

    auto rb_x = factors.data(8);

    auto rb_y = factors.data(9);

    auto rb_z = factors.data(10);

    // Set up components of auxiliary buffer : SS

    auto tg_0_0 = pbuffer.data(idx_ss);

    // Set up components of targeted buffer : SP

    auto tg_0_x = pbuffer.data(idx_sp);

    auto tg_0_y = pbuffer.data(idx_sp + 1);

    auto tg_0_z = pbuffer.data(idx_sp + 2);

    #pragma omp simd aligned(rb_x, rb_y, rb_z, tg_0_0, tg_0_x, tg_0_y, tg_0_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_0_x[i] = tg_0_0[i] * rb_x[i];

        tg_0_y[i] = tg_0_0[i] * rb_y[i];

        tg_0_z[i] = tg_0_0[i] * rb_z[i];
    }
}

} // t2lecp namespace

