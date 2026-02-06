#include "LocalCorePotentialPrimRecPS.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_ps(CSimdArray<double>& pbuffer, 
                                  const size_t idx_ps,
                                  const size_t idx_ss,
                                  const CSimdArray<double>& factors) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(RA) distances

    auto ra_x = factors.data(8);

    auto ra_y = factors.data(9);

    auto ra_z = factors.data(10);

    // Set up components of auxiliary buffer : SS

    auto tg_0_0 = pbuffer.data(idx_ss);

    // Set up components of targeted buffer : PS

    auto tg_x_0 = pbuffer.data(idx_ps);

    auto tg_y_0 = pbuffer.data(idx_ps + 1);

    auto tg_z_0 = pbuffer.data(idx_ps + 2);

    #pragma omp simd aligned(ra_x, ra_y, ra_z, tg_0_0, tg_x_0, tg_y_0, tg_z_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_x_0[i] = tg_0_0[i] * ra_x[i];

        tg_y_0[i] = tg_0_0[i] * ra_y[i];

        tg_z_0[i] = tg_0_0[i] * ra_z[i];
    }
}

} // t2lecp namespace

