#include "LocalCorePotentialPrimRecPP.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_pp(CSimdArray<double>& pbuffer, 
                                  const size_t idx_pp,
                                  const size_t idx_ss,
                                  const size_t idx_sp,
                                  const CSimdArray<double>& factors,
                                  const size_t idx_ra,
                                  const size_t idx_zeta) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(RA) distances

    auto ra_x = factors.data(idx_ra);

    auto ra_y = factors.data(idx_ra + 1);

    auto ra_z = factors.data(idx_ra + 2);

    // Set up inverted 1/2xi

    auto fxi = factors.data(idx_zeta);

    // Set up components of auxiliary buffer : SS

    auto tg_0_0 = pbuffer.data(idx_ss);

    // Set up components of auxiliary buffer : SP

    auto tg_0_x = pbuffer.data(idx_sp);

    auto tg_0_y = pbuffer.data(idx_sp + 1);

    auto tg_0_z = pbuffer.data(idx_sp + 2);

    // Set up components of targeted buffer : PP

    auto tg_x_x = pbuffer.data(idx_pp);

    auto tg_x_y = pbuffer.data(idx_pp + 1);

    auto tg_x_z = pbuffer.data(idx_pp + 2);

    auto tg_y_x = pbuffer.data(idx_pp + 3);

    auto tg_y_y = pbuffer.data(idx_pp + 4);

    auto tg_y_z = pbuffer.data(idx_pp + 5);

    auto tg_z_x = pbuffer.data(idx_pp + 6);

    auto tg_z_y = pbuffer.data(idx_pp + 7);

    auto tg_z_z = pbuffer.data(idx_pp + 8);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_0_0, tg_0_x, tg_0_y, tg_0_z, tg_x_x, tg_x_y, tg_x_z, tg_y_x, tg_y_y, tg_y_z, tg_z_x, tg_z_y, tg_z_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_x_x[i] = tg_0_0[i] * fxi[i] + tg_0_x[i] * ra_x[i];

        tg_x_y[i] = tg_0_y[i] * ra_x[i];

        tg_x_z[i] = tg_0_z[i] * ra_x[i];

        tg_y_x[i] = tg_0_x[i] * ra_y[i];

        tg_y_y[i] = tg_0_0[i] * fxi[i] + tg_0_y[i] * ra_y[i];

        tg_y_z[i] = tg_0_z[i] * ra_y[i];

        tg_z_x[i] = tg_0_x[i] * ra_z[i];

        tg_z_y[i] = tg_0_y[i] * ra_z[i];

        tg_z_z[i] = tg_0_0[i] * fxi[i] + tg_0_z[i] * ra_z[i];
    }
}

} // t2lecp namespace

