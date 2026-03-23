#include "T2CHrrABRecPP.hpp"

namespace t2chrr { // t2chrr namespace

auto
comp_hrr_pp(CSimdArray<double>& cbuffer, 
            const size_t idx_pp,
            const size_t idx_sp,
            const size_t idx_sd,
            const CSimdArray<double>& factors) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    // Set up R(AB) distances

    auto ab_x = factors.data(3);

    auto ab_y = factors.data(4);

    auto ab_z = factors.data(5);

    // Set up components of auxiliary buffer : SP

    auto t_0_x = cbuffer.data(idx_sp);

    auto t_0_y = cbuffer.data(idx_sp + 1);

    auto t_0_z = cbuffer.data(idx_sp + 2);

    // Set up components of auxiliary buffer : SD

    auto t_0_xx = cbuffer.data(idx_sd);

    auto t_0_xy = cbuffer.data(idx_sd + 1);

    auto t_0_xz = cbuffer.data(idx_sd + 2);

    auto t_0_yy = cbuffer.data(idx_sd + 3);

    auto t_0_yz = cbuffer.data(idx_sd + 4);

    auto t_0_zz = cbuffer.data(idx_sd + 5);

    // Set up components of targeted buffer : PP

    auto t_x_x = cbuffer.data(idx_pp);

    auto t_x_y = cbuffer.data(idx_pp + 1);

    auto t_x_z = cbuffer.data(idx_pp + 2);

    auto t_y_x = cbuffer.data(idx_pp + 3);

    auto t_y_y = cbuffer.data(idx_pp + 4);

    auto t_y_z = cbuffer.data(idx_pp + 5);

    auto t_z_x = cbuffer.data(idx_pp + 6);

    auto t_z_y = cbuffer.data(idx_pp + 7);

    auto t_z_z = cbuffer.data(idx_pp + 8);

    #pragma omp simd aligned(ab_x, ab_y, ab_z, t_0_x, t_0_xx, t_0_xy, t_0_xz, t_0_y, t_0_yy, t_0_yz, t_0_z, t_0_zz, t_x_x, t_x_y, t_x_z, t_y_x, t_y_y, t_y_z, t_z_x, t_z_y, t_z_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        t_x_x[i] = -t_0_x[i] * ab_x[i] + t_0_xx[i];

        t_x_y[i] = -t_0_y[i] * ab_x[i] + t_0_xy[i];

        t_x_z[i] = -t_0_z[i] * ab_x[i] + t_0_xz[i];

        t_y_x[i] = -t_0_x[i] * ab_y[i] + t_0_xy[i];

        t_y_y[i] = -t_0_y[i] * ab_y[i] + t_0_yy[i];

        t_y_z[i] = -t_0_z[i] * ab_y[i] + t_0_yz[i];

        t_z_x[i] = -t_0_x[i] * ab_z[i] + t_0_xz[i];

        t_z_y[i] = -t_0_y[i] * ab_z[i] + t_0_yz[i];

        t_z_z[i] = -t_0_z[i] * ab_z[i] + t_0_zz[i];
    }
}

} // t2chrr namespace

