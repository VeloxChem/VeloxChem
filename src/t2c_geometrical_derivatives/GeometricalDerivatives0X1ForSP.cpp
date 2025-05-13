#include "GeometricalDerivatives0X1ForSP.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_geom_deriv_0x1_sp(CSimdArray<double>& pbuffer,
                       const int idx_op_geom_001_sp,
                       const int idx_op_ss,
                       const int idx_op_sd,
                       const CSimdArray<double>& factors,
                       const int op_comps) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : SS

        auto to_0_0 = pbuffer.data(idx_op_ss + i * 1 + 0);

        // Set up components of auxiliary buffer : SD

        auto to_0_xx = pbuffer.data(idx_op_sd + i * 6 + 0);

        auto to_0_xy = pbuffer.data(idx_op_sd + i * 6 + 1);

        auto to_0_xz = pbuffer.data(idx_op_sd + i * 6 + 2);

        auto to_0_yy = pbuffer.data(idx_op_sd + i * 6 + 3);

        auto to_0_yz = pbuffer.data(idx_op_sd + i * 6 + 4);

        auto to_0_zz = pbuffer.data(idx_op_sd + i * 6 + 5);

        // Set up 0-3 components of targeted buffer : SP

        auto to_0_x_0_x = pbuffer.data(idx_op_geom_001_sp + 0 * op_comps * 3 + i * 3 + 0);

        auto to_0_x_0_y = pbuffer.data(idx_op_geom_001_sp + 0 * op_comps * 3 + i * 3 + 1);

        auto to_0_x_0_z = pbuffer.data(idx_op_geom_001_sp + 0 * op_comps * 3 + i * 3 + 2);

        #pragma omp simd aligned(to_0_0, to_0_x_0_x, to_0_x_0_y, to_0_x_0_z, to_0_xx, to_0_xy, to_0_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_0_x[k] = -to_0_0[k] + 2.0 * to_0_xx[k] * tke_0;

            to_0_x_0_y[k] = 2.0 * to_0_xy[k] * tke_0;

            to_0_x_0_z[k] = 2.0 * to_0_xz[k] * tke_0;
        }

        // Set up 3-6 components of targeted buffer : SP

        auto to_0_y_0_x = pbuffer.data(idx_op_geom_001_sp + 1 * op_comps * 3 + i * 3 + 0);

        auto to_0_y_0_y = pbuffer.data(idx_op_geom_001_sp + 1 * op_comps * 3 + i * 3 + 1);

        auto to_0_y_0_z = pbuffer.data(idx_op_geom_001_sp + 1 * op_comps * 3 + i * 3 + 2);

        #pragma omp simd aligned(to_0_0, to_0_xy, to_0_y_0_x, to_0_y_0_y, to_0_y_0_z, to_0_yy, to_0_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_0_x[k] = 2.0 * to_0_xy[k] * tke_0;

            to_0_y_0_y[k] = -to_0_0[k] + 2.0 * to_0_yy[k] * tke_0;

            to_0_y_0_z[k] = 2.0 * to_0_yz[k] * tke_0;
        }

        // Set up 6-9 components of targeted buffer : SP

        auto to_0_z_0_x = pbuffer.data(idx_op_geom_001_sp + 2 * op_comps * 3 + i * 3 + 0);

        auto to_0_z_0_y = pbuffer.data(idx_op_geom_001_sp + 2 * op_comps * 3 + i * 3 + 1);

        auto to_0_z_0_z = pbuffer.data(idx_op_geom_001_sp + 2 * op_comps * 3 + i * 3 + 2);

        #pragma omp simd aligned(to_0_0, to_0_xz, to_0_yz, to_0_z_0_x, to_0_z_0_y, to_0_z_0_z, to_0_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_0_x[k] = 2.0 * to_0_xz[k] * tke_0;

            to_0_z_0_y[k] = 2.0 * to_0_yz[k] * tke_0;

            to_0_z_0_z[k] = -to_0_0[k] + 2.0 * to_0_zz[k] * tke_0;
        }

    }

}

} // t2cgeom namespace

