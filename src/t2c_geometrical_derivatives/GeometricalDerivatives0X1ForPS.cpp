#include "GeometricalDerivatives0X1ForPS.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_geom_deriv_0x1_ps(CSimdArray<double>& pbuffer,
                       const int idx_op_geom_001_ps,
                       const int idx_op_pp,
                       const CSimdArray<double>& factors,
                       const int op_comps) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : PP

        auto to_x_x = pbuffer.data(idx_op_pp + i * 9 + 0);

        auto to_x_y = pbuffer.data(idx_op_pp + i * 9 + 1);

        auto to_x_z = pbuffer.data(idx_op_pp + i * 9 + 2);

        auto to_y_x = pbuffer.data(idx_op_pp + i * 9 + 3);

        auto to_y_y = pbuffer.data(idx_op_pp + i * 9 + 4);

        auto to_y_z = pbuffer.data(idx_op_pp + i * 9 + 5);

        auto to_z_x = pbuffer.data(idx_op_pp + i * 9 + 6);

        auto to_z_y = pbuffer.data(idx_op_pp + i * 9 + 7);

        auto to_z_z = pbuffer.data(idx_op_pp + i * 9 + 8);

        // Set up 0-3 components of targeted buffer : PS

        auto to_0_x_x_0 = pbuffer.data(idx_op_geom_001_ps + 0 * op_comps * 3 + i * 3 + 0);

        auto to_0_x_y_0 = pbuffer.data(idx_op_geom_001_ps + 0 * op_comps * 3 + i * 3 + 1);

        auto to_0_x_z_0 = pbuffer.data(idx_op_geom_001_ps + 0 * op_comps * 3 + i * 3 + 2);

        #pragma omp simd aligned(to_0_x_x_0, to_0_x_y_0, to_0_x_z_0, to_x_x, to_y_x, to_z_x, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_x_0[k] = 2.0 * to_x_x[k] * tke_0;

            to_0_x_y_0[k] = 2.0 * to_y_x[k] * tke_0;

            to_0_x_z_0[k] = 2.0 * to_z_x[k] * tke_0;
        }

        // Set up 3-6 components of targeted buffer : PS

        auto to_0_y_x_0 = pbuffer.data(idx_op_geom_001_ps + 1 * op_comps * 3 + i * 3 + 0);

        auto to_0_y_y_0 = pbuffer.data(idx_op_geom_001_ps + 1 * op_comps * 3 + i * 3 + 1);

        auto to_0_y_z_0 = pbuffer.data(idx_op_geom_001_ps + 1 * op_comps * 3 + i * 3 + 2);

        #pragma omp simd aligned(to_0_y_x_0, to_0_y_y_0, to_0_y_z_0, to_x_y, to_y_y, to_z_y, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_x_0[k] = 2.0 * to_x_y[k] * tke_0;

            to_0_y_y_0[k] = 2.0 * to_y_y[k] * tke_0;

            to_0_y_z_0[k] = 2.0 * to_z_y[k] * tke_0;
        }

        // Set up 6-9 components of targeted buffer : PS

        auto to_0_z_x_0 = pbuffer.data(idx_op_geom_001_ps + 2 * op_comps * 3 + i * 3 + 0);

        auto to_0_z_y_0 = pbuffer.data(idx_op_geom_001_ps + 2 * op_comps * 3 + i * 3 + 1);

        auto to_0_z_z_0 = pbuffer.data(idx_op_geom_001_ps + 2 * op_comps * 3 + i * 3 + 2);

        #pragma omp simd aligned(to_0_z_x_0, to_0_z_y_0, to_0_z_z_0, to_x_z, to_y_z, to_z_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_x_0[k] = 2.0 * to_x_z[k] * tke_0;

            to_0_z_y_0[k] = 2.0 * to_y_z[k] * tke_0;

            to_0_z_z_0[k] = 2.0 * to_z_z[k] * tke_0;
        }

    }

}

} // t2cgeom namespace

