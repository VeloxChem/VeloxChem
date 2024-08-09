#include "GeometricalDerivatives1X1ForSS.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_11_ss(CSimdArray<double>& pbuffer,
                        const size_t idx_op_geom_101_ss,
                        const size_t idx_op_pp,
                        const size_t op_comps,
                        const CSimdArray<double>& factors,
                        const double a_exp) -> void
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

        // Set up components of targeted buffer : SS

        auto to_x_x_0_0 = pbuffer.data(idx_op_geom_101_ss + 0 * op_comps * 1 + i * 1 + 0);

        auto to_x_y_0_0 = pbuffer.data(idx_op_geom_101_ss + 1 * op_comps * 1 + i * 1 + 0);

        auto to_x_z_0_0 = pbuffer.data(idx_op_geom_101_ss + 2 * op_comps * 1 + i * 1 + 0);

        auto to_y_x_0_0 = pbuffer.data(idx_op_geom_101_ss + 3 * op_comps * 1 + i * 1 + 0);

        auto to_y_y_0_0 = pbuffer.data(idx_op_geom_101_ss + 4 * op_comps * 1 + i * 1 + 0);

        auto to_y_z_0_0 = pbuffer.data(idx_op_geom_101_ss + 5 * op_comps * 1 + i * 1 + 0);

        auto to_z_x_0_0 = pbuffer.data(idx_op_geom_101_ss + 6 * op_comps * 1 + i * 1 + 0);

        auto to_z_y_0_0 = pbuffer.data(idx_op_geom_101_ss + 7 * op_comps * 1 + i * 1 + 0);

        auto to_z_z_0_0 = pbuffer.data(idx_op_geom_101_ss + 8 * op_comps * 1 + i * 1 + 0);

        #pragma omp simd aligned(to_x_x, to_x_x_0_0, to_x_y, to_x_y_0_0, to_x_z, to_x_z_0_0, to_y_x, to_y_x_0_0, to_y_y, to_y_y_0_0, to_y_z, to_y_z_0_0, to_z_x, to_z_x_0_0, to_z_y, to_z_y_0_0, to_z_z, to_z_z_0_0, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_0_0[k] = 4.0 * to_x_x[k] * tbe_0 * tke_0;

            to_x_y_0_0[k] = 4.0 * to_x_y[k] * tbe_0 * tke_0;

            to_x_z_0_0[k] = 4.0 * to_x_z[k] * tbe_0 * tke_0;

            to_y_x_0_0[k] = 4.0 * to_y_x[k] * tbe_0 * tke_0;

            to_y_y_0_0[k] = 4.0 * to_y_y[k] * tbe_0 * tke_0;

            to_y_z_0_0[k] = 4.0 * to_y_z[k] * tbe_0 * tke_0;

            to_z_x_0_0[k] = 4.0 * to_z_x[k] * tbe_0 * tke_0;

            to_z_y_0_0[k] = 4.0 * to_z_y[k] * tbe_0 * tke_0;

            to_z_z_0_0[k] = 4.0 * to_z_z[k] * tbe_0 * tke_0;
        }
    }

}

} // t2cgeom namespace

