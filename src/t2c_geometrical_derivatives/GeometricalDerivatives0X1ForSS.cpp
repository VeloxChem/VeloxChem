#include "GeometricalDerivatives0X1ForSS.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_geom_deriv_0x1_ss(CSimdArray<double>& pbuffer,
                       const int idx_op_geom_001_ss,
                       const int idx_op_sp,
                       const CSimdArray<double>& factors,
                       const int op_comps) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : SP

        auto to_0_x = pbuffer.data(idx_op_sp + i * 3 + 0);

        auto to_0_y = pbuffer.data(idx_op_sp + i * 3 + 1);

        auto to_0_z = pbuffer.data(idx_op_sp + i * 3 + 2);

        // Set up components of targeted buffer : SS

        auto to_0_x_0_0 = pbuffer.data(idx_op_geom_001_ss + 0 * op_comps * 1 + i * 1 + 0);

        auto to_0_y_0_0 = pbuffer.data(idx_op_geom_001_ss + 1 * op_comps * 1 + i * 1 + 0);

        auto to_0_z_0_0 = pbuffer.data(idx_op_geom_001_ss + 2 * op_comps * 1 + i * 1 + 0);

        #pragma omp simd aligned(to_0_x, to_0_x_0_0, to_0_y, to_0_y_0_0, to_0_z, to_0_z_0_0, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_0_0[k] = 2.0 * to_0_x[k] * tke_0;

            to_0_y_0_0[k] = 2.0 * to_0_y[k] * tke_0;

            to_0_z_0_0[k] = 2.0 * to_0_z[k] * tke_0;
        }
    }

}

} // t2cgeom namespace

