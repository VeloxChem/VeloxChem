#include "GeometricalDerivatives010ForSS.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_010_ss(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_ss,
                         const int idx_op_sp,
                         const int idx_op_ps,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up components of auxiliary buffer : SP

    auto tr_0_x = pbuffer.data(idx_op_sp);

    auto tr_0_y = pbuffer.data(idx_op_sp + 1);

    auto tr_0_z = pbuffer.data(idx_op_sp + 2);

    // Set up components of auxiliary buffer : PS

    auto tr_x_0 = pbuffer.data(idx_op_ps);

    auto tr_y_0 = pbuffer.data(idx_op_ps + 1);

    auto tr_z_0 = pbuffer.data(idx_op_ps + 2);

    // Set up components of targeted buffer : SS

    auto tr_0_0_x_0_0 = pbuffer.data(idx_op_geom_010_ss);

    auto tr_0_0_y_0_0 = pbuffer.data(idx_op_geom_010_ss + 1);

    auto tr_0_0_z_0_0 = pbuffer.data(idx_op_geom_010_ss + 2);

    #pragma omp simd aligned(tr_0_0_x_0_0, tr_0_0_y_0_0, tr_0_0_z_0_0, tr_0_x, tr_0_y, tr_0_z, tr_x_0, tr_y_0, tr_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_0_0[i] = 2.0 * tr_x_0[i] * tbe_0 + 2.0 * tr_0_x[i] * tke_0;

        tr_0_0_y_0_0[i] = 2.0 * tr_y_0[i] * tbe_0 + 2.0 * tr_0_y[i] * tke_0;

        tr_0_0_z_0_0[i] = 2.0 * tr_z_0[i] * tbe_0 + 2.0 * tr_0_z[i] * tke_0;
    }
}

} // t2cgeom namespace

