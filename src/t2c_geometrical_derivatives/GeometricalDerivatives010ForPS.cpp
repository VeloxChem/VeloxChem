#include "GeometricalDerivatives010ForPS.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_010_ps(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_ps,
                         const int idx_op_ss,
                         const int idx_op_pp,
                         const int idx_op_ds,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up components of auxiliary buffer : SS

    auto tr_0_0 = pbuffer.data(idx_op_ss);

    // Set up components of auxiliary buffer : PP

    auto tr_x_x = pbuffer.data(idx_op_pp);

    auto tr_x_y = pbuffer.data(idx_op_pp + 1);

    auto tr_x_z = pbuffer.data(idx_op_pp + 2);

    auto tr_y_x = pbuffer.data(idx_op_pp + 3);

    auto tr_y_y = pbuffer.data(idx_op_pp + 4);

    auto tr_y_z = pbuffer.data(idx_op_pp + 5);

    auto tr_z_x = pbuffer.data(idx_op_pp + 6);

    auto tr_z_y = pbuffer.data(idx_op_pp + 7);

    auto tr_z_z = pbuffer.data(idx_op_pp + 8);

    // Set up components of auxiliary buffer : DS

    auto tr_xx_0 = pbuffer.data(idx_op_ds);

    auto tr_xy_0 = pbuffer.data(idx_op_ds + 1);

    auto tr_xz_0 = pbuffer.data(idx_op_ds + 2);

    auto tr_yy_0 = pbuffer.data(idx_op_ds + 3);

    auto tr_yz_0 = pbuffer.data(idx_op_ds + 4);

    auto tr_zz_0 = pbuffer.data(idx_op_ds + 5);

    // Set up components of targeted buffer : PS

    auto tr_0_0_x_x_0 = pbuffer.data(idx_op_geom_010_ps);

    auto tr_0_0_x_y_0 = pbuffer.data(idx_op_geom_010_ps + 1);

    auto tr_0_0_x_z_0 = pbuffer.data(idx_op_geom_010_ps + 2);

    auto tr_0_0_y_x_0 = pbuffer.data(idx_op_geom_010_ps + 3);

    auto tr_0_0_y_y_0 = pbuffer.data(idx_op_geom_010_ps + 4);

    auto tr_0_0_y_z_0 = pbuffer.data(idx_op_geom_010_ps + 5);

    auto tr_0_0_z_x_0 = pbuffer.data(idx_op_geom_010_ps + 6);

    auto tr_0_0_z_y_0 = pbuffer.data(idx_op_geom_010_ps + 7);

    auto tr_0_0_z_z_0 = pbuffer.data(idx_op_geom_010_ps + 8);

    #pragma omp simd aligned(tr_0_0, tr_0_0_x_x_0, tr_0_0_x_y_0, tr_0_0_x_z_0, tr_0_0_y_x_0, tr_0_0_y_y_0, tr_0_0_y_z_0, tr_0_0_z_x_0, tr_0_0_z_y_0, tr_0_0_z_z_0, tr_x_x, tr_x_y, tr_x_z, tr_xx_0, tr_xy_0, tr_xz_0, tr_y_x, tr_y_y, tr_y_z, tr_yy_0, tr_yz_0, tr_z_x, tr_z_y, tr_z_z, tr_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_x_0[i] = 2.0 * tr_xx_0[i] * tbe_0 + 2.0 * tr_x_x[i] * tke_0 - tr_0_0[i];

        tr_0_0_x_y_0[i] = 2.0 * tr_xy_0[i] * tbe_0 + 2.0 * tr_y_x[i] * tke_0;

        tr_0_0_x_z_0[i] = 2.0 * tr_xz_0[i] * tbe_0 + 2.0 * tr_z_x[i] * tke_0;

        tr_0_0_y_x_0[i] = 2.0 * tr_xy_0[i] * tbe_0 + 2.0 * tr_x_y[i] * tke_0;

        tr_0_0_y_y_0[i] = 2.0 * tr_yy_0[i] * tbe_0 + 2.0 * tr_y_y[i] * tke_0 - tr_0_0[i];

        tr_0_0_y_z_0[i] = 2.0 * tr_yz_0[i] * tbe_0 + 2.0 * tr_z_y[i] * tke_0;

        tr_0_0_z_x_0[i] = 2.0 * tr_xz_0[i] * tbe_0 + 2.0 * tr_x_z[i] * tke_0;

        tr_0_0_z_y_0[i] = 2.0 * tr_yz_0[i] * tbe_0 + 2.0 * tr_y_z[i] * tke_0;

        tr_0_0_z_z_0[i] = 2.0 * tr_zz_0[i] * tbe_0 + 2.0 * tr_z_z[i] * tke_0 - tr_0_0[i];
    }
}

} // t2cgeom namespace

