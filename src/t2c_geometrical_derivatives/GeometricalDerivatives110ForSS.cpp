#include "GeometricalDerivatives110ForSS.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_110_ss(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_110_ss,
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

    // Set up components of targeted buffer : SS

    auto tr_x_0_x_0_0 = pbuffer.data(idx_op_geom_110_ss);

    auto tr_x_0_y_0_0 = pbuffer.data(idx_op_geom_110_ss + 1);

    auto tr_x_0_z_0_0 = pbuffer.data(idx_op_geom_110_ss + 2);

    auto tr_y_0_x_0_0 = pbuffer.data(idx_op_geom_110_ss + 3);

    auto tr_y_0_y_0_0 = pbuffer.data(idx_op_geom_110_ss + 4);

    auto tr_y_0_z_0_0 = pbuffer.data(idx_op_geom_110_ss + 5);

    auto tr_z_0_x_0_0 = pbuffer.data(idx_op_geom_110_ss + 6);

    auto tr_z_0_y_0_0 = pbuffer.data(idx_op_geom_110_ss + 7);

    auto tr_z_0_z_0_0 = pbuffer.data(idx_op_geom_110_ss + 8);

    #pragma omp simd aligned(tr_0_0, tr_x_0_x_0_0, tr_x_0_y_0_0, tr_x_0_z_0_0, tr_x_x, tr_x_y, tr_x_z, tr_xx_0, tr_xy_0, tr_xz_0, tr_y_0_x_0_0, tr_y_0_y_0_0, tr_y_0_z_0_0, tr_y_x, tr_y_y, tr_y_z, tr_yy_0, tr_yz_0, tr_z_0_x_0_0, tr_z_0_y_0_0, tr_z_0_z_0_0, tr_z_x, tr_z_y, tr_z_z, tr_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_0_0[i] = -2.0 * tr_0_0[i] * tbe_0 + 4.0 * tr_x_x[i] * tbe_0 * tke_0 + 4.0 * tr_xx_0[i] * tbe_0 * tbe_0;

        tr_x_0_y_0_0[i] = 4.0 * tr_x_y[i] * tbe_0 * tke_0 + 4.0 * tr_xy_0[i] * tbe_0 * tbe_0;

        tr_x_0_z_0_0[i] = 4.0 * tr_x_z[i] * tbe_0 * tke_0 + 4.0 * tr_xz_0[i] * tbe_0 * tbe_0;

        tr_y_0_x_0_0[i] = 4.0 * tr_y_x[i] * tbe_0 * tke_0 + 4.0 * tr_xy_0[i] * tbe_0 * tbe_0;

        tr_y_0_y_0_0[i] = -2.0 * tr_0_0[i] * tbe_0 + 4.0 * tr_y_y[i] * tbe_0 * tke_0 + 4.0 * tr_yy_0[i] * tbe_0 * tbe_0;

        tr_y_0_z_0_0[i] = 4.0 * tr_y_z[i] * tbe_0 * tke_0 + 4.0 * tr_yz_0[i] * tbe_0 * tbe_0;

        tr_z_0_x_0_0[i] = 4.0 * tr_z_x[i] * tbe_0 * tke_0 + 4.0 * tr_xz_0[i] * tbe_0 * tbe_0;

        tr_z_0_y_0_0[i] = 4.0 * tr_z_y[i] * tbe_0 * tke_0 + 4.0 * tr_yz_0[i] * tbe_0 * tbe_0;

        tr_z_0_z_0_0[i] = -2.0 * tr_0_0[i] * tbe_0 + 4.0 * tr_z_z[i] * tbe_0 * tke_0 + 4.0 * tr_zz_0[i] * tbe_0 * tbe_0;
    }
}

} // t2cgeom namespace

