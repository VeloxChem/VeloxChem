#include "GeometricalDerivatives010ForSP.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_010_sp(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_sp,
                         const int idx_op_ss,
                         const int idx_op_sd,
                         const int idx_op_pp,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up components of auxiliary buffer : SS

    auto tr_0_0 = pbuffer.data(idx_op_ss);

    // Set up components of auxiliary buffer : SD

    auto tr_0_xx = pbuffer.data(idx_op_sd);

    auto tr_0_xy = pbuffer.data(idx_op_sd + 1);

    auto tr_0_xz = pbuffer.data(idx_op_sd + 2);

    auto tr_0_yy = pbuffer.data(idx_op_sd + 3);

    auto tr_0_yz = pbuffer.data(idx_op_sd + 4);

    auto tr_0_zz = pbuffer.data(idx_op_sd + 5);

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

    // Set up components of targeted buffer : SP

    auto tr_0_0_x_0_x = pbuffer.data(idx_op_geom_010_sp);

    auto tr_0_0_x_0_y = pbuffer.data(idx_op_geom_010_sp + 1);

    auto tr_0_0_x_0_z = pbuffer.data(idx_op_geom_010_sp + 2);

    auto tr_0_0_y_0_x = pbuffer.data(idx_op_geom_010_sp + 3);

    auto tr_0_0_y_0_y = pbuffer.data(idx_op_geom_010_sp + 4);

    auto tr_0_0_y_0_z = pbuffer.data(idx_op_geom_010_sp + 5);

    auto tr_0_0_z_0_x = pbuffer.data(idx_op_geom_010_sp + 6);

    auto tr_0_0_z_0_y = pbuffer.data(idx_op_geom_010_sp + 7);

    auto tr_0_0_z_0_z = pbuffer.data(idx_op_geom_010_sp + 8);

    #pragma omp simd aligned(tr_0_0, tr_0_0_x_0_x, tr_0_0_x_0_y, tr_0_0_x_0_z, tr_0_0_y_0_x, tr_0_0_y_0_y, tr_0_0_y_0_z, tr_0_0_z_0_x, tr_0_0_z_0_y, tr_0_0_z_0_z, tr_0_xx, tr_0_xy, tr_0_xz, tr_0_yy, tr_0_yz, tr_0_zz, tr_x_x, tr_x_y, tr_x_z, tr_y_x, tr_y_y, tr_y_z, tr_z_x, tr_z_y, tr_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_0_x[i] = 2.0 * tr_x_x[i] * tbe_0 + 2.0 * tr_0_xx[i] * tke_0 - tr_0_0[i];

        tr_0_0_x_0_y[i] = 2.0 * tr_x_y[i] * tbe_0 + 2.0 * tr_0_xy[i] * tke_0;

        tr_0_0_x_0_z[i] = 2.0 * tr_x_z[i] * tbe_0 + 2.0 * tr_0_xz[i] * tke_0;

        tr_0_0_y_0_x[i] = 2.0 * tr_y_x[i] * tbe_0 + 2.0 * tr_0_xy[i] * tke_0;

        tr_0_0_y_0_y[i] = 2.0 * tr_y_y[i] * tbe_0 + 2.0 * tr_0_yy[i] * tke_0 - tr_0_0[i];

        tr_0_0_y_0_z[i] = 2.0 * tr_y_z[i] * tbe_0 + 2.0 * tr_0_yz[i] * tke_0;

        tr_0_0_z_0_x[i] = 2.0 * tr_z_x[i] * tbe_0 + 2.0 * tr_0_xz[i] * tke_0;

        tr_0_0_z_0_y[i] = 2.0 * tr_z_y[i] * tbe_0 + 2.0 * tr_0_yz[i] * tke_0;

        tr_0_0_z_0_z[i] = 2.0 * tr_z_z[i] * tbe_0 + 2.0 * tr_0_zz[i] * tke_0 - tr_0_0[i];
    }
}

} // t2cgeom namespace

