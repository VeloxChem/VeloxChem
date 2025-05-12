#include "GeometricalDerivatives0X1ForDS.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_geom_deriv_0x1_ds(CSimdArray<double>& pbuffer,
                       const int idx_op_geom_001_ds,
                       const int idx_op_dp,
                       const CSimdArray<double>& factors,
                       const int op_comps) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : DP

        auto to_xx_x = pbuffer.data(idx_op_dp + i * 18 + 0);

        auto to_xx_y = pbuffer.data(idx_op_dp + i * 18 + 1);

        auto to_xx_z = pbuffer.data(idx_op_dp + i * 18 + 2);

        auto to_xy_x = pbuffer.data(idx_op_dp + i * 18 + 3);

        auto to_xy_y = pbuffer.data(idx_op_dp + i * 18 + 4);

        auto to_xy_z = pbuffer.data(idx_op_dp + i * 18 + 5);

        auto to_xz_x = pbuffer.data(idx_op_dp + i * 18 + 6);

        auto to_xz_y = pbuffer.data(idx_op_dp + i * 18 + 7);

        auto to_xz_z = pbuffer.data(idx_op_dp + i * 18 + 8);

        auto to_yy_x = pbuffer.data(idx_op_dp + i * 18 + 9);

        auto to_yy_y = pbuffer.data(idx_op_dp + i * 18 + 10);

        auto to_yy_z = pbuffer.data(idx_op_dp + i * 18 + 11);

        auto to_yz_x = pbuffer.data(idx_op_dp + i * 18 + 12);

        auto to_yz_y = pbuffer.data(idx_op_dp + i * 18 + 13);

        auto to_yz_z = pbuffer.data(idx_op_dp + i * 18 + 14);

        auto to_zz_x = pbuffer.data(idx_op_dp + i * 18 + 15);

        auto to_zz_y = pbuffer.data(idx_op_dp + i * 18 + 16);

        auto to_zz_z = pbuffer.data(idx_op_dp + i * 18 + 17);

        // Set up 0-6 components of targeted buffer : DS

        auto to_0_x_xx_0 = pbuffer.data(idx_op_geom_001_ds + 0 * op_comps * 6 + i * 6 + 0);

        auto to_0_x_xy_0 = pbuffer.data(idx_op_geom_001_ds + 0 * op_comps * 6 + i * 6 + 1);

        auto to_0_x_xz_0 = pbuffer.data(idx_op_geom_001_ds + 0 * op_comps * 6 + i * 6 + 2);

        auto to_0_x_yy_0 = pbuffer.data(idx_op_geom_001_ds + 0 * op_comps * 6 + i * 6 + 3);

        auto to_0_x_yz_0 = pbuffer.data(idx_op_geom_001_ds + 0 * op_comps * 6 + i * 6 + 4);

        auto to_0_x_zz_0 = pbuffer.data(idx_op_geom_001_ds + 0 * op_comps * 6 + i * 6 + 5);

        #pragma omp simd aligned(to_0_x_xx_0, to_0_x_xy_0, to_0_x_xz_0, to_0_x_yy_0, to_0_x_yz_0, to_0_x_zz_0, to_xx_x, to_xy_x, to_xz_x, to_yy_x, to_yz_x, to_zz_x, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xx_0[k] = 2.0 * to_xx_x[k] * tke_0;

            to_0_x_xy_0[k] = 2.0 * to_xy_x[k] * tke_0;

            to_0_x_xz_0[k] = 2.0 * to_xz_x[k] * tke_0;

            to_0_x_yy_0[k] = 2.0 * to_yy_x[k] * tke_0;

            to_0_x_yz_0[k] = 2.0 * to_yz_x[k] * tke_0;

            to_0_x_zz_0[k] = 2.0 * to_zz_x[k] * tke_0;
        }

        // Set up 6-12 components of targeted buffer : DS

        auto to_0_y_xx_0 = pbuffer.data(idx_op_geom_001_ds + 1 * op_comps * 6 + i * 6 + 0);

        auto to_0_y_xy_0 = pbuffer.data(idx_op_geom_001_ds + 1 * op_comps * 6 + i * 6 + 1);

        auto to_0_y_xz_0 = pbuffer.data(idx_op_geom_001_ds + 1 * op_comps * 6 + i * 6 + 2);

        auto to_0_y_yy_0 = pbuffer.data(idx_op_geom_001_ds + 1 * op_comps * 6 + i * 6 + 3);

        auto to_0_y_yz_0 = pbuffer.data(idx_op_geom_001_ds + 1 * op_comps * 6 + i * 6 + 4);

        auto to_0_y_zz_0 = pbuffer.data(idx_op_geom_001_ds + 1 * op_comps * 6 + i * 6 + 5);

        #pragma omp simd aligned(to_0_y_xx_0, to_0_y_xy_0, to_0_y_xz_0, to_0_y_yy_0, to_0_y_yz_0, to_0_y_zz_0, to_xx_y, to_xy_y, to_xz_y, to_yy_y, to_yz_y, to_zz_y, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xx_0[k] = 2.0 * to_xx_y[k] * tke_0;

            to_0_y_xy_0[k] = 2.0 * to_xy_y[k] * tke_0;

            to_0_y_xz_0[k] = 2.0 * to_xz_y[k] * tke_0;

            to_0_y_yy_0[k] = 2.0 * to_yy_y[k] * tke_0;

            to_0_y_yz_0[k] = 2.0 * to_yz_y[k] * tke_0;

            to_0_y_zz_0[k] = 2.0 * to_zz_y[k] * tke_0;
        }

        // Set up 12-18 components of targeted buffer : DS

        auto to_0_z_xx_0 = pbuffer.data(idx_op_geom_001_ds + 2 * op_comps * 6 + i * 6 + 0);

        auto to_0_z_xy_0 = pbuffer.data(idx_op_geom_001_ds + 2 * op_comps * 6 + i * 6 + 1);

        auto to_0_z_xz_0 = pbuffer.data(idx_op_geom_001_ds + 2 * op_comps * 6 + i * 6 + 2);

        auto to_0_z_yy_0 = pbuffer.data(idx_op_geom_001_ds + 2 * op_comps * 6 + i * 6 + 3);

        auto to_0_z_yz_0 = pbuffer.data(idx_op_geom_001_ds + 2 * op_comps * 6 + i * 6 + 4);

        auto to_0_z_zz_0 = pbuffer.data(idx_op_geom_001_ds + 2 * op_comps * 6 + i * 6 + 5);

        #pragma omp simd aligned(to_0_z_xx_0, to_0_z_xy_0, to_0_z_xz_0, to_0_z_yy_0, to_0_z_yz_0, to_0_z_zz_0, to_xx_z, to_xy_z, to_xz_z, to_yy_z, to_yz_z, to_zz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xx_0[k] = 2.0 * to_xx_z[k] * tke_0;

            to_0_z_xy_0[k] = 2.0 * to_xy_z[k] * tke_0;

            to_0_z_xz_0[k] = 2.0 * to_xz_z[k] * tke_0;

            to_0_z_yy_0[k] = 2.0 * to_yy_z[k] * tke_0;

            to_0_z_yz_0[k] = 2.0 * to_yz_z[k] * tke_0;

            to_0_z_zz_0[k] = 2.0 * to_zz_z[k] * tke_0;
        }

    }

}

} // t2cgeom namespace

