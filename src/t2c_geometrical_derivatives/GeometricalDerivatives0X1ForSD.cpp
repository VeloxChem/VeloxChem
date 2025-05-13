#include "GeometricalDerivatives0X1ForSD.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_geom_deriv_0x1_sd(CSimdArray<double>& pbuffer,
                       const int idx_op_geom_001_sd,
                       const int idx_op_sp,
                       const int idx_op_sf,
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

        // Set up components of auxiliary buffer : SF

        auto to_0_xxx = pbuffer.data(idx_op_sf + i * 10 + 0);

        auto to_0_xxy = pbuffer.data(idx_op_sf + i * 10 + 1);

        auto to_0_xxz = pbuffer.data(idx_op_sf + i * 10 + 2);

        auto to_0_xyy = pbuffer.data(idx_op_sf + i * 10 + 3);

        auto to_0_xyz = pbuffer.data(idx_op_sf + i * 10 + 4);

        auto to_0_xzz = pbuffer.data(idx_op_sf + i * 10 + 5);

        auto to_0_yyy = pbuffer.data(idx_op_sf + i * 10 + 6);

        auto to_0_yyz = pbuffer.data(idx_op_sf + i * 10 + 7);

        auto to_0_yzz = pbuffer.data(idx_op_sf + i * 10 + 8);

        auto to_0_zzz = pbuffer.data(idx_op_sf + i * 10 + 9);

        // Set up 0-6 components of targeted buffer : SD

        auto to_0_x_0_xx = pbuffer.data(idx_op_geom_001_sd + 0 * op_comps * 6 + i * 6 + 0);

        auto to_0_x_0_xy = pbuffer.data(idx_op_geom_001_sd + 0 * op_comps * 6 + i * 6 + 1);

        auto to_0_x_0_xz = pbuffer.data(idx_op_geom_001_sd + 0 * op_comps * 6 + i * 6 + 2);

        auto to_0_x_0_yy = pbuffer.data(idx_op_geom_001_sd + 0 * op_comps * 6 + i * 6 + 3);

        auto to_0_x_0_yz = pbuffer.data(idx_op_geom_001_sd + 0 * op_comps * 6 + i * 6 + 4);

        auto to_0_x_0_zz = pbuffer.data(idx_op_geom_001_sd + 0 * op_comps * 6 + i * 6 + 5);

        #pragma omp simd aligned(to_0_x, to_0_x_0_xx, to_0_x_0_xy, to_0_x_0_xz, to_0_x_0_yy, to_0_x_0_yz, to_0_x_0_zz, to_0_xxx, to_0_xxy, to_0_xxz, to_0_xyy, to_0_xyz, to_0_xzz, to_0_y, to_0_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_0_xx[k] = -2.0 * to_0_x[k] + 2.0 * to_0_xxx[k] * tke_0;

            to_0_x_0_xy[k] = -to_0_y[k] + 2.0 * to_0_xxy[k] * tke_0;

            to_0_x_0_xz[k] = -to_0_z[k] + 2.0 * to_0_xxz[k] * tke_0;

            to_0_x_0_yy[k] = 2.0 * to_0_xyy[k] * tke_0;

            to_0_x_0_yz[k] = 2.0 * to_0_xyz[k] * tke_0;

            to_0_x_0_zz[k] = 2.0 * to_0_xzz[k] * tke_0;
        }

        // Set up 6-12 components of targeted buffer : SD

        auto to_0_y_0_xx = pbuffer.data(idx_op_geom_001_sd + 1 * op_comps * 6 + i * 6 + 0);

        auto to_0_y_0_xy = pbuffer.data(idx_op_geom_001_sd + 1 * op_comps * 6 + i * 6 + 1);

        auto to_0_y_0_xz = pbuffer.data(idx_op_geom_001_sd + 1 * op_comps * 6 + i * 6 + 2);

        auto to_0_y_0_yy = pbuffer.data(idx_op_geom_001_sd + 1 * op_comps * 6 + i * 6 + 3);

        auto to_0_y_0_yz = pbuffer.data(idx_op_geom_001_sd + 1 * op_comps * 6 + i * 6 + 4);

        auto to_0_y_0_zz = pbuffer.data(idx_op_geom_001_sd + 1 * op_comps * 6 + i * 6 + 5);

        #pragma omp simd aligned(to_0_x, to_0_xxy, to_0_xyy, to_0_xyz, to_0_y, to_0_y_0_xx, to_0_y_0_xy, to_0_y_0_xz, to_0_y_0_yy, to_0_y_0_yz, to_0_y_0_zz, to_0_yyy, to_0_yyz, to_0_yzz, to_0_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_0_xx[k] = 2.0 * to_0_xxy[k] * tke_0;

            to_0_y_0_xy[k] = -to_0_x[k] + 2.0 * to_0_xyy[k] * tke_0;

            to_0_y_0_xz[k] = 2.0 * to_0_xyz[k] * tke_0;

            to_0_y_0_yy[k] = -2.0 * to_0_y[k] + 2.0 * to_0_yyy[k] * tke_0;

            to_0_y_0_yz[k] = -to_0_z[k] + 2.0 * to_0_yyz[k] * tke_0;

            to_0_y_0_zz[k] = 2.0 * to_0_yzz[k] * tke_0;
        }

        // Set up 12-18 components of targeted buffer : SD

        auto to_0_z_0_xx = pbuffer.data(idx_op_geom_001_sd + 2 * op_comps * 6 + i * 6 + 0);

        auto to_0_z_0_xy = pbuffer.data(idx_op_geom_001_sd + 2 * op_comps * 6 + i * 6 + 1);

        auto to_0_z_0_xz = pbuffer.data(idx_op_geom_001_sd + 2 * op_comps * 6 + i * 6 + 2);

        auto to_0_z_0_yy = pbuffer.data(idx_op_geom_001_sd + 2 * op_comps * 6 + i * 6 + 3);

        auto to_0_z_0_yz = pbuffer.data(idx_op_geom_001_sd + 2 * op_comps * 6 + i * 6 + 4);

        auto to_0_z_0_zz = pbuffer.data(idx_op_geom_001_sd + 2 * op_comps * 6 + i * 6 + 5);

        #pragma omp simd aligned(to_0_x, to_0_xxz, to_0_xyz, to_0_xzz, to_0_y, to_0_yyz, to_0_yzz, to_0_z, to_0_z_0_xx, to_0_z_0_xy, to_0_z_0_xz, to_0_z_0_yy, to_0_z_0_yz, to_0_z_0_zz, to_0_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_0_xx[k] = 2.0 * to_0_xxz[k] * tke_0;

            to_0_z_0_xy[k] = 2.0 * to_0_xyz[k] * tke_0;

            to_0_z_0_xz[k] = -to_0_x[k] + 2.0 * to_0_xzz[k] * tke_0;

            to_0_z_0_yy[k] = 2.0 * to_0_yyz[k] * tke_0;

            to_0_z_0_yz[k] = -to_0_y[k] + 2.0 * to_0_yzz[k] * tke_0;

            to_0_z_0_zz[k] = -2.0 * to_0_z[k] + 2.0 * to_0_zzz[k] * tke_0;
        }

    }

}

} // t2cgeom namespace

