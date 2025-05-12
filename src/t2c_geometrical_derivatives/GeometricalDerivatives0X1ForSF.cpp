#include "GeometricalDerivatives0X1ForSF.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_geom_deriv_0x1_sf(CSimdArray<double>& pbuffer,
                       const int idx_op_geom_001_sf,
                       const int idx_op_sd,
                       const int idx_op_sg,
                       const CSimdArray<double>& factors,
                       const int op_comps) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : SD

        auto to_0_xx = pbuffer.data(idx_op_sd + i * 6 + 0);

        auto to_0_xy = pbuffer.data(idx_op_sd + i * 6 + 1);

        auto to_0_xz = pbuffer.data(idx_op_sd + i * 6 + 2);

        auto to_0_yy = pbuffer.data(idx_op_sd + i * 6 + 3);

        auto to_0_yz = pbuffer.data(idx_op_sd + i * 6 + 4);

        auto to_0_zz = pbuffer.data(idx_op_sd + i * 6 + 5);

        // Set up components of auxiliary buffer : SG

        auto to_0_xxxx = pbuffer.data(idx_op_sg + i * 15 + 0);

        auto to_0_xxxy = pbuffer.data(idx_op_sg + i * 15 + 1);

        auto to_0_xxxz = pbuffer.data(idx_op_sg + i * 15 + 2);

        auto to_0_xxyy = pbuffer.data(idx_op_sg + i * 15 + 3);

        auto to_0_xxyz = pbuffer.data(idx_op_sg + i * 15 + 4);

        auto to_0_xxzz = pbuffer.data(idx_op_sg + i * 15 + 5);

        auto to_0_xyyy = pbuffer.data(idx_op_sg + i * 15 + 6);

        auto to_0_xyyz = pbuffer.data(idx_op_sg + i * 15 + 7);

        auto to_0_xyzz = pbuffer.data(idx_op_sg + i * 15 + 8);

        auto to_0_xzzz = pbuffer.data(idx_op_sg + i * 15 + 9);

        auto to_0_yyyy = pbuffer.data(idx_op_sg + i * 15 + 10);

        auto to_0_yyyz = pbuffer.data(idx_op_sg + i * 15 + 11);

        auto to_0_yyzz = pbuffer.data(idx_op_sg + i * 15 + 12);

        auto to_0_yzzz = pbuffer.data(idx_op_sg + i * 15 + 13);

        auto to_0_zzzz = pbuffer.data(idx_op_sg + i * 15 + 14);

        // Set up 0-10 components of targeted buffer : SF

        auto to_0_x_0_xxx = pbuffer.data(idx_op_geom_001_sf + 0 * op_comps * 10 + i * 10 + 0);

        auto to_0_x_0_xxy = pbuffer.data(idx_op_geom_001_sf + 0 * op_comps * 10 + i * 10 + 1);

        auto to_0_x_0_xxz = pbuffer.data(idx_op_geom_001_sf + 0 * op_comps * 10 + i * 10 + 2);

        auto to_0_x_0_xyy = pbuffer.data(idx_op_geom_001_sf + 0 * op_comps * 10 + i * 10 + 3);

        auto to_0_x_0_xyz = pbuffer.data(idx_op_geom_001_sf + 0 * op_comps * 10 + i * 10 + 4);

        auto to_0_x_0_xzz = pbuffer.data(idx_op_geom_001_sf + 0 * op_comps * 10 + i * 10 + 5);

        auto to_0_x_0_yyy = pbuffer.data(idx_op_geom_001_sf + 0 * op_comps * 10 + i * 10 + 6);

        auto to_0_x_0_yyz = pbuffer.data(idx_op_geom_001_sf + 0 * op_comps * 10 + i * 10 + 7);

        auto to_0_x_0_yzz = pbuffer.data(idx_op_geom_001_sf + 0 * op_comps * 10 + i * 10 + 8);

        auto to_0_x_0_zzz = pbuffer.data(idx_op_geom_001_sf + 0 * op_comps * 10 + i * 10 + 9);

        #pragma omp simd aligned(to_0_x_0_xxx, to_0_x_0_xxy, to_0_x_0_xxz, to_0_x_0_xyy, to_0_x_0_xyz, to_0_x_0_xzz, to_0_x_0_yyy, to_0_x_0_yyz, to_0_x_0_yzz, to_0_x_0_zzz, to_0_xx, to_0_xxxx, to_0_xxxy, to_0_xxxz, to_0_xxyy, to_0_xxyz, to_0_xxzz, to_0_xy, to_0_xyyy, to_0_xyyz, to_0_xyzz, to_0_xz, to_0_xzzz, to_0_yy, to_0_yz, to_0_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_0_xxx[k] = -3.0 * to_0_xx[k] + 2.0 * to_0_xxxx[k] * tke_0;

            to_0_x_0_xxy[k] = -2.0 * to_0_xy[k] + 2.0 * to_0_xxxy[k] * tke_0;

            to_0_x_0_xxz[k] = -2.0 * to_0_xz[k] + 2.0 * to_0_xxxz[k] * tke_0;

            to_0_x_0_xyy[k] = -to_0_yy[k] + 2.0 * to_0_xxyy[k] * tke_0;

            to_0_x_0_xyz[k] = -to_0_yz[k] + 2.0 * to_0_xxyz[k] * tke_0;

            to_0_x_0_xzz[k] = -to_0_zz[k] + 2.0 * to_0_xxzz[k] * tke_0;

            to_0_x_0_yyy[k] = 2.0 * to_0_xyyy[k] * tke_0;

            to_0_x_0_yyz[k] = 2.0 * to_0_xyyz[k] * tke_0;

            to_0_x_0_yzz[k] = 2.0 * to_0_xyzz[k] * tke_0;

            to_0_x_0_zzz[k] = 2.0 * to_0_xzzz[k] * tke_0;
        }

        // Set up 10-20 components of targeted buffer : SF

        auto to_0_y_0_xxx = pbuffer.data(idx_op_geom_001_sf + 1 * op_comps * 10 + i * 10 + 0);

        auto to_0_y_0_xxy = pbuffer.data(idx_op_geom_001_sf + 1 * op_comps * 10 + i * 10 + 1);

        auto to_0_y_0_xxz = pbuffer.data(idx_op_geom_001_sf + 1 * op_comps * 10 + i * 10 + 2);

        auto to_0_y_0_xyy = pbuffer.data(idx_op_geom_001_sf + 1 * op_comps * 10 + i * 10 + 3);

        auto to_0_y_0_xyz = pbuffer.data(idx_op_geom_001_sf + 1 * op_comps * 10 + i * 10 + 4);

        auto to_0_y_0_xzz = pbuffer.data(idx_op_geom_001_sf + 1 * op_comps * 10 + i * 10 + 5);

        auto to_0_y_0_yyy = pbuffer.data(idx_op_geom_001_sf + 1 * op_comps * 10 + i * 10 + 6);

        auto to_0_y_0_yyz = pbuffer.data(idx_op_geom_001_sf + 1 * op_comps * 10 + i * 10 + 7);

        auto to_0_y_0_yzz = pbuffer.data(idx_op_geom_001_sf + 1 * op_comps * 10 + i * 10 + 8);

        auto to_0_y_0_zzz = pbuffer.data(idx_op_geom_001_sf + 1 * op_comps * 10 + i * 10 + 9);

        #pragma omp simd aligned(to_0_xx, to_0_xxxy, to_0_xxyy, to_0_xxyz, to_0_xy, to_0_xyyy, to_0_xyyz, to_0_xyzz, to_0_xz, to_0_y_0_xxx, to_0_y_0_xxy, to_0_y_0_xxz, to_0_y_0_xyy, to_0_y_0_xyz, to_0_y_0_xzz, to_0_y_0_yyy, to_0_y_0_yyz, to_0_y_0_yzz, to_0_y_0_zzz, to_0_yy, to_0_yyyy, to_0_yyyz, to_0_yyzz, to_0_yz, to_0_yzzz, to_0_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_0_xxx[k] = 2.0 * to_0_xxxy[k] * tke_0;

            to_0_y_0_xxy[k] = -to_0_xx[k] + 2.0 * to_0_xxyy[k] * tke_0;

            to_0_y_0_xxz[k] = 2.0 * to_0_xxyz[k] * tke_0;

            to_0_y_0_xyy[k] = -2.0 * to_0_xy[k] + 2.0 * to_0_xyyy[k] * tke_0;

            to_0_y_0_xyz[k] = -to_0_xz[k] + 2.0 * to_0_xyyz[k] * tke_0;

            to_0_y_0_xzz[k] = 2.0 * to_0_xyzz[k] * tke_0;

            to_0_y_0_yyy[k] = -3.0 * to_0_yy[k] + 2.0 * to_0_yyyy[k] * tke_0;

            to_0_y_0_yyz[k] = -2.0 * to_0_yz[k] + 2.0 * to_0_yyyz[k] * tke_0;

            to_0_y_0_yzz[k] = -to_0_zz[k] + 2.0 * to_0_yyzz[k] * tke_0;

            to_0_y_0_zzz[k] = 2.0 * to_0_yzzz[k] * tke_0;
        }

        // Set up 20-30 components of targeted buffer : SF

        auto to_0_z_0_xxx = pbuffer.data(idx_op_geom_001_sf + 2 * op_comps * 10 + i * 10 + 0);

        auto to_0_z_0_xxy = pbuffer.data(idx_op_geom_001_sf + 2 * op_comps * 10 + i * 10 + 1);

        auto to_0_z_0_xxz = pbuffer.data(idx_op_geom_001_sf + 2 * op_comps * 10 + i * 10 + 2);

        auto to_0_z_0_xyy = pbuffer.data(idx_op_geom_001_sf + 2 * op_comps * 10 + i * 10 + 3);

        auto to_0_z_0_xyz = pbuffer.data(idx_op_geom_001_sf + 2 * op_comps * 10 + i * 10 + 4);

        auto to_0_z_0_xzz = pbuffer.data(idx_op_geom_001_sf + 2 * op_comps * 10 + i * 10 + 5);

        auto to_0_z_0_yyy = pbuffer.data(idx_op_geom_001_sf + 2 * op_comps * 10 + i * 10 + 6);

        auto to_0_z_0_yyz = pbuffer.data(idx_op_geom_001_sf + 2 * op_comps * 10 + i * 10 + 7);

        auto to_0_z_0_yzz = pbuffer.data(idx_op_geom_001_sf + 2 * op_comps * 10 + i * 10 + 8);

        auto to_0_z_0_zzz = pbuffer.data(idx_op_geom_001_sf + 2 * op_comps * 10 + i * 10 + 9);

        #pragma omp simd aligned(to_0_xx, to_0_xxxz, to_0_xxyz, to_0_xxzz, to_0_xy, to_0_xyyz, to_0_xyzz, to_0_xz, to_0_xzzz, to_0_yy, to_0_yyyz, to_0_yyzz, to_0_yz, to_0_yzzz, to_0_z_0_xxx, to_0_z_0_xxy, to_0_z_0_xxz, to_0_z_0_xyy, to_0_z_0_xyz, to_0_z_0_xzz, to_0_z_0_yyy, to_0_z_0_yyz, to_0_z_0_yzz, to_0_z_0_zzz, to_0_zz, to_0_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_0_xxx[k] = 2.0 * to_0_xxxz[k] * tke_0;

            to_0_z_0_xxy[k] = 2.0 * to_0_xxyz[k] * tke_0;

            to_0_z_0_xxz[k] = -to_0_xx[k] + 2.0 * to_0_xxzz[k] * tke_0;

            to_0_z_0_xyy[k] = 2.0 * to_0_xyyz[k] * tke_0;

            to_0_z_0_xyz[k] = -to_0_xy[k] + 2.0 * to_0_xyzz[k] * tke_0;

            to_0_z_0_xzz[k] = -2.0 * to_0_xz[k] + 2.0 * to_0_xzzz[k] * tke_0;

            to_0_z_0_yyy[k] = 2.0 * to_0_yyyz[k] * tke_0;

            to_0_z_0_yyz[k] = -to_0_yy[k] + 2.0 * to_0_yyzz[k] * tke_0;

            to_0_z_0_yzz[k] = -2.0 * to_0_yz[k] + 2.0 * to_0_yzzz[k] * tke_0;

            to_0_z_0_zzz[k] = -3.0 * to_0_zz[k] + 2.0 * to_0_zzzz[k] * tke_0;
        }

    }

}

} // t2cgeom namespace

