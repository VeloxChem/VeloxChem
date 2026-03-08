#include "GeometricalDerivatives010ForSF.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_010_sf(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_sf,
                         const int idx_op_sd,
                         const int idx_op_sg,
                         const int idx_op_pf,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up components of auxiliary buffer : SD

    auto tr_0_xx = pbuffer.data(idx_op_sd);

    auto tr_0_xy = pbuffer.data(idx_op_sd + 1);

    auto tr_0_xz = pbuffer.data(idx_op_sd + 2);

    auto tr_0_yy = pbuffer.data(idx_op_sd + 3);

    auto tr_0_yz = pbuffer.data(idx_op_sd + 4);

    auto tr_0_zz = pbuffer.data(idx_op_sd + 5);

    // Set up components of auxiliary buffer : SG

    auto tr_0_xxxx = pbuffer.data(idx_op_sg);

    auto tr_0_xxxy = pbuffer.data(idx_op_sg + 1);

    auto tr_0_xxxz = pbuffer.data(idx_op_sg + 2);

    auto tr_0_xxyy = pbuffer.data(idx_op_sg + 3);

    auto tr_0_xxyz = pbuffer.data(idx_op_sg + 4);

    auto tr_0_xxzz = pbuffer.data(idx_op_sg + 5);

    auto tr_0_xyyy = pbuffer.data(idx_op_sg + 6);

    auto tr_0_xyyz = pbuffer.data(idx_op_sg + 7);

    auto tr_0_xyzz = pbuffer.data(idx_op_sg + 8);

    auto tr_0_xzzz = pbuffer.data(idx_op_sg + 9);

    auto tr_0_yyyy = pbuffer.data(idx_op_sg + 10);

    auto tr_0_yyyz = pbuffer.data(idx_op_sg + 11);

    auto tr_0_yyzz = pbuffer.data(idx_op_sg + 12);

    auto tr_0_yzzz = pbuffer.data(idx_op_sg + 13);

    auto tr_0_zzzz = pbuffer.data(idx_op_sg + 14);

    // Set up components of auxiliary buffer : PF

    auto tr_x_xxx = pbuffer.data(idx_op_pf);

    auto tr_x_xxy = pbuffer.data(idx_op_pf + 1);

    auto tr_x_xxz = pbuffer.data(idx_op_pf + 2);

    auto tr_x_xyy = pbuffer.data(idx_op_pf + 3);

    auto tr_x_xyz = pbuffer.data(idx_op_pf + 4);

    auto tr_x_xzz = pbuffer.data(idx_op_pf + 5);

    auto tr_x_yyy = pbuffer.data(idx_op_pf + 6);

    auto tr_x_yyz = pbuffer.data(idx_op_pf + 7);

    auto tr_x_yzz = pbuffer.data(idx_op_pf + 8);

    auto tr_x_zzz = pbuffer.data(idx_op_pf + 9);

    auto tr_y_xxx = pbuffer.data(idx_op_pf + 10);

    auto tr_y_xxy = pbuffer.data(idx_op_pf + 11);

    auto tr_y_xxz = pbuffer.data(idx_op_pf + 12);

    auto tr_y_xyy = pbuffer.data(idx_op_pf + 13);

    auto tr_y_xyz = pbuffer.data(idx_op_pf + 14);

    auto tr_y_xzz = pbuffer.data(idx_op_pf + 15);

    auto tr_y_yyy = pbuffer.data(idx_op_pf + 16);

    auto tr_y_yyz = pbuffer.data(idx_op_pf + 17);

    auto tr_y_yzz = pbuffer.data(idx_op_pf + 18);

    auto tr_y_zzz = pbuffer.data(idx_op_pf + 19);

    auto tr_z_xxx = pbuffer.data(idx_op_pf + 20);

    auto tr_z_xxy = pbuffer.data(idx_op_pf + 21);

    auto tr_z_xxz = pbuffer.data(idx_op_pf + 22);

    auto tr_z_xyy = pbuffer.data(idx_op_pf + 23);

    auto tr_z_xyz = pbuffer.data(idx_op_pf + 24);

    auto tr_z_xzz = pbuffer.data(idx_op_pf + 25);

    auto tr_z_yyy = pbuffer.data(idx_op_pf + 26);

    auto tr_z_yyz = pbuffer.data(idx_op_pf + 27);

    auto tr_z_yzz = pbuffer.data(idx_op_pf + 28);

    auto tr_z_zzz = pbuffer.data(idx_op_pf + 29);

    // Set up components of targeted buffer : SF

    auto tr_0_0_x_0_xxx = pbuffer.data(idx_op_geom_010_sf);

    auto tr_0_0_x_0_xxy = pbuffer.data(idx_op_geom_010_sf + 1);

    auto tr_0_0_x_0_xxz = pbuffer.data(idx_op_geom_010_sf + 2);

    auto tr_0_0_x_0_xyy = pbuffer.data(idx_op_geom_010_sf + 3);

    auto tr_0_0_x_0_xyz = pbuffer.data(idx_op_geom_010_sf + 4);

    auto tr_0_0_x_0_xzz = pbuffer.data(idx_op_geom_010_sf + 5);

    auto tr_0_0_x_0_yyy = pbuffer.data(idx_op_geom_010_sf + 6);

    auto tr_0_0_x_0_yyz = pbuffer.data(idx_op_geom_010_sf + 7);

    auto tr_0_0_x_0_yzz = pbuffer.data(idx_op_geom_010_sf + 8);

    auto tr_0_0_x_0_zzz = pbuffer.data(idx_op_geom_010_sf + 9);

    auto tr_0_0_y_0_xxx = pbuffer.data(idx_op_geom_010_sf + 10);

    auto tr_0_0_y_0_xxy = pbuffer.data(idx_op_geom_010_sf + 11);

    auto tr_0_0_y_0_xxz = pbuffer.data(idx_op_geom_010_sf + 12);

    auto tr_0_0_y_0_xyy = pbuffer.data(idx_op_geom_010_sf + 13);

    auto tr_0_0_y_0_xyz = pbuffer.data(idx_op_geom_010_sf + 14);

    auto tr_0_0_y_0_xzz = pbuffer.data(idx_op_geom_010_sf + 15);

    auto tr_0_0_y_0_yyy = pbuffer.data(idx_op_geom_010_sf + 16);

    auto tr_0_0_y_0_yyz = pbuffer.data(idx_op_geom_010_sf + 17);

    auto tr_0_0_y_0_yzz = pbuffer.data(idx_op_geom_010_sf + 18);

    auto tr_0_0_y_0_zzz = pbuffer.data(idx_op_geom_010_sf + 19);

    auto tr_0_0_z_0_xxx = pbuffer.data(idx_op_geom_010_sf + 20);

    auto tr_0_0_z_0_xxy = pbuffer.data(idx_op_geom_010_sf + 21);

    auto tr_0_0_z_0_xxz = pbuffer.data(idx_op_geom_010_sf + 22);

    auto tr_0_0_z_0_xyy = pbuffer.data(idx_op_geom_010_sf + 23);

    auto tr_0_0_z_0_xyz = pbuffer.data(idx_op_geom_010_sf + 24);

    auto tr_0_0_z_0_xzz = pbuffer.data(idx_op_geom_010_sf + 25);

    auto tr_0_0_z_0_yyy = pbuffer.data(idx_op_geom_010_sf + 26);

    auto tr_0_0_z_0_yyz = pbuffer.data(idx_op_geom_010_sf + 27);

    auto tr_0_0_z_0_yzz = pbuffer.data(idx_op_geom_010_sf + 28);

    auto tr_0_0_z_0_zzz = pbuffer.data(idx_op_geom_010_sf + 29);

    #pragma omp simd aligned(tr_0_0_x_0_xxx, tr_0_0_x_0_xxy, tr_0_0_x_0_xxz, tr_0_0_x_0_xyy, tr_0_0_x_0_xyz, tr_0_0_x_0_xzz, tr_0_0_x_0_yyy, tr_0_0_x_0_yyz, tr_0_0_x_0_yzz, tr_0_0_x_0_zzz, tr_0_0_y_0_xxx, tr_0_0_y_0_xxy, tr_0_0_y_0_xxz, tr_0_0_y_0_xyy, tr_0_0_y_0_xyz, tr_0_0_y_0_xzz, tr_0_0_y_0_yyy, tr_0_0_y_0_yyz, tr_0_0_y_0_yzz, tr_0_0_y_0_zzz, tr_0_0_z_0_xxx, tr_0_0_z_0_xxy, tr_0_0_z_0_xxz, tr_0_0_z_0_xyy, tr_0_0_z_0_xyz, tr_0_0_z_0_xzz, tr_0_0_z_0_yyy, tr_0_0_z_0_yyz, tr_0_0_z_0_yzz, tr_0_0_z_0_zzz, tr_0_xx, tr_0_xxxx, tr_0_xxxy, tr_0_xxxz, tr_0_xxyy, tr_0_xxyz, tr_0_xxzz, tr_0_xy, tr_0_xyyy, tr_0_xyyz, tr_0_xyzz, tr_0_xz, tr_0_xzzz, tr_0_yy, tr_0_yyyy, tr_0_yyyz, tr_0_yyzz, tr_0_yz, tr_0_yzzz, tr_0_zz, tr_0_zzzz, tr_x_xxx, tr_x_xxy, tr_x_xxz, tr_x_xyy, tr_x_xyz, tr_x_xzz, tr_x_yyy, tr_x_yyz, tr_x_yzz, tr_x_zzz, tr_y_xxx, tr_y_xxy, tr_y_xxz, tr_y_xyy, tr_y_xyz, tr_y_xzz, tr_y_yyy, tr_y_yyz, tr_y_yzz, tr_y_zzz, tr_z_xxx, tr_z_xxy, tr_z_xxz, tr_z_xyy, tr_z_xyz, tr_z_xzz, tr_z_yyy, tr_z_yyz, tr_z_yzz, tr_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_0_xxx[i] = 2.0 * tr_x_xxx[i] * tbe_0 + 2.0 * tr_0_xxxx[i] * tke_0 - 3.0 * tr_0_xx[i];

        tr_0_0_x_0_xxy[i] = 2.0 * tr_x_xxy[i] * tbe_0 + 2.0 * tr_0_xxxy[i] * tke_0 - 2.0 * tr_0_xy[i];

        tr_0_0_x_0_xxz[i] = 2.0 * tr_x_xxz[i] * tbe_0 + 2.0 * tr_0_xxxz[i] * tke_0 - 2.0 * tr_0_xz[i];

        tr_0_0_x_0_xyy[i] = 2.0 * tr_x_xyy[i] * tbe_0 + 2.0 * tr_0_xxyy[i] * tke_0 - tr_0_yy[i];

        tr_0_0_x_0_xyz[i] = 2.0 * tr_x_xyz[i] * tbe_0 + 2.0 * tr_0_xxyz[i] * tke_0 - tr_0_yz[i];

        tr_0_0_x_0_xzz[i] = 2.0 * tr_x_xzz[i] * tbe_0 + 2.0 * tr_0_xxzz[i] * tke_0 - tr_0_zz[i];

        tr_0_0_x_0_yyy[i] = 2.0 * tr_x_yyy[i] * tbe_0 + 2.0 * tr_0_xyyy[i] * tke_0;

        tr_0_0_x_0_yyz[i] = 2.0 * tr_x_yyz[i] * tbe_0 + 2.0 * tr_0_xyyz[i] * tke_0;

        tr_0_0_x_0_yzz[i] = 2.0 * tr_x_yzz[i] * tbe_0 + 2.0 * tr_0_xyzz[i] * tke_0;

        tr_0_0_x_0_zzz[i] = 2.0 * tr_x_zzz[i] * tbe_0 + 2.0 * tr_0_xzzz[i] * tke_0;

        tr_0_0_y_0_xxx[i] = 2.0 * tr_y_xxx[i] * tbe_0 + 2.0 * tr_0_xxxy[i] * tke_0;

        tr_0_0_y_0_xxy[i] = 2.0 * tr_y_xxy[i] * tbe_0 + 2.0 * tr_0_xxyy[i] * tke_0 - tr_0_xx[i];

        tr_0_0_y_0_xxz[i] = 2.0 * tr_y_xxz[i] * tbe_0 + 2.0 * tr_0_xxyz[i] * tke_0;

        tr_0_0_y_0_xyy[i] = 2.0 * tr_y_xyy[i] * tbe_0 + 2.0 * tr_0_xyyy[i] * tke_0 - 2.0 * tr_0_xy[i];

        tr_0_0_y_0_xyz[i] = 2.0 * tr_y_xyz[i] * tbe_0 + 2.0 * tr_0_xyyz[i] * tke_0 - tr_0_xz[i];

        tr_0_0_y_0_xzz[i] = 2.0 * tr_y_xzz[i] * tbe_0 + 2.0 * tr_0_xyzz[i] * tke_0;

        tr_0_0_y_0_yyy[i] = 2.0 * tr_y_yyy[i] * tbe_0 + 2.0 * tr_0_yyyy[i] * tke_0 - 3.0 * tr_0_yy[i];

        tr_0_0_y_0_yyz[i] = 2.0 * tr_y_yyz[i] * tbe_0 + 2.0 * tr_0_yyyz[i] * tke_0 - 2.0 * tr_0_yz[i];

        tr_0_0_y_0_yzz[i] = 2.0 * tr_y_yzz[i] * tbe_0 + 2.0 * tr_0_yyzz[i] * tke_0 - tr_0_zz[i];

        tr_0_0_y_0_zzz[i] = 2.0 * tr_y_zzz[i] * tbe_0 + 2.0 * tr_0_yzzz[i] * tke_0;

        tr_0_0_z_0_xxx[i] = 2.0 * tr_z_xxx[i] * tbe_0 + 2.0 * tr_0_xxxz[i] * tke_0;

        tr_0_0_z_0_xxy[i] = 2.0 * tr_z_xxy[i] * tbe_0 + 2.0 * tr_0_xxyz[i] * tke_0;

        tr_0_0_z_0_xxz[i] = 2.0 * tr_z_xxz[i] * tbe_0 + 2.0 * tr_0_xxzz[i] * tke_0 - tr_0_xx[i];

        tr_0_0_z_0_xyy[i] = 2.0 * tr_z_xyy[i] * tbe_0 + 2.0 * tr_0_xyyz[i] * tke_0;

        tr_0_0_z_0_xyz[i] = 2.0 * tr_z_xyz[i] * tbe_0 + 2.0 * tr_0_xyzz[i] * tke_0 - tr_0_xy[i];

        tr_0_0_z_0_xzz[i] = 2.0 * tr_z_xzz[i] * tbe_0 + 2.0 * tr_0_xzzz[i] * tke_0 - 2.0 * tr_0_xz[i];

        tr_0_0_z_0_yyy[i] = 2.0 * tr_z_yyy[i] * tbe_0 + 2.0 * tr_0_yyyz[i] * tke_0;

        tr_0_0_z_0_yyz[i] = 2.0 * tr_z_yyz[i] * tbe_0 + 2.0 * tr_0_yyzz[i] * tke_0 - tr_0_yy[i];

        tr_0_0_z_0_yzz[i] = 2.0 * tr_z_yzz[i] * tbe_0 + 2.0 * tr_0_yzzz[i] * tke_0 - 2.0 * tr_0_yz[i];

        tr_0_0_z_0_zzz[i] = 2.0 * tr_z_zzz[i] * tbe_0 + 2.0 * tr_0_zzzz[i] * tke_0 - 3.0 * tr_0_zz[i];
    }
}

} // t2cgeom namespace

