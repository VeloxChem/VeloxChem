#include "ThreeCenterElectronRepulsionGeom010ContrRecXDD.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_ket_geom010_electron_repulsion_xdd(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xdd,
                                        const size_t idx_xpd,
                                        const size_t idx_geom_10_xpd,
                                        const size_t idx_geom_10_xpf,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_cd,
                                        const int a_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom,});

    // Set up R(CD) distances

    auto cd_x = factors.data(idx_cd);

    auto cd_y = factors.data(idx_cd + 1);

    auto cd_z = factors.data(idx_cd + 2);

    for (int i = 0; i < acomps; i++)
    {
        /// Set up components of auxilary buffer : SPD

        const auto pd_off = idx_xpd + i * 18;

        auto g_x_xx = cbuffer.data(pd_off + 0);

        auto g_x_xy = cbuffer.data(pd_off + 1);

        auto g_x_xz = cbuffer.data(pd_off + 2);

        auto g_x_yy = cbuffer.data(pd_off + 3);

        auto g_x_yz = cbuffer.data(pd_off + 4);

        auto g_x_zz = cbuffer.data(pd_off + 5);

        auto g_y_xx = cbuffer.data(pd_off + 6);

        auto g_y_xy = cbuffer.data(pd_off + 7);

        auto g_y_xz = cbuffer.data(pd_off + 8);

        auto g_y_yy = cbuffer.data(pd_off + 9);

        auto g_y_yz = cbuffer.data(pd_off + 10);

        auto g_y_zz = cbuffer.data(pd_off + 11);

        auto g_z_xx = cbuffer.data(pd_off + 12);

        auto g_z_xy = cbuffer.data(pd_off + 13);

        auto g_z_xz = cbuffer.data(pd_off + 14);

        auto g_z_yy = cbuffer.data(pd_off + 15);

        auto g_z_yz = cbuffer.data(pd_off + 16);

        auto g_z_zz = cbuffer.data(pd_off + 17);

        /// Set up components of auxilary buffer : SPD

        const auto pd_geom_10_off = idx_geom_10_xpd + i * 18;

        auto g_x_0_x_xx = cbuffer.data(pd_geom_10_off + 0 * acomps + 0);

        auto g_x_0_x_xy = cbuffer.data(pd_geom_10_off + 0 * acomps + 1);

        auto g_x_0_x_xz = cbuffer.data(pd_geom_10_off + 0 * acomps + 2);

        auto g_x_0_x_yy = cbuffer.data(pd_geom_10_off + 0 * acomps + 3);

        auto g_x_0_x_yz = cbuffer.data(pd_geom_10_off + 0 * acomps + 4);

        auto g_x_0_x_zz = cbuffer.data(pd_geom_10_off + 0 * acomps + 5);

        auto g_x_0_y_xx = cbuffer.data(pd_geom_10_off + 6 * acomps + 0);

        auto g_x_0_y_xy = cbuffer.data(pd_geom_10_off + 6 * acomps + 1);

        auto g_x_0_y_xz = cbuffer.data(pd_geom_10_off + 6 * acomps + 2);

        auto g_x_0_y_yy = cbuffer.data(pd_geom_10_off + 6 * acomps + 3);

        auto g_x_0_y_yz = cbuffer.data(pd_geom_10_off + 6 * acomps + 4);

        auto g_x_0_y_zz = cbuffer.data(pd_geom_10_off + 6 * acomps + 5);

        auto g_x_0_z_xx = cbuffer.data(pd_geom_10_off + 12 * acomps + 0);

        auto g_x_0_z_xy = cbuffer.data(pd_geom_10_off + 12 * acomps + 1);

        auto g_x_0_z_xz = cbuffer.data(pd_geom_10_off + 12 * acomps + 2);

        auto g_x_0_z_yy = cbuffer.data(pd_geom_10_off + 12 * acomps + 3);

        auto g_x_0_z_yz = cbuffer.data(pd_geom_10_off + 12 * acomps + 4);

        auto g_x_0_z_zz = cbuffer.data(pd_geom_10_off + 12 * acomps + 5);

        auto g_y_0_x_xx = cbuffer.data(pd_geom_10_off + 18 * acomps + 0);

        auto g_y_0_x_xy = cbuffer.data(pd_geom_10_off + 18 * acomps + 1);

        auto g_y_0_x_xz = cbuffer.data(pd_geom_10_off + 18 * acomps + 2);

        auto g_y_0_x_yy = cbuffer.data(pd_geom_10_off + 18 * acomps + 3);

        auto g_y_0_x_yz = cbuffer.data(pd_geom_10_off + 18 * acomps + 4);

        auto g_y_0_x_zz = cbuffer.data(pd_geom_10_off + 18 * acomps + 5);

        auto g_y_0_y_xx = cbuffer.data(pd_geom_10_off + 24 * acomps + 0);

        auto g_y_0_y_xy = cbuffer.data(pd_geom_10_off + 24 * acomps + 1);

        auto g_y_0_y_xz = cbuffer.data(pd_geom_10_off + 24 * acomps + 2);

        auto g_y_0_y_yy = cbuffer.data(pd_geom_10_off + 24 * acomps + 3);

        auto g_y_0_y_yz = cbuffer.data(pd_geom_10_off + 24 * acomps + 4);

        auto g_y_0_y_zz = cbuffer.data(pd_geom_10_off + 24 * acomps + 5);

        auto g_y_0_z_xx = cbuffer.data(pd_geom_10_off + 30 * acomps + 0);

        auto g_y_0_z_xy = cbuffer.data(pd_geom_10_off + 30 * acomps + 1);

        auto g_y_0_z_xz = cbuffer.data(pd_geom_10_off + 30 * acomps + 2);

        auto g_y_0_z_yy = cbuffer.data(pd_geom_10_off + 30 * acomps + 3);

        auto g_y_0_z_yz = cbuffer.data(pd_geom_10_off + 30 * acomps + 4);

        auto g_y_0_z_zz = cbuffer.data(pd_geom_10_off + 30 * acomps + 5);

        auto g_z_0_x_xx = cbuffer.data(pd_geom_10_off + 36 * acomps + 0);

        auto g_z_0_x_xy = cbuffer.data(pd_geom_10_off + 36 * acomps + 1);

        auto g_z_0_x_xz = cbuffer.data(pd_geom_10_off + 36 * acomps + 2);

        auto g_z_0_x_yy = cbuffer.data(pd_geom_10_off + 36 * acomps + 3);

        auto g_z_0_x_yz = cbuffer.data(pd_geom_10_off + 36 * acomps + 4);

        auto g_z_0_x_zz = cbuffer.data(pd_geom_10_off + 36 * acomps + 5);

        auto g_z_0_y_xx = cbuffer.data(pd_geom_10_off + 42 * acomps + 0);

        auto g_z_0_y_xy = cbuffer.data(pd_geom_10_off + 42 * acomps + 1);

        auto g_z_0_y_xz = cbuffer.data(pd_geom_10_off + 42 * acomps + 2);

        auto g_z_0_y_yy = cbuffer.data(pd_geom_10_off + 42 * acomps + 3);

        auto g_z_0_y_yz = cbuffer.data(pd_geom_10_off + 42 * acomps + 4);

        auto g_z_0_y_zz = cbuffer.data(pd_geom_10_off + 42 * acomps + 5);

        auto g_z_0_z_xx = cbuffer.data(pd_geom_10_off + 48 * acomps + 0);

        auto g_z_0_z_xy = cbuffer.data(pd_geom_10_off + 48 * acomps + 1);

        auto g_z_0_z_xz = cbuffer.data(pd_geom_10_off + 48 * acomps + 2);

        auto g_z_0_z_yy = cbuffer.data(pd_geom_10_off + 48 * acomps + 3);

        auto g_z_0_z_yz = cbuffer.data(pd_geom_10_off + 48 * acomps + 4);

        auto g_z_0_z_zz = cbuffer.data(pd_geom_10_off + 48 * acomps + 5);

        /// Set up components of auxilary buffer : SPF

        const auto pf_geom_10_off = idx_geom_10_xpf + i * 30;

        auto g_x_0_x_xxx = cbuffer.data(pf_geom_10_off + 0 * acomps + 0);

        auto g_x_0_x_xxy = cbuffer.data(pf_geom_10_off + 0 * acomps + 1);

        auto g_x_0_x_xxz = cbuffer.data(pf_geom_10_off + 0 * acomps + 2);

        auto g_x_0_x_xyy = cbuffer.data(pf_geom_10_off + 0 * acomps + 3);

        auto g_x_0_x_xyz = cbuffer.data(pf_geom_10_off + 0 * acomps + 4);

        auto g_x_0_x_xzz = cbuffer.data(pf_geom_10_off + 0 * acomps + 5);

        auto g_x_0_x_yyy = cbuffer.data(pf_geom_10_off + 0 * acomps + 6);

        auto g_x_0_x_yyz = cbuffer.data(pf_geom_10_off + 0 * acomps + 7);

        auto g_x_0_x_yzz = cbuffer.data(pf_geom_10_off + 0 * acomps + 8);

        auto g_x_0_x_zzz = cbuffer.data(pf_geom_10_off + 0 * acomps + 9);

        auto g_x_0_y_xxx = cbuffer.data(pf_geom_10_off + 10 * acomps + 0);

        auto g_x_0_y_xxy = cbuffer.data(pf_geom_10_off + 10 * acomps + 1);

        auto g_x_0_y_xxz = cbuffer.data(pf_geom_10_off + 10 * acomps + 2);

        auto g_x_0_y_xyy = cbuffer.data(pf_geom_10_off + 10 * acomps + 3);

        auto g_x_0_y_xyz = cbuffer.data(pf_geom_10_off + 10 * acomps + 4);

        auto g_x_0_y_xzz = cbuffer.data(pf_geom_10_off + 10 * acomps + 5);

        auto g_x_0_y_yyy = cbuffer.data(pf_geom_10_off + 10 * acomps + 6);

        auto g_x_0_y_yyz = cbuffer.data(pf_geom_10_off + 10 * acomps + 7);

        auto g_x_0_y_yzz = cbuffer.data(pf_geom_10_off + 10 * acomps + 8);

        auto g_x_0_y_zzz = cbuffer.data(pf_geom_10_off + 10 * acomps + 9);

        auto g_x_0_z_xxx = cbuffer.data(pf_geom_10_off + 20 * acomps + 0);

        auto g_x_0_z_xxy = cbuffer.data(pf_geom_10_off + 20 * acomps + 1);

        auto g_x_0_z_xxz = cbuffer.data(pf_geom_10_off + 20 * acomps + 2);

        auto g_x_0_z_xyy = cbuffer.data(pf_geom_10_off + 20 * acomps + 3);

        auto g_x_0_z_xyz = cbuffer.data(pf_geom_10_off + 20 * acomps + 4);

        auto g_x_0_z_xzz = cbuffer.data(pf_geom_10_off + 20 * acomps + 5);

        auto g_x_0_z_yyy = cbuffer.data(pf_geom_10_off + 20 * acomps + 6);

        auto g_x_0_z_yyz = cbuffer.data(pf_geom_10_off + 20 * acomps + 7);

        auto g_x_0_z_yzz = cbuffer.data(pf_geom_10_off + 20 * acomps + 8);

        auto g_x_0_z_zzz = cbuffer.data(pf_geom_10_off + 20 * acomps + 9);

        auto g_y_0_x_xxx = cbuffer.data(pf_geom_10_off + 30 * acomps + 0);

        auto g_y_0_x_xxy = cbuffer.data(pf_geom_10_off + 30 * acomps + 1);

        auto g_y_0_x_xxz = cbuffer.data(pf_geom_10_off + 30 * acomps + 2);

        auto g_y_0_x_xyy = cbuffer.data(pf_geom_10_off + 30 * acomps + 3);

        auto g_y_0_x_xyz = cbuffer.data(pf_geom_10_off + 30 * acomps + 4);

        auto g_y_0_x_xzz = cbuffer.data(pf_geom_10_off + 30 * acomps + 5);

        auto g_y_0_x_yyy = cbuffer.data(pf_geom_10_off + 30 * acomps + 6);

        auto g_y_0_x_yyz = cbuffer.data(pf_geom_10_off + 30 * acomps + 7);

        auto g_y_0_x_yzz = cbuffer.data(pf_geom_10_off + 30 * acomps + 8);

        auto g_y_0_x_zzz = cbuffer.data(pf_geom_10_off + 30 * acomps + 9);

        auto g_y_0_y_xxx = cbuffer.data(pf_geom_10_off + 40 * acomps + 0);

        auto g_y_0_y_xxy = cbuffer.data(pf_geom_10_off + 40 * acomps + 1);

        auto g_y_0_y_xxz = cbuffer.data(pf_geom_10_off + 40 * acomps + 2);

        auto g_y_0_y_xyy = cbuffer.data(pf_geom_10_off + 40 * acomps + 3);

        auto g_y_0_y_xyz = cbuffer.data(pf_geom_10_off + 40 * acomps + 4);

        auto g_y_0_y_xzz = cbuffer.data(pf_geom_10_off + 40 * acomps + 5);

        auto g_y_0_y_yyy = cbuffer.data(pf_geom_10_off + 40 * acomps + 6);

        auto g_y_0_y_yyz = cbuffer.data(pf_geom_10_off + 40 * acomps + 7);

        auto g_y_0_y_yzz = cbuffer.data(pf_geom_10_off + 40 * acomps + 8);

        auto g_y_0_y_zzz = cbuffer.data(pf_geom_10_off + 40 * acomps + 9);

        auto g_y_0_z_xxx = cbuffer.data(pf_geom_10_off + 50 * acomps + 0);

        auto g_y_0_z_xxy = cbuffer.data(pf_geom_10_off + 50 * acomps + 1);

        auto g_y_0_z_xxz = cbuffer.data(pf_geom_10_off + 50 * acomps + 2);

        auto g_y_0_z_xyy = cbuffer.data(pf_geom_10_off + 50 * acomps + 3);

        auto g_y_0_z_xyz = cbuffer.data(pf_geom_10_off + 50 * acomps + 4);

        auto g_y_0_z_xzz = cbuffer.data(pf_geom_10_off + 50 * acomps + 5);

        auto g_y_0_z_yyy = cbuffer.data(pf_geom_10_off + 50 * acomps + 6);

        auto g_y_0_z_yyz = cbuffer.data(pf_geom_10_off + 50 * acomps + 7);

        auto g_y_0_z_yzz = cbuffer.data(pf_geom_10_off + 50 * acomps + 8);

        auto g_y_0_z_zzz = cbuffer.data(pf_geom_10_off + 50 * acomps + 9);

        auto g_z_0_x_xxx = cbuffer.data(pf_geom_10_off + 60 * acomps + 0);

        auto g_z_0_x_xxy = cbuffer.data(pf_geom_10_off + 60 * acomps + 1);

        auto g_z_0_x_xxz = cbuffer.data(pf_geom_10_off + 60 * acomps + 2);

        auto g_z_0_x_xyy = cbuffer.data(pf_geom_10_off + 60 * acomps + 3);

        auto g_z_0_x_xyz = cbuffer.data(pf_geom_10_off + 60 * acomps + 4);

        auto g_z_0_x_xzz = cbuffer.data(pf_geom_10_off + 60 * acomps + 5);

        auto g_z_0_x_yyy = cbuffer.data(pf_geom_10_off + 60 * acomps + 6);

        auto g_z_0_x_yyz = cbuffer.data(pf_geom_10_off + 60 * acomps + 7);

        auto g_z_0_x_yzz = cbuffer.data(pf_geom_10_off + 60 * acomps + 8);

        auto g_z_0_x_zzz = cbuffer.data(pf_geom_10_off + 60 * acomps + 9);

        auto g_z_0_y_xxx = cbuffer.data(pf_geom_10_off + 70 * acomps + 0);

        auto g_z_0_y_xxy = cbuffer.data(pf_geom_10_off + 70 * acomps + 1);

        auto g_z_0_y_xxz = cbuffer.data(pf_geom_10_off + 70 * acomps + 2);

        auto g_z_0_y_xyy = cbuffer.data(pf_geom_10_off + 70 * acomps + 3);

        auto g_z_0_y_xyz = cbuffer.data(pf_geom_10_off + 70 * acomps + 4);

        auto g_z_0_y_xzz = cbuffer.data(pf_geom_10_off + 70 * acomps + 5);

        auto g_z_0_y_yyy = cbuffer.data(pf_geom_10_off + 70 * acomps + 6);

        auto g_z_0_y_yyz = cbuffer.data(pf_geom_10_off + 70 * acomps + 7);

        auto g_z_0_y_yzz = cbuffer.data(pf_geom_10_off + 70 * acomps + 8);

        auto g_z_0_y_zzz = cbuffer.data(pf_geom_10_off + 70 * acomps + 9);

        auto g_z_0_z_xxx = cbuffer.data(pf_geom_10_off + 80 * acomps + 0);

        auto g_z_0_z_xxy = cbuffer.data(pf_geom_10_off + 80 * acomps + 1);

        auto g_z_0_z_xxz = cbuffer.data(pf_geom_10_off + 80 * acomps + 2);

        auto g_z_0_z_xyy = cbuffer.data(pf_geom_10_off + 80 * acomps + 3);

        auto g_z_0_z_xyz = cbuffer.data(pf_geom_10_off + 80 * acomps + 4);

        auto g_z_0_z_xzz = cbuffer.data(pf_geom_10_off + 80 * acomps + 5);

        auto g_z_0_z_yyy = cbuffer.data(pf_geom_10_off + 80 * acomps + 6);

        auto g_z_0_z_yyz = cbuffer.data(pf_geom_10_off + 80 * acomps + 7);

        auto g_z_0_z_yzz = cbuffer.data(pf_geom_10_off + 80 * acomps + 8);

        auto g_z_0_z_zzz = cbuffer.data(pf_geom_10_off + 80 * acomps + 9);

        /// set up bra offset for contr_buffer_xxdd

        const auto dd_geom_10_off = idx_geom_10_xdd + i * 36;

        /// Set up 0-6 components of targeted buffer : cbuffer.data(

        auto g_x_0_xx_xx = cbuffer.data(dd_geom_10_off + 0 * acomps  + 0);

        auto g_x_0_xx_xy = cbuffer.data(dd_geom_10_off + 0 * acomps  + 1);

        auto g_x_0_xx_xz = cbuffer.data(dd_geom_10_off + 0 * acomps  + 2);

        auto g_x_0_xx_yy = cbuffer.data(dd_geom_10_off + 0 * acomps  + 3);

        auto g_x_0_xx_yz = cbuffer.data(dd_geom_10_off + 0 * acomps  + 4);

        auto g_x_0_xx_zz = cbuffer.data(dd_geom_10_off + 0 * acomps  + 5);

        #pragma omp simd aligned(cd_x, g_x_0_x_xx, g_x_0_x_xxx, g_x_0_x_xxy, g_x_0_x_xxz, g_x_0_x_xy, g_x_0_x_xyy, g_x_0_x_xyz, g_x_0_x_xz, g_x_0_x_xzz, g_x_0_x_yy, g_x_0_x_yz, g_x_0_x_zz, g_x_0_xx_xx, g_x_0_xx_xy, g_x_0_xx_xz, g_x_0_xx_yy, g_x_0_xx_yz, g_x_0_xx_zz, g_x_xx, g_x_xy, g_x_xz, g_x_yy, g_x_yz, g_x_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xx_xx[k] = -g_x_xx[k] - g_x_0_x_xx[k] * cd_x[k] + g_x_0_x_xxx[k];

            g_x_0_xx_xy[k] = -g_x_xy[k] - g_x_0_x_xy[k] * cd_x[k] + g_x_0_x_xxy[k];

            g_x_0_xx_xz[k] = -g_x_xz[k] - g_x_0_x_xz[k] * cd_x[k] + g_x_0_x_xxz[k];

            g_x_0_xx_yy[k] = -g_x_yy[k] - g_x_0_x_yy[k] * cd_x[k] + g_x_0_x_xyy[k];

            g_x_0_xx_yz[k] = -g_x_yz[k] - g_x_0_x_yz[k] * cd_x[k] + g_x_0_x_xyz[k];

            g_x_0_xx_zz[k] = -g_x_zz[k] - g_x_0_x_zz[k] * cd_x[k] + g_x_0_x_xzz[k];
        }

        /// Set up 6-12 components of targeted buffer : cbuffer.data(

        auto g_x_0_xy_xx = cbuffer.data(dd_geom_10_off + 0 * acomps  + 6);

        auto g_x_0_xy_xy = cbuffer.data(dd_geom_10_off + 0 * acomps  + 7);

        auto g_x_0_xy_xz = cbuffer.data(dd_geom_10_off + 0 * acomps  + 8);

        auto g_x_0_xy_yy = cbuffer.data(dd_geom_10_off + 0 * acomps  + 9);

        auto g_x_0_xy_yz = cbuffer.data(dd_geom_10_off + 0 * acomps  + 10);

        auto g_x_0_xy_zz = cbuffer.data(dd_geom_10_off + 0 * acomps  + 11);

        #pragma omp simd aligned(cd_y, g_x_0_x_xx, g_x_0_x_xxy, g_x_0_x_xy, g_x_0_x_xyy, g_x_0_x_xyz, g_x_0_x_xz, g_x_0_x_yy, g_x_0_x_yyy, g_x_0_x_yyz, g_x_0_x_yz, g_x_0_x_yzz, g_x_0_x_zz, g_x_0_xy_xx, g_x_0_xy_xy, g_x_0_xy_xz, g_x_0_xy_yy, g_x_0_xy_yz, g_x_0_xy_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xy_xx[k] = -g_x_0_x_xx[k] * cd_y[k] + g_x_0_x_xxy[k];

            g_x_0_xy_xy[k] = -g_x_0_x_xy[k] * cd_y[k] + g_x_0_x_xyy[k];

            g_x_0_xy_xz[k] = -g_x_0_x_xz[k] * cd_y[k] + g_x_0_x_xyz[k];

            g_x_0_xy_yy[k] = -g_x_0_x_yy[k] * cd_y[k] + g_x_0_x_yyy[k];

            g_x_0_xy_yz[k] = -g_x_0_x_yz[k] * cd_y[k] + g_x_0_x_yyz[k];

            g_x_0_xy_zz[k] = -g_x_0_x_zz[k] * cd_y[k] + g_x_0_x_yzz[k];
        }

        /// Set up 12-18 components of targeted buffer : cbuffer.data(

        auto g_x_0_xz_xx = cbuffer.data(dd_geom_10_off + 0 * acomps  + 12);

        auto g_x_0_xz_xy = cbuffer.data(dd_geom_10_off + 0 * acomps  + 13);

        auto g_x_0_xz_xz = cbuffer.data(dd_geom_10_off + 0 * acomps  + 14);

        auto g_x_0_xz_yy = cbuffer.data(dd_geom_10_off + 0 * acomps  + 15);

        auto g_x_0_xz_yz = cbuffer.data(dd_geom_10_off + 0 * acomps  + 16);

        auto g_x_0_xz_zz = cbuffer.data(dd_geom_10_off + 0 * acomps  + 17);

        #pragma omp simd aligned(cd_z, g_x_0_x_xx, g_x_0_x_xxz, g_x_0_x_xy, g_x_0_x_xyz, g_x_0_x_xz, g_x_0_x_xzz, g_x_0_x_yy, g_x_0_x_yyz, g_x_0_x_yz, g_x_0_x_yzz, g_x_0_x_zz, g_x_0_x_zzz, g_x_0_xz_xx, g_x_0_xz_xy, g_x_0_xz_xz, g_x_0_xz_yy, g_x_0_xz_yz, g_x_0_xz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xz_xx[k] = -g_x_0_x_xx[k] * cd_z[k] + g_x_0_x_xxz[k];

            g_x_0_xz_xy[k] = -g_x_0_x_xy[k] * cd_z[k] + g_x_0_x_xyz[k];

            g_x_0_xz_xz[k] = -g_x_0_x_xz[k] * cd_z[k] + g_x_0_x_xzz[k];

            g_x_0_xz_yy[k] = -g_x_0_x_yy[k] * cd_z[k] + g_x_0_x_yyz[k];

            g_x_0_xz_yz[k] = -g_x_0_x_yz[k] * cd_z[k] + g_x_0_x_yzz[k];

            g_x_0_xz_zz[k] = -g_x_0_x_zz[k] * cd_z[k] + g_x_0_x_zzz[k];
        }

        /// Set up 18-24 components of targeted buffer : cbuffer.data(

        auto g_x_0_yy_xx = cbuffer.data(dd_geom_10_off + 0 * acomps  + 18);

        auto g_x_0_yy_xy = cbuffer.data(dd_geom_10_off + 0 * acomps  + 19);

        auto g_x_0_yy_xz = cbuffer.data(dd_geom_10_off + 0 * acomps  + 20);

        auto g_x_0_yy_yy = cbuffer.data(dd_geom_10_off + 0 * acomps  + 21);

        auto g_x_0_yy_yz = cbuffer.data(dd_geom_10_off + 0 * acomps  + 22);

        auto g_x_0_yy_zz = cbuffer.data(dd_geom_10_off + 0 * acomps  + 23);

        #pragma omp simd aligned(cd_y, g_x_0_y_xx, g_x_0_y_xxy, g_x_0_y_xy, g_x_0_y_xyy, g_x_0_y_xyz, g_x_0_y_xz, g_x_0_y_yy, g_x_0_y_yyy, g_x_0_y_yyz, g_x_0_y_yz, g_x_0_y_yzz, g_x_0_y_zz, g_x_0_yy_xx, g_x_0_yy_xy, g_x_0_yy_xz, g_x_0_yy_yy, g_x_0_yy_yz, g_x_0_yy_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yy_xx[k] = -g_x_0_y_xx[k] * cd_y[k] + g_x_0_y_xxy[k];

            g_x_0_yy_xy[k] = -g_x_0_y_xy[k] * cd_y[k] + g_x_0_y_xyy[k];

            g_x_0_yy_xz[k] = -g_x_0_y_xz[k] * cd_y[k] + g_x_0_y_xyz[k];

            g_x_0_yy_yy[k] = -g_x_0_y_yy[k] * cd_y[k] + g_x_0_y_yyy[k];

            g_x_0_yy_yz[k] = -g_x_0_y_yz[k] * cd_y[k] + g_x_0_y_yyz[k];

            g_x_0_yy_zz[k] = -g_x_0_y_zz[k] * cd_y[k] + g_x_0_y_yzz[k];
        }

        /// Set up 24-30 components of targeted buffer : cbuffer.data(

        auto g_x_0_yz_xx = cbuffer.data(dd_geom_10_off + 0 * acomps  + 24);

        auto g_x_0_yz_xy = cbuffer.data(dd_geom_10_off + 0 * acomps  + 25);

        auto g_x_0_yz_xz = cbuffer.data(dd_geom_10_off + 0 * acomps  + 26);

        auto g_x_0_yz_yy = cbuffer.data(dd_geom_10_off + 0 * acomps  + 27);

        auto g_x_0_yz_yz = cbuffer.data(dd_geom_10_off + 0 * acomps  + 28);

        auto g_x_0_yz_zz = cbuffer.data(dd_geom_10_off + 0 * acomps  + 29);

        #pragma omp simd aligned(cd_y, g_x_0_yz_xx, g_x_0_yz_xy, g_x_0_yz_xz, g_x_0_yz_yy, g_x_0_yz_yz, g_x_0_yz_zz, g_x_0_z_xx, g_x_0_z_xxy, g_x_0_z_xy, g_x_0_z_xyy, g_x_0_z_xyz, g_x_0_z_xz, g_x_0_z_yy, g_x_0_z_yyy, g_x_0_z_yyz, g_x_0_z_yz, g_x_0_z_yzz, g_x_0_z_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yz_xx[k] = -g_x_0_z_xx[k] * cd_y[k] + g_x_0_z_xxy[k];

            g_x_0_yz_xy[k] = -g_x_0_z_xy[k] * cd_y[k] + g_x_0_z_xyy[k];

            g_x_0_yz_xz[k] = -g_x_0_z_xz[k] * cd_y[k] + g_x_0_z_xyz[k];

            g_x_0_yz_yy[k] = -g_x_0_z_yy[k] * cd_y[k] + g_x_0_z_yyy[k];

            g_x_0_yz_yz[k] = -g_x_0_z_yz[k] * cd_y[k] + g_x_0_z_yyz[k];

            g_x_0_yz_zz[k] = -g_x_0_z_zz[k] * cd_y[k] + g_x_0_z_yzz[k];
        }

        /// Set up 30-36 components of targeted buffer : cbuffer.data(

        auto g_x_0_zz_xx = cbuffer.data(dd_geom_10_off + 0 * acomps  + 30);

        auto g_x_0_zz_xy = cbuffer.data(dd_geom_10_off + 0 * acomps  + 31);

        auto g_x_0_zz_xz = cbuffer.data(dd_geom_10_off + 0 * acomps  + 32);

        auto g_x_0_zz_yy = cbuffer.data(dd_geom_10_off + 0 * acomps  + 33);

        auto g_x_0_zz_yz = cbuffer.data(dd_geom_10_off + 0 * acomps  + 34);

        auto g_x_0_zz_zz = cbuffer.data(dd_geom_10_off + 0 * acomps  + 35);

        #pragma omp simd aligned(cd_z, g_x_0_z_xx, g_x_0_z_xxz, g_x_0_z_xy, g_x_0_z_xyz, g_x_0_z_xz, g_x_0_z_xzz, g_x_0_z_yy, g_x_0_z_yyz, g_x_0_z_yz, g_x_0_z_yzz, g_x_0_z_zz, g_x_0_z_zzz, g_x_0_zz_xx, g_x_0_zz_xy, g_x_0_zz_xz, g_x_0_zz_yy, g_x_0_zz_yz, g_x_0_zz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_zz_xx[k] = -g_x_0_z_xx[k] * cd_z[k] + g_x_0_z_xxz[k];

            g_x_0_zz_xy[k] = -g_x_0_z_xy[k] * cd_z[k] + g_x_0_z_xyz[k];

            g_x_0_zz_xz[k] = -g_x_0_z_xz[k] * cd_z[k] + g_x_0_z_xzz[k];

            g_x_0_zz_yy[k] = -g_x_0_z_yy[k] * cd_z[k] + g_x_0_z_yyz[k];

            g_x_0_zz_yz[k] = -g_x_0_z_yz[k] * cd_z[k] + g_x_0_z_yzz[k];

            g_x_0_zz_zz[k] = -g_x_0_z_zz[k] * cd_z[k] + g_x_0_z_zzz[k];
        }
        /// Set up 0-6 components of targeted buffer : cbuffer.data(

        auto g_y_0_xx_xx = cbuffer.data(dd_geom_10_off + 36 * acomps  + 0);

        auto g_y_0_xx_xy = cbuffer.data(dd_geom_10_off + 36 * acomps  + 1);

        auto g_y_0_xx_xz = cbuffer.data(dd_geom_10_off + 36 * acomps  + 2);

        auto g_y_0_xx_yy = cbuffer.data(dd_geom_10_off + 36 * acomps  + 3);

        auto g_y_0_xx_yz = cbuffer.data(dd_geom_10_off + 36 * acomps  + 4);

        auto g_y_0_xx_zz = cbuffer.data(dd_geom_10_off + 36 * acomps  + 5);

        #pragma omp simd aligned(cd_x, g_y_0_x_xx, g_y_0_x_xxx, g_y_0_x_xxy, g_y_0_x_xxz, g_y_0_x_xy, g_y_0_x_xyy, g_y_0_x_xyz, g_y_0_x_xz, g_y_0_x_xzz, g_y_0_x_yy, g_y_0_x_yz, g_y_0_x_zz, g_y_0_xx_xx, g_y_0_xx_xy, g_y_0_xx_xz, g_y_0_xx_yy, g_y_0_xx_yz, g_y_0_xx_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xx_xx[k] = -g_y_0_x_xx[k] * cd_x[k] + g_y_0_x_xxx[k];

            g_y_0_xx_xy[k] = -g_y_0_x_xy[k] * cd_x[k] + g_y_0_x_xxy[k];

            g_y_0_xx_xz[k] = -g_y_0_x_xz[k] * cd_x[k] + g_y_0_x_xxz[k];

            g_y_0_xx_yy[k] = -g_y_0_x_yy[k] * cd_x[k] + g_y_0_x_xyy[k];

            g_y_0_xx_yz[k] = -g_y_0_x_yz[k] * cd_x[k] + g_y_0_x_xyz[k];

            g_y_0_xx_zz[k] = -g_y_0_x_zz[k] * cd_x[k] + g_y_0_x_xzz[k];
        }

        /// Set up 6-12 components of targeted buffer : cbuffer.data(

        auto g_y_0_xy_xx = cbuffer.data(dd_geom_10_off + 36 * acomps  + 6);

        auto g_y_0_xy_xy = cbuffer.data(dd_geom_10_off + 36 * acomps  + 7);

        auto g_y_0_xy_xz = cbuffer.data(dd_geom_10_off + 36 * acomps  + 8);

        auto g_y_0_xy_yy = cbuffer.data(dd_geom_10_off + 36 * acomps  + 9);

        auto g_y_0_xy_yz = cbuffer.data(dd_geom_10_off + 36 * acomps  + 10);

        auto g_y_0_xy_zz = cbuffer.data(dd_geom_10_off + 36 * acomps  + 11);

        #pragma omp simd aligned(cd_x, g_y_0_xy_xx, g_y_0_xy_xy, g_y_0_xy_xz, g_y_0_xy_yy, g_y_0_xy_yz, g_y_0_xy_zz, g_y_0_y_xx, g_y_0_y_xxx, g_y_0_y_xxy, g_y_0_y_xxz, g_y_0_y_xy, g_y_0_y_xyy, g_y_0_y_xyz, g_y_0_y_xz, g_y_0_y_xzz, g_y_0_y_yy, g_y_0_y_yz, g_y_0_y_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xy_xx[k] = -g_y_0_y_xx[k] * cd_x[k] + g_y_0_y_xxx[k];

            g_y_0_xy_xy[k] = -g_y_0_y_xy[k] * cd_x[k] + g_y_0_y_xxy[k];

            g_y_0_xy_xz[k] = -g_y_0_y_xz[k] * cd_x[k] + g_y_0_y_xxz[k];

            g_y_0_xy_yy[k] = -g_y_0_y_yy[k] * cd_x[k] + g_y_0_y_xyy[k];

            g_y_0_xy_yz[k] = -g_y_0_y_yz[k] * cd_x[k] + g_y_0_y_xyz[k];

            g_y_0_xy_zz[k] = -g_y_0_y_zz[k] * cd_x[k] + g_y_0_y_xzz[k];
        }

        /// Set up 12-18 components of targeted buffer : cbuffer.data(

        auto g_y_0_xz_xx = cbuffer.data(dd_geom_10_off + 36 * acomps  + 12);

        auto g_y_0_xz_xy = cbuffer.data(dd_geom_10_off + 36 * acomps  + 13);

        auto g_y_0_xz_xz = cbuffer.data(dd_geom_10_off + 36 * acomps  + 14);

        auto g_y_0_xz_yy = cbuffer.data(dd_geom_10_off + 36 * acomps  + 15);

        auto g_y_0_xz_yz = cbuffer.data(dd_geom_10_off + 36 * acomps  + 16);

        auto g_y_0_xz_zz = cbuffer.data(dd_geom_10_off + 36 * acomps  + 17);

        #pragma omp simd aligned(cd_x, g_y_0_xz_xx, g_y_0_xz_xy, g_y_0_xz_xz, g_y_0_xz_yy, g_y_0_xz_yz, g_y_0_xz_zz, g_y_0_z_xx, g_y_0_z_xxx, g_y_0_z_xxy, g_y_0_z_xxz, g_y_0_z_xy, g_y_0_z_xyy, g_y_0_z_xyz, g_y_0_z_xz, g_y_0_z_xzz, g_y_0_z_yy, g_y_0_z_yz, g_y_0_z_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xz_xx[k] = -g_y_0_z_xx[k] * cd_x[k] + g_y_0_z_xxx[k];

            g_y_0_xz_xy[k] = -g_y_0_z_xy[k] * cd_x[k] + g_y_0_z_xxy[k];

            g_y_0_xz_xz[k] = -g_y_0_z_xz[k] * cd_x[k] + g_y_0_z_xxz[k];

            g_y_0_xz_yy[k] = -g_y_0_z_yy[k] * cd_x[k] + g_y_0_z_xyy[k];

            g_y_0_xz_yz[k] = -g_y_0_z_yz[k] * cd_x[k] + g_y_0_z_xyz[k];

            g_y_0_xz_zz[k] = -g_y_0_z_zz[k] * cd_x[k] + g_y_0_z_xzz[k];
        }

        /// Set up 18-24 components of targeted buffer : cbuffer.data(

        auto g_y_0_yy_xx = cbuffer.data(dd_geom_10_off + 36 * acomps  + 18);

        auto g_y_0_yy_xy = cbuffer.data(dd_geom_10_off + 36 * acomps  + 19);

        auto g_y_0_yy_xz = cbuffer.data(dd_geom_10_off + 36 * acomps  + 20);

        auto g_y_0_yy_yy = cbuffer.data(dd_geom_10_off + 36 * acomps  + 21);

        auto g_y_0_yy_yz = cbuffer.data(dd_geom_10_off + 36 * acomps  + 22);

        auto g_y_0_yy_zz = cbuffer.data(dd_geom_10_off + 36 * acomps  + 23);

        #pragma omp simd aligned(cd_y, g_y_0_y_xx, g_y_0_y_xxy, g_y_0_y_xy, g_y_0_y_xyy, g_y_0_y_xyz, g_y_0_y_xz, g_y_0_y_yy, g_y_0_y_yyy, g_y_0_y_yyz, g_y_0_y_yz, g_y_0_y_yzz, g_y_0_y_zz, g_y_0_yy_xx, g_y_0_yy_xy, g_y_0_yy_xz, g_y_0_yy_yy, g_y_0_yy_yz, g_y_0_yy_zz, g_y_xx, g_y_xy, g_y_xz, g_y_yy, g_y_yz, g_y_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yy_xx[k] = -g_y_xx[k] - g_y_0_y_xx[k] * cd_y[k] + g_y_0_y_xxy[k];

            g_y_0_yy_xy[k] = -g_y_xy[k] - g_y_0_y_xy[k] * cd_y[k] + g_y_0_y_xyy[k];

            g_y_0_yy_xz[k] = -g_y_xz[k] - g_y_0_y_xz[k] * cd_y[k] + g_y_0_y_xyz[k];

            g_y_0_yy_yy[k] = -g_y_yy[k] - g_y_0_y_yy[k] * cd_y[k] + g_y_0_y_yyy[k];

            g_y_0_yy_yz[k] = -g_y_yz[k] - g_y_0_y_yz[k] * cd_y[k] + g_y_0_y_yyz[k];

            g_y_0_yy_zz[k] = -g_y_zz[k] - g_y_0_y_zz[k] * cd_y[k] + g_y_0_y_yzz[k];
        }

        /// Set up 24-30 components of targeted buffer : cbuffer.data(

        auto g_y_0_yz_xx = cbuffer.data(dd_geom_10_off + 36 * acomps  + 24);

        auto g_y_0_yz_xy = cbuffer.data(dd_geom_10_off + 36 * acomps  + 25);

        auto g_y_0_yz_xz = cbuffer.data(dd_geom_10_off + 36 * acomps  + 26);

        auto g_y_0_yz_yy = cbuffer.data(dd_geom_10_off + 36 * acomps  + 27);

        auto g_y_0_yz_yz = cbuffer.data(dd_geom_10_off + 36 * acomps  + 28);

        auto g_y_0_yz_zz = cbuffer.data(dd_geom_10_off + 36 * acomps  + 29);

        #pragma omp simd aligned(cd_z, g_y_0_y_xx, g_y_0_y_xxz, g_y_0_y_xy, g_y_0_y_xyz, g_y_0_y_xz, g_y_0_y_xzz, g_y_0_y_yy, g_y_0_y_yyz, g_y_0_y_yz, g_y_0_y_yzz, g_y_0_y_zz, g_y_0_y_zzz, g_y_0_yz_xx, g_y_0_yz_xy, g_y_0_yz_xz, g_y_0_yz_yy, g_y_0_yz_yz, g_y_0_yz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yz_xx[k] = -g_y_0_y_xx[k] * cd_z[k] + g_y_0_y_xxz[k];

            g_y_0_yz_xy[k] = -g_y_0_y_xy[k] * cd_z[k] + g_y_0_y_xyz[k];

            g_y_0_yz_xz[k] = -g_y_0_y_xz[k] * cd_z[k] + g_y_0_y_xzz[k];

            g_y_0_yz_yy[k] = -g_y_0_y_yy[k] * cd_z[k] + g_y_0_y_yyz[k];

            g_y_0_yz_yz[k] = -g_y_0_y_yz[k] * cd_z[k] + g_y_0_y_yzz[k];

            g_y_0_yz_zz[k] = -g_y_0_y_zz[k] * cd_z[k] + g_y_0_y_zzz[k];
        }

        /// Set up 30-36 components of targeted buffer : cbuffer.data(

        auto g_y_0_zz_xx = cbuffer.data(dd_geom_10_off + 36 * acomps  + 30);

        auto g_y_0_zz_xy = cbuffer.data(dd_geom_10_off + 36 * acomps  + 31);

        auto g_y_0_zz_xz = cbuffer.data(dd_geom_10_off + 36 * acomps  + 32);

        auto g_y_0_zz_yy = cbuffer.data(dd_geom_10_off + 36 * acomps  + 33);

        auto g_y_0_zz_yz = cbuffer.data(dd_geom_10_off + 36 * acomps  + 34);

        auto g_y_0_zz_zz = cbuffer.data(dd_geom_10_off + 36 * acomps  + 35);

        #pragma omp simd aligned(cd_z, g_y_0_z_xx, g_y_0_z_xxz, g_y_0_z_xy, g_y_0_z_xyz, g_y_0_z_xz, g_y_0_z_xzz, g_y_0_z_yy, g_y_0_z_yyz, g_y_0_z_yz, g_y_0_z_yzz, g_y_0_z_zz, g_y_0_z_zzz, g_y_0_zz_xx, g_y_0_zz_xy, g_y_0_zz_xz, g_y_0_zz_yy, g_y_0_zz_yz, g_y_0_zz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_zz_xx[k] = -g_y_0_z_xx[k] * cd_z[k] + g_y_0_z_xxz[k];

            g_y_0_zz_xy[k] = -g_y_0_z_xy[k] * cd_z[k] + g_y_0_z_xyz[k];

            g_y_0_zz_xz[k] = -g_y_0_z_xz[k] * cd_z[k] + g_y_0_z_xzz[k];

            g_y_0_zz_yy[k] = -g_y_0_z_yy[k] * cd_z[k] + g_y_0_z_yyz[k];

            g_y_0_zz_yz[k] = -g_y_0_z_yz[k] * cd_z[k] + g_y_0_z_yzz[k];

            g_y_0_zz_zz[k] = -g_y_0_z_zz[k] * cd_z[k] + g_y_0_z_zzz[k];
        }
        /// Set up 0-6 components of targeted buffer : cbuffer.data(

        auto g_z_0_xx_xx = cbuffer.data(dd_geom_10_off + 72 * acomps  + 0);

        auto g_z_0_xx_xy = cbuffer.data(dd_geom_10_off + 72 * acomps  + 1);

        auto g_z_0_xx_xz = cbuffer.data(dd_geom_10_off + 72 * acomps  + 2);

        auto g_z_0_xx_yy = cbuffer.data(dd_geom_10_off + 72 * acomps  + 3);

        auto g_z_0_xx_yz = cbuffer.data(dd_geom_10_off + 72 * acomps  + 4);

        auto g_z_0_xx_zz = cbuffer.data(dd_geom_10_off + 72 * acomps  + 5);

        #pragma omp simd aligned(cd_x, g_z_0_x_xx, g_z_0_x_xxx, g_z_0_x_xxy, g_z_0_x_xxz, g_z_0_x_xy, g_z_0_x_xyy, g_z_0_x_xyz, g_z_0_x_xz, g_z_0_x_xzz, g_z_0_x_yy, g_z_0_x_yz, g_z_0_x_zz, g_z_0_xx_xx, g_z_0_xx_xy, g_z_0_xx_xz, g_z_0_xx_yy, g_z_0_xx_yz, g_z_0_xx_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xx_xx[k] = -g_z_0_x_xx[k] * cd_x[k] + g_z_0_x_xxx[k];

            g_z_0_xx_xy[k] = -g_z_0_x_xy[k] * cd_x[k] + g_z_0_x_xxy[k];

            g_z_0_xx_xz[k] = -g_z_0_x_xz[k] * cd_x[k] + g_z_0_x_xxz[k];

            g_z_0_xx_yy[k] = -g_z_0_x_yy[k] * cd_x[k] + g_z_0_x_xyy[k];

            g_z_0_xx_yz[k] = -g_z_0_x_yz[k] * cd_x[k] + g_z_0_x_xyz[k];

            g_z_0_xx_zz[k] = -g_z_0_x_zz[k] * cd_x[k] + g_z_0_x_xzz[k];
        }

        /// Set up 6-12 components of targeted buffer : cbuffer.data(

        auto g_z_0_xy_xx = cbuffer.data(dd_geom_10_off + 72 * acomps  + 6);

        auto g_z_0_xy_xy = cbuffer.data(dd_geom_10_off + 72 * acomps  + 7);

        auto g_z_0_xy_xz = cbuffer.data(dd_geom_10_off + 72 * acomps  + 8);

        auto g_z_0_xy_yy = cbuffer.data(dd_geom_10_off + 72 * acomps  + 9);

        auto g_z_0_xy_yz = cbuffer.data(dd_geom_10_off + 72 * acomps  + 10);

        auto g_z_0_xy_zz = cbuffer.data(dd_geom_10_off + 72 * acomps  + 11);

        #pragma omp simd aligned(cd_x, g_z_0_xy_xx, g_z_0_xy_xy, g_z_0_xy_xz, g_z_0_xy_yy, g_z_0_xy_yz, g_z_0_xy_zz, g_z_0_y_xx, g_z_0_y_xxx, g_z_0_y_xxy, g_z_0_y_xxz, g_z_0_y_xy, g_z_0_y_xyy, g_z_0_y_xyz, g_z_0_y_xz, g_z_0_y_xzz, g_z_0_y_yy, g_z_0_y_yz, g_z_0_y_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xy_xx[k] = -g_z_0_y_xx[k] * cd_x[k] + g_z_0_y_xxx[k];

            g_z_0_xy_xy[k] = -g_z_0_y_xy[k] * cd_x[k] + g_z_0_y_xxy[k];

            g_z_0_xy_xz[k] = -g_z_0_y_xz[k] * cd_x[k] + g_z_0_y_xxz[k];

            g_z_0_xy_yy[k] = -g_z_0_y_yy[k] * cd_x[k] + g_z_0_y_xyy[k];

            g_z_0_xy_yz[k] = -g_z_0_y_yz[k] * cd_x[k] + g_z_0_y_xyz[k];

            g_z_0_xy_zz[k] = -g_z_0_y_zz[k] * cd_x[k] + g_z_0_y_xzz[k];
        }

        /// Set up 12-18 components of targeted buffer : cbuffer.data(

        auto g_z_0_xz_xx = cbuffer.data(dd_geom_10_off + 72 * acomps  + 12);

        auto g_z_0_xz_xy = cbuffer.data(dd_geom_10_off + 72 * acomps  + 13);

        auto g_z_0_xz_xz = cbuffer.data(dd_geom_10_off + 72 * acomps  + 14);

        auto g_z_0_xz_yy = cbuffer.data(dd_geom_10_off + 72 * acomps  + 15);

        auto g_z_0_xz_yz = cbuffer.data(dd_geom_10_off + 72 * acomps  + 16);

        auto g_z_0_xz_zz = cbuffer.data(dd_geom_10_off + 72 * acomps  + 17);

        #pragma omp simd aligned(cd_x, g_z_0_xz_xx, g_z_0_xz_xy, g_z_0_xz_xz, g_z_0_xz_yy, g_z_0_xz_yz, g_z_0_xz_zz, g_z_0_z_xx, g_z_0_z_xxx, g_z_0_z_xxy, g_z_0_z_xxz, g_z_0_z_xy, g_z_0_z_xyy, g_z_0_z_xyz, g_z_0_z_xz, g_z_0_z_xzz, g_z_0_z_yy, g_z_0_z_yz, g_z_0_z_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xz_xx[k] = -g_z_0_z_xx[k] * cd_x[k] + g_z_0_z_xxx[k];

            g_z_0_xz_xy[k] = -g_z_0_z_xy[k] * cd_x[k] + g_z_0_z_xxy[k];

            g_z_0_xz_xz[k] = -g_z_0_z_xz[k] * cd_x[k] + g_z_0_z_xxz[k];

            g_z_0_xz_yy[k] = -g_z_0_z_yy[k] * cd_x[k] + g_z_0_z_xyy[k];

            g_z_0_xz_yz[k] = -g_z_0_z_yz[k] * cd_x[k] + g_z_0_z_xyz[k];

            g_z_0_xz_zz[k] = -g_z_0_z_zz[k] * cd_x[k] + g_z_0_z_xzz[k];
        }

        /// Set up 18-24 components of targeted buffer : cbuffer.data(

        auto g_z_0_yy_xx = cbuffer.data(dd_geom_10_off + 72 * acomps  + 18);

        auto g_z_0_yy_xy = cbuffer.data(dd_geom_10_off + 72 * acomps  + 19);

        auto g_z_0_yy_xz = cbuffer.data(dd_geom_10_off + 72 * acomps  + 20);

        auto g_z_0_yy_yy = cbuffer.data(dd_geom_10_off + 72 * acomps  + 21);

        auto g_z_0_yy_yz = cbuffer.data(dd_geom_10_off + 72 * acomps  + 22);

        auto g_z_0_yy_zz = cbuffer.data(dd_geom_10_off + 72 * acomps  + 23);

        #pragma omp simd aligned(cd_y, g_z_0_y_xx, g_z_0_y_xxy, g_z_0_y_xy, g_z_0_y_xyy, g_z_0_y_xyz, g_z_0_y_xz, g_z_0_y_yy, g_z_0_y_yyy, g_z_0_y_yyz, g_z_0_y_yz, g_z_0_y_yzz, g_z_0_y_zz, g_z_0_yy_xx, g_z_0_yy_xy, g_z_0_yy_xz, g_z_0_yy_yy, g_z_0_yy_yz, g_z_0_yy_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yy_xx[k] = -g_z_0_y_xx[k] * cd_y[k] + g_z_0_y_xxy[k];

            g_z_0_yy_xy[k] = -g_z_0_y_xy[k] * cd_y[k] + g_z_0_y_xyy[k];

            g_z_0_yy_xz[k] = -g_z_0_y_xz[k] * cd_y[k] + g_z_0_y_xyz[k];

            g_z_0_yy_yy[k] = -g_z_0_y_yy[k] * cd_y[k] + g_z_0_y_yyy[k];

            g_z_0_yy_yz[k] = -g_z_0_y_yz[k] * cd_y[k] + g_z_0_y_yyz[k];

            g_z_0_yy_zz[k] = -g_z_0_y_zz[k] * cd_y[k] + g_z_0_y_yzz[k];
        }

        /// Set up 24-30 components of targeted buffer : cbuffer.data(

        auto g_z_0_yz_xx = cbuffer.data(dd_geom_10_off + 72 * acomps  + 24);

        auto g_z_0_yz_xy = cbuffer.data(dd_geom_10_off + 72 * acomps  + 25);

        auto g_z_0_yz_xz = cbuffer.data(dd_geom_10_off + 72 * acomps  + 26);

        auto g_z_0_yz_yy = cbuffer.data(dd_geom_10_off + 72 * acomps  + 27);

        auto g_z_0_yz_yz = cbuffer.data(dd_geom_10_off + 72 * acomps  + 28);

        auto g_z_0_yz_zz = cbuffer.data(dd_geom_10_off + 72 * acomps  + 29);

        #pragma omp simd aligned(cd_y, g_z_0_yz_xx, g_z_0_yz_xy, g_z_0_yz_xz, g_z_0_yz_yy, g_z_0_yz_yz, g_z_0_yz_zz, g_z_0_z_xx, g_z_0_z_xxy, g_z_0_z_xy, g_z_0_z_xyy, g_z_0_z_xyz, g_z_0_z_xz, g_z_0_z_yy, g_z_0_z_yyy, g_z_0_z_yyz, g_z_0_z_yz, g_z_0_z_yzz, g_z_0_z_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yz_xx[k] = -g_z_0_z_xx[k] * cd_y[k] + g_z_0_z_xxy[k];

            g_z_0_yz_xy[k] = -g_z_0_z_xy[k] * cd_y[k] + g_z_0_z_xyy[k];

            g_z_0_yz_xz[k] = -g_z_0_z_xz[k] * cd_y[k] + g_z_0_z_xyz[k];

            g_z_0_yz_yy[k] = -g_z_0_z_yy[k] * cd_y[k] + g_z_0_z_yyy[k];

            g_z_0_yz_yz[k] = -g_z_0_z_yz[k] * cd_y[k] + g_z_0_z_yyz[k];

            g_z_0_yz_zz[k] = -g_z_0_z_zz[k] * cd_y[k] + g_z_0_z_yzz[k];
        }

        /// Set up 30-36 components of targeted buffer : cbuffer.data(

        auto g_z_0_zz_xx = cbuffer.data(dd_geom_10_off + 72 * acomps  + 30);

        auto g_z_0_zz_xy = cbuffer.data(dd_geom_10_off + 72 * acomps  + 31);

        auto g_z_0_zz_xz = cbuffer.data(dd_geom_10_off + 72 * acomps  + 32);

        auto g_z_0_zz_yy = cbuffer.data(dd_geom_10_off + 72 * acomps  + 33);

        auto g_z_0_zz_yz = cbuffer.data(dd_geom_10_off + 72 * acomps  + 34);

        auto g_z_0_zz_zz = cbuffer.data(dd_geom_10_off + 72 * acomps  + 35);

        #pragma omp simd aligned(cd_z, g_z_0_z_xx, g_z_0_z_xxz, g_z_0_z_xy, g_z_0_z_xyz, g_z_0_z_xz, g_z_0_z_xzz, g_z_0_z_yy, g_z_0_z_yyz, g_z_0_z_yz, g_z_0_z_yzz, g_z_0_z_zz, g_z_0_z_zzz, g_z_0_zz_xx, g_z_0_zz_xy, g_z_0_zz_xz, g_z_0_zz_yy, g_z_0_zz_yz, g_z_0_zz_zz, g_z_xx, g_z_xy, g_z_xz, g_z_yy, g_z_yz, g_z_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_zz_xx[k] = -g_z_xx[k] - g_z_0_z_xx[k] * cd_z[k] + g_z_0_z_xxz[k];

            g_z_0_zz_xy[k] = -g_z_xy[k] - g_z_0_z_xy[k] * cd_z[k] + g_z_0_z_xyz[k];

            g_z_0_zz_xz[k] = -g_z_xz[k] - g_z_0_z_xz[k] * cd_z[k] + g_z_0_z_xzz[k];

            g_z_0_zz_yy[k] = -g_z_yy[k] - g_z_0_z_yy[k] * cd_z[k] + g_z_0_z_yyz[k];

            g_z_0_zz_yz[k] = -g_z_yz[k] - g_z_0_z_yz[k] * cd_z[k] + g_z_0_z_yzz[k];

            g_z_0_zz_zz[k] = -g_z_zz[k] - g_z_0_z_zz[k] * cd_z[k] + g_z_0_z_zzz[k];
        }
    }
}

} // t3ceri namespace

