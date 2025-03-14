#include "ThreeCenterElectronRepulsionContrRecXPD.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_hrr_electron_repulsion_xpd(CSimdArray<double>& cbuffer,
                                const size_t idx_xpd,
                                const size_t idx_xsd,
                                const size_t idx_xsf,
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
        /// Set up components of auxilary buffer : SSD

        const auto sd_off = idx_xsd + i * 6;

        auto g_0_xx = cbuffer.data(sd_off + 0);

        auto g_0_xy = cbuffer.data(sd_off + 1);

        auto g_0_xz = cbuffer.data(sd_off + 2);

        auto g_0_yy = cbuffer.data(sd_off + 3);

        auto g_0_yz = cbuffer.data(sd_off + 4);

        auto g_0_zz = cbuffer.data(sd_off + 5);

        /// Set up components of auxilary buffer : SSF

        const auto sf_off = idx_xsf + i * 10;

        auto g_0_xxx = cbuffer.data(sf_off + 0);

        auto g_0_xxy = cbuffer.data(sf_off + 1);

        auto g_0_xxz = cbuffer.data(sf_off + 2);

        auto g_0_xyy = cbuffer.data(sf_off + 3);

        auto g_0_xyz = cbuffer.data(sf_off + 4);

        auto g_0_xzz = cbuffer.data(sf_off + 5);

        auto g_0_yyy = cbuffer.data(sf_off + 6);

        auto g_0_yyz = cbuffer.data(sf_off + 7);

        auto g_0_yzz = cbuffer.data(sf_off + 8);

        auto g_0_zzz = cbuffer.data(sf_off + 9);

        /// set up bra offset for contr_buffer_xpd

        const auto pd_off = idx_xpd + i * 18;

        /// Set up 0-6 components of targeted buffer : cbuffer.data(

        auto g_x_xx = cbuffer.data(pd_off + 0);

        auto g_x_xy = cbuffer.data(pd_off + 1);

        auto g_x_xz = cbuffer.data(pd_off + 2);

        auto g_x_yy = cbuffer.data(pd_off + 3);

        auto g_x_yz = cbuffer.data(pd_off + 4);

        auto g_x_zz = cbuffer.data(pd_off + 5);

        #pragma omp simd aligned(cd_x, g_0_xx, g_0_xxx, g_0_xxy, g_0_xxz, g_0_xy, g_0_xyy, g_0_xyz, g_0_xz, g_0_xzz, g_0_yy, g_0_yz, g_0_zz, g_x_xx, g_x_xy, g_x_xz, g_x_yy, g_x_yz, g_x_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_xx[k] = -g_0_xx[k] * cd_x[k] + g_0_xxx[k];

            g_x_xy[k] = -g_0_xy[k] * cd_x[k] + g_0_xxy[k];

            g_x_xz[k] = -g_0_xz[k] * cd_x[k] + g_0_xxz[k];

            g_x_yy[k] = -g_0_yy[k] * cd_x[k] + g_0_xyy[k];

            g_x_yz[k] = -g_0_yz[k] * cd_x[k] + g_0_xyz[k];

            g_x_zz[k] = -g_0_zz[k] * cd_x[k] + g_0_xzz[k];
        }

        /// Set up 6-12 components of targeted buffer : cbuffer.data(

        auto g_y_xx = cbuffer.data(pd_off + 6);

        auto g_y_xy = cbuffer.data(pd_off + 7);

        auto g_y_xz = cbuffer.data(pd_off + 8);

        auto g_y_yy = cbuffer.data(pd_off + 9);

        auto g_y_yz = cbuffer.data(pd_off + 10);

        auto g_y_zz = cbuffer.data(pd_off + 11);

        #pragma omp simd aligned(cd_y, g_0_xx, g_0_xxy, g_0_xy, g_0_xyy, g_0_xyz, g_0_xz, g_0_yy, g_0_yyy, g_0_yyz, g_0_yz, g_0_yzz, g_0_zz, g_y_xx, g_y_xy, g_y_xz, g_y_yy, g_y_yz, g_y_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_xx[k] = -g_0_xx[k] * cd_y[k] + g_0_xxy[k];

            g_y_xy[k] = -g_0_xy[k] * cd_y[k] + g_0_xyy[k];

            g_y_xz[k] = -g_0_xz[k] * cd_y[k] + g_0_xyz[k];

            g_y_yy[k] = -g_0_yy[k] * cd_y[k] + g_0_yyy[k];

            g_y_yz[k] = -g_0_yz[k] * cd_y[k] + g_0_yyz[k];

            g_y_zz[k] = -g_0_zz[k] * cd_y[k] + g_0_yzz[k];
        }

        /// Set up 12-18 components of targeted buffer : cbuffer.data(

        auto g_z_xx = cbuffer.data(pd_off + 12);

        auto g_z_xy = cbuffer.data(pd_off + 13);

        auto g_z_xz = cbuffer.data(pd_off + 14);

        auto g_z_yy = cbuffer.data(pd_off + 15);

        auto g_z_yz = cbuffer.data(pd_off + 16);

        auto g_z_zz = cbuffer.data(pd_off + 17);

        #pragma omp simd aligned(cd_z, g_0_xx, g_0_xxz, g_0_xy, g_0_xyz, g_0_xz, g_0_xzz, g_0_yy, g_0_yyz, g_0_yz, g_0_yzz, g_0_zz, g_0_zzz, g_z_xx, g_z_xy, g_z_xz, g_z_yy, g_z_yz, g_z_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_xx[k] = -g_0_xx[k] * cd_z[k] + g_0_xxz[k];

            g_z_xy[k] = -g_0_xy[k] * cd_z[k] + g_0_xyz[k];

            g_z_xz[k] = -g_0_xz[k] * cd_z[k] + g_0_xzz[k];

            g_z_yy[k] = -g_0_yy[k] * cd_z[k] + g_0_yyz[k];

            g_z_yz[k] = -g_0_yz[k] * cd_z[k] + g_0_yzz[k];

            g_z_zz[k] = -g_0_zz[k] * cd_z[k] + g_0_zzz[k];
        }
    }
}

} // t3ceri namespace

