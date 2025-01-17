#include "ElectronRepulsionGeom0100ContrRecSDXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_sdxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_sdxx,
                                            const size_t idx_spxx,
                                            const size_t idx_sfxx,
                                            const int c_angmom,
                                            const int d_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom,});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom,});

    for (int i = 0; i < ccomps; i++)
    {
        for (int j = 0; j < dcomps; j++)
        {
            /// Set up components of auxilary buffer : SPSS

            const auto sp_off = idx_spxx + i * dcomps + j;

            auto g_0_x = cbuffer.data(sp_off + 0 * ccomps * dcomps);

            auto g_0_y = cbuffer.data(sp_off + 1 * ccomps * dcomps);

            auto g_0_z = cbuffer.data(sp_off + 2 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SFSS

            const auto sf_off = idx_sfxx + i * dcomps + j;

            auto g_0_xxx = cbuffer.data(sf_off + 0 * ccomps * dcomps);

            auto g_0_xxy = cbuffer.data(sf_off + 1 * ccomps * dcomps);

            auto g_0_xxz = cbuffer.data(sf_off + 2 * ccomps * dcomps);

            auto g_0_xyy = cbuffer.data(sf_off + 3 * ccomps * dcomps);

            auto g_0_xyz = cbuffer.data(sf_off + 4 * ccomps * dcomps);

            auto g_0_xzz = cbuffer.data(sf_off + 5 * ccomps * dcomps);

            auto g_0_yyy = cbuffer.data(sf_off + 6 * ccomps * dcomps);

            auto g_0_yyz = cbuffer.data(sf_off + 7 * ccomps * dcomps);

            auto g_0_yzz = cbuffer.data(sf_off + 8 * ccomps * dcomps);

            auto g_0_zzz = cbuffer.data(sf_off + 9 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_sdxx

            const auto sd_geom_01_off = idx_geom_01_sdxx + i * dcomps + j;

            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_0_x_0_xx = cbuffer.data(sd_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_0_xy = cbuffer.data(sd_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_0_xz = cbuffer.data(sd_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_0_yy = cbuffer.data(sd_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_0_yz = cbuffer.data(sd_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_0_zz = cbuffer.data(sd_geom_01_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x, g_0_x_0_xx, g_0_x_0_xy, g_0_x_0_xz, g_0_x_0_yy, g_0_x_0_yz, g_0_x_0_zz, g_0_xxx, g_0_xxy, g_0_xxz, g_0_xyy, g_0_xyz, g_0_xzz, g_0_y, g_0_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_0_xx[k] = -2.0 * g_0_x[k] + g_0_xxx[k];

                g_0_x_0_xy[k] = -g_0_y[k] + g_0_xxy[k];

                g_0_x_0_xz[k] = -g_0_z[k] + g_0_xxz[k];

                g_0_x_0_yy[k] = g_0_xyy[k];

                g_0_x_0_yz[k] = g_0_xyz[k];

                g_0_x_0_zz[k] = g_0_xzz[k];
            }

            /// Set up 6-12 components of targeted buffer : cbuffer.data(

            auto g_0_y_0_xx = cbuffer.data(sd_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_y_0_xy = cbuffer.data(sd_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_y_0_xz = cbuffer.data(sd_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_y_0_yy = cbuffer.data(sd_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_y_0_yz = cbuffer.data(sd_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_y_0_zz = cbuffer.data(sd_geom_01_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x, g_0_xxy, g_0_xyy, g_0_xyz, g_0_y, g_0_y_0_xx, g_0_y_0_xy, g_0_y_0_xz, g_0_y_0_yy, g_0_y_0_yz, g_0_y_0_zz, g_0_yyy, g_0_yyz, g_0_yzz, g_0_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_0_xx[k] = g_0_xxy[k];

                g_0_y_0_xy[k] = -g_0_x[k] + g_0_xyy[k];

                g_0_y_0_xz[k] = g_0_xyz[k];

                g_0_y_0_yy[k] = -2.0 * g_0_y[k] + g_0_yyy[k];

                g_0_y_0_yz[k] = -g_0_z[k] + g_0_yyz[k];

                g_0_y_0_zz[k] = g_0_yzz[k];
            }

            /// Set up 12-18 components of targeted buffer : cbuffer.data(

            auto g_0_z_0_xx = cbuffer.data(sd_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_z_0_xy = cbuffer.data(sd_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_z_0_xz = cbuffer.data(sd_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_z_0_yy = cbuffer.data(sd_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_z_0_yz = cbuffer.data(sd_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_z_0_zz = cbuffer.data(sd_geom_01_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x, g_0_xxz, g_0_xyz, g_0_xzz, g_0_y, g_0_yyz, g_0_yzz, g_0_z, g_0_z_0_xx, g_0_z_0_xy, g_0_z_0_xz, g_0_z_0_yy, g_0_z_0_yz, g_0_z_0_zz, g_0_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_0_xx[k] = g_0_xxz[k];

                g_0_z_0_xy[k] = g_0_xyz[k];

                g_0_z_0_xz[k] = -g_0_x[k] + g_0_xzz[k];

                g_0_z_0_yy[k] = g_0_yyz[k];

                g_0_z_0_yz[k] = -g_0_y[k] + g_0_yzz[k];

                g_0_z_0_zz[k] = -2.0 * g_0_z[k] + g_0_zzz[k];
            }
        }
    }
}

} // erirec namespace

