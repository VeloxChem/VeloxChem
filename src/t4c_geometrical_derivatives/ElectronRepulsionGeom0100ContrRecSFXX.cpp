#include "ElectronRepulsionGeom0100ContrRecSFXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_sfxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_sfxx,
                                            const size_t idx_sdxx,
                                            const size_t idx_sgxx,
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
            /// Set up components of auxilary buffer : SDSS

            const auto sd_off = idx_sdxx + i * dcomps + j;

            auto g_0_xx = cbuffer.data(sd_off + 0 * ccomps * dcomps);

            auto g_0_xy = cbuffer.data(sd_off + 1 * ccomps * dcomps);

            auto g_0_xz = cbuffer.data(sd_off + 2 * ccomps * dcomps);

            auto g_0_yy = cbuffer.data(sd_off + 3 * ccomps * dcomps);

            auto g_0_yz = cbuffer.data(sd_off + 4 * ccomps * dcomps);

            auto g_0_zz = cbuffer.data(sd_off + 5 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SGSS

            const auto sg_off = idx_sgxx + i * dcomps + j;

            auto g_0_xxxx = cbuffer.data(sg_off + 0 * ccomps * dcomps);

            auto g_0_xxxy = cbuffer.data(sg_off + 1 * ccomps * dcomps);

            auto g_0_xxxz = cbuffer.data(sg_off + 2 * ccomps * dcomps);

            auto g_0_xxyy = cbuffer.data(sg_off + 3 * ccomps * dcomps);

            auto g_0_xxyz = cbuffer.data(sg_off + 4 * ccomps * dcomps);

            auto g_0_xxzz = cbuffer.data(sg_off + 5 * ccomps * dcomps);

            auto g_0_xyyy = cbuffer.data(sg_off + 6 * ccomps * dcomps);

            auto g_0_xyyz = cbuffer.data(sg_off + 7 * ccomps * dcomps);

            auto g_0_xyzz = cbuffer.data(sg_off + 8 * ccomps * dcomps);

            auto g_0_xzzz = cbuffer.data(sg_off + 9 * ccomps * dcomps);

            auto g_0_yyyy = cbuffer.data(sg_off + 10 * ccomps * dcomps);

            auto g_0_yyyz = cbuffer.data(sg_off + 11 * ccomps * dcomps);

            auto g_0_yyzz = cbuffer.data(sg_off + 12 * ccomps * dcomps);

            auto g_0_yzzz = cbuffer.data(sg_off + 13 * ccomps * dcomps);

            auto g_0_zzzz = cbuffer.data(sg_off + 14 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_sfxx

            const auto sf_geom_01_off = idx_geom_01_sfxx + i * dcomps + j;

            /// Set up 0-10 components of targeted buffer : cbuffer.data(

            auto g_0_x_0_xxx = cbuffer.data(sf_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_0_xxy = cbuffer.data(sf_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_0_xxz = cbuffer.data(sf_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_0_xyy = cbuffer.data(sf_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_0_xyz = cbuffer.data(sf_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_0_xzz = cbuffer.data(sf_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_0_yyy = cbuffer.data(sf_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_0_yyz = cbuffer.data(sf_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_0_yzz = cbuffer.data(sf_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_0_zzz = cbuffer.data(sf_geom_01_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xxx, g_0_x_0_xxy, g_0_x_0_xxz, g_0_x_0_xyy, g_0_x_0_xyz, g_0_x_0_xzz, g_0_x_0_yyy, g_0_x_0_yyz, g_0_x_0_yzz, g_0_x_0_zzz, g_0_xx, g_0_xxxx, g_0_xxxy, g_0_xxxz, g_0_xxyy, g_0_xxyz, g_0_xxzz, g_0_xy, g_0_xyyy, g_0_xyyz, g_0_xyzz, g_0_xz, g_0_xzzz, g_0_yy, g_0_yz, g_0_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_0_xxx[k] = -3.0 * g_0_xx[k] + g_0_xxxx[k];

                g_0_x_0_xxy[k] = -2.0 * g_0_xy[k] + g_0_xxxy[k];

                g_0_x_0_xxz[k] = -2.0 * g_0_xz[k] + g_0_xxxz[k];

                g_0_x_0_xyy[k] = -g_0_yy[k] + g_0_xxyy[k];

                g_0_x_0_xyz[k] = -g_0_yz[k] + g_0_xxyz[k];

                g_0_x_0_xzz[k] = -g_0_zz[k] + g_0_xxzz[k];

                g_0_x_0_yyy[k] = g_0_xyyy[k];

                g_0_x_0_yyz[k] = g_0_xyyz[k];

                g_0_x_0_yzz[k] = g_0_xyzz[k];

                g_0_x_0_zzz[k] = g_0_xzzz[k];
            }

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

            auto g_0_y_0_xxx = cbuffer.data(sf_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_y_0_xxy = cbuffer.data(sf_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_y_0_xxz = cbuffer.data(sf_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_y_0_xyy = cbuffer.data(sf_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_y_0_xyz = cbuffer.data(sf_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_y_0_xzz = cbuffer.data(sf_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_y_0_yyy = cbuffer.data(sf_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_y_0_yyz = cbuffer.data(sf_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_y_0_yzz = cbuffer.data(sf_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_y_0_zzz = cbuffer.data(sf_geom_01_off + 19 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_xx, g_0_xxxy, g_0_xxyy, g_0_xxyz, g_0_xy, g_0_xyyy, g_0_xyyz, g_0_xyzz, g_0_xz, g_0_y_0_xxx, g_0_y_0_xxy, g_0_y_0_xxz, g_0_y_0_xyy, g_0_y_0_xyz, g_0_y_0_xzz, g_0_y_0_yyy, g_0_y_0_yyz, g_0_y_0_yzz, g_0_y_0_zzz, g_0_yy, g_0_yyyy, g_0_yyyz, g_0_yyzz, g_0_yz, g_0_yzzz, g_0_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_0_xxx[k] = g_0_xxxy[k];

                g_0_y_0_xxy[k] = -g_0_xx[k] + g_0_xxyy[k];

                g_0_y_0_xxz[k] = g_0_xxyz[k];

                g_0_y_0_xyy[k] = -2.0 * g_0_xy[k] + g_0_xyyy[k];

                g_0_y_0_xyz[k] = -g_0_xz[k] + g_0_xyyz[k];

                g_0_y_0_xzz[k] = g_0_xyzz[k];

                g_0_y_0_yyy[k] = -3.0 * g_0_yy[k] + g_0_yyyy[k];

                g_0_y_0_yyz[k] = -2.0 * g_0_yz[k] + g_0_yyyz[k];

                g_0_y_0_yzz[k] = -g_0_zz[k] + g_0_yyzz[k];

                g_0_y_0_zzz[k] = g_0_yzzz[k];
            }

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

            auto g_0_z_0_xxx = cbuffer.data(sf_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_z_0_xxy = cbuffer.data(sf_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_z_0_xxz = cbuffer.data(sf_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_z_0_xyy = cbuffer.data(sf_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_z_0_xyz = cbuffer.data(sf_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_z_0_xzz = cbuffer.data(sf_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_z_0_yyy = cbuffer.data(sf_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_z_0_yyz = cbuffer.data(sf_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_z_0_yzz = cbuffer.data(sf_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_z_0_zzz = cbuffer.data(sf_geom_01_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_xx, g_0_xxxz, g_0_xxyz, g_0_xxzz, g_0_xy, g_0_xyyz, g_0_xyzz, g_0_xz, g_0_xzzz, g_0_yy, g_0_yyyz, g_0_yyzz, g_0_yz, g_0_yzzz, g_0_z_0_xxx, g_0_z_0_xxy, g_0_z_0_xxz, g_0_z_0_xyy, g_0_z_0_xyz, g_0_z_0_xzz, g_0_z_0_yyy, g_0_z_0_yyz, g_0_z_0_yzz, g_0_z_0_zzz, g_0_zz, g_0_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_0_xxx[k] = g_0_xxxz[k];

                g_0_z_0_xxy[k] = g_0_xxyz[k];

                g_0_z_0_xxz[k] = -g_0_xx[k] + g_0_xxzz[k];

                g_0_z_0_xyy[k] = g_0_xyyz[k];

                g_0_z_0_xyz[k] = -g_0_xy[k] + g_0_xyzz[k];

                g_0_z_0_xzz[k] = -2.0 * g_0_xz[k] + g_0_xzzz[k];

                g_0_z_0_yyy[k] = g_0_yyyz[k];

                g_0_z_0_yyz[k] = -g_0_yy[k] + g_0_yyzz[k];

                g_0_z_0_yzz[k] = -2.0 * g_0_yz[k] + g_0_yzzz[k];

                g_0_z_0_zzz[k] = -3.0 * g_0_zz[k] + g_0_zzzz[k];
            }
        }
    }
}

} // erirec namespace

