#include "ThreeCenterElectronRepulsionGeom100ContrRecFXX.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_bra_geom1_electron_repulsion_fxx(CSimdArray<double>& cbuffer,
                                      const size_t idx_geom_100_fxx,
                                      const size_t idx_dxx,
                                      const size_t idx_gxx,
                                      const int c_angmom,
                                      const int d_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto ccomps = tensor::number_of_cartesian_components(std::array<int, 1>{c_angmom,});

    const auto dcomps = tensor::number_of_cartesian_components(std::array<int, 1>{d_angmom,});

    for (int i = 0; i < ccomps; i++)
    {
        for (int j = 0; j < dcomps; j++)
        {
            /// Set up components of auxilary buffer : DSS

            const auto d_off = idx_dxx + i * dcomps + j;

            auto g_xx_0 = cbuffer.data(d_off + 0 * ccomps * dcomps);

            auto g_xy_0 = cbuffer.data(d_off + 1 * ccomps * dcomps);

            auto g_xz_0 = cbuffer.data(d_off + 2 * ccomps * dcomps);

            auto g_yy_0 = cbuffer.data(d_off + 3 * ccomps * dcomps);

            auto g_yz_0 = cbuffer.data(d_off + 4 * ccomps * dcomps);

            auto g_zz_0 = cbuffer.data(d_off + 5 * ccomps * dcomps);

            /// Set up components of auxilary buffer : GSS

            const auto g_off = idx_gxx + i * dcomps + j;

            auto g_xxxx_0 = cbuffer.data(g_off + 0 * ccomps * dcomps);

            auto g_xxxy_0 = cbuffer.data(g_off + 1 * ccomps * dcomps);

            auto g_xxxz_0 = cbuffer.data(g_off + 2 * ccomps * dcomps);

            auto g_xxyy_0 = cbuffer.data(g_off + 3 * ccomps * dcomps);

            auto g_xxyz_0 = cbuffer.data(g_off + 4 * ccomps * dcomps);

            auto g_xxzz_0 = cbuffer.data(g_off + 5 * ccomps * dcomps);

            auto g_xyyy_0 = cbuffer.data(g_off + 6 * ccomps * dcomps);

            auto g_xyyz_0 = cbuffer.data(g_off + 7 * ccomps * dcomps);

            auto g_xyzz_0 = cbuffer.data(g_off + 8 * ccomps * dcomps);

            auto g_xzzz_0 = cbuffer.data(g_off + 9 * ccomps * dcomps);

            auto g_yyyy_0 = cbuffer.data(g_off + 10 * ccomps * dcomps);

            auto g_yyyz_0 = cbuffer.data(g_off + 11 * ccomps * dcomps);

            auto g_yyzz_0 = cbuffer.data(g_off + 12 * ccomps * dcomps);

            auto g_yzzz_0 = cbuffer.data(g_off + 13 * ccomps * dcomps);

            auto g_zzzz_0 = cbuffer.data(g_off + 14 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_fxx

            const auto f_geom_100_off = idx_geom_100_fxx + i * dcomps + j;

            /// Set up 0-10 components of targeted buffer : cbuffer.data(

            auto g_x_0_0_xxx_0 = cbuffer.data(f_geom_100_off + 0 * ccomps * dcomps);

            auto g_x_0_0_xxy_0 = cbuffer.data(f_geom_100_off + 1 * ccomps * dcomps);

            auto g_x_0_0_xxz_0 = cbuffer.data(f_geom_100_off + 2 * ccomps * dcomps);

            auto g_x_0_0_xyy_0 = cbuffer.data(f_geom_100_off + 3 * ccomps * dcomps);

            auto g_x_0_0_xyz_0 = cbuffer.data(f_geom_100_off + 4 * ccomps * dcomps);

            auto g_x_0_0_xzz_0 = cbuffer.data(f_geom_100_off + 5 * ccomps * dcomps);

            auto g_x_0_0_yyy_0 = cbuffer.data(f_geom_100_off + 6 * ccomps * dcomps);

            auto g_x_0_0_yyz_0 = cbuffer.data(f_geom_100_off + 7 * ccomps * dcomps);

            auto g_x_0_0_yzz_0 = cbuffer.data(f_geom_100_off + 8 * ccomps * dcomps);

            auto g_x_0_0_zzz_0 = cbuffer.data(f_geom_100_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_0_xxx_0, g_x_0_0_xxy_0, g_x_0_0_xxz_0, g_x_0_0_xyy_0, g_x_0_0_xyz_0, g_x_0_0_xzz_0, g_x_0_0_yyy_0, g_x_0_0_yyz_0, g_x_0_0_yzz_0, g_x_0_0_zzz_0, g_xx_0, g_xxxx_0, g_xxxy_0, g_xxxz_0, g_xxyy_0, g_xxyz_0, g_xxzz_0, g_xy_0, g_xyyy_0, g_xyyz_0, g_xyzz_0, g_xz_0, g_xzzz_0, g_yy_0, g_yz_0, g_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_0_xxx_0[k] = -3.0 * g_xx_0[k] + g_xxxx_0[k];

                g_x_0_0_xxy_0[k] = -2.0 * g_xy_0[k] + g_xxxy_0[k];

                g_x_0_0_xxz_0[k] = -2.0 * g_xz_0[k] + g_xxxz_0[k];

                g_x_0_0_xyy_0[k] = -g_yy_0[k] + g_xxyy_0[k];

                g_x_0_0_xyz_0[k] = -g_yz_0[k] + g_xxyz_0[k];

                g_x_0_0_xzz_0[k] = -g_zz_0[k] + g_xxzz_0[k];

                g_x_0_0_yyy_0[k] = g_xyyy_0[k];

                g_x_0_0_yyz_0[k] = g_xyyz_0[k];

                g_x_0_0_yzz_0[k] = g_xyzz_0[k];

                g_x_0_0_zzz_0[k] = g_xzzz_0[k];
            }

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

            auto g_y_0_0_xxx_0 = cbuffer.data(f_geom_100_off + 10 * ccomps * dcomps);

            auto g_y_0_0_xxy_0 = cbuffer.data(f_geom_100_off + 11 * ccomps * dcomps);

            auto g_y_0_0_xxz_0 = cbuffer.data(f_geom_100_off + 12 * ccomps * dcomps);

            auto g_y_0_0_xyy_0 = cbuffer.data(f_geom_100_off + 13 * ccomps * dcomps);

            auto g_y_0_0_xyz_0 = cbuffer.data(f_geom_100_off + 14 * ccomps * dcomps);

            auto g_y_0_0_xzz_0 = cbuffer.data(f_geom_100_off + 15 * ccomps * dcomps);

            auto g_y_0_0_yyy_0 = cbuffer.data(f_geom_100_off + 16 * ccomps * dcomps);

            auto g_y_0_0_yyz_0 = cbuffer.data(f_geom_100_off + 17 * ccomps * dcomps);

            auto g_y_0_0_yzz_0 = cbuffer.data(f_geom_100_off + 18 * ccomps * dcomps);

            auto g_y_0_0_zzz_0 = cbuffer.data(f_geom_100_off + 19 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0, g_xxxy_0, g_xxyy_0, g_xxyz_0, g_xy_0, g_xyyy_0, g_xyyz_0, g_xyzz_0, g_xz_0, g_y_0_0_xxx_0, g_y_0_0_xxy_0, g_y_0_0_xxz_0, g_y_0_0_xyy_0, g_y_0_0_xyz_0, g_y_0_0_xzz_0, g_y_0_0_yyy_0, g_y_0_0_yyz_0, g_y_0_0_yzz_0, g_y_0_0_zzz_0, g_yy_0, g_yyyy_0, g_yyyz_0, g_yyzz_0, g_yz_0, g_yzzz_0, g_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_0_xxx_0[k] = g_xxxy_0[k];

                g_y_0_0_xxy_0[k] = -g_xx_0[k] + g_xxyy_0[k];

                g_y_0_0_xxz_0[k] = g_xxyz_0[k];

                g_y_0_0_xyy_0[k] = -2.0 * g_xy_0[k] + g_xyyy_0[k];

                g_y_0_0_xyz_0[k] = -g_xz_0[k] + g_xyyz_0[k];

                g_y_0_0_xzz_0[k] = g_xyzz_0[k];

                g_y_0_0_yyy_0[k] = -3.0 * g_yy_0[k] + g_yyyy_0[k];

                g_y_0_0_yyz_0[k] = -2.0 * g_yz_0[k] + g_yyyz_0[k];

                g_y_0_0_yzz_0[k] = -g_zz_0[k] + g_yyzz_0[k];

                g_y_0_0_zzz_0[k] = g_yzzz_0[k];
            }

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

            auto g_z_0_0_xxx_0 = cbuffer.data(f_geom_100_off + 20 * ccomps * dcomps);

            auto g_z_0_0_xxy_0 = cbuffer.data(f_geom_100_off + 21 * ccomps * dcomps);

            auto g_z_0_0_xxz_0 = cbuffer.data(f_geom_100_off + 22 * ccomps * dcomps);

            auto g_z_0_0_xyy_0 = cbuffer.data(f_geom_100_off + 23 * ccomps * dcomps);

            auto g_z_0_0_xyz_0 = cbuffer.data(f_geom_100_off + 24 * ccomps * dcomps);

            auto g_z_0_0_xzz_0 = cbuffer.data(f_geom_100_off + 25 * ccomps * dcomps);

            auto g_z_0_0_yyy_0 = cbuffer.data(f_geom_100_off + 26 * ccomps * dcomps);

            auto g_z_0_0_yyz_0 = cbuffer.data(f_geom_100_off + 27 * ccomps * dcomps);

            auto g_z_0_0_yzz_0 = cbuffer.data(f_geom_100_off + 28 * ccomps * dcomps);

            auto g_z_0_0_zzz_0 = cbuffer.data(f_geom_100_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0, g_xxxz_0, g_xxyz_0, g_xxzz_0, g_xy_0, g_xyyz_0, g_xyzz_0, g_xz_0, g_xzzz_0, g_yy_0, g_yyyz_0, g_yyzz_0, g_yz_0, g_yzzz_0, g_z_0_0_xxx_0, g_z_0_0_xxy_0, g_z_0_0_xxz_0, g_z_0_0_xyy_0, g_z_0_0_xyz_0, g_z_0_0_xzz_0, g_z_0_0_yyy_0, g_z_0_0_yyz_0, g_z_0_0_yzz_0, g_z_0_0_zzz_0, g_zz_0, g_zzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_0_xxx_0[k] = g_xxxz_0[k];

                g_z_0_0_xxy_0[k] = g_xxyz_0[k];

                g_z_0_0_xxz_0[k] = -g_xx_0[k] + g_xxzz_0[k];

                g_z_0_0_xyy_0[k] = g_xyyz_0[k];

                g_z_0_0_xyz_0[k] = -g_xy_0[k] + g_xyzz_0[k];

                g_z_0_0_xzz_0[k] = -2.0 * g_xz_0[k] + g_xzzz_0[k];

                g_z_0_0_yyy_0[k] = g_yyyz_0[k];

                g_z_0_0_yyz_0[k] = -g_yy_0[k] + g_yyzz_0[k];

                g_z_0_0_yzz_0[k] = -2.0 * g_yz_0[k] + g_yzzz_0[k];

                g_z_0_0_zzz_0[k] = -3.0 * g_zz_0[k] + g_zzzz_0[k];
            }
        }
    }
}

} // t3ceri namespace

