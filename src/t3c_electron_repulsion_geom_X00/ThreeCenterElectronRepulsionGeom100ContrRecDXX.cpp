#include "ThreeCenterElectronRepulsionGeom100ContrRecDXX.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_bra_geom1_electron_repulsion_dxx(CSimdArray<double>& cbuffer,
                                      const size_t idx_geom_100_dxx,
                                      const size_t idx_pxx,
                                      const size_t idx_fxx,
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
            /// Set up components of auxilary buffer : PSS

            const auto p_off = idx_pxx + i * dcomps + j;

            auto g_x_0 = cbuffer.data(p_off + 0 * ccomps * dcomps);

            auto g_y_0 = cbuffer.data(p_off + 1 * ccomps * dcomps);

            auto g_z_0 = cbuffer.data(p_off + 2 * ccomps * dcomps);

            /// Set up components of auxilary buffer : FSS

            const auto f_off = idx_fxx + i * dcomps + j;

            auto g_xxx_0 = cbuffer.data(f_off + 0 * ccomps * dcomps);

            auto g_xxy_0 = cbuffer.data(f_off + 1 * ccomps * dcomps);

            auto g_xxz_0 = cbuffer.data(f_off + 2 * ccomps * dcomps);

            auto g_xyy_0 = cbuffer.data(f_off + 3 * ccomps * dcomps);

            auto g_xyz_0 = cbuffer.data(f_off + 4 * ccomps * dcomps);

            auto g_xzz_0 = cbuffer.data(f_off + 5 * ccomps * dcomps);

            auto g_yyy_0 = cbuffer.data(f_off + 6 * ccomps * dcomps);

            auto g_yyz_0 = cbuffer.data(f_off + 7 * ccomps * dcomps);

            auto g_yzz_0 = cbuffer.data(f_off + 8 * ccomps * dcomps);

            auto g_zzz_0 = cbuffer.data(f_off + 9 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_dxx

            const auto d_geom_100_off = idx_geom_100_dxx + i * dcomps + j;

            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_x_0_0_xx_0 = cbuffer.data(d_geom_100_off + 0 * ccomps * dcomps);

            auto g_x_0_0_xy_0 = cbuffer.data(d_geom_100_off + 1 * ccomps * dcomps);

            auto g_x_0_0_xz_0 = cbuffer.data(d_geom_100_off + 2 * ccomps * dcomps);

            auto g_x_0_0_yy_0 = cbuffer.data(d_geom_100_off + 3 * ccomps * dcomps);

            auto g_x_0_0_yz_0 = cbuffer.data(d_geom_100_off + 4 * ccomps * dcomps);

            auto g_x_0_0_zz_0 = cbuffer.data(d_geom_100_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0, g_x_0_0_xx_0, g_x_0_0_xy_0, g_x_0_0_xz_0, g_x_0_0_yy_0, g_x_0_0_yz_0, g_x_0_0_zz_0, g_xxx_0, g_xxy_0, g_xxz_0, g_xyy_0, g_xyz_0, g_xzz_0, g_y_0, g_z_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_0_xx_0[k] = -2.0 * g_x_0[k] + g_xxx_0[k];

                g_x_0_0_xy_0[k] = -g_y_0[k] + g_xxy_0[k];

                g_x_0_0_xz_0[k] = -g_z_0[k] + g_xxz_0[k];

                g_x_0_0_yy_0[k] = g_xyy_0[k];

                g_x_0_0_yz_0[k] = g_xyz_0[k];

                g_x_0_0_zz_0[k] = g_xzz_0[k];
            }

            /// Set up 6-12 components of targeted buffer : cbuffer.data(

            auto g_y_0_0_xx_0 = cbuffer.data(d_geom_100_off + 6 * ccomps * dcomps);

            auto g_y_0_0_xy_0 = cbuffer.data(d_geom_100_off + 7 * ccomps * dcomps);

            auto g_y_0_0_xz_0 = cbuffer.data(d_geom_100_off + 8 * ccomps * dcomps);

            auto g_y_0_0_yy_0 = cbuffer.data(d_geom_100_off + 9 * ccomps * dcomps);

            auto g_y_0_0_yz_0 = cbuffer.data(d_geom_100_off + 10 * ccomps * dcomps);

            auto g_y_0_0_zz_0 = cbuffer.data(d_geom_100_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0, g_xxy_0, g_xyy_0, g_xyz_0, g_y_0, g_y_0_0_xx_0, g_y_0_0_xy_0, g_y_0_0_xz_0, g_y_0_0_yy_0, g_y_0_0_yz_0, g_y_0_0_zz_0, g_yyy_0, g_yyz_0, g_yzz_0, g_z_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_0_xx_0[k] = g_xxy_0[k];

                g_y_0_0_xy_0[k] = -g_x_0[k] + g_xyy_0[k];

                g_y_0_0_xz_0[k] = g_xyz_0[k];

                g_y_0_0_yy_0[k] = -2.0 * g_y_0[k] + g_yyy_0[k];

                g_y_0_0_yz_0[k] = -g_z_0[k] + g_yyz_0[k];

                g_y_0_0_zz_0[k] = g_yzz_0[k];
            }

            /// Set up 12-18 components of targeted buffer : cbuffer.data(

            auto g_z_0_0_xx_0 = cbuffer.data(d_geom_100_off + 12 * ccomps * dcomps);

            auto g_z_0_0_xy_0 = cbuffer.data(d_geom_100_off + 13 * ccomps * dcomps);

            auto g_z_0_0_xz_0 = cbuffer.data(d_geom_100_off + 14 * ccomps * dcomps);

            auto g_z_0_0_yy_0 = cbuffer.data(d_geom_100_off + 15 * ccomps * dcomps);

            auto g_z_0_0_yz_0 = cbuffer.data(d_geom_100_off + 16 * ccomps * dcomps);

            auto g_z_0_0_zz_0 = cbuffer.data(d_geom_100_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0, g_xxz_0, g_xyz_0, g_xzz_0, g_y_0, g_yyz_0, g_yzz_0, g_z_0, g_z_0_0_xx_0, g_z_0_0_xy_0, g_z_0_0_xz_0, g_z_0_0_yy_0, g_z_0_0_yz_0, g_z_0_0_zz_0, g_zzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_0_xx_0[k] = g_xxz_0[k];

                g_z_0_0_xy_0[k] = g_xyz_0[k];

                g_z_0_0_xz_0[k] = -g_x_0[k] + g_xzz_0[k];

                g_z_0_0_yy_0[k] = g_yyz_0[k];

                g_z_0_0_yz_0[k] = -g_y_0[k] + g_yzz_0[k];

                g_z_0_0_zz_0[k] = -2.0 * g_z_0[k] + g_zzz_0[k];
            }
        }
    }
}

} // t3ceri namespace

