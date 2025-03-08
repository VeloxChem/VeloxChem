#include "ThreeCenterElectronRepulsionGeom100ContrRecGXX.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_bra_geom1_electron_repulsion_gxx(CSimdArray<double>& cbuffer,
                                      const size_t idx_geom_100_gxx,
                                      const size_t idx_fxx,
                                      const size_t idx_hxx,
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

            /// Set up components of auxilary buffer : HSS

            const auto h_off = idx_hxx + i * dcomps + j;

            auto g_xxxxx_0 = cbuffer.data(h_off + 0 * ccomps * dcomps);

            auto g_xxxxy_0 = cbuffer.data(h_off + 1 * ccomps * dcomps);

            auto g_xxxxz_0 = cbuffer.data(h_off + 2 * ccomps * dcomps);

            auto g_xxxyy_0 = cbuffer.data(h_off + 3 * ccomps * dcomps);

            auto g_xxxyz_0 = cbuffer.data(h_off + 4 * ccomps * dcomps);

            auto g_xxxzz_0 = cbuffer.data(h_off + 5 * ccomps * dcomps);

            auto g_xxyyy_0 = cbuffer.data(h_off + 6 * ccomps * dcomps);

            auto g_xxyyz_0 = cbuffer.data(h_off + 7 * ccomps * dcomps);

            auto g_xxyzz_0 = cbuffer.data(h_off + 8 * ccomps * dcomps);

            auto g_xxzzz_0 = cbuffer.data(h_off + 9 * ccomps * dcomps);

            auto g_xyyyy_0 = cbuffer.data(h_off + 10 * ccomps * dcomps);

            auto g_xyyyz_0 = cbuffer.data(h_off + 11 * ccomps * dcomps);

            auto g_xyyzz_0 = cbuffer.data(h_off + 12 * ccomps * dcomps);

            auto g_xyzzz_0 = cbuffer.data(h_off + 13 * ccomps * dcomps);

            auto g_xzzzz_0 = cbuffer.data(h_off + 14 * ccomps * dcomps);

            auto g_yyyyy_0 = cbuffer.data(h_off + 15 * ccomps * dcomps);

            auto g_yyyyz_0 = cbuffer.data(h_off + 16 * ccomps * dcomps);

            auto g_yyyzz_0 = cbuffer.data(h_off + 17 * ccomps * dcomps);

            auto g_yyzzz_0 = cbuffer.data(h_off + 18 * ccomps * dcomps);

            auto g_yzzzz_0 = cbuffer.data(h_off + 19 * ccomps * dcomps);

            auto g_zzzzz_0 = cbuffer.data(h_off + 20 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_gxx

            const auto g_geom_100_off = idx_geom_100_gxx + i * dcomps + j;

            /// Set up 0-15 components of targeted buffer : cbuffer.data(

            auto g_x_0_0_xxxx_0 = cbuffer.data(g_geom_100_off + 0 * ccomps * dcomps);

            auto g_x_0_0_xxxy_0 = cbuffer.data(g_geom_100_off + 1 * ccomps * dcomps);

            auto g_x_0_0_xxxz_0 = cbuffer.data(g_geom_100_off + 2 * ccomps * dcomps);

            auto g_x_0_0_xxyy_0 = cbuffer.data(g_geom_100_off + 3 * ccomps * dcomps);

            auto g_x_0_0_xxyz_0 = cbuffer.data(g_geom_100_off + 4 * ccomps * dcomps);

            auto g_x_0_0_xxzz_0 = cbuffer.data(g_geom_100_off + 5 * ccomps * dcomps);

            auto g_x_0_0_xyyy_0 = cbuffer.data(g_geom_100_off + 6 * ccomps * dcomps);

            auto g_x_0_0_xyyz_0 = cbuffer.data(g_geom_100_off + 7 * ccomps * dcomps);

            auto g_x_0_0_xyzz_0 = cbuffer.data(g_geom_100_off + 8 * ccomps * dcomps);

            auto g_x_0_0_xzzz_0 = cbuffer.data(g_geom_100_off + 9 * ccomps * dcomps);

            auto g_x_0_0_yyyy_0 = cbuffer.data(g_geom_100_off + 10 * ccomps * dcomps);

            auto g_x_0_0_yyyz_0 = cbuffer.data(g_geom_100_off + 11 * ccomps * dcomps);

            auto g_x_0_0_yyzz_0 = cbuffer.data(g_geom_100_off + 12 * ccomps * dcomps);

            auto g_x_0_0_yzzz_0 = cbuffer.data(g_geom_100_off + 13 * ccomps * dcomps);

            auto g_x_0_0_zzzz_0 = cbuffer.data(g_geom_100_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_0_xxxx_0, g_x_0_0_xxxy_0, g_x_0_0_xxxz_0, g_x_0_0_xxyy_0, g_x_0_0_xxyz_0, g_x_0_0_xxzz_0, g_x_0_0_xyyy_0, g_x_0_0_xyyz_0, g_x_0_0_xyzz_0, g_x_0_0_xzzz_0, g_x_0_0_yyyy_0, g_x_0_0_yyyz_0, g_x_0_0_yyzz_0, g_x_0_0_yzzz_0, g_x_0_0_zzzz_0, g_xxx_0, g_xxxxx_0, g_xxxxy_0, g_xxxxz_0, g_xxxyy_0, g_xxxyz_0, g_xxxzz_0, g_xxy_0, g_xxyyy_0, g_xxyyz_0, g_xxyzz_0, g_xxz_0, g_xxzzz_0, g_xyy_0, g_xyyyy_0, g_xyyyz_0, g_xyyzz_0, g_xyz_0, g_xyzzz_0, g_xzz_0, g_xzzzz_0, g_yyy_0, g_yyz_0, g_yzz_0, g_zzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_0_xxxx_0[k] = -4.0 * g_xxx_0[k] + g_xxxxx_0[k];

                g_x_0_0_xxxy_0[k] = -3.0 * g_xxy_0[k] + g_xxxxy_0[k];

                g_x_0_0_xxxz_0[k] = -3.0 * g_xxz_0[k] + g_xxxxz_0[k];

                g_x_0_0_xxyy_0[k] = -2.0 * g_xyy_0[k] + g_xxxyy_0[k];

                g_x_0_0_xxyz_0[k] = -2.0 * g_xyz_0[k] + g_xxxyz_0[k];

                g_x_0_0_xxzz_0[k] = -2.0 * g_xzz_0[k] + g_xxxzz_0[k];

                g_x_0_0_xyyy_0[k] = -g_yyy_0[k] + g_xxyyy_0[k];

                g_x_0_0_xyyz_0[k] = -g_yyz_0[k] + g_xxyyz_0[k];

                g_x_0_0_xyzz_0[k] = -g_yzz_0[k] + g_xxyzz_0[k];

                g_x_0_0_xzzz_0[k] = -g_zzz_0[k] + g_xxzzz_0[k];

                g_x_0_0_yyyy_0[k] = g_xyyyy_0[k];

                g_x_0_0_yyyz_0[k] = g_xyyyz_0[k];

                g_x_0_0_yyzz_0[k] = g_xyyzz_0[k];

                g_x_0_0_yzzz_0[k] = g_xyzzz_0[k];

                g_x_0_0_zzzz_0[k] = g_xzzzz_0[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

            auto g_y_0_0_xxxx_0 = cbuffer.data(g_geom_100_off + 15 * ccomps * dcomps);

            auto g_y_0_0_xxxy_0 = cbuffer.data(g_geom_100_off + 16 * ccomps * dcomps);

            auto g_y_0_0_xxxz_0 = cbuffer.data(g_geom_100_off + 17 * ccomps * dcomps);

            auto g_y_0_0_xxyy_0 = cbuffer.data(g_geom_100_off + 18 * ccomps * dcomps);

            auto g_y_0_0_xxyz_0 = cbuffer.data(g_geom_100_off + 19 * ccomps * dcomps);

            auto g_y_0_0_xxzz_0 = cbuffer.data(g_geom_100_off + 20 * ccomps * dcomps);

            auto g_y_0_0_xyyy_0 = cbuffer.data(g_geom_100_off + 21 * ccomps * dcomps);

            auto g_y_0_0_xyyz_0 = cbuffer.data(g_geom_100_off + 22 * ccomps * dcomps);

            auto g_y_0_0_xyzz_0 = cbuffer.data(g_geom_100_off + 23 * ccomps * dcomps);

            auto g_y_0_0_xzzz_0 = cbuffer.data(g_geom_100_off + 24 * ccomps * dcomps);

            auto g_y_0_0_yyyy_0 = cbuffer.data(g_geom_100_off + 25 * ccomps * dcomps);

            auto g_y_0_0_yyyz_0 = cbuffer.data(g_geom_100_off + 26 * ccomps * dcomps);

            auto g_y_0_0_yyzz_0 = cbuffer.data(g_geom_100_off + 27 * ccomps * dcomps);

            auto g_y_0_0_yzzz_0 = cbuffer.data(g_geom_100_off + 28 * ccomps * dcomps);

            auto g_y_0_0_zzzz_0 = cbuffer.data(g_geom_100_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_xxx_0, g_xxxxy_0, g_xxxyy_0, g_xxxyz_0, g_xxy_0, g_xxyyy_0, g_xxyyz_0, g_xxyzz_0, g_xxz_0, g_xyy_0, g_xyyyy_0, g_xyyyz_0, g_xyyzz_0, g_xyz_0, g_xyzzz_0, g_xzz_0, g_y_0_0_xxxx_0, g_y_0_0_xxxy_0, g_y_0_0_xxxz_0, g_y_0_0_xxyy_0, g_y_0_0_xxyz_0, g_y_0_0_xxzz_0, g_y_0_0_xyyy_0, g_y_0_0_xyyz_0, g_y_0_0_xyzz_0, g_y_0_0_xzzz_0, g_y_0_0_yyyy_0, g_y_0_0_yyyz_0, g_y_0_0_yyzz_0, g_y_0_0_yzzz_0, g_y_0_0_zzzz_0, g_yyy_0, g_yyyyy_0, g_yyyyz_0, g_yyyzz_0, g_yyz_0, g_yyzzz_0, g_yzz_0, g_yzzzz_0, g_zzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_0_xxxx_0[k] = g_xxxxy_0[k];

                g_y_0_0_xxxy_0[k] = -g_xxx_0[k] + g_xxxyy_0[k];

                g_y_0_0_xxxz_0[k] = g_xxxyz_0[k];

                g_y_0_0_xxyy_0[k] = -2.0 * g_xxy_0[k] + g_xxyyy_0[k];

                g_y_0_0_xxyz_0[k] = -g_xxz_0[k] + g_xxyyz_0[k];

                g_y_0_0_xxzz_0[k] = g_xxyzz_0[k];

                g_y_0_0_xyyy_0[k] = -3.0 * g_xyy_0[k] + g_xyyyy_0[k];

                g_y_0_0_xyyz_0[k] = -2.0 * g_xyz_0[k] + g_xyyyz_0[k];

                g_y_0_0_xyzz_0[k] = -g_xzz_0[k] + g_xyyzz_0[k];

                g_y_0_0_xzzz_0[k] = g_xyzzz_0[k];

                g_y_0_0_yyyy_0[k] = -4.0 * g_yyy_0[k] + g_yyyyy_0[k];

                g_y_0_0_yyyz_0[k] = -3.0 * g_yyz_0[k] + g_yyyyz_0[k];

                g_y_0_0_yyzz_0[k] = -2.0 * g_yzz_0[k] + g_yyyzz_0[k];

                g_y_0_0_yzzz_0[k] = -g_zzz_0[k] + g_yyzzz_0[k];

                g_y_0_0_zzzz_0[k] = g_yzzzz_0[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

            auto g_z_0_0_xxxx_0 = cbuffer.data(g_geom_100_off + 30 * ccomps * dcomps);

            auto g_z_0_0_xxxy_0 = cbuffer.data(g_geom_100_off + 31 * ccomps * dcomps);

            auto g_z_0_0_xxxz_0 = cbuffer.data(g_geom_100_off + 32 * ccomps * dcomps);

            auto g_z_0_0_xxyy_0 = cbuffer.data(g_geom_100_off + 33 * ccomps * dcomps);

            auto g_z_0_0_xxyz_0 = cbuffer.data(g_geom_100_off + 34 * ccomps * dcomps);

            auto g_z_0_0_xxzz_0 = cbuffer.data(g_geom_100_off + 35 * ccomps * dcomps);

            auto g_z_0_0_xyyy_0 = cbuffer.data(g_geom_100_off + 36 * ccomps * dcomps);

            auto g_z_0_0_xyyz_0 = cbuffer.data(g_geom_100_off + 37 * ccomps * dcomps);

            auto g_z_0_0_xyzz_0 = cbuffer.data(g_geom_100_off + 38 * ccomps * dcomps);

            auto g_z_0_0_xzzz_0 = cbuffer.data(g_geom_100_off + 39 * ccomps * dcomps);

            auto g_z_0_0_yyyy_0 = cbuffer.data(g_geom_100_off + 40 * ccomps * dcomps);

            auto g_z_0_0_yyyz_0 = cbuffer.data(g_geom_100_off + 41 * ccomps * dcomps);

            auto g_z_0_0_yyzz_0 = cbuffer.data(g_geom_100_off + 42 * ccomps * dcomps);

            auto g_z_0_0_yzzz_0 = cbuffer.data(g_geom_100_off + 43 * ccomps * dcomps);

            auto g_z_0_0_zzzz_0 = cbuffer.data(g_geom_100_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_xxx_0, g_xxxxz_0, g_xxxyz_0, g_xxxzz_0, g_xxy_0, g_xxyyz_0, g_xxyzz_0, g_xxz_0, g_xxzzz_0, g_xyy_0, g_xyyyz_0, g_xyyzz_0, g_xyz_0, g_xyzzz_0, g_xzz_0, g_xzzzz_0, g_yyy_0, g_yyyyz_0, g_yyyzz_0, g_yyz_0, g_yyzzz_0, g_yzz_0, g_yzzzz_0, g_z_0_0_xxxx_0, g_z_0_0_xxxy_0, g_z_0_0_xxxz_0, g_z_0_0_xxyy_0, g_z_0_0_xxyz_0, g_z_0_0_xxzz_0, g_z_0_0_xyyy_0, g_z_0_0_xyyz_0, g_z_0_0_xyzz_0, g_z_0_0_xzzz_0, g_z_0_0_yyyy_0, g_z_0_0_yyyz_0, g_z_0_0_yyzz_0, g_z_0_0_yzzz_0, g_z_0_0_zzzz_0, g_zzz_0, g_zzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_0_xxxx_0[k] = g_xxxxz_0[k];

                g_z_0_0_xxxy_0[k] = g_xxxyz_0[k];

                g_z_0_0_xxxz_0[k] = -g_xxx_0[k] + g_xxxzz_0[k];

                g_z_0_0_xxyy_0[k] = g_xxyyz_0[k];

                g_z_0_0_xxyz_0[k] = -g_xxy_0[k] + g_xxyzz_0[k];

                g_z_0_0_xxzz_0[k] = -2.0 * g_xxz_0[k] + g_xxzzz_0[k];

                g_z_0_0_xyyy_0[k] = g_xyyyz_0[k];

                g_z_0_0_xyyz_0[k] = -g_xyy_0[k] + g_xyyzz_0[k];

                g_z_0_0_xyzz_0[k] = -2.0 * g_xyz_0[k] + g_xyzzz_0[k];

                g_z_0_0_xzzz_0[k] = -3.0 * g_xzz_0[k] + g_xzzzz_0[k];

                g_z_0_0_yyyy_0[k] = g_yyyyz_0[k];

                g_z_0_0_yyyz_0[k] = -g_yyy_0[k] + g_yyyzz_0[k];

                g_z_0_0_yyzz_0[k] = -2.0 * g_yyz_0[k] + g_yyzzz_0[k];

                g_z_0_0_yzzz_0[k] = -3.0 * g_yzz_0[k] + g_yzzzz_0[k];

                g_z_0_0_zzzz_0[k] = -4.0 * g_zzz_0[k] + g_zzzzz_0[k];
            }
        }
    }
}

} // t3ceri namespace

