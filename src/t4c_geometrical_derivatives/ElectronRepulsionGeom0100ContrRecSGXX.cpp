#include "ElectronRepulsionGeom0100ContrRecSGXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_sgxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_sgxx,
                                            const size_t idx_sfxx,
                                            const size_t idx_shxx,
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

            /// Set up components of auxilary buffer : SHSS

            const auto sh_off = idx_shxx + i * dcomps + j;

            auto g_0_xxxxx = cbuffer.data(sh_off + 0 * ccomps * dcomps);

            auto g_0_xxxxy = cbuffer.data(sh_off + 1 * ccomps * dcomps);

            auto g_0_xxxxz = cbuffer.data(sh_off + 2 * ccomps * dcomps);

            auto g_0_xxxyy = cbuffer.data(sh_off + 3 * ccomps * dcomps);

            auto g_0_xxxyz = cbuffer.data(sh_off + 4 * ccomps * dcomps);

            auto g_0_xxxzz = cbuffer.data(sh_off + 5 * ccomps * dcomps);

            auto g_0_xxyyy = cbuffer.data(sh_off + 6 * ccomps * dcomps);

            auto g_0_xxyyz = cbuffer.data(sh_off + 7 * ccomps * dcomps);

            auto g_0_xxyzz = cbuffer.data(sh_off + 8 * ccomps * dcomps);

            auto g_0_xxzzz = cbuffer.data(sh_off + 9 * ccomps * dcomps);

            auto g_0_xyyyy = cbuffer.data(sh_off + 10 * ccomps * dcomps);

            auto g_0_xyyyz = cbuffer.data(sh_off + 11 * ccomps * dcomps);

            auto g_0_xyyzz = cbuffer.data(sh_off + 12 * ccomps * dcomps);

            auto g_0_xyzzz = cbuffer.data(sh_off + 13 * ccomps * dcomps);

            auto g_0_xzzzz = cbuffer.data(sh_off + 14 * ccomps * dcomps);

            auto g_0_yyyyy = cbuffer.data(sh_off + 15 * ccomps * dcomps);

            auto g_0_yyyyz = cbuffer.data(sh_off + 16 * ccomps * dcomps);

            auto g_0_yyyzz = cbuffer.data(sh_off + 17 * ccomps * dcomps);

            auto g_0_yyzzz = cbuffer.data(sh_off + 18 * ccomps * dcomps);

            auto g_0_yzzzz = cbuffer.data(sh_off + 19 * ccomps * dcomps);

            auto g_0_zzzzz = cbuffer.data(sh_off + 20 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_sgxx

            const auto sg_geom_01_off = idx_geom_01_sgxx + i * dcomps + j;

            /// Set up 0-15 components of targeted buffer : cbuffer.data(

            auto g_0_x_0_xxxx = cbuffer.data(sg_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_0_xxxy = cbuffer.data(sg_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_0_xxxz = cbuffer.data(sg_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_0_xxyy = cbuffer.data(sg_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_0_xxyz = cbuffer.data(sg_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_0_xxzz = cbuffer.data(sg_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_0_xyyy = cbuffer.data(sg_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_0_xyyz = cbuffer.data(sg_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_0_xyzz = cbuffer.data(sg_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_0_xzzz = cbuffer.data(sg_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_0_yyyy = cbuffer.data(sg_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_0_yyyz = cbuffer.data(sg_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_0_yyzz = cbuffer.data(sg_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_0_yzzz = cbuffer.data(sg_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_0_zzzz = cbuffer.data(sg_geom_01_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xxxx, g_0_x_0_xxxy, g_0_x_0_xxxz, g_0_x_0_xxyy, g_0_x_0_xxyz, g_0_x_0_xxzz, g_0_x_0_xyyy, g_0_x_0_xyyz, g_0_x_0_xyzz, g_0_x_0_xzzz, g_0_x_0_yyyy, g_0_x_0_yyyz, g_0_x_0_yyzz, g_0_x_0_yzzz, g_0_x_0_zzzz, g_0_xxx, g_0_xxxxx, g_0_xxxxy, g_0_xxxxz, g_0_xxxyy, g_0_xxxyz, g_0_xxxzz, g_0_xxy, g_0_xxyyy, g_0_xxyyz, g_0_xxyzz, g_0_xxz, g_0_xxzzz, g_0_xyy, g_0_xyyyy, g_0_xyyyz, g_0_xyyzz, g_0_xyz, g_0_xyzzz, g_0_xzz, g_0_xzzzz, g_0_yyy, g_0_yyz, g_0_yzz, g_0_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_0_xxxx[k] = -4.0 * g_0_xxx[k] + g_0_xxxxx[k];

                g_0_x_0_xxxy[k] = -3.0 * g_0_xxy[k] + g_0_xxxxy[k];

                g_0_x_0_xxxz[k] = -3.0 * g_0_xxz[k] + g_0_xxxxz[k];

                g_0_x_0_xxyy[k] = -2.0 * g_0_xyy[k] + g_0_xxxyy[k];

                g_0_x_0_xxyz[k] = -2.0 * g_0_xyz[k] + g_0_xxxyz[k];

                g_0_x_0_xxzz[k] = -2.0 * g_0_xzz[k] + g_0_xxxzz[k];

                g_0_x_0_xyyy[k] = -g_0_yyy[k] + g_0_xxyyy[k];

                g_0_x_0_xyyz[k] = -g_0_yyz[k] + g_0_xxyyz[k];

                g_0_x_0_xyzz[k] = -g_0_yzz[k] + g_0_xxyzz[k];

                g_0_x_0_xzzz[k] = -g_0_zzz[k] + g_0_xxzzz[k];

                g_0_x_0_yyyy[k] = g_0_xyyyy[k];

                g_0_x_0_yyyz[k] = g_0_xyyyz[k];

                g_0_x_0_yyzz[k] = g_0_xyyzz[k];

                g_0_x_0_yzzz[k] = g_0_xyzzz[k];

                g_0_x_0_zzzz[k] = g_0_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

            auto g_0_y_0_xxxx = cbuffer.data(sg_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_y_0_xxxy = cbuffer.data(sg_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_y_0_xxxz = cbuffer.data(sg_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_y_0_xxyy = cbuffer.data(sg_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_y_0_xxyz = cbuffer.data(sg_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_y_0_xxzz = cbuffer.data(sg_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_y_0_xyyy = cbuffer.data(sg_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_y_0_xyyz = cbuffer.data(sg_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_y_0_xyzz = cbuffer.data(sg_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_y_0_xzzz = cbuffer.data(sg_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_y_0_yyyy = cbuffer.data(sg_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_y_0_yyyz = cbuffer.data(sg_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_y_0_yyzz = cbuffer.data(sg_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_y_0_yzzz = cbuffer.data(sg_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_y_0_zzzz = cbuffer.data(sg_geom_01_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_xxx, g_0_xxxxy, g_0_xxxyy, g_0_xxxyz, g_0_xxy, g_0_xxyyy, g_0_xxyyz, g_0_xxyzz, g_0_xxz, g_0_xyy, g_0_xyyyy, g_0_xyyyz, g_0_xyyzz, g_0_xyz, g_0_xyzzz, g_0_xzz, g_0_y_0_xxxx, g_0_y_0_xxxy, g_0_y_0_xxxz, g_0_y_0_xxyy, g_0_y_0_xxyz, g_0_y_0_xxzz, g_0_y_0_xyyy, g_0_y_0_xyyz, g_0_y_0_xyzz, g_0_y_0_xzzz, g_0_y_0_yyyy, g_0_y_0_yyyz, g_0_y_0_yyzz, g_0_y_0_yzzz, g_0_y_0_zzzz, g_0_yyy, g_0_yyyyy, g_0_yyyyz, g_0_yyyzz, g_0_yyz, g_0_yyzzz, g_0_yzz, g_0_yzzzz, g_0_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_0_xxxx[k] = g_0_xxxxy[k];

                g_0_y_0_xxxy[k] = -g_0_xxx[k] + g_0_xxxyy[k];

                g_0_y_0_xxxz[k] = g_0_xxxyz[k];

                g_0_y_0_xxyy[k] = -2.0 * g_0_xxy[k] + g_0_xxyyy[k];

                g_0_y_0_xxyz[k] = -g_0_xxz[k] + g_0_xxyyz[k];

                g_0_y_0_xxzz[k] = g_0_xxyzz[k];

                g_0_y_0_xyyy[k] = -3.0 * g_0_xyy[k] + g_0_xyyyy[k];

                g_0_y_0_xyyz[k] = -2.0 * g_0_xyz[k] + g_0_xyyyz[k];

                g_0_y_0_xyzz[k] = -g_0_xzz[k] + g_0_xyyzz[k];

                g_0_y_0_xzzz[k] = g_0_xyzzz[k];

                g_0_y_0_yyyy[k] = -4.0 * g_0_yyy[k] + g_0_yyyyy[k];

                g_0_y_0_yyyz[k] = -3.0 * g_0_yyz[k] + g_0_yyyyz[k];

                g_0_y_0_yyzz[k] = -2.0 * g_0_yzz[k] + g_0_yyyzz[k];

                g_0_y_0_yzzz[k] = -g_0_zzz[k] + g_0_yyzzz[k];

                g_0_y_0_zzzz[k] = g_0_yzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

            auto g_0_z_0_xxxx = cbuffer.data(sg_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_z_0_xxxy = cbuffer.data(sg_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_z_0_xxxz = cbuffer.data(sg_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_z_0_xxyy = cbuffer.data(sg_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_z_0_xxyz = cbuffer.data(sg_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_z_0_xxzz = cbuffer.data(sg_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_z_0_xyyy = cbuffer.data(sg_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_z_0_xyyz = cbuffer.data(sg_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_z_0_xyzz = cbuffer.data(sg_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_z_0_xzzz = cbuffer.data(sg_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_z_0_yyyy = cbuffer.data(sg_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_z_0_yyyz = cbuffer.data(sg_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_z_0_yyzz = cbuffer.data(sg_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_z_0_yzzz = cbuffer.data(sg_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_z_0_zzzz = cbuffer.data(sg_geom_01_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_xxx, g_0_xxxxz, g_0_xxxyz, g_0_xxxzz, g_0_xxy, g_0_xxyyz, g_0_xxyzz, g_0_xxz, g_0_xxzzz, g_0_xyy, g_0_xyyyz, g_0_xyyzz, g_0_xyz, g_0_xyzzz, g_0_xzz, g_0_xzzzz, g_0_yyy, g_0_yyyyz, g_0_yyyzz, g_0_yyz, g_0_yyzzz, g_0_yzz, g_0_yzzzz, g_0_z_0_xxxx, g_0_z_0_xxxy, g_0_z_0_xxxz, g_0_z_0_xxyy, g_0_z_0_xxyz, g_0_z_0_xxzz, g_0_z_0_xyyy, g_0_z_0_xyyz, g_0_z_0_xyzz, g_0_z_0_xzzz, g_0_z_0_yyyy, g_0_z_0_yyyz, g_0_z_0_yyzz, g_0_z_0_yzzz, g_0_z_0_zzzz, g_0_zzz, g_0_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_0_xxxx[k] = g_0_xxxxz[k];

                g_0_z_0_xxxy[k] = g_0_xxxyz[k];

                g_0_z_0_xxxz[k] = -g_0_xxx[k] + g_0_xxxzz[k];

                g_0_z_0_xxyy[k] = g_0_xxyyz[k];

                g_0_z_0_xxyz[k] = -g_0_xxy[k] + g_0_xxyzz[k];

                g_0_z_0_xxzz[k] = -2.0 * g_0_xxz[k] + g_0_xxzzz[k];

                g_0_z_0_xyyy[k] = g_0_xyyyz[k];

                g_0_z_0_xyyz[k] = -g_0_xyy[k] + g_0_xyyzz[k];

                g_0_z_0_xyzz[k] = -2.0 * g_0_xyz[k] + g_0_xyzzz[k];

                g_0_z_0_xzzz[k] = -3.0 * g_0_xzz[k] + g_0_xzzzz[k];

                g_0_z_0_yyyy[k] = g_0_yyyyz[k];

                g_0_z_0_yyyz[k] = -g_0_yyy[k] + g_0_yyyzz[k];

                g_0_z_0_yyzz[k] = -2.0 * g_0_yyz[k] + g_0_yyzzz[k];

                g_0_z_0_yzzz[k] = -3.0 * g_0_yzz[k] + g_0_yzzzz[k];

                g_0_z_0_zzzz[k] = -4.0 * g_0_zzz[k] + g_0_zzzzz[k];
            }
        }
    }
}

} // erirec namespace

