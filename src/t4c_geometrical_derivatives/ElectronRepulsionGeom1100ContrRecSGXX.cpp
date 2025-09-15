#include "ElectronRepulsionGeom1100ContrRecSGXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom11_hrr_electron_repulsion_sgxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_11_sgxx,
                                            const size_t idx_sgxx,
                                            const size_t idx_geom_01_sgxx,
                                            const size_t idx_geom_01_shxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom,});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom,});

    // set up R(AB) distances

    const auto xyz = r_ab.coordinates();

    const auto ab_x = xyz[0];

    const auto ab_y = xyz[1];

    const auto ab_z = xyz[2];

    for (int i = 0; i < ccomps; i++)
    {
        for (int j = 0; j < dcomps; j++)
        {
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

            /// Set up components of auxilary buffer : SGSS

            const auto sg_geom_01_off = idx_geom_01_sgxx + i * dcomps + j;

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

            /// Set up components of auxilary buffer : SHSS

            const auto sh_geom_01_off = idx_geom_01_shxx + i * dcomps + j;

            auto g_0_x_0_xxxxx = cbuffer.data(sh_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_0_xxxxy = cbuffer.data(sh_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_0_xxxxz = cbuffer.data(sh_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_0_xxxyy = cbuffer.data(sh_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_0_xxxyz = cbuffer.data(sh_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_0_xxxzz = cbuffer.data(sh_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_0_xxyyy = cbuffer.data(sh_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_0_xxyyz = cbuffer.data(sh_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_0_xxyzz = cbuffer.data(sh_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_0_xxzzz = cbuffer.data(sh_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_0_xyyyy = cbuffer.data(sh_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_0_xyyyz = cbuffer.data(sh_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_0_xyyzz = cbuffer.data(sh_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_0_xyzzz = cbuffer.data(sh_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_0_xzzzz = cbuffer.data(sh_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_0_yyyyy = cbuffer.data(sh_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_0_yyyyz = cbuffer.data(sh_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_0_yyyzz = cbuffer.data(sh_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_0_yyzzz = cbuffer.data(sh_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_0_yzzzz = cbuffer.data(sh_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_0_zzzzz = cbuffer.data(sh_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_y_0_xxxxx = cbuffer.data(sh_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_y_0_xxxxy = cbuffer.data(sh_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_y_0_xxxxz = cbuffer.data(sh_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_y_0_xxxyy = cbuffer.data(sh_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_y_0_xxxyz = cbuffer.data(sh_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_y_0_xxxzz = cbuffer.data(sh_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_y_0_xxyyy = cbuffer.data(sh_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_y_0_xxyyz = cbuffer.data(sh_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_y_0_xxyzz = cbuffer.data(sh_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_y_0_xxzzz = cbuffer.data(sh_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_y_0_xyyyy = cbuffer.data(sh_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_y_0_xyyyz = cbuffer.data(sh_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_y_0_xyyzz = cbuffer.data(sh_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_y_0_xyzzz = cbuffer.data(sh_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_y_0_xzzzz = cbuffer.data(sh_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_y_0_yyyyy = cbuffer.data(sh_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_y_0_yyyyz = cbuffer.data(sh_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_y_0_yyyzz = cbuffer.data(sh_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_y_0_yyzzz = cbuffer.data(sh_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_y_0_yzzzz = cbuffer.data(sh_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_y_0_zzzzz = cbuffer.data(sh_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_z_0_xxxxx = cbuffer.data(sh_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_z_0_xxxxy = cbuffer.data(sh_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_z_0_xxxxz = cbuffer.data(sh_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_z_0_xxxyy = cbuffer.data(sh_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_z_0_xxxyz = cbuffer.data(sh_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_z_0_xxxzz = cbuffer.data(sh_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_z_0_xxyyy = cbuffer.data(sh_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_z_0_xxyyz = cbuffer.data(sh_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_z_0_xxyzz = cbuffer.data(sh_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_z_0_xxzzz = cbuffer.data(sh_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_z_0_xyyyy = cbuffer.data(sh_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_z_0_xyyyz = cbuffer.data(sh_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_z_0_xyyzz = cbuffer.data(sh_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_z_0_xyzzz = cbuffer.data(sh_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_z_0_xzzzz = cbuffer.data(sh_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_z_0_yyyyy = cbuffer.data(sh_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_z_0_yyyyz = cbuffer.data(sh_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_z_0_yyyzz = cbuffer.data(sh_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_z_0_yyzzz = cbuffer.data(sh_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_z_0_yzzzz = cbuffer.data(sh_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_z_0_zzzzz = cbuffer.data(sh_geom_01_off + 62 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_sgxx

            const auto sg_geom_11_off = idx_geom_11_sgxx + i * dcomps + j;

            /// Set up 0-15 components of targeted buffer : cbuffer.data(

            auto g_x_x_0_xxxx = cbuffer.data(sg_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_0_xxxy = cbuffer.data(sg_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_0_xxxz = cbuffer.data(sg_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_0_xxyy = cbuffer.data(sg_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_0_xxyz = cbuffer.data(sg_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_0_xxzz = cbuffer.data(sg_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_0_xyyy = cbuffer.data(sg_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_0_xyyz = cbuffer.data(sg_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_0_xyzz = cbuffer.data(sg_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_0_xzzz = cbuffer.data(sg_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_0_yyyy = cbuffer.data(sg_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_0_yyyz = cbuffer.data(sg_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_0_yyzz = cbuffer.data(sg_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_0_yzzz = cbuffer.data(sg_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_0_zzzz = cbuffer.data(sg_geom_11_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xxxx, g_0_x_0_xxxxx, g_0_x_0_xxxxy, g_0_x_0_xxxxz, g_0_x_0_xxxy, g_0_x_0_xxxyy, g_0_x_0_xxxyz, g_0_x_0_xxxz, g_0_x_0_xxxzz, g_0_x_0_xxyy, g_0_x_0_xxyyy, g_0_x_0_xxyyz, g_0_x_0_xxyz, g_0_x_0_xxyzz, g_0_x_0_xxzz, g_0_x_0_xxzzz, g_0_x_0_xyyy, g_0_x_0_xyyyy, g_0_x_0_xyyyz, g_0_x_0_xyyz, g_0_x_0_xyyzz, g_0_x_0_xyzz, g_0_x_0_xyzzz, g_0_x_0_xzzz, g_0_x_0_xzzzz, g_0_x_0_yyyy, g_0_x_0_yyyz, g_0_x_0_yyzz, g_0_x_0_yzzz, g_0_x_0_zzzz, g_0_xxxx, g_0_xxxy, g_0_xxxz, g_0_xxyy, g_0_xxyz, g_0_xxzz, g_0_xyyy, g_0_xyyz, g_0_xyzz, g_0_xzzz, g_0_yyyy, g_0_yyyz, g_0_yyzz, g_0_yzzz, g_0_zzzz, g_x_x_0_xxxx, g_x_x_0_xxxy, g_x_x_0_xxxz, g_x_x_0_xxyy, g_x_x_0_xxyz, g_x_x_0_xxzz, g_x_x_0_xyyy, g_x_x_0_xyyz, g_x_x_0_xyzz, g_x_x_0_xzzz, g_x_x_0_yyyy, g_x_x_0_yyyz, g_x_x_0_yyzz, g_x_x_0_yzzz, g_x_x_0_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_0_xxxx[k] = g_0_xxxx[k] - g_0_x_0_xxxx[k] * ab_x + g_0_x_0_xxxxx[k];

                g_x_x_0_xxxy[k] = g_0_xxxy[k] - g_0_x_0_xxxy[k] * ab_x + g_0_x_0_xxxxy[k];

                g_x_x_0_xxxz[k] = g_0_xxxz[k] - g_0_x_0_xxxz[k] * ab_x + g_0_x_0_xxxxz[k];

                g_x_x_0_xxyy[k] = g_0_xxyy[k] - g_0_x_0_xxyy[k] * ab_x + g_0_x_0_xxxyy[k];

                g_x_x_0_xxyz[k] = g_0_xxyz[k] - g_0_x_0_xxyz[k] * ab_x + g_0_x_0_xxxyz[k];

                g_x_x_0_xxzz[k] = g_0_xxzz[k] - g_0_x_0_xxzz[k] * ab_x + g_0_x_0_xxxzz[k];

                g_x_x_0_xyyy[k] = g_0_xyyy[k] - g_0_x_0_xyyy[k] * ab_x + g_0_x_0_xxyyy[k];

                g_x_x_0_xyyz[k] = g_0_xyyz[k] - g_0_x_0_xyyz[k] * ab_x + g_0_x_0_xxyyz[k];

                g_x_x_0_xyzz[k] = g_0_xyzz[k] - g_0_x_0_xyzz[k] * ab_x + g_0_x_0_xxyzz[k];

                g_x_x_0_xzzz[k] = g_0_xzzz[k] - g_0_x_0_xzzz[k] * ab_x + g_0_x_0_xxzzz[k];

                g_x_x_0_yyyy[k] = g_0_yyyy[k] - g_0_x_0_yyyy[k] * ab_x + g_0_x_0_xyyyy[k];

                g_x_x_0_yyyz[k] = g_0_yyyz[k] - g_0_x_0_yyyz[k] * ab_x + g_0_x_0_xyyyz[k];

                g_x_x_0_yyzz[k] = g_0_yyzz[k] - g_0_x_0_yyzz[k] * ab_x + g_0_x_0_xyyzz[k];

                g_x_x_0_yzzz[k] = g_0_yzzz[k] - g_0_x_0_yzzz[k] * ab_x + g_0_x_0_xyzzz[k];

                g_x_x_0_zzzz[k] = g_0_zzzz[k] - g_0_x_0_zzzz[k] * ab_x + g_0_x_0_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

            auto g_x_y_0_xxxx = cbuffer.data(sg_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_y_0_xxxy = cbuffer.data(sg_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_y_0_xxxz = cbuffer.data(sg_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_y_0_xxyy = cbuffer.data(sg_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_y_0_xxyz = cbuffer.data(sg_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_y_0_xxzz = cbuffer.data(sg_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_y_0_xyyy = cbuffer.data(sg_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_y_0_xyyz = cbuffer.data(sg_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_y_0_xyzz = cbuffer.data(sg_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_y_0_xzzz = cbuffer.data(sg_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_y_0_yyyy = cbuffer.data(sg_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_y_0_yyyz = cbuffer.data(sg_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_y_0_yyzz = cbuffer.data(sg_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_y_0_yzzz = cbuffer.data(sg_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_y_0_zzzz = cbuffer.data(sg_geom_11_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_0_xxxx, g_0_y_0_xxxxx, g_0_y_0_xxxxy, g_0_y_0_xxxxz, g_0_y_0_xxxy, g_0_y_0_xxxyy, g_0_y_0_xxxyz, g_0_y_0_xxxz, g_0_y_0_xxxzz, g_0_y_0_xxyy, g_0_y_0_xxyyy, g_0_y_0_xxyyz, g_0_y_0_xxyz, g_0_y_0_xxyzz, g_0_y_0_xxzz, g_0_y_0_xxzzz, g_0_y_0_xyyy, g_0_y_0_xyyyy, g_0_y_0_xyyyz, g_0_y_0_xyyz, g_0_y_0_xyyzz, g_0_y_0_xyzz, g_0_y_0_xyzzz, g_0_y_0_xzzz, g_0_y_0_xzzzz, g_0_y_0_yyyy, g_0_y_0_yyyz, g_0_y_0_yyzz, g_0_y_0_yzzz, g_0_y_0_zzzz, g_x_y_0_xxxx, g_x_y_0_xxxy, g_x_y_0_xxxz, g_x_y_0_xxyy, g_x_y_0_xxyz, g_x_y_0_xxzz, g_x_y_0_xyyy, g_x_y_0_xyyz, g_x_y_0_xyzz, g_x_y_0_xzzz, g_x_y_0_yyyy, g_x_y_0_yyyz, g_x_y_0_yyzz, g_x_y_0_yzzz, g_x_y_0_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_0_xxxx[k] = -g_0_y_0_xxxx[k] * ab_x + g_0_y_0_xxxxx[k];

                g_x_y_0_xxxy[k] = -g_0_y_0_xxxy[k] * ab_x + g_0_y_0_xxxxy[k];

                g_x_y_0_xxxz[k] = -g_0_y_0_xxxz[k] * ab_x + g_0_y_0_xxxxz[k];

                g_x_y_0_xxyy[k] = -g_0_y_0_xxyy[k] * ab_x + g_0_y_0_xxxyy[k];

                g_x_y_0_xxyz[k] = -g_0_y_0_xxyz[k] * ab_x + g_0_y_0_xxxyz[k];

                g_x_y_0_xxzz[k] = -g_0_y_0_xxzz[k] * ab_x + g_0_y_0_xxxzz[k];

                g_x_y_0_xyyy[k] = -g_0_y_0_xyyy[k] * ab_x + g_0_y_0_xxyyy[k];

                g_x_y_0_xyyz[k] = -g_0_y_0_xyyz[k] * ab_x + g_0_y_0_xxyyz[k];

                g_x_y_0_xyzz[k] = -g_0_y_0_xyzz[k] * ab_x + g_0_y_0_xxyzz[k];

                g_x_y_0_xzzz[k] = -g_0_y_0_xzzz[k] * ab_x + g_0_y_0_xxzzz[k];

                g_x_y_0_yyyy[k] = -g_0_y_0_yyyy[k] * ab_x + g_0_y_0_xyyyy[k];

                g_x_y_0_yyyz[k] = -g_0_y_0_yyyz[k] * ab_x + g_0_y_0_xyyyz[k];

                g_x_y_0_yyzz[k] = -g_0_y_0_yyzz[k] * ab_x + g_0_y_0_xyyzz[k];

                g_x_y_0_yzzz[k] = -g_0_y_0_yzzz[k] * ab_x + g_0_y_0_xyzzz[k];

                g_x_y_0_zzzz[k] = -g_0_y_0_zzzz[k] * ab_x + g_0_y_0_xzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

            auto g_x_z_0_xxxx = cbuffer.data(sg_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_z_0_xxxy = cbuffer.data(sg_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_z_0_xxxz = cbuffer.data(sg_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_z_0_xxyy = cbuffer.data(sg_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_z_0_xxyz = cbuffer.data(sg_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_z_0_xxzz = cbuffer.data(sg_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_z_0_xyyy = cbuffer.data(sg_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_z_0_xyyz = cbuffer.data(sg_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_z_0_xyzz = cbuffer.data(sg_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_z_0_xzzz = cbuffer.data(sg_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_z_0_yyyy = cbuffer.data(sg_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_z_0_yyyz = cbuffer.data(sg_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_z_0_yyzz = cbuffer.data(sg_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_z_0_yzzz = cbuffer.data(sg_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_z_0_zzzz = cbuffer.data(sg_geom_11_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_0_xxxx, g_0_z_0_xxxxx, g_0_z_0_xxxxy, g_0_z_0_xxxxz, g_0_z_0_xxxy, g_0_z_0_xxxyy, g_0_z_0_xxxyz, g_0_z_0_xxxz, g_0_z_0_xxxzz, g_0_z_0_xxyy, g_0_z_0_xxyyy, g_0_z_0_xxyyz, g_0_z_0_xxyz, g_0_z_0_xxyzz, g_0_z_0_xxzz, g_0_z_0_xxzzz, g_0_z_0_xyyy, g_0_z_0_xyyyy, g_0_z_0_xyyyz, g_0_z_0_xyyz, g_0_z_0_xyyzz, g_0_z_0_xyzz, g_0_z_0_xyzzz, g_0_z_0_xzzz, g_0_z_0_xzzzz, g_0_z_0_yyyy, g_0_z_0_yyyz, g_0_z_0_yyzz, g_0_z_0_yzzz, g_0_z_0_zzzz, g_x_z_0_xxxx, g_x_z_0_xxxy, g_x_z_0_xxxz, g_x_z_0_xxyy, g_x_z_0_xxyz, g_x_z_0_xxzz, g_x_z_0_xyyy, g_x_z_0_xyyz, g_x_z_0_xyzz, g_x_z_0_xzzz, g_x_z_0_yyyy, g_x_z_0_yyyz, g_x_z_0_yyzz, g_x_z_0_yzzz, g_x_z_0_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_0_xxxx[k] = -g_0_z_0_xxxx[k] * ab_x + g_0_z_0_xxxxx[k];

                g_x_z_0_xxxy[k] = -g_0_z_0_xxxy[k] * ab_x + g_0_z_0_xxxxy[k];

                g_x_z_0_xxxz[k] = -g_0_z_0_xxxz[k] * ab_x + g_0_z_0_xxxxz[k];

                g_x_z_0_xxyy[k] = -g_0_z_0_xxyy[k] * ab_x + g_0_z_0_xxxyy[k];

                g_x_z_0_xxyz[k] = -g_0_z_0_xxyz[k] * ab_x + g_0_z_0_xxxyz[k];

                g_x_z_0_xxzz[k] = -g_0_z_0_xxzz[k] * ab_x + g_0_z_0_xxxzz[k];

                g_x_z_0_xyyy[k] = -g_0_z_0_xyyy[k] * ab_x + g_0_z_0_xxyyy[k];

                g_x_z_0_xyyz[k] = -g_0_z_0_xyyz[k] * ab_x + g_0_z_0_xxyyz[k];

                g_x_z_0_xyzz[k] = -g_0_z_0_xyzz[k] * ab_x + g_0_z_0_xxyzz[k];

                g_x_z_0_xzzz[k] = -g_0_z_0_xzzz[k] * ab_x + g_0_z_0_xxzzz[k];

                g_x_z_0_yyyy[k] = -g_0_z_0_yyyy[k] * ab_x + g_0_z_0_xyyyy[k];

                g_x_z_0_yyyz[k] = -g_0_z_0_yyyz[k] * ab_x + g_0_z_0_xyyyz[k];

                g_x_z_0_yyzz[k] = -g_0_z_0_yyzz[k] * ab_x + g_0_z_0_xyyzz[k];

                g_x_z_0_yzzz[k] = -g_0_z_0_yzzz[k] * ab_x + g_0_z_0_xyzzz[k];

                g_x_z_0_zzzz[k] = -g_0_z_0_zzzz[k] * ab_x + g_0_z_0_xzzzz[k];
            }

            /// Set up 45-60 components of targeted buffer : cbuffer.data(

            auto g_y_x_0_xxxx = cbuffer.data(sg_geom_11_off + 45 * ccomps * dcomps);

            auto g_y_x_0_xxxy = cbuffer.data(sg_geom_11_off + 46 * ccomps * dcomps);

            auto g_y_x_0_xxxz = cbuffer.data(sg_geom_11_off + 47 * ccomps * dcomps);

            auto g_y_x_0_xxyy = cbuffer.data(sg_geom_11_off + 48 * ccomps * dcomps);

            auto g_y_x_0_xxyz = cbuffer.data(sg_geom_11_off + 49 * ccomps * dcomps);

            auto g_y_x_0_xxzz = cbuffer.data(sg_geom_11_off + 50 * ccomps * dcomps);

            auto g_y_x_0_xyyy = cbuffer.data(sg_geom_11_off + 51 * ccomps * dcomps);

            auto g_y_x_0_xyyz = cbuffer.data(sg_geom_11_off + 52 * ccomps * dcomps);

            auto g_y_x_0_xyzz = cbuffer.data(sg_geom_11_off + 53 * ccomps * dcomps);

            auto g_y_x_0_xzzz = cbuffer.data(sg_geom_11_off + 54 * ccomps * dcomps);

            auto g_y_x_0_yyyy = cbuffer.data(sg_geom_11_off + 55 * ccomps * dcomps);

            auto g_y_x_0_yyyz = cbuffer.data(sg_geom_11_off + 56 * ccomps * dcomps);

            auto g_y_x_0_yyzz = cbuffer.data(sg_geom_11_off + 57 * ccomps * dcomps);

            auto g_y_x_0_yzzz = cbuffer.data(sg_geom_11_off + 58 * ccomps * dcomps);

            auto g_y_x_0_zzzz = cbuffer.data(sg_geom_11_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xxxx, g_0_x_0_xxxxy, g_0_x_0_xxxy, g_0_x_0_xxxyy, g_0_x_0_xxxyz, g_0_x_0_xxxz, g_0_x_0_xxyy, g_0_x_0_xxyyy, g_0_x_0_xxyyz, g_0_x_0_xxyz, g_0_x_0_xxyzz, g_0_x_0_xxzz, g_0_x_0_xyyy, g_0_x_0_xyyyy, g_0_x_0_xyyyz, g_0_x_0_xyyz, g_0_x_0_xyyzz, g_0_x_0_xyzz, g_0_x_0_xyzzz, g_0_x_0_xzzz, g_0_x_0_yyyy, g_0_x_0_yyyyy, g_0_x_0_yyyyz, g_0_x_0_yyyz, g_0_x_0_yyyzz, g_0_x_0_yyzz, g_0_x_0_yyzzz, g_0_x_0_yzzz, g_0_x_0_yzzzz, g_0_x_0_zzzz, g_y_x_0_xxxx, g_y_x_0_xxxy, g_y_x_0_xxxz, g_y_x_0_xxyy, g_y_x_0_xxyz, g_y_x_0_xxzz, g_y_x_0_xyyy, g_y_x_0_xyyz, g_y_x_0_xyzz, g_y_x_0_xzzz, g_y_x_0_yyyy, g_y_x_0_yyyz, g_y_x_0_yyzz, g_y_x_0_yzzz, g_y_x_0_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_0_xxxx[k] = -g_0_x_0_xxxx[k] * ab_y + g_0_x_0_xxxxy[k];

                g_y_x_0_xxxy[k] = -g_0_x_0_xxxy[k] * ab_y + g_0_x_0_xxxyy[k];

                g_y_x_0_xxxz[k] = -g_0_x_0_xxxz[k] * ab_y + g_0_x_0_xxxyz[k];

                g_y_x_0_xxyy[k] = -g_0_x_0_xxyy[k] * ab_y + g_0_x_0_xxyyy[k];

                g_y_x_0_xxyz[k] = -g_0_x_0_xxyz[k] * ab_y + g_0_x_0_xxyyz[k];

                g_y_x_0_xxzz[k] = -g_0_x_0_xxzz[k] * ab_y + g_0_x_0_xxyzz[k];

                g_y_x_0_xyyy[k] = -g_0_x_0_xyyy[k] * ab_y + g_0_x_0_xyyyy[k];

                g_y_x_0_xyyz[k] = -g_0_x_0_xyyz[k] * ab_y + g_0_x_0_xyyyz[k];

                g_y_x_0_xyzz[k] = -g_0_x_0_xyzz[k] * ab_y + g_0_x_0_xyyzz[k];

                g_y_x_0_xzzz[k] = -g_0_x_0_xzzz[k] * ab_y + g_0_x_0_xyzzz[k];

                g_y_x_0_yyyy[k] = -g_0_x_0_yyyy[k] * ab_y + g_0_x_0_yyyyy[k];

                g_y_x_0_yyyz[k] = -g_0_x_0_yyyz[k] * ab_y + g_0_x_0_yyyyz[k];

                g_y_x_0_yyzz[k] = -g_0_x_0_yyzz[k] * ab_y + g_0_x_0_yyyzz[k];

                g_y_x_0_yzzz[k] = -g_0_x_0_yzzz[k] * ab_y + g_0_x_0_yyzzz[k];

                g_y_x_0_zzzz[k] = -g_0_x_0_zzzz[k] * ab_y + g_0_x_0_yzzzz[k];
            }

            /// Set up 60-75 components of targeted buffer : cbuffer.data(

            auto g_y_y_0_xxxx = cbuffer.data(sg_geom_11_off + 60 * ccomps * dcomps);

            auto g_y_y_0_xxxy = cbuffer.data(sg_geom_11_off + 61 * ccomps * dcomps);

            auto g_y_y_0_xxxz = cbuffer.data(sg_geom_11_off + 62 * ccomps * dcomps);

            auto g_y_y_0_xxyy = cbuffer.data(sg_geom_11_off + 63 * ccomps * dcomps);

            auto g_y_y_0_xxyz = cbuffer.data(sg_geom_11_off + 64 * ccomps * dcomps);

            auto g_y_y_0_xxzz = cbuffer.data(sg_geom_11_off + 65 * ccomps * dcomps);

            auto g_y_y_0_xyyy = cbuffer.data(sg_geom_11_off + 66 * ccomps * dcomps);

            auto g_y_y_0_xyyz = cbuffer.data(sg_geom_11_off + 67 * ccomps * dcomps);

            auto g_y_y_0_xyzz = cbuffer.data(sg_geom_11_off + 68 * ccomps * dcomps);

            auto g_y_y_0_xzzz = cbuffer.data(sg_geom_11_off + 69 * ccomps * dcomps);

            auto g_y_y_0_yyyy = cbuffer.data(sg_geom_11_off + 70 * ccomps * dcomps);

            auto g_y_y_0_yyyz = cbuffer.data(sg_geom_11_off + 71 * ccomps * dcomps);

            auto g_y_y_0_yyzz = cbuffer.data(sg_geom_11_off + 72 * ccomps * dcomps);

            auto g_y_y_0_yzzz = cbuffer.data(sg_geom_11_off + 73 * ccomps * dcomps);

            auto g_y_y_0_zzzz = cbuffer.data(sg_geom_11_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_xxxx, g_0_xxxy, g_0_xxxz, g_0_xxyy, g_0_xxyz, g_0_xxzz, g_0_xyyy, g_0_xyyz, g_0_xyzz, g_0_xzzz, g_0_y_0_xxxx, g_0_y_0_xxxxy, g_0_y_0_xxxy, g_0_y_0_xxxyy, g_0_y_0_xxxyz, g_0_y_0_xxxz, g_0_y_0_xxyy, g_0_y_0_xxyyy, g_0_y_0_xxyyz, g_0_y_0_xxyz, g_0_y_0_xxyzz, g_0_y_0_xxzz, g_0_y_0_xyyy, g_0_y_0_xyyyy, g_0_y_0_xyyyz, g_0_y_0_xyyz, g_0_y_0_xyyzz, g_0_y_0_xyzz, g_0_y_0_xyzzz, g_0_y_0_xzzz, g_0_y_0_yyyy, g_0_y_0_yyyyy, g_0_y_0_yyyyz, g_0_y_0_yyyz, g_0_y_0_yyyzz, g_0_y_0_yyzz, g_0_y_0_yyzzz, g_0_y_0_yzzz, g_0_y_0_yzzzz, g_0_y_0_zzzz, g_0_yyyy, g_0_yyyz, g_0_yyzz, g_0_yzzz, g_0_zzzz, g_y_y_0_xxxx, g_y_y_0_xxxy, g_y_y_0_xxxz, g_y_y_0_xxyy, g_y_y_0_xxyz, g_y_y_0_xxzz, g_y_y_0_xyyy, g_y_y_0_xyyz, g_y_y_0_xyzz, g_y_y_0_xzzz, g_y_y_0_yyyy, g_y_y_0_yyyz, g_y_y_0_yyzz, g_y_y_0_yzzz, g_y_y_0_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_0_xxxx[k] = g_0_xxxx[k] - g_0_y_0_xxxx[k] * ab_y + g_0_y_0_xxxxy[k];

                g_y_y_0_xxxy[k] = g_0_xxxy[k] - g_0_y_0_xxxy[k] * ab_y + g_0_y_0_xxxyy[k];

                g_y_y_0_xxxz[k] = g_0_xxxz[k] - g_0_y_0_xxxz[k] * ab_y + g_0_y_0_xxxyz[k];

                g_y_y_0_xxyy[k] = g_0_xxyy[k] - g_0_y_0_xxyy[k] * ab_y + g_0_y_0_xxyyy[k];

                g_y_y_0_xxyz[k] = g_0_xxyz[k] - g_0_y_0_xxyz[k] * ab_y + g_0_y_0_xxyyz[k];

                g_y_y_0_xxzz[k] = g_0_xxzz[k] - g_0_y_0_xxzz[k] * ab_y + g_0_y_0_xxyzz[k];

                g_y_y_0_xyyy[k] = g_0_xyyy[k] - g_0_y_0_xyyy[k] * ab_y + g_0_y_0_xyyyy[k];

                g_y_y_0_xyyz[k] = g_0_xyyz[k] - g_0_y_0_xyyz[k] * ab_y + g_0_y_0_xyyyz[k];

                g_y_y_0_xyzz[k] = g_0_xyzz[k] - g_0_y_0_xyzz[k] * ab_y + g_0_y_0_xyyzz[k];

                g_y_y_0_xzzz[k] = g_0_xzzz[k] - g_0_y_0_xzzz[k] * ab_y + g_0_y_0_xyzzz[k];

                g_y_y_0_yyyy[k] = g_0_yyyy[k] - g_0_y_0_yyyy[k] * ab_y + g_0_y_0_yyyyy[k];

                g_y_y_0_yyyz[k] = g_0_yyyz[k] - g_0_y_0_yyyz[k] * ab_y + g_0_y_0_yyyyz[k];

                g_y_y_0_yyzz[k] = g_0_yyzz[k] - g_0_y_0_yyzz[k] * ab_y + g_0_y_0_yyyzz[k];

                g_y_y_0_yzzz[k] = g_0_yzzz[k] - g_0_y_0_yzzz[k] * ab_y + g_0_y_0_yyzzz[k];

                g_y_y_0_zzzz[k] = g_0_zzzz[k] - g_0_y_0_zzzz[k] * ab_y + g_0_y_0_yzzzz[k];
            }

            /// Set up 75-90 components of targeted buffer : cbuffer.data(

            auto g_y_z_0_xxxx = cbuffer.data(sg_geom_11_off + 75 * ccomps * dcomps);

            auto g_y_z_0_xxxy = cbuffer.data(sg_geom_11_off + 76 * ccomps * dcomps);

            auto g_y_z_0_xxxz = cbuffer.data(sg_geom_11_off + 77 * ccomps * dcomps);

            auto g_y_z_0_xxyy = cbuffer.data(sg_geom_11_off + 78 * ccomps * dcomps);

            auto g_y_z_0_xxyz = cbuffer.data(sg_geom_11_off + 79 * ccomps * dcomps);

            auto g_y_z_0_xxzz = cbuffer.data(sg_geom_11_off + 80 * ccomps * dcomps);

            auto g_y_z_0_xyyy = cbuffer.data(sg_geom_11_off + 81 * ccomps * dcomps);

            auto g_y_z_0_xyyz = cbuffer.data(sg_geom_11_off + 82 * ccomps * dcomps);

            auto g_y_z_0_xyzz = cbuffer.data(sg_geom_11_off + 83 * ccomps * dcomps);

            auto g_y_z_0_xzzz = cbuffer.data(sg_geom_11_off + 84 * ccomps * dcomps);

            auto g_y_z_0_yyyy = cbuffer.data(sg_geom_11_off + 85 * ccomps * dcomps);

            auto g_y_z_0_yyyz = cbuffer.data(sg_geom_11_off + 86 * ccomps * dcomps);

            auto g_y_z_0_yyzz = cbuffer.data(sg_geom_11_off + 87 * ccomps * dcomps);

            auto g_y_z_0_yzzz = cbuffer.data(sg_geom_11_off + 88 * ccomps * dcomps);

            auto g_y_z_0_zzzz = cbuffer.data(sg_geom_11_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_0_xxxx, g_0_z_0_xxxxy, g_0_z_0_xxxy, g_0_z_0_xxxyy, g_0_z_0_xxxyz, g_0_z_0_xxxz, g_0_z_0_xxyy, g_0_z_0_xxyyy, g_0_z_0_xxyyz, g_0_z_0_xxyz, g_0_z_0_xxyzz, g_0_z_0_xxzz, g_0_z_0_xyyy, g_0_z_0_xyyyy, g_0_z_0_xyyyz, g_0_z_0_xyyz, g_0_z_0_xyyzz, g_0_z_0_xyzz, g_0_z_0_xyzzz, g_0_z_0_xzzz, g_0_z_0_yyyy, g_0_z_0_yyyyy, g_0_z_0_yyyyz, g_0_z_0_yyyz, g_0_z_0_yyyzz, g_0_z_0_yyzz, g_0_z_0_yyzzz, g_0_z_0_yzzz, g_0_z_0_yzzzz, g_0_z_0_zzzz, g_y_z_0_xxxx, g_y_z_0_xxxy, g_y_z_0_xxxz, g_y_z_0_xxyy, g_y_z_0_xxyz, g_y_z_0_xxzz, g_y_z_0_xyyy, g_y_z_0_xyyz, g_y_z_0_xyzz, g_y_z_0_xzzz, g_y_z_0_yyyy, g_y_z_0_yyyz, g_y_z_0_yyzz, g_y_z_0_yzzz, g_y_z_0_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_0_xxxx[k] = -g_0_z_0_xxxx[k] * ab_y + g_0_z_0_xxxxy[k];

                g_y_z_0_xxxy[k] = -g_0_z_0_xxxy[k] * ab_y + g_0_z_0_xxxyy[k];

                g_y_z_0_xxxz[k] = -g_0_z_0_xxxz[k] * ab_y + g_0_z_0_xxxyz[k];

                g_y_z_0_xxyy[k] = -g_0_z_0_xxyy[k] * ab_y + g_0_z_0_xxyyy[k];

                g_y_z_0_xxyz[k] = -g_0_z_0_xxyz[k] * ab_y + g_0_z_0_xxyyz[k];

                g_y_z_0_xxzz[k] = -g_0_z_0_xxzz[k] * ab_y + g_0_z_0_xxyzz[k];

                g_y_z_0_xyyy[k] = -g_0_z_0_xyyy[k] * ab_y + g_0_z_0_xyyyy[k];

                g_y_z_0_xyyz[k] = -g_0_z_0_xyyz[k] * ab_y + g_0_z_0_xyyyz[k];

                g_y_z_0_xyzz[k] = -g_0_z_0_xyzz[k] * ab_y + g_0_z_0_xyyzz[k];

                g_y_z_0_xzzz[k] = -g_0_z_0_xzzz[k] * ab_y + g_0_z_0_xyzzz[k];

                g_y_z_0_yyyy[k] = -g_0_z_0_yyyy[k] * ab_y + g_0_z_0_yyyyy[k];

                g_y_z_0_yyyz[k] = -g_0_z_0_yyyz[k] * ab_y + g_0_z_0_yyyyz[k];

                g_y_z_0_yyzz[k] = -g_0_z_0_yyzz[k] * ab_y + g_0_z_0_yyyzz[k];

                g_y_z_0_yzzz[k] = -g_0_z_0_yzzz[k] * ab_y + g_0_z_0_yyzzz[k];

                g_y_z_0_zzzz[k] = -g_0_z_0_zzzz[k] * ab_y + g_0_z_0_yzzzz[k];
            }

            /// Set up 90-105 components of targeted buffer : cbuffer.data(

            auto g_z_x_0_xxxx = cbuffer.data(sg_geom_11_off + 90 * ccomps * dcomps);

            auto g_z_x_0_xxxy = cbuffer.data(sg_geom_11_off + 91 * ccomps * dcomps);

            auto g_z_x_0_xxxz = cbuffer.data(sg_geom_11_off + 92 * ccomps * dcomps);

            auto g_z_x_0_xxyy = cbuffer.data(sg_geom_11_off + 93 * ccomps * dcomps);

            auto g_z_x_0_xxyz = cbuffer.data(sg_geom_11_off + 94 * ccomps * dcomps);

            auto g_z_x_0_xxzz = cbuffer.data(sg_geom_11_off + 95 * ccomps * dcomps);

            auto g_z_x_0_xyyy = cbuffer.data(sg_geom_11_off + 96 * ccomps * dcomps);

            auto g_z_x_0_xyyz = cbuffer.data(sg_geom_11_off + 97 * ccomps * dcomps);

            auto g_z_x_0_xyzz = cbuffer.data(sg_geom_11_off + 98 * ccomps * dcomps);

            auto g_z_x_0_xzzz = cbuffer.data(sg_geom_11_off + 99 * ccomps * dcomps);

            auto g_z_x_0_yyyy = cbuffer.data(sg_geom_11_off + 100 * ccomps * dcomps);

            auto g_z_x_0_yyyz = cbuffer.data(sg_geom_11_off + 101 * ccomps * dcomps);

            auto g_z_x_0_yyzz = cbuffer.data(sg_geom_11_off + 102 * ccomps * dcomps);

            auto g_z_x_0_yzzz = cbuffer.data(sg_geom_11_off + 103 * ccomps * dcomps);

            auto g_z_x_0_zzzz = cbuffer.data(sg_geom_11_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xxxx, g_0_x_0_xxxxz, g_0_x_0_xxxy, g_0_x_0_xxxyz, g_0_x_0_xxxz, g_0_x_0_xxxzz, g_0_x_0_xxyy, g_0_x_0_xxyyz, g_0_x_0_xxyz, g_0_x_0_xxyzz, g_0_x_0_xxzz, g_0_x_0_xxzzz, g_0_x_0_xyyy, g_0_x_0_xyyyz, g_0_x_0_xyyz, g_0_x_0_xyyzz, g_0_x_0_xyzz, g_0_x_0_xyzzz, g_0_x_0_xzzz, g_0_x_0_xzzzz, g_0_x_0_yyyy, g_0_x_0_yyyyz, g_0_x_0_yyyz, g_0_x_0_yyyzz, g_0_x_0_yyzz, g_0_x_0_yyzzz, g_0_x_0_yzzz, g_0_x_0_yzzzz, g_0_x_0_zzzz, g_0_x_0_zzzzz, g_z_x_0_xxxx, g_z_x_0_xxxy, g_z_x_0_xxxz, g_z_x_0_xxyy, g_z_x_0_xxyz, g_z_x_0_xxzz, g_z_x_0_xyyy, g_z_x_0_xyyz, g_z_x_0_xyzz, g_z_x_0_xzzz, g_z_x_0_yyyy, g_z_x_0_yyyz, g_z_x_0_yyzz, g_z_x_0_yzzz, g_z_x_0_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_0_xxxx[k] = -g_0_x_0_xxxx[k] * ab_z + g_0_x_0_xxxxz[k];

                g_z_x_0_xxxy[k] = -g_0_x_0_xxxy[k] * ab_z + g_0_x_0_xxxyz[k];

                g_z_x_0_xxxz[k] = -g_0_x_0_xxxz[k] * ab_z + g_0_x_0_xxxzz[k];

                g_z_x_0_xxyy[k] = -g_0_x_0_xxyy[k] * ab_z + g_0_x_0_xxyyz[k];

                g_z_x_0_xxyz[k] = -g_0_x_0_xxyz[k] * ab_z + g_0_x_0_xxyzz[k];

                g_z_x_0_xxzz[k] = -g_0_x_0_xxzz[k] * ab_z + g_0_x_0_xxzzz[k];

                g_z_x_0_xyyy[k] = -g_0_x_0_xyyy[k] * ab_z + g_0_x_0_xyyyz[k];

                g_z_x_0_xyyz[k] = -g_0_x_0_xyyz[k] * ab_z + g_0_x_0_xyyzz[k];

                g_z_x_0_xyzz[k] = -g_0_x_0_xyzz[k] * ab_z + g_0_x_0_xyzzz[k];

                g_z_x_0_xzzz[k] = -g_0_x_0_xzzz[k] * ab_z + g_0_x_0_xzzzz[k];

                g_z_x_0_yyyy[k] = -g_0_x_0_yyyy[k] * ab_z + g_0_x_0_yyyyz[k];

                g_z_x_0_yyyz[k] = -g_0_x_0_yyyz[k] * ab_z + g_0_x_0_yyyzz[k];

                g_z_x_0_yyzz[k] = -g_0_x_0_yyzz[k] * ab_z + g_0_x_0_yyzzz[k];

                g_z_x_0_yzzz[k] = -g_0_x_0_yzzz[k] * ab_z + g_0_x_0_yzzzz[k];

                g_z_x_0_zzzz[k] = -g_0_x_0_zzzz[k] * ab_z + g_0_x_0_zzzzz[k];
            }

            /// Set up 105-120 components of targeted buffer : cbuffer.data(

            auto g_z_y_0_xxxx = cbuffer.data(sg_geom_11_off + 105 * ccomps * dcomps);

            auto g_z_y_0_xxxy = cbuffer.data(sg_geom_11_off + 106 * ccomps * dcomps);

            auto g_z_y_0_xxxz = cbuffer.data(sg_geom_11_off + 107 * ccomps * dcomps);

            auto g_z_y_0_xxyy = cbuffer.data(sg_geom_11_off + 108 * ccomps * dcomps);

            auto g_z_y_0_xxyz = cbuffer.data(sg_geom_11_off + 109 * ccomps * dcomps);

            auto g_z_y_0_xxzz = cbuffer.data(sg_geom_11_off + 110 * ccomps * dcomps);

            auto g_z_y_0_xyyy = cbuffer.data(sg_geom_11_off + 111 * ccomps * dcomps);

            auto g_z_y_0_xyyz = cbuffer.data(sg_geom_11_off + 112 * ccomps * dcomps);

            auto g_z_y_0_xyzz = cbuffer.data(sg_geom_11_off + 113 * ccomps * dcomps);

            auto g_z_y_0_xzzz = cbuffer.data(sg_geom_11_off + 114 * ccomps * dcomps);

            auto g_z_y_0_yyyy = cbuffer.data(sg_geom_11_off + 115 * ccomps * dcomps);

            auto g_z_y_0_yyyz = cbuffer.data(sg_geom_11_off + 116 * ccomps * dcomps);

            auto g_z_y_0_yyzz = cbuffer.data(sg_geom_11_off + 117 * ccomps * dcomps);

            auto g_z_y_0_yzzz = cbuffer.data(sg_geom_11_off + 118 * ccomps * dcomps);

            auto g_z_y_0_zzzz = cbuffer.data(sg_geom_11_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_0_xxxx, g_0_y_0_xxxxz, g_0_y_0_xxxy, g_0_y_0_xxxyz, g_0_y_0_xxxz, g_0_y_0_xxxzz, g_0_y_0_xxyy, g_0_y_0_xxyyz, g_0_y_0_xxyz, g_0_y_0_xxyzz, g_0_y_0_xxzz, g_0_y_0_xxzzz, g_0_y_0_xyyy, g_0_y_0_xyyyz, g_0_y_0_xyyz, g_0_y_0_xyyzz, g_0_y_0_xyzz, g_0_y_0_xyzzz, g_0_y_0_xzzz, g_0_y_0_xzzzz, g_0_y_0_yyyy, g_0_y_0_yyyyz, g_0_y_0_yyyz, g_0_y_0_yyyzz, g_0_y_0_yyzz, g_0_y_0_yyzzz, g_0_y_0_yzzz, g_0_y_0_yzzzz, g_0_y_0_zzzz, g_0_y_0_zzzzz, g_z_y_0_xxxx, g_z_y_0_xxxy, g_z_y_0_xxxz, g_z_y_0_xxyy, g_z_y_0_xxyz, g_z_y_0_xxzz, g_z_y_0_xyyy, g_z_y_0_xyyz, g_z_y_0_xyzz, g_z_y_0_xzzz, g_z_y_0_yyyy, g_z_y_0_yyyz, g_z_y_0_yyzz, g_z_y_0_yzzz, g_z_y_0_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_0_xxxx[k] = -g_0_y_0_xxxx[k] * ab_z + g_0_y_0_xxxxz[k];

                g_z_y_0_xxxy[k] = -g_0_y_0_xxxy[k] * ab_z + g_0_y_0_xxxyz[k];

                g_z_y_0_xxxz[k] = -g_0_y_0_xxxz[k] * ab_z + g_0_y_0_xxxzz[k];

                g_z_y_0_xxyy[k] = -g_0_y_0_xxyy[k] * ab_z + g_0_y_0_xxyyz[k];

                g_z_y_0_xxyz[k] = -g_0_y_0_xxyz[k] * ab_z + g_0_y_0_xxyzz[k];

                g_z_y_0_xxzz[k] = -g_0_y_0_xxzz[k] * ab_z + g_0_y_0_xxzzz[k];

                g_z_y_0_xyyy[k] = -g_0_y_0_xyyy[k] * ab_z + g_0_y_0_xyyyz[k];

                g_z_y_0_xyyz[k] = -g_0_y_0_xyyz[k] * ab_z + g_0_y_0_xyyzz[k];

                g_z_y_0_xyzz[k] = -g_0_y_0_xyzz[k] * ab_z + g_0_y_0_xyzzz[k];

                g_z_y_0_xzzz[k] = -g_0_y_0_xzzz[k] * ab_z + g_0_y_0_xzzzz[k];

                g_z_y_0_yyyy[k] = -g_0_y_0_yyyy[k] * ab_z + g_0_y_0_yyyyz[k];

                g_z_y_0_yyyz[k] = -g_0_y_0_yyyz[k] * ab_z + g_0_y_0_yyyzz[k];

                g_z_y_0_yyzz[k] = -g_0_y_0_yyzz[k] * ab_z + g_0_y_0_yyzzz[k];

                g_z_y_0_yzzz[k] = -g_0_y_0_yzzz[k] * ab_z + g_0_y_0_yzzzz[k];

                g_z_y_0_zzzz[k] = -g_0_y_0_zzzz[k] * ab_z + g_0_y_0_zzzzz[k];
            }

            /// Set up 120-135 components of targeted buffer : cbuffer.data(

            auto g_z_z_0_xxxx = cbuffer.data(sg_geom_11_off + 120 * ccomps * dcomps);

            auto g_z_z_0_xxxy = cbuffer.data(sg_geom_11_off + 121 * ccomps * dcomps);

            auto g_z_z_0_xxxz = cbuffer.data(sg_geom_11_off + 122 * ccomps * dcomps);

            auto g_z_z_0_xxyy = cbuffer.data(sg_geom_11_off + 123 * ccomps * dcomps);

            auto g_z_z_0_xxyz = cbuffer.data(sg_geom_11_off + 124 * ccomps * dcomps);

            auto g_z_z_0_xxzz = cbuffer.data(sg_geom_11_off + 125 * ccomps * dcomps);

            auto g_z_z_0_xyyy = cbuffer.data(sg_geom_11_off + 126 * ccomps * dcomps);

            auto g_z_z_0_xyyz = cbuffer.data(sg_geom_11_off + 127 * ccomps * dcomps);

            auto g_z_z_0_xyzz = cbuffer.data(sg_geom_11_off + 128 * ccomps * dcomps);

            auto g_z_z_0_xzzz = cbuffer.data(sg_geom_11_off + 129 * ccomps * dcomps);

            auto g_z_z_0_yyyy = cbuffer.data(sg_geom_11_off + 130 * ccomps * dcomps);

            auto g_z_z_0_yyyz = cbuffer.data(sg_geom_11_off + 131 * ccomps * dcomps);

            auto g_z_z_0_yyzz = cbuffer.data(sg_geom_11_off + 132 * ccomps * dcomps);

            auto g_z_z_0_yzzz = cbuffer.data(sg_geom_11_off + 133 * ccomps * dcomps);

            auto g_z_z_0_zzzz = cbuffer.data(sg_geom_11_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_xxxx, g_0_xxxy, g_0_xxxz, g_0_xxyy, g_0_xxyz, g_0_xxzz, g_0_xyyy, g_0_xyyz, g_0_xyzz, g_0_xzzz, g_0_yyyy, g_0_yyyz, g_0_yyzz, g_0_yzzz, g_0_z_0_xxxx, g_0_z_0_xxxxz, g_0_z_0_xxxy, g_0_z_0_xxxyz, g_0_z_0_xxxz, g_0_z_0_xxxzz, g_0_z_0_xxyy, g_0_z_0_xxyyz, g_0_z_0_xxyz, g_0_z_0_xxyzz, g_0_z_0_xxzz, g_0_z_0_xxzzz, g_0_z_0_xyyy, g_0_z_0_xyyyz, g_0_z_0_xyyz, g_0_z_0_xyyzz, g_0_z_0_xyzz, g_0_z_0_xyzzz, g_0_z_0_xzzz, g_0_z_0_xzzzz, g_0_z_0_yyyy, g_0_z_0_yyyyz, g_0_z_0_yyyz, g_0_z_0_yyyzz, g_0_z_0_yyzz, g_0_z_0_yyzzz, g_0_z_0_yzzz, g_0_z_0_yzzzz, g_0_z_0_zzzz, g_0_z_0_zzzzz, g_0_zzzz, g_z_z_0_xxxx, g_z_z_0_xxxy, g_z_z_0_xxxz, g_z_z_0_xxyy, g_z_z_0_xxyz, g_z_z_0_xxzz, g_z_z_0_xyyy, g_z_z_0_xyyz, g_z_z_0_xyzz, g_z_z_0_xzzz, g_z_z_0_yyyy, g_z_z_0_yyyz, g_z_z_0_yyzz, g_z_z_0_yzzz, g_z_z_0_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_0_xxxx[k] = g_0_xxxx[k] - g_0_z_0_xxxx[k] * ab_z + g_0_z_0_xxxxz[k];

                g_z_z_0_xxxy[k] = g_0_xxxy[k] - g_0_z_0_xxxy[k] * ab_z + g_0_z_0_xxxyz[k];

                g_z_z_0_xxxz[k] = g_0_xxxz[k] - g_0_z_0_xxxz[k] * ab_z + g_0_z_0_xxxzz[k];

                g_z_z_0_xxyy[k] = g_0_xxyy[k] - g_0_z_0_xxyy[k] * ab_z + g_0_z_0_xxyyz[k];

                g_z_z_0_xxyz[k] = g_0_xxyz[k] - g_0_z_0_xxyz[k] * ab_z + g_0_z_0_xxyzz[k];

                g_z_z_0_xxzz[k] = g_0_xxzz[k] - g_0_z_0_xxzz[k] * ab_z + g_0_z_0_xxzzz[k];

                g_z_z_0_xyyy[k] = g_0_xyyy[k] - g_0_z_0_xyyy[k] * ab_z + g_0_z_0_xyyyz[k];

                g_z_z_0_xyyz[k] = g_0_xyyz[k] - g_0_z_0_xyyz[k] * ab_z + g_0_z_0_xyyzz[k];

                g_z_z_0_xyzz[k] = g_0_xyzz[k] - g_0_z_0_xyzz[k] * ab_z + g_0_z_0_xyzzz[k];

                g_z_z_0_xzzz[k] = g_0_xzzz[k] - g_0_z_0_xzzz[k] * ab_z + g_0_z_0_xzzzz[k];

                g_z_z_0_yyyy[k] = g_0_yyyy[k] - g_0_z_0_yyyy[k] * ab_z + g_0_z_0_yyyyz[k];

                g_z_z_0_yyyz[k] = g_0_yyyz[k] - g_0_z_0_yyyz[k] * ab_z + g_0_z_0_yyyzz[k];

                g_z_z_0_yyzz[k] = g_0_yyzz[k] - g_0_z_0_yyzz[k] * ab_z + g_0_z_0_yyzzz[k];

                g_z_z_0_yzzz[k] = g_0_yzzz[k] - g_0_z_0_yzzz[k] * ab_z + g_0_z_0_yzzzz[k];

                g_z_z_0_zzzz[k] = g_0_zzzz[k] - g_0_z_0_zzzz[k] * ab_z + g_0_z_0_zzzzz[k];
            }
        }
    }
}

} // erirec namespace

