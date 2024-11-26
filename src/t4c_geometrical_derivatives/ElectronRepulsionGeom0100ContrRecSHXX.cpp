#include "ElectronRepulsionGeom0100ContrRecSHXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_shxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_shxx,
                                            const size_t idx_sgxx,
                                            const size_t idx_sixx,
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

            /// Set up components of auxilary buffer : SISS

            const auto si_off = idx_sixx + i * dcomps + j;

            auto g_0_xxxxxx = cbuffer.data(si_off + 0 * ccomps * dcomps);

            auto g_0_xxxxxy = cbuffer.data(si_off + 1 * ccomps * dcomps);

            auto g_0_xxxxxz = cbuffer.data(si_off + 2 * ccomps * dcomps);

            auto g_0_xxxxyy = cbuffer.data(si_off + 3 * ccomps * dcomps);

            auto g_0_xxxxyz = cbuffer.data(si_off + 4 * ccomps * dcomps);

            auto g_0_xxxxzz = cbuffer.data(si_off + 5 * ccomps * dcomps);

            auto g_0_xxxyyy = cbuffer.data(si_off + 6 * ccomps * dcomps);

            auto g_0_xxxyyz = cbuffer.data(si_off + 7 * ccomps * dcomps);

            auto g_0_xxxyzz = cbuffer.data(si_off + 8 * ccomps * dcomps);

            auto g_0_xxxzzz = cbuffer.data(si_off + 9 * ccomps * dcomps);

            auto g_0_xxyyyy = cbuffer.data(si_off + 10 * ccomps * dcomps);

            auto g_0_xxyyyz = cbuffer.data(si_off + 11 * ccomps * dcomps);

            auto g_0_xxyyzz = cbuffer.data(si_off + 12 * ccomps * dcomps);

            auto g_0_xxyzzz = cbuffer.data(si_off + 13 * ccomps * dcomps);

            auto g_0_xxzzzz = cbuffer.data(si_off + 14 * ccomps * dcomps);

            auto g_0_xyyyyy = cbuffer.data(si_off + 15 * ccomps * dcomps);

            auto g_0_xyyyyz = cbuffer.data(si_off + 16 * ccomps * dcomps);

            auto g_0_xyyyzz = cbuffer.data(si_off + 17 * ccomps * dcomps);

            auto g_0_xyyzzz = cbuffer.data(si_off + 18 * ccomps * dcomps);

            auto g_0_xyzzzz = cbuffer.data(si_off + 19 * ccomps * dcomps);

            auto g_0_xzzzzz = cbuffer.data(si_off + 20 * ccomps * dcomps);

            auto g_0_yyyyyy = cbuffer.data(si_off + 21 * ccomps * dcomps);

            auto g_0_yyyyyz = cbuffer.data(si_off + 22 * ccomps * dcomps);

            auto g_0_yyyyzz = cbuffer.data(si_off + 23 * ccomps * dcomps);

            auto g_0_yyyzzz = cbuffer.data(si_off + 24 * ccomps * dcomps);

            auto g_0_yyzzzz = cbuffer.data(si_off + 25 * ccomps * dcomps);

            auto g_0_yzzzzz = cbuffer.data(si_off + 26 * ccomps * dcomps);

            auto g_0_zzzzzz = cbuffer.data(si_off + 27 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_shxx

            const auto sh_geom_01_off = idx_geom_01_shxx + i * dcomps + j;

            /// Set up 0-21 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_x_0_xxxxx, g_0_x_0_xxxxy, g_0_x_0_xxxxz, g_0_x_0_xxxyy, g_0_x_0_xxxyz, g_0_x_0_xxxzz, g_0_x_0_xxyyy, g_0_x_0_xxyyz, g_0_x_0_xxyzz, g_0_x_0_xxzzz, g_0_x_0_xyyyy, g_0_x_0_xyyyz, g_0_x_0_xyyzz, g_0_x_0_xyzzz, g_0_x_0_xzzzz, g_0_x_0_yyyyy, g_0_x_0_yyyyz, g_0_x_0_yyyzz, g_0_x_0_yyzzz, g_0_x_0_yzzzz, g_0_x_0_zzzzz, g_0_xxxx, g_0_xxxxxx, g_0_xxxxxy, g_0_xxxxxz, g_0_xxxxyy, g_0_xxxxyz, g_0_xxxxzz, g_0_xxxy, g_0_xxxyyy, g_0_xxxyyz, g_0_xxxyzz, g_0_xxxz, g_0_xxxzzz, g_0_xxyy, g_0_xxyyyy, g_0_xxyyyz, g_0_xxyyzz, g_0_xxyz, g_0_xxyzzz, g_0_xxzz, g_0_xxzzzz, g_0_xyyy, g_0_xyyyyy, g_0_xyyyyz, g_0_xyyyzz, g_0_xyyz, g_0_xyyzzz, g_0_xyzz, g_0_xyzzzz, g_0_xzzz, g_0_xzzzzz, g_0_yyyy, g_0_yyyz, g_0_yyzz, g_0_yzzz, g_0_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_0_xxxxx[k] = -5.0 * g_0_xxxx[k] + g_0_xxxxxx[k];

                g_0_x_0_xxxxy[k] = -4.0 * g_0_xxxy[k] + g_0_xxxxxy[k];

                g_0_x_0_xxxxz[k] = -4.0 * g_0_xxxz[k] + g_0_xxxxxz[k];

                g_0_x_0_xxxyy[k] = -3.0 * g_0_xxyy[k] + g_0_xxxxyy[k];

                g_0_x_0_xxxyz[k] = -3.0 * g_0_xxyz[k] + g_0_xxxxyz[k];

                g_0_x_0_xxxzz[k] = -3.0 * g_0_xxzz[k] + g_0_xxxxzz[k];

                g_0_x_0_xxyyy[k] = -2.0 * g_0_xyyy[k] + g_0_xxxyyy[k];

                g_0_x_0_xxyyz[k] = -2.0 * g_0_xyyz[k] + g_0_xxxyyz[k];

                g_0_x_0_xxyzz[k] = -2.0 * g_0_xyzz[k] + g_0_xxxyzz[k];

                g_0_x_0_xxzzz[k] = -2.0 * g_0_xzzz[k] + g_0_xxxzzz[k];

                g_0_x_0_xyyyy[k] = -g_0_yyyy[k] + g_0_xxyyyy[k];

                g_0_x_0_xyyyz[k] = -g_0_yyyz[k] + g_0_xxyyyz[k];

                g_0_x_0_xyyzz[k] = -g_0_yyzz[k] + g_0_xxyyzz[k];

                g_0_x_0_xyzzz[k] = -g_0_yzzz[k] + g_0_xxyzzz[k];

                g_0_x_0_xzzzz[k] = -g_0_zzzz[k] + g_0_xxzzzz[k];

                g_0_x_0_yyyyy[k] = g_0_xyyyyy[k];

                g_0_x_0_yyyyz[k] = g_0_xyyyyz[k];

                g_0_x_0_yyyzz[k] = g_0_xyyyzz[k];

                g_0_x_0_yyzzz[k] = g_0_xyyzzz[k];

                g_0_x_0_yzzzz[k] = g_0_xyzzzz[k];

                g_0_x_0_zzzzz[k] = g_0_xzzzzz[k];
            }

            /// Set up 21-42 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_xxxx, g_0_xxxxxy, g_0_xxxxyy, g_0_xxxxyz, g_0_xxxy, g_0_xxxyyy, g_0_xxxyyz, g_0_xxxyzz, g_0_xxxz, g_0_xxyy, g_0_xxyyyy, g_0_xxyyyz, g_0_xxyyzz, g_0_xxyz, g_0_xxyzzz, g_0_xxzz, g_0_xyyy, g_0_xyyyyy, g_0_xyyyyz, g_0_xyyyzz, g_0_xyyz, g_0_xyyzzz, g_0_xyzz, g_0_xyzzzz, g_0_xzzz, g_0_y_0_xxxxx, g_0_y_0_xxxxy, g_0_y_0_xxxxz, g_0_y_0_xxxyy, g_0_y_0_xxxyz, g_0_y_0_xxxzz, g_0_y_0_xxyyy, g_0_y_0_xxyyz, g_0_y_0_xxyzz, g_0_y_0_xxzzz, g_0_y_0_xyyyy, g_0_y_0_xyyyz, g_0_y_0_xyyzz, g_0_y_0_xyzzz, g_0_y_0_xzzzz, g_0_y_0_yyyyy, g_0_y_0_yyyyz, g_0_y_0_yyyzz, g_0_y_0_yyzzz, g_0_y_0_yzzzz, g_0_y_0_zzzzz, g_0_yyyy, g_0_yyyyyy, g_0_yyyyyz, g_0_yyyyzz, g_0_yyyz, g_0_yyyzzz, g_0_yyzz, g_0_yyzzzz, g_0_yzzz, g_0_yzzzzz, g_0_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_0_xxxxx[k] = g_0_xxxxxy[k];

                g_0_y_0_xxxxy[k] = -g_0_xxxx[k] + g_0_xxxxyy[k];

                g_0_y_0_xxxxz[k] = g_0_xxxxyz[k];

                g_0_y_0_xxxyy[k] = -2.0 * g_0_xxxy[k] + g_0_xxxyyy[k];

                g_0_y_0_xxxyz[k] = -g_0_xxxz[k] + g_0_xxxyyz[k];

                g_0_y_0_xxxzz[k] = g_0_xxxyzz[k];

                g_0_y_0_xxyyy[k] = -3.0 * g_0_xxyy[k] + g_0_xxyyyy[k];

                g_0_y_0_xxyyz[k] = -2.0 * g_0_xxyz[k] + g_0_xxyyyz[k];

                g_0_y_0_xxyzz[k] = -g_0_xxzz[k] + g_0_xxyyzz[k];

                g_0_y_0_xxzzz[k] = g_0_xxyzzz[k];

                g_0_y_0_xyyyy[k] = -4.0 * g_0_xyyy[k] + g_0_xyyyyy[k];

                g_0_y_0_xyyyz[k] = -3.0 * g_0_xyyz[k] + g_0_xyyyyz[k];

                g_0_y_0_xyyzz[k] = -2.0 * g_0_xyzz[k] + g_0_xyyyzz[k];

                g_0_y_0_xyzzz[k] = -g_0_xzzz[k] + g_0_xyyzzz[k];

                g_0_y_0_xzzzz[k] = g_0_xyzzzz[k];

                g_0_y_0_yyyyy[k] = -5.0 * g_0_yyyy[k] + g_0_yyyyyy[k];

                g_0_y_0_yyyyz[k] = -4.0 * g_0_yyyz[k] + g_0_yyyyyz[k];

                g_0_y_0_yyyzz[k] = -3.0 * g_0_yyzz[k] + g_0_yyyyzz[k];

                g_0_y_0_yyzzz[k] = -2.0 * g_0_yzzz[k] + g_0_yyyzzz[k];

                g_0_y_0_yzzzz[k] = -g_0_zzzz[k] + g_0_yyzzzz[k];

                g_0_y_0_zzzzz[k] = g_0_yzzzzz[k];
            }

            /// Set up 42-63 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(g_0_xxxx, g_0_xxxxxz, g_0_xxxxyz, g_0_xxxxzz, g_0_xxxy, g_0_xxxyyz, g_0_xxxyzz, g_0_xxxz, g_0_xxxzzz, g_0_xxyy, g_0_xxyyyz, g_0_xxyyzz, g_0_xxyz, g_0_xxyzzz, g_0_xxzz, g_0_xxzzzz, g_0_xyyy, g_0_xyyyyz, g_0_xyyyzz, g_0_xyyz, g_0_xyyzzz, g_0_xyzz, g_0_xyzzzz, g_0_xzzz, g_0_xzzzzz, g_0_yyyy, g_0_yyyyyz, g_0_yyyyzz, g_0_yyyz, g_0_yyyzzz, g_0_yyzz, g_0_yyzzzz, g_0_yzzz, g_0_yzzzzz, g_0_z_0_xxxxx, g_0_z_0_xxxxy, g_0_z_0_xxxxz, g_0_z_0_xxxyy, g_0_z_0_xxxyz, g_0_z_0_xxxzz, g_0_z_0_xxyyy, g_0_z_0_xxyyz, g_0_z_0_xxyzz, g_0_z_0_xxzzz, g_0_z_0_xyyyy, g_0_z_0_xyyyz, g_0_z_0_xyyzz, g_0_z_0_xyzzz, g_0_z_0_xzzzz, g_0_z_0_yyyyy, g_0_z_0_yyyyz, g_0_z_0_yyyzz, g_0_z_0_yyzzz, g_0_z_0_yzzzz, g_0_z_0_zzzzz, g_0_zzzz, g_0_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_0_xxxxx[k] = g_0_xxxxxz[k];

                g_0_z_0_xxxxy[k] = g_0_xxxxyz[k];

                g_0_z_0_xxxxz[k] = -g_0_xxxx[k] + g_0_xxxxzz[k];

                g_0_z_0_xxxyy[k] = g_0_xxxyyz[k];

                g_0_z_0_xxxyz[k] = -g_0_xxxy[k] + g_0_xxxyzz[k];

                g_0_z_0_xxxzz[k] = -2.0 * g_0_xxxz[k] + g_0_xxxzzz[k];

                g_0_z_0_xxyyy[k] = g_0_xxyyyz[k];

                g_0_z_0_xxyyz[k] = -g_0_xxyy[k] + g_0_xxyyzz[k];

                g_0_z_0_xxyzz[k] = -2.0 * g_0_xxyz[k] + g_0_xxyzzz[k];

                g_0_z_0_xxzzz[k] = -3.0 * g_0_xxzz[k] + g_0_xxzzzz[k];

                g_0_z_0_xyyyy[k] = g_0_xyyyyz[k];

                g_0_z_0_xyyyz[k] = -g_0_xyyy[k] + g_0_xyyyzz[k];

                g_0_z_0_xyyzz[k] = -2.0 * g_0_xyyz[k] + g_0_xyyzzz[k];

                g_0_z_0_xyzzz[k] = -3.0 * g_0_xyzz[k] + g_0_xyzzzz[k];

                g_0_z_0_xzzzz[k] = -4.0 * g_0_xzzz[k] + g_0_xzzzzz[k];

                g_0_z_0_yyyyy[k] = g_0_yyyyyz[k];

                g_0_z_0_yyyyz[k] = -g_0_yyyy[k] + g_0_yyyyzz[k];

                g_0_z_0_yyyzz[k] = -2.0 * g_0_yyyz[k] + g_0_yyyzzz[k];

                g_0_z_0_yyzzz[k] = -3.0 * g_0_yyzz[k] + g_0_yyzzzz[k];

                g_0_z_0_yzzzz[k] = -4.0 * g_0_yzzz[k] + g_0_yzzzzz[k];

                g_0_z_0_zzzzz[k] = -5.0 * g_0_zzzz[k] + g_0_zzzzzz[k];
            }
        }
    }
}

} // erirec namespace

