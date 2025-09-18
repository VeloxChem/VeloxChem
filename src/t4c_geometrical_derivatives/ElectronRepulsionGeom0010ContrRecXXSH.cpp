#include "ElectronRepulsionGeom0010ContrRecXXSH.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxsh(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxsh,
                                            CSimdArray<double>& pbuffer,
                                            const size_t idx_xxsh,
                                            const size_t idx_xxsi,
                                            const CSimdArray<double>& factors,
                                            const size_t idx_cd,
                                            const int a_angmom,
                                            const int b_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto acomps = tensor::number_of_cartesian_components(std::array<int, 1>{a_angmom,});

    const auto bcomps = tensor::number_of_cartesian_components(std::array<int, 1>{b_angmom,});

    // Set up R(CD) distances

    auto cd_x = factors.data(idx_cd);

    auto cd_y = factors.data(idx_cd + 1);

    auto cd_z = factors.data(idx_cd + 2);

    for (int i = 0; i < acomps; i++)
    {
        for (int j = 0; j < bcomps; j++)
        {
            /// Set up components of auxilary buffer : SSSH

            const auto sh_off = idx_xxsh + (i * bcomps + j) * 21;

            auto g_0_xxxxx = pbuffer.data(sh_off + 0);

            auto g_0_xxxxy = pbuffer.data(sh_off + 1);

            auto g_0_xxxxz = pbuffer.data(sh_off + 2);

            auto g_0_xxxyy = pbuffer.data(sh_off + 3);

            auto g_0_xxxyz = pbuffer.data(sh_off + 4);

            auto g_0_xxxzz = pbuffer.data(sh_off + 5);

            auto g_0_xxyyy = pbuffer.data(sh_off + 6);

            auto g_0_xxyyz = pbuffer.data(sh_off + 7);

            auto g_0_xxyzz = pbuffer.data(sh_off + 8);

            auto g_0_xxzzz = pbuffer.data(sh_off + 9);

            auto g_0_xyyyy = pbuffer.data(sh_off + 10);

            auto g_0_xyyyz = pbuffer.data(sh_off + 11);

            auto g_0_xyyzz = pbuffer.data(sh_off + 12);

            auto g_0_xyzzz = pbuffer.data(sh_off + 13);

            auto g_0_xzzzz = pbuffer.data(sh_off + 14);

            auto g_0_yyyyy = pbuffer.data(sh_off + 15);

            auto g_0_yyyyz = pbuffer.data(sh_off + 16);

            auto g_0_yyyzz = pbuffer.data(sh_off + 17);

            auto g_0_yyzzz = pbuffer.data(sh_off + 18);

            auto g_0_yzzzz = pbuffer.data(sh_off + 19);

            auto g_0_zzzzz = pbuffer.data(sh_off + 20);

            /// Set up components of auxilary buffer : SSSI

            const auto si_off = idx_xxsi + (i * bcomps + j) * 28;

            auto g_0_xxxxxx = pbuffer.data(si_off + 0);

            auto g_0_xxxxxy = pbuffer.data(si_off + 1);

            auto g_0_xxxxxz = pbuffer.data(si_off + 2);

            auto g_0_xxxxyy = pbuffer.data(si_off + 3);

            auto g_0_xxxxyz = pbuffer.data(si_off + 4);

            auto g_0_xxxxzz = pbuffer.data(si_off + 5);

            auto g_0_xxxyyy = pbuffer.data(si_off + 6);

            auto g_0_xxxyyz = pbuffer.data(si_off + 7);

            auto g_0_xxxyzz = pbuffer.data(si_off + 8);

            auto g_0_xxxzzz = pbuffer.data(si_off + 9);

            auto g_0_xxyyyy = pbuffer.data(si_off + 10);

            auto g_0_xxyyyz = pbuffer.data(si_off + 11);

            auto g_0_xxyyzz = pbuffer.data(si_off + 12);

            auto g_0_xxyzzz = pbuffer.data(si_off + 13);

            auto g_0_xxzzzz = pbuffer.data(si_off + 14);

            auto g_0_xyyyyy = pbuffer.data(si_off + 15);

            auto g_0_xyyyyz = pbuffer.data(si_off + 16);

            auto g_0_xyyyzz = pbuffer.data(si_off + 17);

            auto g_0_xyyzzz = pbuffer.data(si_off + 18);

            auto g_0_xyzzzz = pbuffer.data(si_off + 19);

            auto g_0_xzzzzz = pbuffer.data(si_off + 20);

            auto g_0_yyyyyy = pbuffer.data(si_off + 21);

            auto g_0_yyyyyz = pbuffer.data(si_off + 22);

            auto g_0_yyyyzz = pbuffer.data(si_off + 23);

            auto g_0_yyyzzz = pbuffer.data(si_off + 24);

            auto g_0_yyzzzz = pbuffer.data(si_off + 25);

            auto g_0_yzzzzz = pbuffer.data(si_off + 26);

            auto g_0_zzzzzz = pbuffer.data(si_off + 27);

            /// set up bra offset for contr_buffer_xxsh

            const auto sh_geom_10_off = idx_geom_10_xxsh + (i * bcomps + j) * 21;

            /// Set up 0-21 components of targeted buffer : cbuffer.data(

            auto g_x_0_0_xxxxx = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_0_xxxxy = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_0_xxxxz = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_0_xxxyy = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_0_xxxyz = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_0_xxxzz = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_0_xxyyy = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_0_xxyyz = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_0_xxyzz = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_0_xxzzz = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_0_xyyyy = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_0_xyyyz = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_0_xyyzz = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_0_xyzzz = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_0_xzzzz = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_0_yyyyy = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_0_yyyyz = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_0_yyyzz = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_0_yyzzz = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_0_yzzzz = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_0_zzzzz = cbuffer.data(sh_geom_10_off + 0 * acomps * bcomps + 20);

            #pragma omp simd aligned(cd_x, g_0_xxxxx, g_0_xxxxxx, g_0_xxxxxy, g_0_xxxxxz, g_0_xxxxy, g_0_xxxxyy, g_0_xxxxyz, g_0_xxxxz, g_0_xxxxzz, g_0_xxxyy, g_0_xxxyyy, g_0_xxxyyz, g_0_xxxyz, g_0_xxxyzz, g_0_xxxzz, g_0_xxxzzz, g_0_xxyyy, g_0_xxyyyy, g_0_xxyyyz, g_0_xxyyz, g_0_xxyyzz, g_0_xxyzz, g_0_xxyzzz, g_0_xxzzz, g_0_xxzzzz, g_0_xyyyy, g_0_xyyyyy, g_0_xyyyyz, g_0_xyyyz, g_0_xyyyzz, g_0_xyyzz, g_0_xyyzzz, g_0_xyzzz, g_0_xyzzzz, g_0_xzzzz, g_0_xzzzzz, g_0_yyyyy, g_0_yyyyz, g_0_yyyzz, g_0_yyzzz, g_0_yzzzz, g_0_zzzzz, g_x_0_0_xxxxx, g_x_0_0_xxxxy, g_x_0_0_xxxxz, g_x_0_0_xxxyy, g_x_0_0_xxxyz, g_x_0_0_xxxzz, g_x_0_0_xxyyy, g_x_0_0_xxyyz, g_x_0_0_xxyzz, g_x_0_0_xxzzz, g_x_0_0_xyyyy, g_x_0_0_xyyyz, g_x_0_0_xyyzz, g_x_0_0_xyzzz, g_x_0_0_xzzzz, g_x_0_0_yyyyy, g_x_0_0_yyyyz, g_x_0_0_yyyzz, g_x_0_0_yyzzz, g_x_0_0_yzzzz, g_x_0_0_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_0_xxxxx[k] = -g_0_xxxxx[k] * cd_x[k] + g_0_xxxxxx[k];

                g_x_0_0_xxxxy[k] = -g_0_xxxxy[k] * cd_x[k] + g_0_xxxxxy[k];

                g_x_0_0_xxxxz[k] = -g_0_xxxxz[k] * cd_x[k] + g_0_xxxxxz[k];

                g_x_0_0_xxxyy[k] = -g_0_xxxyy[k] * cd_x[k] + g_0_xxxxyy[k];

                g_x_0_0_xxxyz[k] = -g_0_xxxyz[k] * cd_x[k] + g_0_xxxxyz[k];

                g_x_0_0_xxxzz[k] = -g_0_xxxzz[k] * cd_x[k] + g_0_xxxxzz[k];

                g_x_0_0_xxyyy[k] = -g_0_xxyyy[k] * cd_x[k] + g_0_xxxyyy[k];

                g_x_0_0_xxyyz[k] = -g_0_xxyyz[k] * cd_x[k] + g_0_xxxyyz[k];

                g_x_0_0_xxyzz[k] = -g_0_xxyzz[k] * cd_x[k] + g_0_xxxyzz[k];

                g_x_0_0_xxzzz[k] = -g_0_xxzzz[k] * cd_x[k] + g_0_xxxzzz[k];

                g_x_0_0_xyyyy[k] = -g_0_xyyyy[k] * cd_x[k] + g_0_xxyyyy[k];

                g_x_0_0_xyyyz[k] = -g_0_xyyyz[k] * cd_x[k] + g_0_xxyyyz[k];

                g_x_0_0_xyyzz[k] = -g_0_xyyzz[k] * cd_x[k] + g_0_xxyyzz[k];

                g_x_0_0_xyzzz[k] = -g_0_xyzzz[k] * cd_x[k] + g_0_xxyzzz[k];

                g_x_0_0_xzzzz[k] = -g_0_xzzzz[k] * cd_x[k] + g_0_xxzzzz[k];

                g_x_0_0_yyyyy[k] = -g_0_yyyyy[k] * cd_x[k] + g_0_xyyyyy[k];

                g_x_0_0_yyyyz[k] = -g_0_yyyyz[k] * cd_x[k] + g_0_xyyyyz[k];

                g_x_0_0_yyyzz[k] = -g_0_yyyzz[k] * cd_x[k] + g_0_xyyyzz[k];

                g_x_0_0_yyzzz[k] = -g_0_yyzzz[k] * cd_x[k] + g_0_xyyzzz[k];

                g_x_0_0_yzzzz[k] = -g_0_yzzzz[k] * cd_x[k] + g_0_xyzzzz[k];

                g_x_0_0_zzzzz[k] = -g_0_zzzzz[k] * cd_x[k] + g_0_xzzzzz[k];
            }
            /// Set up 0-21 components of targeted buffer : cbuffer.data(

            auto g_y_0_0_xxxxx = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 0);

            auto g_y_0_0_xxxxy = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 1);

            auto g_y_0_0_xxxxz = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 2);

            auto g_y_0_0_xxxyy = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 3);

            auto g_y_0_0_xxxyz = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 4);

            auto g_y_0_0_xxxzz = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 5);

            auto g_y_0_0_xxyyy = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 6);

            auto g_y_0_0_xxyyz = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 7);

            auto g_y_0_0_xxyzz = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 8);

            auto g_y_0_0_xxzzz = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 9);

            auto g_y_0_0_xyyyy = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 10);

            auto g_y_0_0_xyyyz = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 11);

            auto g_y_0_0_xyyzz = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 12);

            auto g_y_0_0_xyzzz = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 13);

            auto g_y_0_0_xzzzz = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 14);

            auto g_y_0_0_yyyyy = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 15);

            auto g_y_0_0_yyyyz = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 16);

            auto g_y_0_0_yyyzz = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 17);

            auto g_y_0_0_yyzzz = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 18);

            auto g_y_0_0_yzzzz = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 19);

            auto g_y_0_0_zzzzz = cbuffer.data(sh_geom_10_off + 21 * acomps * bcomps + 20);

            #pragma omp simd aligned(cd_y, g_0_xxxxx, g_0_xxxxxy, g_0_xxxxy, g_0_xxxxyy, g_0_xxxxyz, g_0_xxxxz, g_0_xxxyy, g_0_xxxyyy, g_0_xxxyyz, g_0_xxxyz, g_0_xxxyzz, g_0_xxxzz, g_0_xxyyy, g_0_xxyyyy, g_0_xxyyyz, g_0_xxyyz, g_0_xxyyzz, g_0_xxyzz, g_0_xxyzzz, g_0_xxzzz, g_0_xyyyy, g_0_xyyyyy, g_0_xyyyyz, g_0_xyyyz, g_0_xyyyzz, g_0_xyyzz, g_0_xyyzzz, g_0_xyzzz, g_0_xyzzzz, g_0_xzzzz, g_0_yyyyy, g_0_yyyyyy, g_0_yyyyyz, g_0_yyyyz, g_0_yyyyzz, g_0_yyyzz, g_0_yyyzzz, g_0_yyzzz, g_0_yyzzzz, g_0_yzzzz, g_0_yzzzzz, g_0_zzzzz, g_y_0_0_xxxxx, g_y_0_0_xxxxy, g_y_0_0_xxxxz, g_y_0_0_xxxyy, g_y_0_0_xxxyz, g_y_0_0_xxxzz, g_y_0_0_xxyyy, g_y_0_0_xxyyz, g_y_0_0_xxyzz, g_y_0_0_xxzzz, g_y_0_0_xyyyy, g_y_0_0_xyyyz, g_y_0_0_xyyzz, g_y_0_0_xyzzz, g_y_0_0_xzzzz, g_y_0_0_yyyyy, g_y_0_0_yyyyz, g_y_0_0_yyyzz, g_y_0_0_yyzzz, g_y_0_0_yzzzz, g_y_0_0_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_0_xxxxx[k] = -g_0_xxxxx[k] * cd_y[k] + g_0_xxxxxy[k];

                g_y_0_0_xxxxy[k] = -g_0_xxxxy[k] * cd_y[k] + g_0_xxxxyy[k];

                g_y_0_0_xxxxz[k] = -g_0_xxxxz[k] * cd_y[k] + g_0_xxxxyz[k];

                g_y_0_0_xxxyy[k] = -g_0_xxxyy[k] * cd_y[k] + g_0_xxxyyy[k];

                g_y_0_0_xxxyz[k] = -g_0_xxxyz[k] * cd_y[k] + g_0_xxxyyz[k];

                g_y_0_0_xxxzz[k] = -g_0_xxxzz[k] * cd_y[k] + g_0_xxxyzz[k];

                g_y_0_0_xxyyy[k] = -g_0_xxyyy[k] * cd_y[k] + g_0_xxyyyy[k];

                g_y_0_0_xxyyz[k] = -g_0_xxyyz[k] * cd_y[k] + g_0_xxyyyz[k];

                g_y_0_0_xxyzz[k] = -g_0_xxyzz[k] * cd_y[k] + g_0_xxyyzz[k];

                g_y_0_0_xxzzz[k] = -g_0_xxzzz[k] * cd_y[k] + g_0_xxyzzz[k];

                g_y_0_0_xyyyy[k] = -g_0_xyyyy[k] * cd_y[k] + g_0_xyyyyy[k];

                g_y_0_0_xyyyz[k] = -g_0_xyyyz[k] * cd_y[k] + g_0_xyyyyz[k];

                g_y_0_0_xyyzz[k] = -g_0_xyyzz[k] * cd_y[k] + g_0_xyyyzz[k];

                g_y_0_0_xyzzz[k] = -g_0_xyzzz[k] * cd_y[k] + g_0_xyyzzz[k];

                g_y_0_0_xzzzz[k] = -g_0_xzzzz[k] * cd_y[k] + g_0_xyzzzz[k];

                g_y_0_0_yyyyy[k] = -g_0_yyyyy[k] * cd_y[k] + g_0_yyyyyy[k];

                g_y_0_0_yyyyz[k] = -g_0_yyyyz[k] * cd_y[k] + g_0_yyyyyz[k];

                g_y_0_0_yyyzz[k] = -g_0_yyyzz[k] * cd_y[k] + g_0_yyyyzz[k];

                g_y_0_0_yyzzz[k] = -g_0_yyzzz[k] * cd_y[k] + g_0_yyyzzz[k];

                g_y_0_0_yzzzz[k] = -g_0_yzzzz[k] * cd_y[k] + g_0_yyzzzz[k];

                g_y_0_0_zzzzz[k] = -g_0_zzzzz[k] * cd_y[k] + g_0_yzzzzz[k];
            }
            /// Set up 0-21 components of targeted buffer : cbuffer.data(

            auto g_z_0_0_xxxxx = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 0);

            auto g_z_0_0_xxxxy = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 1);

            auto g_z_0_0_xxxxz = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 2);

            auto g_z_0_0_xxxyy = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 3);

            auto g_z_0_0_xxxyz = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 4);

            auto g_z_0_0_xxxzz = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 5);

            auto g_z_0_0_xxyyy = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 6);

            auto g_z_0_0_xxyyz = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 7);

            auto g_z_0_0_xxyzz = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 8);

            auto g_z_0_0_xxzzz = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 9);

            auto g_z_0_0_xyyyy = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 10);

            auto g_z_0_0_xyyyz = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 11);

            auto g_z_0_0_xyyzz = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 12);

            auto g_z_0_0_xyzzz = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 13);

            auto g_z_0_0_xzzzz = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 14);

            auto g_z_0_0_yyyyy = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 15);

            auto g_z_0_0_yyyyz = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 16);

            auto g_z_0_0_yyyzz = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 17);

            auto g_z_0_0_yyzzz = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 18);

            auto g_z_0_0_yzzzz = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 19);

            auto g_z_0_0_zzzzz = cbuffer.data(sh_geom_10_off + 42 * acomps * bcomps + 20);

            #pragma omp simd aligned(cd_z, g_0_xxxxx, g_0_xxxxxz, g_0_xxxxy, g_0_xxxxyz, g_0_xxxxz, g_0_xxxxzz, g_0_xxxyy, g_0_xxxyyz, g_0_xxxyz, g_0_xxxyzz, g_0_xxxzz, g_0_xxxzzz, g_0_xxyyy, g_0_xxyyyz, g_0_xxyyz, g_0_xxyyzz, g_0_xxyzz, g_0_xxyzzz, g_0_xxzzz, g_0_xxzzzz, g_0_xyyyy, g_0_xyyyyz, g_0_xyyyz, g_0_xyyyzz, g_0_xyyzz, g_0_xyyzzz, g_0_xyzzz, g_0_xyzzzz, g_0_xzzzz, g_0_xzzzzz, g_0_yyyyy, g_0_yyyyyz, g_0_yyyyz, g_0_yyyyzz, g_0_yyyzz, g_0_yyyzzz, g_0_yyzzz, g_0_yyzzzz, g_0_yzzzz, g_0_yzzzzz, g_0_zzzzz, g_0_zzzzzz, g_z_0_0_xxxxx, g_z_0_0_xxxxy, g_z_0_0_xxxxz, g_z_0_0_xxxyy, g_z_0_0_xxxyz, g_z_0_0_xxxzz, g_z_0_0_xxyyy, g_z_0_0_xxyyz, g_z_0_0_xxyzz, g_z_0_0_xxzzz, g_z_0_0_xyyyy, g_z_0_0_xyyyz, g_z_0_0_xyyzz, g_z_0_0_xyzzz, g_z_0_0_xzzzz, g_z_0_0_yyyyy, g_z_0_0_yyyyz, g_z_0_0_yyyzz, g_z_0_0_yyzzz, g_z_0_0_yzzzz, g_z_0_0_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_0_xxxxx[k] = -g_0_xxxxx[k] * cd_z[k] + g_0_xxxxxz[k];

                g_z_0_0_xxxxy[k] = -g_0_xxxxy[k] * cd_z[k] + g_0_xxxxyz[k];

                g_z_0_0_xxxxz[k] = -g_0_xxxxz[k] * cd_z[k] + g_0_xxxxzz[k];

                g_z_0_0_xxxyy[k] = -g_0_xxxyy[k] * cd_z[k] + g_0_xxxyyz[k];

                g_z_0_0_xxxyz[k] = -g_0_xxxyz[k] * cd_z[k] + g_0_xxxyzz[k];

                g_z_0_0_xxxzz[k] = -g_0_xxxzz[k] * cd_z[k] + g_0_xxxzzz[k];

                g_z_0_0_xxyyy[k] = -g_0_xxyyy[k] * cd_z[k] + g_0_xxyyyz[k];

                g_z_0_0_xxyyz[k] = -g_0_xxyyz[k] * cd_z[k] + g_0_xxyyzz[k];

                g_z_0_0_xxyzz[k] = -g_0_xxyzz[k] * cd_z[k] + g_0_xxyzzz[k];

                g_z_0_0_xxzzz[k] = -g_0_xxzzz[k] * cd_z[k] + g_0_xxzzzz[k];

                g_z_0_0_xyyyy[k] = -g_0_xyyyy[k] * cd_z[k] + g_0_xyyyyz[k];

                g_z_0_0_xyyyz[k] = -g_0_xyyyz[k] * cd_z[k] + g_0_xyyyzz[k];

                g_z_0_0_xyyzz[k] = -g_0_xyyzz[k] * cd_z[k] + g_0_xyyzzz[k];

                g_z_0_0_xyzzz[k] = -g_0_xyzzz[k] * cd_z[k] + g_0_xyzzzz[k];

                g_z_0_0_xzzzz[k] = -g_0_xzzzz[k] * cd_z[k] + g_0_xzzzzz[k];

                g_z_0_0_yyyyy[k] = -g_0_yyyyy[k] * cd_z[k] + g_0_yyyyyz[k];

                g_z_0_0_yyyyz[k] = -g_0_yyyyz[k] * cd_z[k] + g_0_yyyyzz[k];

                g_z_0_0_yyyzz[k] = -g_0_yyyzz[k] * cd_z[k] + g_0_yyyzzz[k];

                g_z_0_0_yyzzz[k] = -g_0_yyzzz[k] * cd_z[k] + g_0_yyzzzz[k];

                g_z_0_0_yzzzz[k] = -g_0_yzzzz[k] * cd_z[k] + g_0_yzzzzz[k];

                g_z_0_0_zzzzz[k] = -g_0_zzzzz[k] * cd_z[k] + g_0_zzzzzz[k];
            }
        }
    }
}

} // erirec namespace

