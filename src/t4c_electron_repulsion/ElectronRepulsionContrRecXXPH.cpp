#include "ElectronRepulsionContrRecXXPH.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_hrr_electron_repulsion_xxph(CSimdArray<double>& cbuffer,
                                     const size_t idx_xxph,
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

            /// set up bra offset for contr_buffer_xxph

            const auto ph_off = idx_xxph + (i * bcomps + j) * 63;

            /// Set up 0-21 components of targeted buffer : cbuffer.data(

            auto g_x_xxxxx = cbuffer.data(ph_off + 0);

            auto g_x_xxxxy = cbuffer.data(ph_off + 1);

            auto g_x_xxxxz = cbuffer.data(ph_off + 2);

            auto g_x_xxxyy = cbuffer.data(ph_off + 3);

            auto g_x_xxxyz = cbuffer.data(ph_off + 4);

            auto g_x_xxxzz = cbuffer.data(ph_off + 5);

            auto g_x_xxyyy = cbuffer.data(ph_off + 6);

            auto g_x_xxyyz = cbuffer.data(ph_off + 7);

            auto g_x_xxyzz = cbuffer.data(ph_off + 8);

            auto g_x_xxzzz = cbuffer.data(ph_off + 9);

            auto g_x_xyyyy = cbuffer.data(ph_off + 10);

            auto g_x_xyyyz = cbuffer.data(ph_off + 11);

            auto g_x_xyyzz = cbuffer.data(ph_off + 12);

            auto g_x_xyzzz = cbuffer.data(ph_off + 13);

            auto g_x_xzzzz = cbuffer.data(ph_off + 14);

            auto g_x_yyyyy = cbuffer.data(ph_off + 15);

            auto g_x_yyyyz = cbuffer.data(ph_off + 16);

            auto g_x_yyyzz = cbuffer.data(ph_off + 17);

            auto g_x_yyzzz = cbuffer.data(ph_off + 18);

            auto g_x_yzzzz = cbuffer.data(ph_off + 19);

            auto g_x_zzzzz = cbuffer.data(ph_off + 20);

            #pragma omp simd aligned(cd_x, g_0_xxxxx, g_0_xxxxxx, g_0_xxxxxy, g_0_xxxxxz, g_0_xxxxy, g_0_xxxxyy, g_0_xxxxyz, g_0_xxxxz, g_0_xxxxzz, g_0_xxxyy, g_0_xxxyyy, g_0_xxxyyz, g_0_xxxyz, g_0_xxxyzz, g_0_xxxzz, g_0_xxxzzz, g_0_xxyyy, g_0_xxyyyy, g_0_xxyyyz, g_0_xxyyz, g_0_xxyyzz, g_0_xxyzz, g_0_xxyzzz, g_0_xxzzz, g_0_xxzzzz, g_0_xyyyy, g_0_xyyyyy, g_0_xyyyyz, g_0_xyyyz, g_0_xyyyzz, g_0_xyyzz, g_0_xyyzzz, g_0_xyzzz, g_0_xyzzzz, g_0_xzzzz, g_0_xzzzzz, g_0_yyyyy, g_0_yyyyz, g_0_yyyzz, g_0_yyzzz, g_0_yzzzz, g_0_zzzzz, g_x_xxxxx, g_x_xxxxy, g_x_xxxxz, g_x_xxxyy, g_x_xxxyz, g_x_xxxzz, g_x_xxyyy, g_x_xxyyz, g_x_xxyzz, g_x_xxzzz, g_x_xyyyy, g_x_xyyyz, g_x_xyyzz, g_x_xyzzz, g_x_xzzzz, g_x_yyyyy, g_x_yyyyz, g_x_yyyzz, g_x_yyzzz, g_x_yzzzz, g_x_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_xxxxx[k] = -g_0_xxxxx[k] * cd_x[k] + g_0_xxxxxx[k];

                g_x_xxxxy[k] = -g_0_xxxxy[k] * cd_x[k] + g_0_xxxxxy[k];

                g_x_xxxxz[k] = -g_0_xxxxz[k] * cd_x[k] + g_0_xxxxxz[k];

                g_x_xxxyy[k] = -g_0_xxxyy[k] * cd_x[k] + g_0_xxxxyy[k];

                g_x_xxxyz[k] = -g_0_xxxyz[k] * cd_x[k] + g_0_xxxxyz[k];

                g_x_xxxzz[k] = -g_0_xxxzz[k] * cd_x[k] + g_0_xxxxzz[k];

                g_x_xxyyy[k] = -g_0_xxyyy[k] * cd_x[k] + g_0_xxxyyy[k];

                g_x_xxyyz[k] = -g_0_xxyyz[k] * cd_x[k] + g_0_xxxyyz[k];

                g_x_xxyzz[k] = -g_0_xxyzz[k] * cd_x[k] + g_0_xxxyzz[k];

                g_x_xxzzz[k] = -g_0_xxzzz[k] * cd_x[k] + g_0_xxxzzz[k];

                g_x_xyyyy[k] = -g_0_xyyyy[k] * cd_x[k] + g_0_xxyyyy[k];

                g_x_xyyyz[k] = -g_0_xyyyz[k] * cd_x[k] + g_0_xxyyyz[k];

                g_x_xyyzz[k] = -g_0_xyyzz[k] * cd_x[k] + g_0_xxyyzz[k];

                g_x_xyzzz[k] = -g_0_xyzzz[k] * cd_x[k] + g_0_xxyzzz[k];

                g_x_xzzzz[k] = -g_0_xzzzz[k] * cd_x[k] + g_0_xxzzzz[k];

                g_x_yyyyy[k] = -g_0_yyyyy[k] * cd_x[k] + g_0_xyyyyy[k];

                g_x_yyyyz[k] = -g_0_yyyyz[k] * cd_x[k] + g_0_xyyyyz[k];

                g_x_yyyzz[k] = -g_0_yyyzz[k] * cd_x[k] + g_0_xyyyzz[k];

                g_x_yyzzz[k] = -g_0_yyzzz[k] * cd_x[k] + g_0_xyyzzz[k];

                g_x_yzzzz[k] = -g_0_yzzzz[k] * cd_x[k] + g_0_xyzzzz[k];

                g_x_zzzzz[k] = -g_0_zzzzz[k] * cd_x[k] + g_0_xzzzzz[k];
            }

            /// Set up 21-42 components of targeted buffer : cbuffer.data(

            auto g_y_xxxxx = cbuffer.data(ph_off + 21);

            auto g_y_xxxxy = cbuffer.data(ph_off + 22);

            auto g_y_xxxxz = cbuffer.data(ph_off + 23);

            auto g_y_xxxyy = cbuffer.data(ph_off + 24);

            auto g_y_xxxyz = cbuffer.data(ph_off + 25);

            auto g_y_xxxzz = cbuffer.data(ph_off + 26);

            auto g_y_xxyyy = cbuffer.data(ph_off + 27);

            auto g_y_xxyyz = cbuffer.data(ph_off + 28);

            auto g_y_xxyzz = cbuffer.data(ph_off + 29);

            auto g_y_xxzzz = cbuffer.data(ph_off + 30);

            auto g_y_xyyyy = cbuffer.data(ph_off + 31);

            auto g_y_xyyyz = cbuffer.data(ph_off + 32);

            auto g_y_xyyzz = cbuffer.data(ph_off + 33);

            auto g_y_xyzzz = cbuffer.data(ph_off + 34);

            auto g_y_xzzzz = cbuffer.data(ph_off + 35);

            auto g_y_yyyyy = cbuffer.data(ph_off + 36);

            auto g_y_yyyyz = cbuffer.data(ph_off + 37);

            auto g_y_yyyzz = cbuffer.data(ph_off + 38);

            auto g_y_yyzzz = cbuffer.data(ph_off + 39);

            auto g_y_yzzzz = cbuffer.data(ph_off + 40);

            auto g_y_zzzzz = cbuffer.data(ph_off + 41);

            #pragma omp simd aligned(cd_y, g_0_xxxxx, g_0_xxxxxy, g_0_xxxxy, g_0_xxxxyy, g_0_xxxxyz, g_0_xxxxz, g_0_xxxyy, g_0_xxxyyy, g_0_xxxyyz, g_0_xxxyz, g_0_xxxyzz, g_0_xxxzz, g_0_xxyyy, g_0_xxyyyy, g_0_xxyyyz, g_0_xxyyz, g_0_xxyyzz, g_0_xxyzz, g_0_xxyzzz, g_0_xxzzz, g_0_xyyyy, g_0_xyyyyy, g_0_xyyyyz, g_0_xyyyz, g_0_xyyyzz, g_0_xyyzz, g_0_xyyzzz, g_0_xyzzz, g_0_xyzzzz, g_0_xzzzz, g_0_yyyyy, g_0_yyyyyy, g_0_yyyyyz, g_0_yyyyz, g_0_yyyyzz, g_0_yyyzz, g_0_yyyzzz, g_0_yyzzz, g_0_yyzzzz, g_0_yzzzz, g_0_yzzzzz, g_0_zzzzz, g_y_xxxxx, g_y_xxxxy, g_y_xxxxz, g_y_xxxyy, g_y_xxxyz, g_y_xxxzz, g_y_xxyyy, g_y_xxyyz, g_y_xxyzz, g_y_xxzzz, g_y_xyyyy, g_y_xyyyz, g_y_xyyzz, g_y_xyzzz, g_y_xzzzz, g_y_yyyyy, g_y_yyyyz, g_y_yyyzz, g_y_yyzzz, g_y_yzzzz, g_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_xxxxx[k] = -g_0_xxxxx[k] * cd_y[k] + g_0_xxxxxy[k];

                g_y_xxxxy[k] = -g_0_xxxxy[k] * cd_y[k] + g_0_xxxxyy[k];

                g_y_xxxxz[k] = -g_0_xxxxz[k] * cd_y[k] + g_0_xxxxyz[k];

                g_y_xxxyy[k] = -g_0_xxxyy[k] * cd_y[k] + g_0_xxxyyy[k];

                g_y_xxxyz[k] = -g_0_xxxyz[k] * cd_y[k] + g_0_xxxyyz[k];

                g_y_xxxzz[k] = -g_0_xxxzz[k] * cd_y[k] + g_0_xxxyzz[k];

                g_y_xxyyy[k] = -g_0_xxyyy[k] * cd_y[k] + g_0_xxyyyy[k];

                g_y_xxyyz[k] = -g_0_xxyyz[k] * cd_y[k] + g_0_xxyyyz[k];

                g_y_xxyzz[k] = -g_0_xxyzz[k] * cd_y[k] + g_0_xxyyzz[k];

                g_y_xxzzz[k] = -g_0_xxzzz[k] * cd_y[k] + g_0_xxyzzz[k];

                g_y_xyyyy[k] = -g_0_xyyyy[k] * cd_y[k] + g_0_xyyyyy[k];

                g_y_xyyyz[k] = -g_0_xyyyz[k] * cd_y[k] + g_0_xyyyyz[k];

                g_y_xyyzz[k] = -g_0_xyyzz[k] * cd_y[k] + g_0_xyyyzz[k];

                g_y_xyzzz[k] = -g_0_xyzzz[k] * cd_y[k] + g_0_xyyzzz[k];

                g_y_xzzzz[k] = -g_0_xzzzz[k] * cd_y[k] + g_0_xyzzzz[k];

                g_y_yyyyy[k] = -g_0_yyyyy[k] * cd_y[k] + g_0_yyyyyy[k];

                g_y_yyyyz[k] = -g_0_yyyyz[k] * cd_y[k] + g_0_yyyyyz[k];

                g_y_yyyzz[k] = -g_0_yyyzz[k] * cd_y[k] + g_0_yyyyzz[k];

                g_y_yyzzz[k] = -g_0_yyzzz[k] * cd_y[k] + g_0_yyyzzz[k];

                g_y_yzzzz[k] = -g_0_yzzzz[k] * cd_y[k] + g_0_yyzzzz[k];

                g_y_zzzzz[k] = -g_0_zzzzz[k] * cd_y[k] + g_0_yzzzzz[k];
            }

            /// Set up 42-63 components of targeted buffer : cbuffer.data(

            auto g_z_xxxxx = cbuffer.data(ph_off + 42);

            auto g_z_xxxxy = cbuffer.data(ph_off + 43);

            auto g_z_xxxxz = cbuffer.data(ph_off + 44);

            auto g_z_xxxyy = cbuffer.data(ph_off + 45);

            auto g_z_xxxyz = cbuffer.data(ph_off + 46);

            auto g_z_xxxzz = cbuffer.data(ph_off + 47);

            auto g_z_xxyyy = cbuffer.data(ph_off + 48);

            auto g_z_xxyyz = cbuffer.data(ph_off + 49);

            auto g_z_xxyzz = cbuffer.data(ph_off + 50);

            auto g_z_xxzzz = cbuffer.data(ph_off + 51);

            auto g_z_xyyyy = cbuffer.data(ph_off + 52);

            auto g_z_xyyyz = cbuffer.data(ph_off + 53);

            auto g_z_xyyzz = cbuffer.data(ph_off + 54);

            auto g_z_xyzzz = cbuffer.data(ph_off + 55);

            auto g_z_xzzzz = cbuffer.data(ph_off + 56);

            auto g_z_yyyyy = cbuffer.data(ph_off + 57);

            auto g_z_yyyyz = cbuffer.data(ph_off + 58);

            auto g_z_yyyzz = cbuffer.data(ph_off + 59);

            auto g_z_yyzzz = cbuffer.data(ph_off + 60);

            auto g_z_yzzzz = cbuffer.data(ph_off + 61);

            auto g_z_zzzzz = cbuffer.data(ph_off + 62);

            #pragma omp simd aligned(cd_z, g_0_xxxxx, g_0_xxxxxz, g_0_xxxxy, g_0_xxxxyz, g_0_xxxxz, g_0_xxxxzz, g_0_xxxyy, g_0_xxxyyz, g_0_xxxyz, g_0_xxxyzz, g_0_xxxzz, g_0_xxxzzz, g_0_xxyyy, g_0_xxyyyz, g_0_xxyyz, g_0_xxyyzz, g_0_xxyzz, g_0_xxyzzz, g_0_xxzzz, g_0_xxzzzz, g_0_xyyyy, g_0_xyyyyz, g_0_xyyyz, g_0_xyyyzz, g_0_xyyzz, g_0_xyyzzz, g_0_xyzzz, g_0_xyzzzz, g_0_xzzzz, g_0_xzzzzz, g_0_yyyyy, g_0_yyyyyz, g_0_yyyyz, g_0_yyyyzz, g_0_yyyzz, g_0_yyyzzz, g_0_yyzzz, g_0_yyzzzz, g_0_yzzzz, g_0_yzzzzz, g_0_zzzzz, g_0_zzzzzz, g_z_xxxxx, g_z_xxxxy, g_z_xxxxz, g_z_xxxyy, g_z_xxxyz, g_z_xxxzz, g_z_xxyyy, g_z_xxyyz, g_z_xxyzz, g_z_xxzzz, g_z_xyyyy, g_z_xyyyz, g_z_xyyzz, g_z_xyzzz, g_z_xzzzz, g_z_yyyyy, g_z_yyyyz, g_z_yyyzz, g_z_yyzzz, g_z_yzzzz, g_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_xxxxx[k] = -g_0_xxxxx[k] * cd_z[k] + g_0_xxxxxz[k];

                g_z_xxxxy[k] = -g_0_xxxxy[k] * cd_z[k] + g_0_xxxxyz[k];

                g_z_xxxxz[k] = -g_0_xxxxz[k] * cd_z[k] + g_0_xxxxzz[k];

                g_z_xxxyy[k] = -g_0_xxxyy[k] * cd_z[k] + g_0_xxxyyz[k];

                g_z_xxxyz[k] = -g_0_xxxyz[k] * cd_z[k] + g_0_xxxyzz[k];

                g_z_xxxzz[k] = -g_0_xxxzz[k] * cd_z[k] + g_0_xxxzzz[k];

                g_z_xxyyy[k] = -g_0_xxyyy[k] * cd_z[k] + g_0_xxyyyz[k];

                g_z_xxyyz[k] = -g_0_xxyyz[k] * cd_z[k] + g_0_xxyyzz[k];

                g_z_xxyzz[k] = -g_0_xxyzz[k] * cd_z[k] + g_0_xxyzzz[k];

                g_z_xxzzz[k] = -g_0_xxzzz[k] * cd_z[k] + g_0_xxzzzz[k];

                g_z_xyyyy[k] = -g_0_xyyyy[k] * cd_z[k] + g_0_xyyyyz[k];

                g_z_xyyyz[k] = -g_0_xyyyz[k] * cd_z[k] + g_0_xyyyzz[k];

                g_z_xyyzz[k] = -g_0_xyyzz[k] * cd_z[k] + g_0_xyyzzz[k];

                g_z_xyzzz[k] = -g_0_xyzzz[k] * cd_z[k] + g_0_xyzzzz[k];

                g_z_xzzzz[k] = -g_0_xzzzz[k] * cd_z[k] + g_0_xzzzzz[k];

                g_z_yyyyy[k] = -g_0_yyyyy[k] * cd_z[k] + g_0_yyyyyz[k];

                g_z_yyyyz[k] = -g_0_yyyyz[k] * cd_z[k] + g_0_yyyyzz[k];

                g_z_yyyzz[k] = -g_0_yyyzz[k] * cd_z[k] + g_0_yyyzzz[k];

                g_z_yyzzz[k] = -g_0_yyzzz[k] * cd_z[k] + g_0_yyzzzz[k];

                g_z_yzzzz[k] = -g_0_yzzzz[k] * cd_z[k] + g_0_yzzzzz[k];

                g_z_zzzzz[k] = -g_0_zzzzz[k] * cd_z[k] + g_0_zzzzzz[k];
            }
        }
    }
}

} // erirec namespace

