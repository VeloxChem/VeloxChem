#include "ThreeCenterElectronRepulsionGeom010ContrRecXPH.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_ket_geom010_electron_repulsion_xph(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xph,
                                        const size_t idx_xsh,
                                        const size_t idx_geom_10_xsh,
                                        const size_t idx_geom_10_xsi,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_cd,
                                        const int a_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom,});

    // Set up R(CD) distances

    auto cd_x = factors.data(idx_cd);

    auto cd_y = factors.data(idx_cd + 1);

    auto cd_z = factors.data(idx_cd + 2);

    for (int i = 0; i < acomps; i++)
    {
        /// Set up components of auxilary buffer : SSH

        const auto sh_off = idx_xsh + i * 21;

        auto g_0_xxxxx = cbuffer.data(sh_off + 0);

        auto g_0_xxxxy = cbuffer.data(sh_off + 1);

        auto g_0_xxxxz = cbuffer.data(sh_off + 2);

        auto g_0_xxxyy = cbuffer.data(sh_off + 3);

        auto g_0_xxxyz = cbuffer.data(sh_off + 4);

        auto g_0_xxxzz = cbuffer.data(sh_off + 5);

        auto g_0_xxyyy = cbuffer.data(sh_off + 6);

        auto g_0_xxyyz = cbuffer.data(sh_off + 7);

        auto g_0_xxyzz = cbuffer.data(sh_off + 8);

        auto g_0_xxzzz = cbuffer.data(sh_off + 9);

        auto g_0_xyyyy = cbuffer.data(sh_off + 10);

        auto g_0_xyyyz = cbuffer.data(sh_off + 11);

        auto g_0_xyyzz = cbuffer.data(sh_off + 12);

        auto g_0_xyzzz = cbuffer.data(sh_off + 13);

        auto g_0_xzzzz = cbuffer.data(sh_off + 14);

        auto g_0_yyyyy = cbuffer.data(sh_off + 15);

        auto g_0_yyyyz = cbuffer.data(sh_off + 16);

        auto g_0_yyyzz = cbuffer.data(sh_off + 17);

        auto g_0_yyzzz = cbuffer.data(sh_off + 18);

        auto g_0_yzzzz = cbuffer.data(sh_off + 19);

        auto g_0_zzzzz = cbuffer.data(sh_off + 20);

        /// Set up components of auxilary buffer : SSH

        const auto sh_geom_10_off = idx_geom_10_xsh + i * 21;

        auto g_x_0_0_xxxxx = cbuffer.data(sh_geom_10_off + 0 * acomps + 0);

        auto g_x_0_0_xxxxy = cbuffer.data(sh_geom_10_off + 0 * acomps + 1);

        auto g_x_0_0_xxxxz = cbuffer.data(sh_geom_10_off + 0 * acomps + 2);

        auto g_x_0_0_xxxyy = cbuffer.data(sh_geom_10_off + 0 * acomps + 3);

        auto g_x_0_0_xxxyz = cbuffer.data(sh_geom_10_off + 0 * acomps + 4);

        auto g_x_0_0_xxxzz = cbuffer.data(sh_geom_10_off + 0 * acomps + 5);

        auto g_x_0_0_xxyyy = cbuffer.data(sh_geom_10_off + 0 * acomps + 6);

        auto g_x_0_0_xxyyz = cbuffer.data(sh_geom_10_off + 0 * acomps + 7);

        auto g_x_0_0_xxyzz = cbuffer.data(sh_geom_10_off + 0 * acomps + 8);

        auto g_x_0_0_xxzzz = cbuffer.data(sh_geom_10_off + 0 * acomps + 9);

        auto g_x_0_0_xyyyy = cbuffer.data(sh_geom_10_off + 0 * acomps + 10);

        auto g_x_0_0_xyyyz = cbuffer.data(sh_geom_10_off + 0 * acomps + 11);

        auto g_x_0_0_xyyzz = cbuffer.data(sh_geom_10_off + 0 * acomps + 12);

        auto g_x_0_0_xyzzz = cbuffer.data(sh_geom_10_off + 0 * acomps + 13);

        auto g_x_0_0_xzzzz = cbuffer.data(sh_geom_10_off + 0 * acomps + 14);

        auto g_x_0_0_yyyyy = cbuffer.data(sh_geom_10_off + 0 * acomps + 15);

        auto g_x_0_0_yyyyz = cbuffer.data(sh_geom_10_off + 0 * acomps + 16);

        auto g_x_0_0_yyyzz = cbuffer.data(sh_geom_10_off + 0 * acomps + 17);

        auto g_x_0_0_yyzzz = cbuffer.data(sh_geom_10_off + 0 * acomps + 18);

        auto g_x_0_0_yzzzz = cbuffer.data(sh_geom_10_off + 0 * acomps + 19);

        auto g_x_0_0_zzzzz = cbuffer.data(sh_geom_10_off + 0 * acomps + 20);

        auto g_y_0_0_xxxxx = cbuffer.data(sh_geom_10_off + 21 * acomps + 0);

        auto g_y_0_0_xxxxy = cbuffer.data(sh_geom_10_off + 21 * acomps + 1);

        auto g_y_0_0_xxxxz = cbuffer.data(sh_geom_10_off + 21 * acomps + 2);

        auto g_y_0_0_xxxyy = cbuffer.data(sh_geom_10_off + 21 * acomps + 3);

        auto g_y_0_0_xxxyz = cbuffer.data(sh_geom_10_off + 21 * acomps + 4);

        auto g_y_0_0_xxxzz = cbuffer.data(sh_geom_10_off + 21 * acomps + 5);

        auto g_y_0_0_xxyyy = cbuffer.data(sh_geom_10_off + 21 * acomps + 6);

        auto g_y_0_0_xxyyz = cbuffer.data(sh_geom_10_off + 21 * acomps + 7);

        auto g_y_0_0_xxyzz = cbuffer.data(sh_geom_10_off + 21 * acomps + 8);

        auto g_y_0_0_xxzzz = cbuffer.data(sh_geom_10_off + 21 * acomps + 9);

        auto g_y_0_0_xyyyy = cbuffer.data(sh_geom_10_off + 21 * acomps + 10);

        auto g_y_0_0_xyyyz = cbuffer.data(sh_geom_10_off + 21 * acomps + 11);

        auto g_y_0_0_xyyzz = cbuffer.data(sh_geom_10_off + 21 * acomps + 12);

        auto g_y_0_0_xyzzz = cbuffer.data(sh_geom_10_off + 21 * acomps + 13);

        auto g_y_0_0_xzzzz = cbuffer.data(sh_geom_10_off + 21 * acomps + 14);

        auto g_y_0_0_yyyyy = cbuffer.data(sh_geom_10_off + 21 * acomps + 15);

        auto g_y_0_0_yyyyz = cbuffer.data(sh_geom_10_off + 21 * acomps + 16);

        auto g_y_0_0_yyyzz = cbuffer.data(sh_geom_10_off + 21 * acomps + 17);

        auto g_y_0_0_yyzzz = cbuffer.data(sh_geom_10_off + 21 * acomps + 18);

        auto g_y_0_0_yzzzz = cbuffer.data(sh_geom_10_off + 21 * acomps + 19);

        auto g_y_0_0_zzzzz = cbuffer.data(sh_geom_10_off + 21 * acomps + 20);

        auto g_z_0_0_xxxxx = cbuffer.data(sh_geom_10_off + 42 * acomps + 0);

        auto g_z_0_0_xxxxy = cbuffer.data(sh_geom_10_off + 42 * acomps + 1);

        auto g_z_0_0_xxxxz = cbuffer.data(sh_geom_10_off + 42 * acomps + 2);

        auto g_z_0_0_xxxyy = cbuffer.data(sh_geom_10_off + 42 * acomps + 3);

        auto g_z_0_0_xxxyz = cbuffer.data(sh_geom_10_off + 42 * acomps + 4);

        auto g_z_0_0_xxxzz = cbuffer.data(sh_geom_10_off + 42 * acomps + 5);

        auto g_z_0_0_xxyyy = cbuffer.data(sh_geom_10_off + 42 * acomps + 6);

        auto g_z_0_0_xxyyz = cbuffer.data(sh_geom_10_off + 42 * acomps + 7);

        auto g_z_0_0_xxyzz = cbuffer.data(sh_geom_10_off + 42 * acomps + 8);

        auto g_z_0_0_xxzzz = cbuffer.data(sh_geom_10_off + 42 * acomps + 9);

        auto g_z_0_0_xyyyy = cbuffer.data(sh_geom_10_off + 42 * acomps + 10);

        auto g_z_0_0_xyyyz = cbuffer.data(sh_geom_10_off + 42 * acomps + 11);

        auto g_z_0_0_xyyzz = cbuffer.data(sh_geom_10_off + 42 * acomps + 12);

        auto g_z_0_0_xyzzz = cbuffer.data(sh_geom_10_off + 42 * acomps + 13);

        auto g_z_0_0_xzzzz = cbuffer.data(sh_geom_10_off + 42 * acomps + 14);

        auto g_z_0_0_yyyyy = cbuffer.data(sh_geom_10_off + 42 * acomps + 15);

        auto g_z_0_0_yyyyz = cbuffer.data(sh_geom_10_off + 42 * acomps + 16);

        auto g_z_0_0_yyyzz = cbuffer.data(sh_geom_10_off + 42 * acomps + 17);

        auto g_z_0_0_yyzzz = cbuffer.data(sh_geom_10_off + 42 * acomps + 18);

        auto g_z_0_0_yzzzz = cbuffer.data(sh_geom_10_off + 42 * acomps + 19);

        auto g_z_0_0_zzzzz = cbuffer.data(sh_geom_10_off + 42 * acomps + 20);

        /// Set up components of auxilary buffer : SSI

        const auto si_geom_10_off = idx_geom_10_xsi + i * 28;

        auto g_x_0_0_xxxxxx = cbuffer.data(si_geom_10_off + 0 * acomps + 0);

        auto g_x_0_0_xxxxxy = cbuffer.data(si_geom_10_off + 0 * acomps + 1);

        auto g_x_0_0_xxxxxz = cbuffer.data(si_geom_10_off + 0 * acomps + 2);

        auto g_x_0_0_xxxxyy = cbuffer.data(si_geom_10_off + 0 * acomps + 3);

        auto g_x_0_0_xxxxyz = cbuffer.data(si_geom_10_off + 0 * acomps + 4);

        auto g_x_0_0_xxxxzz = cbuffer.data(si_geom_10_off + 0 * acomps + 5);

        auto g_x_0_0_xxxyyy = cbuffer.data(si_geom_10_off + 0 * acomps + 6);

        auto g_x_0_0_xxxyyz = cbuffer.data(si_geom_10_off + 0 * acomps + 7);

        auto g_x_0_0_xxxyzz = cbuffer.data(si_geom_10_off + 0 * acomps + 8);

        auto g_x_0_0_xxxzzz = cbuffer.data(si_geom_10_off + 0 * acomps + 9);

        auto g_x_0_0_xxyyyy = cbuffer.data(si_geom_10_off + 0 * acomps + 10);

        auto g_x_0_0_xxyyyz = cbuffer.data(si_geom_10_off + 0 * acomps + 11);

        auto g_x_0_0_xxyyzz = cbuffer.data(si_geom_10_off + 0 * acomps + 12);

        auto g_x_0_0_xxyzzz = cbuffer.data(si_geom_10_off + 0 * acomps + 13);

        auto g_x_0_0_xxzzzz = cbuffer.data(si_geom_10_off + 0 * acomps + 14);

        auto g_x_0_0_xyyyyy = cbuffer.data(si_geom_10_off + 0 * acomps + 15);

        auto g_x_0_0_xyyyyz = cbuffer.data(si_geom_10_off + 0 * acomps + 16);

        auto g_x_0_0_xyyyzz = cbuffer.data(si_geom_10_off + 0 * acomps + 17);

        auto g_x_0_0_xyyzzz = cbuffer.data(si_geom_10_off + 0 * acomps + 18);

        auto g_x_0_0_xyzzzz = cbuffer.data(si_geom_10_off + 0 * acomps + 19);

        auto g_x_0_0_xzzzzz = cbuffer.data(si_geom_10_off + 0 * acomps + 20);

        auto g_x_0_0_yyyyyy = cbuffer.data(si_geom_10_off + 0 * acomps + 21);

        auto g_x_0_0_yyyyyz = cbuffer.data(si_geom_10_off + 0 * acomps + 22);

        auto g_x_0_0_yyyyzz = cbuffer.data(si_geom_10_off + 0 * acomps + 23);

        auto g_x_0_0_yyyzzz = cbuffer.data(si_geom_10_off + 0 * acomps + 24);

        auto g_x_0_0_yyzzzz = cbuffer.data(si_geom_10_off + 0 * acomps + 25);

        auto g_x_0_0_yzzzzz = cbuffer.data(si_geom_10_off + 0 * acomps + 26);

        auto g_x_0_0_zzzzzz = cbuffer.data(si_geom_10_off + 0 * acomps + 27);

        auto g_y_0_0_xxxxxx = cbuffer.data(si_geom_10_off + 28 * acomps + 0);

        auto g_y_0_0_xxxxxy = cbuffer.data(si_geom_10_off + 28 * acomps + 1);

        auto g_y_0_0_xxxxxz = cbuffer.data(si_geom_10_off + 28 * acomps + 2);

        auto g_y_0_0_xxxxyy = cbuffer.data(si_geom_10_off + 28 * acomps + 3);

        auto g_y_0_0_xxxxyz = cbuffer.data(si_geom_10_off + 28 * acomps + 4);

        auto g_y_0_0_xxxxzz = cbuffer.data(si_geom_10_off + 28 * acomps + 5);

        auto g_y_0_0_xxxyyy = cbuffer.data(si_geom_10_off + 28 * acomps + 6);

        auto g_y_0_0_xxxyyz = cbuffer.data(si_geom_10_off + 28 * acomps + 7);

        auto g_y_0_0_xxxyzz = cbuffer.data(si_geom_10_off + 28 * acomps + 8);

        auto g_y_0_0_xxxzzz = cbuffer.data(si_geom_10_off + 28 * acomps + 9);

        auto g_y_0_0_xxyyyy = cbuffer.data(si_geom_10_off + 28 * acomps + 10);

        auto g_y_0_0_xxyyyz = cbuffer.data(si_geom_10_off + 28 * acomps + 11);

        auto g_y_0_0_xxyyzz = cbuffer.data(si_geom_10_off + 28 * acomps + 12);

        auto g_y_0_0_xxyzzz = cbuffer.data(si_geom_10_off + 28 * acomps + 13);

        auto g_y_0_0_xxzzzz = cbuffer.data(si_geom_10_off + 28 * acomps + 14);

        auto g_y_0_0_xyyyyy = cbuffer.data(si_geom_10_off + 28 * acomps + 15);

        auto g_y_0_0_xyyyyz = cbuffer.data(si_geom_10_off + 28 * acomps + 16);

        auto g_y_0_0_xyyyzz = cbuffer.data(si_geom_10_off + 28 * acomps + 17);

        auto g_y_0_0_xyyzzz = cbuffer.data(si_geom_10_off + 28 * acomps + 18);

        auto g_y_0_0_xyzzzz = cbuffer.data(si_geom_10_off + 28 * acomps + 19);

        auto g_y_0_0_xzzzzz = cbuffer.data(si_geom_10_off + 28 * acomps + 20);

        auto g_y_0_0_yyyyyy = cbuffer.data(si_geom_10_off + 28 * acomps + 21);

        auto g_y_0_0_yyyyyz = cbuffer.data(si_geom_10_off + 28 * acomps + 22);

        auto g_y_0_0_yyyyzz = cbuffer.data(si_geom_10_off + 28 * acomps + 23);

        auto g_y_0_0_yyyzzz = cbuffer.data(si_geom_10_off + 28 * acomps + 24);

        auto g_y_0_0_yyzzzz = cbuffer.data(si_geom_10_off + 28 * acomps + 25);

        auto g_y_0_0_yzzzzz = cbuffer.data(si_geom_10_off + 28 * acomps + 26);

        auto g_y_0_0_zzzzzz = cbuffer.data(si_geom_10_off + 28 * acomps + 27);

        auto g_z_0_0_xxxxxx = cbuffer.data(si_geom_10_off + 56 * acomps + 0);

        auto g_z_0_0_xxxxxy = cbuffer.data(si_geom_10_off + 56 * acomps + 1);

        auto g_z_0_0_xxxxxz = cbuffer.data(si_geom_10_off + 56 * acomps + 2);

        auto g_z_0_0_xxxxyy = cbuffer.data(si_geom_10_off + 56 * acomps + 3);

        auto g_z_0_0_xxxxyz = cbuffer.data(si_geom_10_off + 56 * acomps + 4);

        auto g_z_0_0_xxxxzz = cbuffer.data(si_geom_10_off + 56 * acomps + 5);

        auto g_z_0_0_xxxyyy = cbuffer.data(si_geom_10_off + 56 * acomps + 6);

        auto g_z_0_0_xxxyyz = cbuffer.data(si_geom_10_off + 56 * acomps + 7);

        auto g_z_0_0_xxxyzz = cbuffer.data(si_geom_10_off + 56 * acomps + 8);

        auto g_z_0_0_xxxzzz = cbuffer.data(si_geom_10_off + 56 * acomps + 9);

        auto g_z_0_0_xxyyyy = cbuffer.data(si_geom_10_off + 56 * acomps + 10);

        auto g_z_0_0_xxyyyz = cbuffer.data(si_geom_10_off + 56 * acomps + 11);

        auto g_z_0_0_xxyyzz = cbuffer.data(si_geom_10_off + 56 * acomps + 12);

        auto g_z_0_0_xxyzzz = cbuffer.data(si_geom_10_off + 56 * acomps + 13);

        auto g_z_0_0_xxzzzz = cbuffer.data(si_geom_10_off + 56 * acomps + 14);

        auto g_z_0_0_xyyyyy = cbuffer.data(si_geom_10_off + 56 * acomps + 15);

        auto g_z_0_0_xyyyyz = cbuffer.data(si_geom_10_off + 56 * acomps + 16);

        auto g_z_0_0_xyyyzz = cbuffer.data(si_geom_10_off + 56 * acomps + 17);

        auto g_z_0_0_xyyzzz = cbuffer.data(si_geom_10_off + 56 * acomps + 18);

        auto g_z_0_0_xyzzzz = cbuffer.data(si_geom_10_off + 56 * acomps + 19);

        auto g_z_0_0_xzzzzz = cbuffer.data(si_geom_10_off + 56 * acomps + 20);

        auto g_z_0_0_yyyyyy = cbuffer.data(si_geom_10_off + 56 * acomps + 21);

        auto g_z_0_0_yyyyyz = cbuffer.data(si_geom_10_off + 56 * acomps + 22);

        auto g_z_0_0_yyyyzz = cbuffer.data(si_geom_10_off + 56 * acomps + 23);

        auto g_z_0_0_yyyzzz = cbuffer.data(si_geom_10_off + 56 * acomps + 24);

        auto g_z_0_0_yyzzzz = cbuffer.data(si_geom_10_off + 56 * acomps + 25);

        auto g_z_0_0_yzzzzz = cbuffer.data(si_geom_10_off + 56 * acomps + 26);

        auto g_z_0_0_zzzzzz = cbuffer.data(si_geom_10_off + 56 * acomps + 27);

        /// set up bra offset for contr_buffer_xxph

        const auto ph_geom_10_off = idx_geom_10_xph + i * 63;

        /// Set up 0-21 components of targeted buffer : cbuffer.data(

        auto g_x_0_x_xxxxx = cbuffer.data(ph_geom_10_off + 0 * acomps  + 0);

        auto g_x_0_x_xxxxy = cbuffer.data(ph_geom_10_off + 0 * acomps  + 1);

        auto g_x_0_x_xxxxz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 2);

        auto g_x_0_x_xxxyy = cbuffer.data(ph_geom_10_off + 0 * acomps  + 3);

        auto g_x_0_x_xxxyz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 4);

        auto g_x_0_x_xxxzz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 5);

        auto g_x_0_x_xxyyy = cbuffer.data(ph_geom_10_off + 0 * acomps  + 6);

        auto g_x_0_x_xxyyz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 7);

        auto g_x_0_x_xxyzz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 8);

        auto g_x_0_x_xxzzz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 9);

        auto g_x_0_x_xyyyy = cbuffer.data(ph_geom_10_off + 0 * acomps  + 10);

        auto g_x_0_x_xyyyz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 11);

        auto g_x_0_x_xyyzz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 12);

        auto g_x_0_x_xyzzz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 13);

        auto g_x_0_x_xzzzz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 14);

        auto g_x_0_x_yyyyy = cbuffer.data(ph_geom_10_off + 0 * acomps  + 15);

        auto g_x_0_x_yyyyz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 16);

        auto g_x_0_x_yyyzz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 17);

        auto g_x_0_x_yyzzz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 18);

        auto g_x_0_x_yzzzz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 19);

        auto g_x_0_x_zzzzz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 20);

        #pragma omp simd aligned(cd_x, g_0_xxxxx, g_0_xxxxy, g_0_xxxxz, g_0_xxxyy, g_0_xxxyz, g_0_xxxzz, g_0_xxyyy, g_0_xxyyz, g_0_xxyzz, g_0_xxzzz, g_0_xyyyy, g_0_xyyyz, g_0_xyyzz, g_0_xyzzz, g_0_xzzzz, g_0_yyyyy, g_0_yyyyz, g_0_yyyzz, g_0_yyzzz, g_0_yzzzz, g_0_zzzzz, g_x_0_0_xxxxx, g_x_0_0_xxxxxx, g_x_0_0_xxxxxy, g_x_0_0_xxxxxz, g_x_0_0_xxxxy, g_x_0_0_xxxxyy, g_x_0_0_xxxxyz, g_x_0_0_xxxxz, g_x_0_0_xxxxzz, g_x_0_0_xxxyy, g_x_0_0_xxxyyy, g_x_0_0_xxxyyz, g_x_0_0_xxxyz, g_x_0_0_xxxyzz, g_x_0_0_xxxzz, g_x_0_0_xxxzzz, g_x_0_0_xxyyy, g_x_0_0_xxyyyy, g_x_0_0_xxyyyz, g_x_0_0_xxyyz, g_x_0_0_xxyyzz, g_x_0_0_xxyzz, g_x_0_0_xxyzzz, g_x_0_0_xxzzz, g_x_0_0_xxzzzz, g_x_0_0_xyyyy, g_x_0_0_xyyyyy, g_x_0_0_xyyyyz, g_x_0_0_xyyyz, g_x_0_0_xyyyzz, g_x_0_0_xyyzz, g_x_0_0_xyyzzz, g_x_0_0_xyzzz, g_x_0_0_xyzzzz, g_x_0_0_xzzzz, g_x_0_0_xzzzzz, g_x_0_0_yyyyy, g_x_0_0_yyyyz, g_x_0_0_yyyzz, g_x_0_0_yyzzz, g_x_0_0_yzzzz, g_x_0_0_zzzzz, g_x_0_x_xxxxx, g_x_0_x_xxxxy, g_x_0_x_xxxxz, g_x_0_x_xxxyy, g_x_0_x_xxxyz, g_x_0_x_xxxzz, g_x_0_x_xxyyy, g_x_0_x_xxyyz, g_x_0_x_xxyzz, g_x_0_x_xxzzz, g_x_0_x_xyyyy, g_x_0_x_xyyyz, g_x_0_x_xyyzz, g_x_0_x_xyzzz, g_x_0_x_xzzzz, g_x_0_x_yyyyy, g_x_0_x_yyyyz, g_x_0_x_yyyzz, g_x_0_x_yyzzz, g_x_0_x_yzzzz, g_x_0_x_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_x_xxxxx[k] = -g_0_xxxxx[k] - g_x_0_0_xxxxx[k] * cd_x[k] + g_x_0_0_xxxxxx[k];

            g_x_0_x_xxxxy[k] = -g_0_xxxxy[k] - g_x_0_0_xxxxy[k] * cd_x[k] + g_x_0_0_xxxxxy[k];

            g_x_0_x_xxxxz[k] = -g_0_xxxxz[k] - g_x_0_0_xxxxz[k] * cd_x[k] + g_x_0_0_xxxxxz[k];

            g_x_0_x_xxxyy[k] = -g_0_xxxyy[k] - g_x_0_0_xxxyy[k] * cd_x[k] + g_x_0_0_xxxxyy[k];

            g_x_0_x_xxxyz[k] = -g_0_xxxyz[k] - g_x_0_0_xxxyz[k] * cd_x[k] + g_x_0_0_xxxxyz[k];

            g_x_0_x_xxxzz[k] = -g_0_xxxzz[k] - g_x_0_0_xxxzz[k] * cd_x[k] + g_x_0_0_xxxxzz[k];

            g_x_0_x_xxyyy[k] = -g_0_xxyyy[k] - g_x_0_0_xxyyy[k] * cd_x[k] + g_x_0_0_xxxyyy[k];

            g_x_0_x_xxyyz[k] = -g_0_xxyyz[k] - g_x_0_0_xxyyz[k] * cd_x[k] + g_x_0_0_xxxyyz[k];

            g_x_0_x_xxyzz[k] = -g_0_xxyzz[k] - g_x_0_0_xxyzz[k] * cd_x[k] + g_x_0_0_xxxyzz[k];

            g_x_0_x_xxzzz[k] = -g_0_xxzzz[k] - g_x_0_0_xxzzz[k] * cd_x[k] + g_x_0_0_xxxzzz[k];

            g_x_0_x_xyyyy[k] = -g_0_xyyyy[k] - g_x_0_0_xyyyy[k] * cd_x[k] + g_x_0_0_xxyyyy[k];

            g_x_0_x_xyyyz[k] = -g_0_xyyyz[k] - g_x_0_0_xyyyz[k] * cd_x[k] + g_x_0_0_xxyyyz[k];

            g_x_0_x_xyyzz[k] = -g_0_xyyzz[k] - g_x_0_0_xyyzz[k] * cd_x[k] + g_x_0_0_xxyyzz[k];

            g_x_0_x_xyzzz[k] = -g_0_xyzzz[k] - g_x_0_0_xyzzz[k] * cd_x[k] + g_x_0_0_xxyzzz[k];

            g_x_0_x_xzzzz[k] = -g_0_xzzzz[k] - g_x_0_0_xzzzz[k] * cd_x[k] + g_x_0_0_xxzzzz[k];

            g_x_0_x_yyyyy[k] = -g_0_yyyyy[k] - g_x_0_0_yyyyy[k] * cd_x[k] + g_x_0_0_xyyyyy[k];

            g_x_0_x_yyyyz[k] = -g_0_yyyyz[k] - g_x_0_0_yyyyz[k] * cd_x[k] + g_x_0_0_xyyyyz[k];

            g_x_0_x_yyyzz[k] = -g_0_yyyzz[k] - g_x_0_0_yyyzz[k] * cd_x[k] + g_x_0_0_xyyyzz[k];

            g_x_0_x_yyzzz[k] = -g_0_yyzzz[k] - g_x_0_0_yyzzz[k] * cd_x[k] + g_x_0_0_xyyzzz[k];

            g_x_0_x_yzzzz[k] = -g_0_yzzzz[k] - g_x_0_0_yzzzz[k] * cd_x[k] + g_x_0_0_xyzzzz[k];

            g_x_0_x_zzzzz[k] = -g_0_zzzzz[k] - g_x_0_0_zzzzz[k] * cd_x[k] + g_x_0_0_xzzzzz[k];
        }

        /// Set up 21-42 components of targeted buffer : cbuffer.data(

        auto g_x_0_y_xxxxx = cbuffer.data(ph_geom_10_off + 0 * acomps  + 21);

        auto g_x_0_y_xxxxy = cbuffer.data(ph_geom_10_off + 0 * acomps  + 22);

        auto g_x_0_y_xxxxz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 23);

        auto g_x_0_y_xxxyy = cbuffer.data(ph_geom_10_off + 0 * acomps  + 24);

        auto g_x_0_y_xxxyz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 25);

        auto g_x_0_y_xxxzz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 26);

        auto g_x_0_y_xxyyy = cbuffer.data(ph_geom_10_off + 0 * acomps  + 27);

        auto g_x_0_y_xxyyz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 28);

        auto g_x_0_y_xxyzz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 29);

        auto g_x_0_y_xxzzz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 30);

        auto g_x_0_y_xyyyy = cbuffer.data(ph_geom_10_off + 0 * acomps  + 31);

        auto g_x_0_y_xyyyz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 32);

        auto g_x_0_y_xyyzz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 33);

        auto g_x_0_y_xyzzz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 34);

        auto g_x_0_y_xzzzz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 35);

        auto g_x_0_y_yyyyy = cbuffer.data(ph_geom_10_off + 0 * acomps  + 36);

        auto g_x_0_y_yyyyz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 37);

        auto g_x_0_y_yyyzz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 38);

        auto g_x_0_y_yyzzz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 39);

        auto g_x_0_y_yzzzz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 40);

        auto g_x_0_y_zzzzz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 41);

        #pragma omp simd aligned(cd_y, g_x_0_0_xxxxx, g_x_0_0_xxxxxy, g_x_0_0_xxxxy, g_x_0_0_xxxxyy, g_x_0_0_xxxxyz, g_x_0_0_xxxxz, g_x_0_0_xxxyy, g_x_0_0_xxxyyy, g_x_0_0_xxxyyz, g_x_0_0_xxxyz, g_x_0_0_xxxyzz, g_x_0_0_xxxzz, g_x_0_0_xxyyy, g_x_0_0_xxyyyy, g_x_0_0_xxyyyz, g_x_0_0_xxyyz, g_x_0_0_xxyyzz, g_x_0_0_xxyzz, g_x_0_0_xxyzzz, g_x_0_0_xxzzz, g_x_0_0_xyyyy, g_x_0_0_xyyyyy, g_x_0_0_xyyyyz, g_x_0_0_xyyyz, g_x_0_0_xyyyzz, g_x_0_0_xyyzz, g_x_0_0_xyyzzz, g_x_0_0_xyzzz, g_x_0_0_xyzzzz, g_x_0_0_xzzzz, g_x_0_0_yyyyy, g_x_0_0_yyyyyy, g_x_0_0_yyyyyz, g_x_0_0_yyyyz, g_x_0_0_yyyyzz, g_x_0_0_yyyzz, g_x_0_0_yyyzzz, g_x_0_0_yyzzz, g_x_0_0_yyzzzz, g_x_0_0_yzzzz, g_x_0_0_yzzzzz, g_x_0_0_zzzzz, g_x_0_y_xxxxx, g_x_0_y_xxxxy, g_x_0_y_xxxxz, g_x_0_y_xxxyy, g_x_0_y_xxxyz, g_x_0_y_xxxzz, g_x_0_y_xxyyy, g_x_0_y_xxyyz, g_x_0_y_xxyzz, g_x_0_y_xxzzz, g_x_0_y_xyyyy, g_x_0_y_xyyyz, g_x_0_y_xyyzz, g_x_0_y_xyzzz, g_x_0_y_xzzzz, g_x_0_y_yyyyy, g_x_0_y_yyyyz, g_x_0_y_yyyzz, g_x_0_y_yyzzz, g_x_0_y_yzzzz, g_x_0_y_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_y_xxxxx[k] = -g_x_0_0_xxxxx[k] * cd_y[k] + g_x_0_0_xxxxxy[k];

            g_x_0_y_xxxxy[k] = -g_x_0_0_xxxxy[k] * cd_y[k] + g_x_0_0_xxxxyy[k];

            g_x_0_y_xxxxz[k] = -g_x_0_0_xxxxz[k] * cd_y[k] + g_x_0_0_xxxxyz[k];

            g_x_0_y_xxxyy[k] = -g_x_0_0_xxxyy[k] * cd_y[k] + g_x_0_0_xxxyyy[k];

            g_x_0_y_xxxyz[k] = -g_x_0_0_xxxyz[k] * cd_y[k] + g_x_0_0_xxxyyz[k];

            g_x_0_y_xxxzz[k] = -g_x_0_0_xxxzz[k] * cd_y[k] + g_x_0_0_xxxyzz[k];

            g_x_0_y_xxyyy[k] = -g_x_0_0_xxyyy[k] * cd_y[k] + g_x_0_0_xxyyyy[k];

            g_x_0_y_xxyyz[k] = -g_x_0_0_xxyyz[k] * cd_y[k] + g_x_0_0_xxyyyz[k];

            g_x_0_y_xxyzz[k] = -g_x_0_0_xxyzz[k] * cd_y[k] + g_x_0_0_xxyyzz[k];

            g_x_0_y_xxzzz[k] = -g_x_0_0_xxzzz[k] * cd_y[k] + g_x_0_0_xxyzzz[k];

            g_x_0_y_xyyyy[k] = -g_x_0_0_xyyyy[k] * cd_y[k] + g_x_0_0_xyyyyy[k];

            g_x_0_y_xyyyz[k] = -g_x_0_0_xyyyz[k] * cd_y[k] + g_x_0_0_xyyyyz[k];

            g_x_0_y_xyyzz[k] = -g_x_0_0_xyyzz[k] * cd_y[k] + g_x_0_0_xyyyzz[k];

            g_x_0_y_xyzzz[k] = -g_x_0_0_xyzzz[k] * cd_y[k] + g_x_0_0_xyyzzz[k];

            g_x_0_y_xzzzz[k] = -g_x_0_0_xzzzz[k] * cd_y[k] + g_x_0_0_xyzzzz[k];

            g_x_0_y_yyyyy[k] = -g_x_0_0_yyyyy[k] * cd_y[k] + g_x_0_0_yyyyyy[k];

            g_x_0_y_yyyyz[k] = -g_x_0_0_yyyyz[k] * cd_y[k] + g_x_0_0_yyyyyz[k];

            g_x_0_y_yyyzz[k] = -g_x_0_0_yyyzz[k] * cd_y[k] + g_x_0_0_yyyyzz[k];

            g_x_0_y_yyzzz[k] = -g_x_0_0_yyzzz[k] * cd_y[k] + g_x_0_0_yyyzzz[k];

            g_x_0_y_yzzzz[k] = -g_x_0_0_yzzzz[k] * cd_y[k] + g_x_0_0_yyzzzz[k];

            g_x_0_y_zzzzz[k] = -g_x_0_0_zzzzz[k] * cd_y[k] + g_x_0_0_yzzzzz[k];
        }

        /// Set up 42-63 components of targeted buffer : cbuffer.data(

        auto g_x_0_z_xxxxx = cbuffer.data(ph_geom_10_off + 0 * acomps  + 42);

        auto g_x_0_z_xxxxy = cbuffer.data(ph_geom_10_off + 0 * acomps  + 43);

        auto g_x_0_z_xxxxz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 44);

        auto g_x_0_z_xxxyy = cbuffer.data(ph_geom_10_off + 0 * acomps  + 45);

        auto g_x_0_z_xxxyz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 46);

        auto g_x_0_z_xxxzz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 47);

        auto g_x_0_z_xxyyy = cbuffer.data(ph_geom_10_off + 0 * acomps  + 48);

        auto g_x_0_z_xxyyz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 49);

        auto g_x_0_z_xxyzz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 50);

        auto g_x_0_z_xxzzz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 51);

        auto g_x_0_z_xyyyy = cbuffer.data(ph_geom_10_off + 0 * acomps  + 52);

        auto g_x_0_z_xyyyz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 53);

        auto g_x_0_z_xyyzz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 54);

        auto g_x_0_z_xyzzz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 55);

        auto g_x_0_z_xzzzz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 56);

        auto g_x_0_z_yyyyy = cbuffer.data(ph_geom_10_off + 0 * acomps  + 57);

        auto g_x_0_z_yyyyz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 58);

        auto g_x_0_z_yyyzz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 59);

        auto g_x_0_z_yyzzz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 60);

        auto g_x_0_z_yzzzz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 61);

        auto g_x_0_z_zzzzz = cbuffer.data(ph_geom_10_off + 0 * acomps  + 62);

        #pragma omp simd aligned(cd_z, g_x_0_0_xxxxx, g_x_0_0_xxxxxz, g_x_0_0_xxxxy, g_x_0_0_xxxxyz, g_x_0_0_xxxxz, g_x_0_0_xxxxzz, g_x_0_0_xxxyy, g_x_0_0_xxxyyz, g_x_0_0_xxxyz, g_x_0_0_xxxyzz, g_x_0_0_xxxzz, g_x_0_0_xxxzzz, g_x_0_0_xxyyy, g_x_0_0_xxyyyz, g_x_0_0_xxyyz, g_x_0_0_xxyyzz, g_x_0_0_xxyzz, g_x_0_0_xxyzzz, g_x_0_0_xxzzz, g_x_0_0_xxzzzz, g_x_0_0_xyyyy, g_x_0_0_xyyyyz, g_x_0_0_xyyyz, g_x_0_0_xyyyzz, g_x_0_0_xyyzz, g_x_0_0_xyyzzz, g_x_0_0_xyzzz, g_x_0_0_xyzzzz, g_x_0_0_xzzzz, g_x_0_0_xzzzzz, g_x_0_0_yyyyy, g_x_0_0_yyyyyz, g_x_0_0_yyyyz, g_x_0_0_yyyyzz, g_x_0_0_yyyzz, g_x_0_0_yyyzzz, g_x_0_0_yyzzz, g_x_0_0_yyzzzz, g_x_0_0_yzzzz, g_x_0_0_yzzzzz, g_x_0_0_zzzzz, g_x_0_0_zzzzzz, g_x_0_z_xxxxx, g_x_0_z_xxxxy, g_x_0_z_xxxxz, g_x_0_z_xxxyy, g_x_0_z_xxxyz, g_x_0_z_xxxzz, g_x_0_z_xxyyy, g_x_0_z_xxyyz, g_x_0_z_xxyzz, g_x_0_z_xxzzz, g_x_0_z_xyyyy, g_x_0_z_xyyyz, g_x_0_z_xyyzz, g_x_0_z_xyzzz, g_x_0_z_xzzzz, g_x_0_z_yyyyy, g_x_0_z_yyyyz, g_x_0_z_yyyzz, g_x_0_z_yyzzz, g_x_0_z_yzzzz, g_x_0_z_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_z_xxxxx[k] = -g_x_0_0_xxxxx[k] * cd_z[k] + g_x_0_0_xxxxxz[k];

            g_x_0_z_xxxxy[k] = -g_x_0_0_xxxxy[k] * cd_z[k] + g_x_0_0_xxxxyz[k];

            g_x_0_z_xxxxz[k] = -g_x_0_0_xxxxz[k] * cd_z[k] + g_x_0_0_xxxxzz[k];

            g_x_0_z_xxxyy[k] = -g_x_0_0_xxxyy[k] * cd_z[k] + g_x_0_0_xxxyyz[k];

            g_x_0_z_xxxyz[k] = -g_x_0_0_xxxyz[k] * cd_z[k] + g_x_0_0_xxxyzz[k];

            g_x_0_z_xxxzz[k] = -g_x_0_0_xxxzz[k] * cd_z[k] + g_x_0_0_xxxzzz[k];

            g_x_0_z_xxyyy[k] = -g_x_0_0_xxyyy[k] * cd_z[k] + g_x_0_0_xxyyyz[k];

            g_x_0_z_xxyyz[k] = -g_x_0_0_xxyyz[k] * cd_z[k] + g_x_0_0_xxyyzz[k];

            g_x_0_z_xxyzz[k] = -g_x_0_0_xxyzz[k] * cd_z[k] + g_x_0_0_xxyzzz[k];

            g_x_0_z_xxzzz[k] = -g_x_0_0_xxzzz[k] * cd_z[k] + g_x_0_0_xxzzzz[k];

            g_x_0_z_xyyyy[k] = -g_x_0_0_xyyyy[k] * cd_z[k] + g_x_0_0_xyyyyz[k];

            g_x_0_z_xyyyz[k] = -g_x_0_0_xyyyz[k] * cd_z[k] + g_x_0_0_xyyyzz[k];

            g_x_0_z_xyyzz[k] = -g_x_0_0_xyyzz[k] * cd_z[k] + g_x_0_0_xyyzzz[k];

            g_x_0_z_xyzzz[k] = -g_x_0_0_xyzzz[k] * cd_z[k] + g_x_0_0_xyzzzz[k];

            g_x_0_z_xzzzz[k] = -g_x_0_0_xzzzz[k] * cd_z[k] + g_x_0_0_xzzzzz[k];

            g_x_0_z_yyyyy[k] = -g_x_0_0_yyyyy[k] * cd_z[k] + g_x_0_0_yyyyyz[k];

            g_x_0_z_yyyyz[k] = -g_x_0_0_yyyyz[k] * cd_z[k] + g_x_0_0_yyyyzz[k];

            g_x_0_z_yyyzz[k] = -g_x_0_0_yyyzz[k] * cd_z[k] + g_x_0_0_yyyzzz[k];

            g_x_0_z_yyzzz[k] = -g_x_0_0_yyzzz[k] * cd_z[k] + g_x_0_0_yyzzzz[k];

            g_x_0_z_yzzzz[k] = -g_x_0_0_yzzzz[k] * cd_z[k] + g_x_0_0_yzzzzz[k];

            g_x_0_z_zzzzz[k] = -g_x_0_0_zzzzz[k] * cd_z[k] + g_x_0_0_zzzzzz[k];
        }
        /// Set up 0-21 components of targeted buffer : cbuffer.data(

        auto g_y_0_x_xxxxx = cbuffer.data(ph_geom_10_off + 63 * acomps  + 0);

        auto g_y_0_x_xxxxy = cbuffer.data(ph_geom_10_off + 63 * acomps  + 1);

        auto g_y_0_x_xxxxz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 2);

        auto g_y_0_x_xxxyy = cbuffer.data(ph_geom_10_off + 63 * acomps  + 3);

        auto g_y_0_x_xxxyz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 4);

        auto g_y_0_x_xxxzz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 5);

        auto g_y_0_x_xxyyy = cbuffer.data(ph_geom_10_off + 63 * acomps  + 6);

        auto g_y_0_x_xxyyz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 7);

        auto g_y_0_x_xxyzz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 8);

        auto g_y_0_x_xxzzz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 9);

        auto g_y_0_x_xyyyy = cbuffer.data(ph_geom_10_off + 63 * acomps  + 10);

        auto g_y_0_x_xyyyz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 11);

        auto g_y_0_x_xyyzz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 12);

        auto g_y_0_x_xyzzz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 13);

        auto g_y_0_x_xzzzz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 14);

        auto g_y_0_x_yyyyy = cbuffer.data(ph_geom_10_off + 63 * acomps  + 15);

        auto g_y_0_x_yyyyz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 16);

        auto g_y_0_x_yyyzz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 17);

        auto g_y_0_x_yyzzz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 18);

        auto g_y_0_x_yzzzz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 19);

        auto g_y_0_x_zzzzz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 20);

        #pragma omp simd aligned(cd_x, g_y_0_0_xxxxx, g_y_0_0_xxxxxx, g_y_0_0_xxxxxy, g_y_0_0_xxxxxz, g_y_0_0_xxxxy, g_y_0_0_xxxxyy, g_y_0_0_xxxxyz, g_y_0_0_xxxxz, g_y_0_0_xxxxzz, g_y_0_0_xxxyy, g_y_0_0_xxxyyy, g_y_0_0_xxxyyz, g_y_0_0_xxxyz, g_y_0_0_xxxyzz, g_y_0_0_xxxzz, g_y_0_0_xxxzzz, g_y_0_0_xxyyy, g_y_0_0_xxyyyy, g_y_0_0_xxyyyz, g_y_0_0_xxyyz, g_y_0_0_xxyyzz, g_y_0_0_xxyzz, g_y_0_0_xxyzzz, g_y_0_0_xxzzz, g_y_0_0_xxzzzz, g_y_0_0_xyyyy, g_y_0_0_xyyyyy, g_y_0_0_xyyyyz, g_y_0_0_xyyyz, g_y_0_0_xyyyzz, g_y_0_0_xyyzz, g_y_0_0_xyyzzz, g_y_0_0_xyzzz, g_y_0_0_xyzzzz, g_y_0_0_xzzzz, g_y_0_0_xzzzzz, g_y_0_0_yyyyy, g_y_0_0_yyyyz, g_y_0_0_yyyzz, g_y_0_0_yyzzz, g_y_0_0_yzzzz, g_y_0_0_zzzzz, g_y_0_x_xxxxx, g_y_0_x_xxxxy, g_y_0_x_xxxxz, g_y_0_x_xxxyy, g_y_0_x_xxxyz, g_y_0_x_xxxzz, g_y_0_x_xxyyy, g_y_0_x_xxyyz, g_y_0_x_xxyzz, g_y_0_x_xxzzz, g_y_0_x_xyyyy, g_y_0_x_xyyyz, g_y_0_x_xyyzz, g_y_0_x_xyzzz, g_y_0_x_xzzzz, g_y_0_x_yyyyy, g_y_0_x_yyyyz, g_y_0_x_yyyzz, g_y_0_x_yyzzz, g_y_0_x_yzzzz, g_y_0_x_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_x_xxxxx[k] = -g_y_0_0_xxxxx[k] * cd_x[k] + g_y_0_0_xxxxxx[k];

            g_y_0_x_xxxxy[k] = -g_y_0_0_xxxxy[k] * cd_x[k] + g_y_0_0_xxxxxy[k];

            g_y_0_x_xxxxz[k] = -g_y_0_0_xxxxz[k] * cd_x[k] + g_y_0_0_xxxxxz[k];

            g_y_0_x_xxxyy[k] = -g_y_0_0_xxxyy[k] * cd_x[k] + g_y_0_0_xxxxyy[k];

            g_y_0_x_xxxyz[k] = -g_y_0_0_xxxyz[k] * cd_x[k] + g_y_0_0_xxxxyz[k];

            g_y_0_x_xxxzz[k] = -g_y_0_0_xxxzz[k] * cd_x[k] + g_y_0_0_xxxxzz[k];

            g_y_0_x_xxyyy[k] = -g_y_0_0_xxyyy[k] * cd_x[k] + g_y_0_0_xxxyyy[k];

            g_y_0_x_xxyyz[k] = -g_y_0_0_xxyyz[k] * cd_x[k] + g_y_0_0_xxxyyz[k];

            g_y_0_x_xxyzz[k] = -g_y_0_0_xxyzz[k] * cd_x[k] + g_y_0_0_xxxyzz[k];

            g_y_0_x_xxzzz[k] = -g_y_0_0_xxzzz[k] * cd_x[k] + g_y_0_0_xxxzzz[k];

            g_y_0_x_xyyyy[k] = -g_y_0_0_xyyyy[k] * cd_x[k] + g_y_0_0_xxyyyy[k];

            g_y_0_x_xyyyz[k] = -g_y_0_0_xyyyz[k] * cd_x[k] + g_y_0_0_xxyyyz[k];

            g_y_0_x_xyyzz[k] = -g_y_0_0_xyyzz[k] * cd_x[k] + g_y_0_0_xxyyzz[k];

            g_y_0_x_xyzzz[k] = -g_y_0_0_xyzzz[k] * cd_x[k] + g_y_0_0_xxyzzz[k];

            g_y_0_x_xzzzz[k] = -g_y_0_0_xzzzz[k] * cd_x[k] + g_y_0_0_xxzzzz[k];

            g_y_0_x_yyyyy[k] = -g_y_0_0_yyyyy[k] * cd_x[k] + g_y_0_0_xyyyyy[k];

            g_y_0_x_yyyyz[k] = -g_y_0_0_yyyyz[k] * cd_x[k] + g_y_0_0_xyyyyz[k];

            g_y_0_x_yyyzz[k] = -g_y_0_0_yyyzz[k] * cd_x[k] + g_y_0_0_xyyyzz[k];

            g_y_0_x_yyzzz[k] = -g_y_0_0_yyzzz[k] * cd_x[k] + g_y_0_0_xyyzzz[k];

            g_y_0_x_yzzzz[k] = -g_y_0_0_yzzzz[k] * cd_x[k] + g_y_0_0_xyzzzz[k];

            g_y_0_x_zzzzz[k] = -g_y_0_0_zzzzz[k] * cd_x[k] + g_y_0_0_xzzzzz[k];
        }

        /// Set up 21-42 components of targeted buffer : cbuffer.data(

        auto g_y_0_y_xxxxx = cbuffer.data(ph_geom_10_off + 63 * acomps  + 21);

        auto g_y_0_y_xxxxy = cbuffer.data(ph_geom_10_off + 63 * acomps  + 22);

        auto g_y_0_y_xxxxz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 23);

        auto g_y_0_y_xxxyy = cbuffer.data(ph_geom_10_off + 63 * acomps  + 24);

        auto g_y_0_y_xxxyz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 25);

        auto g_y_0_y_xxxzz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 26);

        auto g_y_0_y_xxyyy = cbuffer.data(ph_geom_10_off + 63 * acomps  + 27);

        auto g_y_0_y_xxyyz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 28);

        auto g_y_0_y_xxyzz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 29);

        auto g_y_0_y_xxzzz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 30);

        auto g_y_0_y_xyyyy = cbuffer.data(ph_geom_10_off + 63 * acomps  + 31);

        auto g_y_0_y_xyyyz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 32);

        auto g_y_0_y_xyyzz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 33);

        auto g_y_0_y_xyzzz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 34);

        auto g_y_0_y_xzzzz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 35);

        auto g_y_0_y_yyyyy = cbuffer.data(ph_geom_10_off + 63 * acomps  + 36);

        auto g_y_0_y_yyyyz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 37);

        auto g_y_0_y_yyyzz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 38);

        auto g_y_0_y_yyzzz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 39);

        auto g_y_0_y_yzzzz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 40);

        auto g_y_0_y_zzzzz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 41);

        #pragma omp simd aligned(cd_y, g_0_xxxxx, g_0_xxxxy, g_0_xxxxz, g_0_xxxyy, g_0_xxxyz, g_0_xxxzz, g_0_xxyyy, g_0_xxyyz, g_0_xxyzz, g_0_xxzzz, g_0_xyyyy, g_0_xyyyz, g_0_xyyzz, g_0_xyzzz, g_0_xzzzz, g_0_yyyyy, g_0_yyyyz, g_0_yyyzz, g_0_yyzzz, g_0_yzzzz, g_0_zzzzz, g_y_0_0_xxxxx, g_y_0_0_xxxxxy, g_y_0_0_xxxxy, g_y_0_0_xxxxyy, g_y_0_0_xxxxyz, g_y_0_0_xxxxz, g_y_0_0_xxxyy, g_y_0_0_xxxyyy, g_y_0_0_xxxyyz, g_y_0_0_xxxyz, g_y_0_0_xxxyzz, g_y_0_0_xxxzz, g_y_0_0_xxyyy, g_y_0_0_xxyyyy, g_y_0_0_xxyyyz, g_y_0_0_xxyyz, g_y_0_0_xxyyzz, g_y_0_0_xxyzz, g_y_0_0_xxyzzz, g_y_0_0_xxzzz, g_y_0_0_xyyyy, g_y_0_0_xyyyyy, g_y_0_0_xyyyyz, g_y_0_0_xyyyz, g_y_0_0_xyyyzz, g_y_0_0_xyyzz, g_y_0_0_xyyzzz, g_y_0_0_xyzzz, g_y_0_0_xyzzzz, g_y_0_0_xzzzz, g_y_0_0_yyyyy, g_y_0_0_yyyyyy, g_y_0_0_yyyyyz, g_y_0_0_yyyyz, g_y_0_0_yyyyzz, g_y_0_0_yyyzz, g_y_0_0_yyyzzz, g_y_0_0_yyzzz, g_y_0_0_yyzzzz, g_y_0_0_yzzzz, g_y_0_0_yzzzzz, g_y_0_0_zzzzz, g_y_0_y_xxxxx, g_y_0_y_xxxxy, g_y_0_y_xxxxz, g_y_0_y_xxxyy, g_y_0_y_xxxyz, g_y_0_y_xxxzz, g_y_0_y_xxyyy, g_y_0_y_xxyyz, g_y_0_y_xxyzz, g_y_0_y_xxzzz, g_y_0_y_xyyyy, g_y_0_y_xyyyz, g_y_0_y_xyyzz, g_y_0_y_xyzzz, g_y_0_y_xzzzz, g_y_0_y_yyyyy, g_y_0_y_yyyyz, g_y_0_y_yyyzz, g_y_0_y_yyzzz, g_y_0_y_yzzzz, g_y_0_y_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_y_xxxxx[k] = -g_0_xxxxx[k] - g_y_0_0_xxxxx[k] * cd_y[k] + g_y_0_0_xxxxxy[k];

            g_y_0_y_xxxxy[k] = -g_0_xxxxy[k] - g_y_0_0_xxxxy[k] * cd_y[k] + g_y_0_0_xxxxyy[k];

            g_y_0_y_xxxxz[k] = -g_0_xxxxz[k] - g_y_0_0_xxxxz[k] * cd_y[k] + g_y_0_0_xxxxyz[k];

            g_y_0_y_xxxyy[k] = -g_0_xxxyy[k] - g_y_0_0_xxxyy[k] * cd_y[k] + g_y_0_0_xxxyyy[k];

            g_y_0_y_xxxyz[k] = -g_0_xxxyz[k] - g_y_0_0_xxxyz[k] * cd_y[k] + g_y_0_0_xxxyyz[k];

            g_y_0_y_xxxzz[k] = -g_0_xxxzz[k] - g_y_0_0_xxxzz[k] * cd_y[k] + g_y_0_0_xxxyzz[k];

            g_y_0_y_xxyyy[k] = -g_0_xxyyy[k] - g_y_0_0_xxyyy[k] * cd_y[k] + g_y_0_0_xxyyyy[k];

            g_y_0_y_xxyyz[k] = -g_0_xxyyz[k] - g_y_0_0_xxyyz[k] * cd_y[k] + g_y_0_0_xxyyyz[k];

            g_y_0_y_xxyzz[k] = -g_0_xxyzz[k] - g_y_0_0_xxyzz[k] * cd_y[k] + g_y_0_0_xxyyzz[k];

            g_y_0_y_xxzzz[k] = -g_0_xxzzz[k] - g_y_0_0_xxzzz[k] * cd_y[k] + g_y_0_0_xxyzzz[k];

            g_y_0_y_xyyyy[k] = -g_0_xyyyy[k] - g_y_0_0_xyyyy[k] * cd_y[k] + g_y_0_0_xyyyyy[k];

            g_y_0_y_xyyyz[k] = -g_0_xyyyz[k] - g_y_0_0_xyyyz[k] * cd_y[k] + g_y_0_0_xyyyyz[k];

            g_y_0_y_xyyzz[k] = -g_0_xyyzz[k] - g_y_0_0_xyyzz[k] * cd_y[k] + g_y_0_0_xyyyzz[k];

            g_y_0_y_xyzzz[k] = -g_0_xyzzz[k] - g_y_0_0_xyzzz[k] * cd_y[k] + g_y_0_0_xyyzzz[k];

            g_y_0_y_xzzzz[k] = -g_0_xzzzz[k] - g_y_0_0_xzzzz[k] * cd_y[k] + g_y_0_0_xyzzzz[k];

            g_y_0_y_yyyyy[k] = -g_0_yyyyy[k] - g_y_0_0_yyyyy[k] * cd_y[k] + g_y_0_0_yyyyyy[k];

            g_y_0_y_yyyyz[k] = -g_0_yyyyz[k] - g_y_0_0_yyyyz[k] * cd_y[k] + g_y_0_0_yyyyyz[k];

            g_y_0_y_yyyzz[k] = -g_0_yyyzz[k] - g_y_0_0_yyyzz[k] * cd_y[k] + g_y_0_0_yyyyzz[k];

            g_y_0_y_yyzzz[k] = -g_0_yyzzz[k] - g_y_0_0_yyzzz[k] * cd_y[k] + g_y_0_0_yyyzzz[k];

            g_y_0_y_yzzzz[k] = -g_0_yzzzz[k] - g_y_0_0_yzzzz[k] * cd_y[k] + g_y_0_0_yyzzzz[k];

            g_y_0_y_zzzzz[k] = -g_0_zzzzz[k] - g_y_0_0_zzzzz[k] * cd_y[k] + g_y_0_0_yzzzzz[k];
        }

        /// Set up 42-63 components of targeted buffer : cbuffer.data(

        auto g_y_0_z_xxxxx = cbuffer.data(ph_geom_10_off + 63 * acomps  + 42);

        auto g_y_0_z_xxxxy = cbuffer.data(ph_geom_10_off + 63 * acomps  + 43);

        auto g_y_0_z_xxxxz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 44);

        auto g_y_0_z_xxxyy = cbuffer.data(ph_geom_10_off + 63 * acomps  + 45);

        auto g_y_0_z_xxxyz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 46);

        auto g_y_0_z_xxxzz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 47);

        auto g_y_0_z_xxyyy = cbuffer.data(ph_geom_10_off + 63 * acomps  + 48);

        auto g_y_0_z_xxyyz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 49);

        auto g_y_0_z_xxyzz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 50);

        auto g_y_0_z_xxzzz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 51);

        auto g_y_0_z_xyyyy = cbuffer.data(ph_geom_10_off + 63 * acomps  + 52);

        auto g_y_0_z_xyyyz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 53);

        auto g_y_0_z_xyyzz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 54);

        auto g_y_0_z_xyzzz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 55);

        auto g_y_0_z_xzzzz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 56);

        auto g_y_0_z_yyyyy = cbuffer.data(ph_geom_10_off + 63 * acomps  + 57);

        auto g_y_0_z_yyyyz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 58);

        auto g_y_0_z_yyyzz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 59);

        auto g_y_0_z_yyzzz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 60);

        auto g_y_0_z_yzzzz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 61);

        auto g_y_0_z_zzzzz = cbuffer.data(ph_geom_10_off + 63 * acomps  + 62);

        #pragma omp simd aligned(cd_z, g_y_0_0_xxxxx, g_y_0_0_xxxxxz, g_y_0_0_xxxxy, g_y_0_0_xxxxyz, g_y_0_0_xxxxz, g_y_0_0_xxxxzz, g_y_0_0_xxxyy, g_y_0_0_xxxyyz, g_y_0_0_xxxyz, g_y_0_0_xxxyzz, g_y_0_0_xxxzz, g_y_0_0_xxxzzz, g_y_0_0_xxyyy, g_y_0_0_xxyyyz, g_y_0_0_xxyyz, g_y_0_0_xxyyzz, g_y_0_0_xxyzz, g_y_0_0_xxyzzz, g_y_0_0_xxzzz, g_y_0_0_xxzzzz, g_y_0_0_xyyyy, g_y_0_0_xyyyyz, g_y_0_0_xyyyz, g_y_0_0_xyyyzz, g_y_0_0_xyyzz, g_y_0_0_xyyzzz, g_y_0_0_xyzzz, g_y_0_0_xyzzzz, g_y_0_0_xzzzz, g_y_0_0_xzzzzz, g_y_0_0_yyyyy, g_y_0_0_yyyyyz, g_y_0_0_yyyyz, g_y_0_0_yyyyzz, g_y_0_0_yyyzz, g_y_0_0_yyyzzz, g_y_0_0_yyzzz, g_y_0_0_yyzzzz, g_y_0_0_yzzzz, g_y_0_0_yzzzzz, g_y_0_0_zzzzz, g_y_0_0_zzzzzz, g_y_0_z_xxxxx, g_y_0_z_xxxxy, g_y_0_z_xxxxz, g_y_0_z_xxxyy, g_y_0_z_xxxyz, g_y_0_z_xxxzz, g_y_0_z_xxyyy, g_y_0_z_xxyyz, g_y_0_z_xxyzz, g_y_0_z_xxzzz, g_y_0_z_xyyyy, g_y_0_z_xyyyz, g_y_0_z_xyyzz, g_y_0_z_xyzzz, g_y_0_z_xzzzz, g_y_0_z_yyyyy, g_y_0_z_yyyyz, g_y_0_z_yyyzz, g_y_0_z_yyzzz, g_y_0_z_yzzzz, g_y_0_z_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_z_xxxxx[k] = -g_y_0_0_xxxxx[k] * cd_z[k] + g_y_0_0_xxxxxz[k];

            g_y_0_z_xxxxy[k] = -g_y_0_0_xxxxy[k] * cd_z[k] + g_y_0_0_xxxxyz[k];

            g_y_0_z_xxxxz[k] = -g_y_0_0_xxxxz[k] * cd_z[k] + g_y_0_0_xxxxzz[k];

            g_y_0_z_xxxyy[k] = -g_y_0_0_xxxyy[k] * cd_z[k] + g_y_0_0_xxxyyz[k];

            g_y_0_z_xxxyz[k] = -g_y_0_0_xxxyz[k] * cd_z[k] + g_y_0_0_xxxyzz[k];

            g_y_0_z_xxxzz[k] = -g_y_0_0_xxxzz[k] * cd_z[k] + g_y_0_0_xxxzzz[k];

            g_y_0_z_xxyyy[k] = -g_y_0_0_xxyyy[k] * cd_z[k] + g_y_0_0_xxyyyz[k];

            g_y_0_z_xxyyz[k] = -g_y_0_0_xxyyz[k] * cd_z[k] + g_y_0_0_xxyyzz[k];

            g_y_0_z_xxyzz[k] = -g_y_0_0_xxyzz[k] * cd_z[k] + g_y_0_0_xxyzzz[k];

            g_y_0_z_xxzzz[k] = -g_y_0_0_xxzzz[k] * cd_z[k] + g_y_0_0_xxzzzz[k];

            g_y_0_z_xyyyy[k] = -g_y_0_0_xyyyy[k] * cd_z[k] + g_y_0_0_xyyyyz[k];

            g_y_0_z_xyyyz[k] = -g_y_0_0_xyyyz[k] * cd_z[k] + g_y_0_0_xyyyzz[k];

            g_y_0_z_xyyzz[k] = -g_y_0_0_xyyzz[k] * cd_z[k] + g_y_0_0_xyyzzz[k];

            g_y_0_z_xyzzz[k] = -g_y_0_0_xyzzz[k] * cd_z[k] + g_y_0_0_xyzzzz[k];

            g_y_0_z_xzzzz[k] = -g_y_0_0_xzzzz[k] * cd_z[k] + g_y_0_0_xzzzzz[k];

            g_y_0_z_yyyyy[k] = -g_y_0_0_yyyyy[k] * cd_z[k] + g_y_0_0_yyyyyz[k];

            g_y_0_z_yyyyz[k] = -g_y_0_0_yyyyz[k] * cd_z[k] + g_y_0_0_yyyyzz[k];

            g_y_0_z_yyyzz[k] = -g_y_0_0_yyyzz[k] * cd_z[k] + g_y_0_0_yyyzzz[k];

            g_y_0_z_yyzzz[k] = -g_y_0_0_yyzzz[k] * cd_z[k] + g_y_0_0_yyzzzz[k];

            g_y_0_z_yzzzz[k] = -g_y_0_0_yzzzz[k] * cd_z[k] + g_y_0_0_yzzzzz[k];

            g_y_0_z_zzzzz[k] = -g_y_0_0_zzzzz[k] * cd_z[k] + g_y_0_0_zzzzzz[k];
        }
        /// Set up 0-21 components of targeted buffer : cbuffer.data(

        auto g_z_0_x_xxxxx = cbuffer.data(ph_geom_10_off + 126 * acomps  + 0);

        auto g_z_0_x_xxxxy = cbuffer.data(ph_geom_10_off + 126 * acomps  + 1);

        auto g_z_0_x_xxxxz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 2);

        auto g_z_0_x_xxxyy = cbuffer.data(ph_geom_10_off + 126 * acomps  + 3);

        auto g_z_0_x_xxxyz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 4);

        auto g_z_0_x_xxxzz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 5);

        auto g_z_0_x_xxyyy = cbuffer.data(ph_geom_10_off + 126 * acomps  + 6);

        auto g_z_0_x_xxyyz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 7);

        auto g_z_0_x_xxyzz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 8);

        auto g_z_0_x_xxzzz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 9);

        auto g_z_0_x_xyyyy = cbuffer.data(ph_geom_10_off + 126 * acomps  + 10);

        auto g_z_0_x_xyyyz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 11);

        auto g_z_0_x_xyyzz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 12);

        auto g_z_0_x_xyzzz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 13);

        auto g_z_0_x_xzzzz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 14);

        auto g_z_0_x_yyyyy = cbuffer.data(ph_geom_10_off + 126 * acomps  + 15);

        auto g_z_0_x_yyyyz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 16);

        auto g_z_0_x_yyyzz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 17);

        auto g_z_0_x_yyzzz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 18);

        auto g_z_0_x_yzzzz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 19);

        auto g_z_0_x_zzzzz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 20);

        #pragma omp simd aligned(cd_x, g_z_0_0_xxxxx, g_z_0_0_xxxxxx, g_z_0_0_xxxxxy, g_z_0_0_xxxxxz, g_z_0_0_xxxxy, g_z_0_0_xxxxyy, g_z_0_0_xxxxyz, g_z_0_0_xxxxz, g_z_0_0_xxxxzz, g_z_0_0_xxxyy, g_z_0_0_xxxyyy, g_z_0_0_xxxyyz, g_z_0_0_xxxyz, g_z_0_0_xxxyzz, g_z_0_0_xxxzz, g_z_0_0_xxxzzz, g_z_0_0_xxyyy, g_z_0_0_xxyyyy, g_z_0_0_xxyyyz, g_z_0_0_xxyyz, g_z_0_0_xxyyzz, g_z_0_0_xxyzz, g_z_0_0_xxyzzz, g_z_0_0_xxzzz, g_z_0_0_xxzzzz, g_z_0_0_xyyyy, g_z_0_0_xyyyyy, g_z_0_0_xyyyyz, g_z_0_0_xyyyz, g_z_0_0_xyyyzz, g_z_0_0_xyyzz, g_z_0_0_xyyzzz, g_z_0_0_xyzzz, g_z_0_0_xyzzzz, g_z_0_0_xzzzz, g_z_0_0_xzzzzz, g_z_0_0_yyyyy, g_z_0_0_yyyyz, g_z_0_0_yyyzz, g_z_0_0_yyzzz, g_z_0_0_yzzzz, g_z_0_0_zzzzz, g_z_0_x_xxxxx, g_z_0_x_xxxxy, g_z_0_x_xxxxz, g_z_0_x_xxxyy, g_z_0_x_xxxyz, g_z_0_x_xxxzz, g_z_0_x_xxyyy, g_z_0_x_xxyyz, g_z_0_x_xxyzz, g_z_0_x_xxzzz, g_z_0_x_xyyyy, g_z_0_x_xyyyz, g_z_0_x_xyyzz, g_z_0_x_xyzzz, g_z_0_x_xzzzz, g_z_0_x_yyyyy, g_z_0_x_yyyyz, g_z_0_x_yyyzz, g_z_0_x_yyzzz, g_z_0_x_yzzzz, g_z_0_x_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_x_xxxxx[k] = -g_z_0_0_xxxxx[k] * cd_x[k] + g_z_0_0_xxxxxx[k];

            g_z_0_x_xxxxy[k] = -g_z_0_0_xxxxy[k] * cd_x[k] + g_z_0_0_xxxxxy[k];

            g_z_0_x_xxxxz[k] = -g_z_0_0_xxxxz[k] * cd_x[k] + g_z_0_0_xxxxxz[k];

            g_z_0_x_xxxyy[k] = -g_z_0_0_xxxyy[k] * cd_x[k] + g_z_0_0_xxxxyy[k];

            g_z_0_x_xxxyz[k] = -g_z_0_0_xxxyz[k] * cd_x[k] + g_z_0_0_xxxxyz[k];

            g_z_0_x_xxxzz[k] = -g_z_0_0_xxxzz[k] * cd_x[k] + g_z_0_0_xxxxzz[k];

            g_z_0_x_xxyyy[k] = -g_z_0_0_xxyyy[k] * cd_x[k] + g_z_0_0_xxxyyy[k];

            g_z_0_x_xxyyz[k] = -g_z_0_0_xxyyz[k] * cd_x[k] + g_z_0_0_xxxyyz[k];

            g_z_0_x_xxyzz[k] = -g_z_0_0_xxyzz[k] * cd_x[k] + g_z_0_0_xxxyzz[k];

            g_z_0_x_xxzzz[k] = -g_z_0_0_xxzzz[k] * cd_x[k] + g_z_0_0_xxxzzz[k];

            g_z_0_x_xyyyy[k] = -g_z_0_0_xyyyy[k] * cd_x[k] + g_z_0_0_xxyyyy[k];

            g_z_0_x_xyyyz[k] = -g_z_0_0_xyyyz[k] * cd_x[k] + g_z_0_0_xxyyyz[k];

            g_z_0_x_xyyzz[k] = -g_z_0_0_xyyzz[k] * cd_x[k] + g_z_0_0_xxyyzz[k];

            g_z_0_x_xyzzz[k] = -g_z_0_0_xyzzz[k] * cd_x[k] + g_z_0_0_xxyzzz[k];

            g_z_0_x_xzzzz[k] = -g_z_0_0_xzzzz[k] * cd_x[k] + g_z_0_0_xxzzzz[k];

            g_z_0_x_yyyyy[k] = -g_z_0_0_yyyyy[k] * cd_x[k] + g_z_0_0_xyyyyy[k];

            g_z_0_x_yyyyz[k] = -g_z_0_0_yyyyz[k] * cd_x[k] + g_z_0_0_xyyyyz[k];

            g_z_0_x_yyyzz[k] = -g_z_0_0_yyyzz[k] * cd_x[k] + g_z_0_0_xyyyzz[k];

            g_z_0_x_yyzzz[k] = -g_z_0_0_yyzzz[k] * cd_x[k] + g_z_0_0_xyyzzz[k];

            g_z_0_x_yzzzz[k] = -g_z_0_0_yzzzz[k] * cd_x[k] + g_z_0_0_xyzzzz[k];

            g_z_0_x_zzzzz[k] = -g_z_0_0_zzzzz[k] * cd_x[k] + g_z_0_0_xzzzzz[k];
        }

        /// Set up 21-42 components of targeted buffer : cbuffer.data(

        auto g_z_0_y_xxxxx = cbuffer.data(ph_geom_10_off + 126 * acomps  + 21);

        auto g_z_0_y_xxxxy = cbuffer.data(ph_geom_10_off + 126 * acomps  + 22);

        auto g_z_0_y_xxxxz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 23);

        auto g_z_0_y_xxxyy = cbuffer.data(ph_geom_10_off + 126 * acomps  + 24);

        auto g_z_0_y_xxxyz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 25);

        auto g_z_0_y_xxxzz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 26);

        auto g_z_0_y_xxyyy = cbuffer.data(ph_geom_10_off + 126 * acomps  + 27);

        auto g_z_0_y_xxyyz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 28);

        auto g_z_0_y_xxyzz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 29);

        auto g_z_0_y_xxzzz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 30);

        auto g_z_0_y_xyyyy = cbuffer.data(ph_geom_10_off + 126 * acomps  + 31);

        auto g_z_0_y_xyyyz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 32);

        auto g_z_0_y_xyyzz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 33);

        auto g_z_0_y_xyzzz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 34);

        auto g_z_0_y_xzzzz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 35);

        auto g_z_0_y_yyyyy = cbuffer.data(ph_geom_10_off + 126 * acomps  + 36);

        auto g_z_0_y_yyyyz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 37);

        auto g_z_0_y_yyyzz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 38);

        auto g_z_0_y_yyzzz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 39);

        auto g_z_0_y_yzzzz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 40);

        auto g_z_0_y_zzzzz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 41);

        #pragma omp simd aligned(cd_y, g_z_0_0_xxxxx, g_z_0_0_xxxxxy, g_z_0_0_xxxxy, g_z_0_0_xxxxyy, g_z_0_0_xxxxyz, g_z_0_0_xxxxz, g_z_0_0_xxxyy, g_z_0_0_xxxyyy, g_z_0_0_xxxyyz, g_z_0_0_xxxyz, g_z_0_0_xxxyzz, g_z_0_0_xxxzz, g_z_0_0_xxyyy, g_z_0_0_xxyyyy, g_z_0_0_xxyyyz, g_z_0_0_xxyyz, g_z_0_0_xxyyzz, g_z_0_0_xxyzz, g_z_0_0_xxyzzz, g_z_0_0_xxzzz, g_z_0_0_xyyyy, g_z_0_0_xyyyyy, g_z_0_0_xyyyyz, g_z_0_0_xyyyz, g_z_0_0_xyyyzz, g_z_0_0_xyyzz, g_z_0_0_xyyzzz, g_z_0_0_xyzzz, g_z_0_0_xyzzzz, g_z_0_0_xzzzz, g_z_0_0_yyyyy, g_z_0_0_yyyyyy, g_z_0_0_yyyyyz, g_z_0_0_yyyyz, g_z_0_0_yyyyzz, g_z_0_0_yyyzz, g_z_0_0_yyyzzz, g_z_0_0_yyzzz, g_z_0_0_yyzzzz, g_z_0_0_yzzzz, g_z_0_0_yzzzzz, g_z_0_0_zzzzz, g_z_0_y_xxxxx, g_z_0_y_xxxxy, g_z_0_y_xxxxz, g_z_0_y_xxxyy, g_z_0_y_xxxyz, g_z_0_y_xxxzz, g_z_0_y_xxyyy, g_z_0_y_xxyyz, g_z_0_y_xxyzz, g_z_0_y_xxzzz, g_z_0_y_xyyyy, g_z_0_y_xyyyz, g_z_0_y_xyyzz, g_z_0_y_xyzzz, g_z_0_y_xzzzz, g_z_0_y_yyyyy, g_z_0_y_yyyyz, g_z_0_y_yyyzz, g_z_0_y_yyzzz, g_z_0_y_yzzzz, g_z_0_y_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_y_xxxxx[k] = -g_z_0_0_xxxxx[k] * cd_y[k] + g_z_0_0_xxxxxy[k];

            g_z_0_y_xxxxy[k] = -g_z_0_0_xxxxy[k] * cd_y[k] + g_z_0_0_xxxxyy[k];

            g_z_0_y_xxxxz[k] = -g_z_0_0_xxxxz[k] * cd_y[k] + g_z_0_0_xxxxyz[k];

            g_z_0_y_xxxyy[k] = -g_z_0_0_xxxyy[k] * cd_y[k] + g_z_0_0_xxxyyy[k];

            g_z_0_y_xxxyz[k] = -g_z_0_0_xxxyz[k] * cd_y[k] + g_z_0_0_xxxyyz[k];

            g_z_0_y_xxxzz[k] = -g_z_0_0_xxxzz[k] * cd_y[k] + g_z_0_0_xxxyzz[k];

            g_z_0_y_xxyyy[k] = -g_z_0_0_xxyyy[k] * cd_y[k] + g_z_0_0_xxyyyy[k];

            g_z_0_y_xxyyz[k] = -g_z_0_0_xxyyz[k] * cd_y[k] + g_z_0_0_xxyyyz[k];

            g_z_0_y_xxyzz[k] = -g_z_0_0_xxyzz[k] * cd_y[k] + g_z_0_0_xxyyzz[k];

            g_z_0_y_xxzzz[k] = -g_z_0_0_xxzzz[k] * cd_y[k] + g_z_0_0_xxyzzz[k];

            g_z_0_y_xyyyy[k] = -g_z_0_0_xyyyy[k] * cd_y[k] + g_z_0_0_xyyyyy[k];

            g_z_0_y_xyyyz[k] = -g_z_0_0_xyyyz[k] * cd_y[k] + g_z_0_0_xyyyyz[k];

            g_z_0_y_xyyzz[k] = -g_z_0_0_xyyzz[k] * cd_y[k] + g_z_0_0_xyyyzz[k];

            g_z_0_y_xyzzz[k] = -g_z_0_0_xyzzz[k] * cd_y[k] + g_z_0_0_xyyzzz[k];

            g_z_0_y_xzzzz[k] = -g_z_0_0_xzzzz[k] * cd_y[k] + g_z_0_0_xyzzzz[k];

            g_z_0_y_yyyyy[k] = -g_z_0_0_yyyyy[k] * cd_y[k] + g_z_0_0_yyyyyy[k];

            g_z_0_y_yyyyz[k] = -g_z_0_0_yyyyz[k] * cd_y[k] + g_z_0_0_yyyyyz[k];

            g_z_0_y_yyyzz[k] = -g_z_0_0_yyyzz[k] * cd_y[k] + g_z_0_0_yyyyzz[k];

            g_z_0_y_yyzzz[k] = -g_z_0_0_yyzzz[k] * cd_y[k] + g_z_0_0_yyyzzz[k];

            g_z_0_y_yzzzz[k] = -g_z_0_0_yzzzz[k] * cd_y[k] + g_z_0_0_yyzzzz[k];

            g_z_0_y_zzzzz[k] = -g_z_0_0_zzzzz[k] * cd_y[k] + g_z_0_0_yzzzzz[k];
        }

        /// Set up 42-63 components of targeted buffer : cbuffer.data(

        auto g_z_0_z_xxxxx = cbuffer.data(ph_geom_10_off + 126 * acomps  + 42);

        auto g_z_0_z_xxxxy = cbuffer.data(ph_geom_10_off + 126 * acomps  + 43);

        auto g_z_0_z_xxxxz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 44);

        auto g_z_0_z_xxxyy = cbuffer.data(ph_geom_10_off + 126 * acomps  + 45);

        auto g_z_0_z_xxxyz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 46);

        auto g_z_0_z_xxxzz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 47);

        auto g_z_0_z_xxyyy = cbuffer.data(ph_geom_10_off + 126 * acomps  + 48);

        auto g_z_0_z_xxyyz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 49);

        auto g_z_0_z_xxyzz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 50);

        auto g_z_0_z_xxzzz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 51);

        auto g_z_0_z_xyyyy = cbuffer.data(ph_geom_10_off + 126 * acomps  + 52);

        auto g_z_0_z_xyyyz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 53);

        auto g_z_0_z_xyyzz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 54);

        auto g_z_0_z_xyzzz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 55);

        auto g_z_0_z_xzzzz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 56);

        auto g_z_0_z_yyyyy = cbuffer.data(ph_geom_10_off + 126 * acomps  + 57);

        auto g_z_0_z_yyyyz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 58);

        auto g_z_0_z_yyyzz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 59);

        auto g_z_0_z_yyzzz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 60);

        auto g_z_0_z_yzzzz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 61);

        auto g_z_0_z_zzzzz = cbuffer.data(ph_geom_10_off + 126 * acomps  + 62);

        #pragma omp simd aligned(cd_z, g_0_xxxxx, g_0_xxxxy, g_0_xxxxz, g_0_xxxyy, g_0_xxxyz, g_0_xxxzz, g_0_xxyyy, g_0_xxyyz, g_0_xxyzz, g_0_xxzzz, g_0_xyyyy, g_0_xyyyz, g_0_xyyzz, g_0_xyzzz, g_0_xzzzz, g_0_yyyyy, g_0_yyyyz, g_0_yyyzz, g_0_yyzzz, g_0_yzzzz, g_0_zzzzz, g_z_0_0_xxxxx, g_z_0_0_xxxxxz, g_z_0_0_xxxxy, g_z_0_0_xxxxyz, g_z_0_0_xxxxz, g_z_0_0_xxxxzz, g_z_0_0_xxxyy, g_z_0_0_xxxyyz, g_z_0_0_xxxyz, g_z_0_0_xxxyzz, g_z_0_0_xxxzz, g_z_0_0_xxxzzz, g_z_0_0_xxyyy, g_z_0_0_xxyyyz, g_z_0_0_xxyyz, g_z_0_0_xxyyzz, g_z_0_0_xxyzz, g_z_0_0_xxyzzz, g_z_0_0_xxzzz, g_z_0_0_xxzzzz, g_z_0_0_xyyyy, g_z_0_0_xyyyyz, g_z_0_0_xyyyz, g_z_0_0_xyyyzz, g_z_0_0_xyyzz, g_z_0_0_xyyzzz, g_z_0_0_xyzzz, g_z_0_0_xyzzzz, g_z_0_0_xzzzz, g_z_0_0_xzzzzz, g_z_0_0_yyyyy, g_z_0_0_yyyyyz, g_z_0_0_yyyyz, g_z_0_0_yyyyzz, g_z_0_0_yyyzz, g_z_0_0_yyyzzz, g_z_0_0_yyzzz, g_z_0_0_yyzzzz, g_z_0_0_yzzzz, g_z_0_0_yzzzzz, g_z_0_0_zzzzz, g_z_0_0_zzzzzz, g_z_0_z_xxxxx, g_z_0_z_xxxxy, g_z_0_z_xxxxz, g_z_0_z_xxxyy, g_z_0_z_xxxyz, g_z_0_z_xxxzz, g_z_0_z_xxyyy, g_z_0_z_xxyyz, g_z_0_z_xxyzz, g_z_0_z_xxzzz, g_z_0_z_xyyyy, g_z_0_z_xyyyz, g_z_0_z_xyyzz, g_z_0_z_xyzzz, g_z_0_z_xzzzz, g_z_0_z_yyyyy, g_z_0_z_yyyyz, g_z_0_z_yyyzz, g_z_0_z_yyzzz, g_z_0_z_yzzzz, g_z_0_z_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_z_xxxxx[k] = -g_0_xxxxx[k] - g_z_0_0_xxxxx[k] * cd_z[k] + g_z_0_0_xxxxxz[k];

            g_z_0_z_xxxxy[k] = -g_0_xxxxy[k] - g_z_0_0_xxxxy[k] * cd_z[k] + g_z_0_0_xxxxyz[k];

            g_z_0_z_xxxxz[k] = -g_0_xxxxz[k] - g_z_0_0_xxxxz[k] * cd_z[k] + g_z_0_0_xxxxzz[k];

            g_z_0_z_xxxyy[k] = -g_0_xxxyy[k] - g_z_0_0_xxxyy[k] * cd_z[k] + g_z_0_0_xxxyyz[k];

            g_z_0_z_xxxyz[k] = -g_0_xxxyz[k] - g_z_0_0_xxxyz[k] * cd_z[k] + g_z_0_0_xxxyzz[k];

            g_z_0_z_xxxzz[k] = -g_0_xxxzz[k] - g_z_0_0_xxxzz[k] * cd_z[k] + g_z_0_0_xxxzzz[k];

            g_z_0_z_xxyyy[k] = -g_0_xxyyy[k] - g_z_0_0_xxyyy[k] * cd_z[k] + g_z_0_0_xxyyyz[k];

            g_z_0_z_xxyyz[k] = -g_0_xxyyz[k] - g_z_0_0_xxyyz[k] * cd_z[k] + g_z_0_0_xxyyzz[k];

            g_z_0_z_xxyzz[k] = -g_0_xxyzz[k] - g_z_0_0_xxyzz[k] * cd_z[k] + g_z_0_0_xxyzzz[k];

            g_z_0_z_xxzzz[k] = -g_0_xxzzz[k] - g_z_0_0_xxzzz[k] * cd_z[k] + g_z_0_0_xxzzzz[k];

            g_z_0_z_xyyyy[k] = -g_0_xyyyy[k] - g_z_0_0_xyyyy[k] * cd_z[k] + g_z_0_0_xyyyyz[k];

            g_z_0_z_xyyyz[k] = -g_0_xyyyz[k] - g_z_0_0_xyyyz[k] * cd_z[k] + g_z_0_0_xyyyzz[k];

            g_z_0_z_xyyzz[k] = -g_0_xyyzz[k] - g_z_0_0_xyyzz[k] * cd_z[k] + g_z_0_0_xyyzzz[k];

            g_z_0_z_xyzzz[k] = -g_0_xyzzz[k] - g_z_0_0_xyzzz[k] * cd_z[k] + g_z_0_0_xyzzzz[k];

            g_z_0_z_xzzzz[k] = -g_0_xzzzz[k] - g_z_0_0_xzzzz[k] * cd_z[k] + g_z_0_0_xzzzzz[k];

            g_z_0_z_yyyyy[k] = -g_0_yyyyy[k] - g_z_0_0_yyyyy[k] * cd_z[k] + g_z_0_0_yyyyyz[k];

            g_z_0_z_yyyyz[k] = -g_0_yyyyz[k] - g_z_0_0_yyyyz[k] * cd_z[k] + g_z_0_0_yyyyzz[k];

            g_z_0_z_yyyzz[k] = -g_0_yyyzz[k] - g_z_0_0_yyyzz[k] * cd_z[k] + g_z_0_0_yyyzzz[k];

            g_z_0_z_yyzzz[k] = -g_0_yyzzz[k] - g_z_0_0_yyzzz[k] * cd_z[k] + g_z_0_0_yyzzzz[k];

            g_z_0_z_yzzzz[k] = -g_0_yzzzz[k] - g_z_0_0_yzzzz[k] * cd_z[k] + g_z_0_0_yzzzzz[k];

            g_z_0_z_zzzzz[k] = -g_0_zzzzz[k] - g_z_0_0_zzzzz[k] * cd_z[k] + g_z_0_0_zzzzzz[k];
        }
    }
}

} // t3ceri namespace

