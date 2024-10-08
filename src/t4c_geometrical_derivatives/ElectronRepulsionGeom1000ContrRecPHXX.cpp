#include "ElectronRepulsionGeom1000ContrRecPHXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom10_hrr_electron_repulsion_phxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_phxx,
                                            const size_t idx_shxx,
                                            const size_t idx_geom_10_shxx,
                                            const size_t idx_geom_10_sixx,
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

            /// Set up components of auxilary buffer : SHSS

            const auto sh_geom_10_off = idx_geom_10_shxx + i * dcomps + j;

            auto g_x_0_0_xxxxx = cbuffer.data(sh_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_0_xxxxy = cbuffer.data(sh_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_0_xxxxz = cbuffer.data(sh_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_0_xxxyy = cbuffer.data(sh_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_0_xxxyz = cbuffer.data(sh_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_0_xxxzz = cbuffer.data(sh_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_0_xxyyy = cbuffer.data(sh_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_0_xxyyz = cbuffer.data(sh_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_0_xxyzz = cbuffer.data(sh_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_0_xxzzz = cbuffer.data(sh_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_0_xyyyy = cbuffer.data(sh_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_0_xyyyz = cbuffer.data(sh_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_0_xyyzz = cbuffer.data(sh_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_0_xyzzz = cbuffer.data(sh_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_0_xzzzz = cbuffer.data(sh_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_0_yyyyy = cbuffer.data(sh_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_0_yyyyz = cbuffer.data(sh_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_0_yyyzz = cbuffer.data(sh_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_0_yyzzz = cbuffer.data(sh_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_0_yzzzz = cbuffer.data(sh_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_0_zzzzz = cbuffer.data(sh_geom_10_off + 20 * ccomps * dcomps);

            auto g_y_0_0_xxxxx = cbuffer.data(sh_geom_10_off + 21 * ccomps * dcomps);

            auto g_y_0_0_xxxxy = cbuffer.data(sh_geom_10_off + 22 * ccomps * dcomps);

            auto g_y_0_0_xxxxz = cbuffer.data(sh_geom_10_off + 23 * ccomps * dcomps);

            auto g_y_0_0_xxxyy = cbuffer.data(sh_geom_10_off + 24 * ccomps * dcomps);

            auto g_y_0_0_xxxyz = cbuffer.data(sh_geom_10_off + 25 * ccomps * dcomps);

            auto g_y_0_0_xxxzz = cbuffer.data(sh_geom_10_off + 26 * ccomps * dcomps);

            auto g_y_0_0_xxyyy = cbuffer.data(sh_geom_10_off + 27 * ccomps * dcomps);

            auto g_y_0_0_xxyyz = cbuffer.data(sh_geom_10_off + 28 * ccomps * dcomps);

            auto g_y_0_0_xxyzz = cbuffer.data(sh_geom_10_off + 29 * ccomps * dcomps);

            auto g_y_0_0_xxzzz = cbuffer.data(sh_geom_10_off + 30 * ccomps * dcomps);

            auto g_y_0_0_xyyyy = cbuffer.data(sh_geom_10_off + 31 * ccomps * dcomps);

            auto g_y_0_0_xyyyz = cbuffer.data(sh_geom_10_off + 32 * ccomps * dcomps);

            auto g_y_0_0_xyyzz = cbuffer.data(sh_geom_10_off + 33 * ccomps * dcomps);

            auto g_y_0_0_xyzzz = cbuffer.data(sh_geom_10_off + 34 * ccomps * dcomps);

            auto g_y_0_0_xzzzz = cbuffer.data(sh_geom_10_off + 35 * ccomps * dcomps);

            auto g_y_0_0_yyyyy = cbuffer.data(sh_geom_10_off + 36 * ccomps * dcomps);

            auto g_y_0_0_yyyyz = cbuffer.data(sh_geom_10_off + 37 * ccomps * dcomps);

            auto g_y_0_0_yyyzz = cbuffer.data(sh_geom_10_off + 38 * ccomps * dcomps);

            auto g_y_0_0_yyzzz = cbuffer.data(sh_geom_10_off + 39 * ccomps * dcomps);

            auto g_y_0_0_yzzzz = cbuffer.data(sh_geom_10_off + 40 * ccomps * dcomps);

            auto g_y_0_0_zzzzz = cbuffer.data(sh_geom_10_off + 41 * ccomps * dcomps);

            auto g_z_0_0_xxxxx = cbuffer.data(sh_geom_10_off + 42 * ccomps * dcomps);

            auto g_z_0_0_xxxxy = cbuffer.data(sh_geom_10_off + 43 * ccomps * dcomps);

            auto g_z_0_0_xxxxz = cbuffer.data(sh_geom_10_off + 44 * ccomps * dcomps);

            auto g_z_0_0_xxxyy = cbuffer.data(sh_geom_10_off + 45 * ccomps * dcomps);

            auto g_z_0_0_xxxyz = cbuffer.data(sh_geom_10_off + 46 * ccomps * dcomps);

            auto g_z_0_0_xxxzz = cbuffer.data(sh_geom_10_off + 47 * ccomps * dcomps);

            auto g_z_0_0_xxyyy = cbuffer.data(sh_geom_10_off + 48 * ccomps * dcomps);

            auto g_z_0_0_xxyyz = cbuffer.data(sh_geom_10_off + 49 * ccomps * dcomps);

            auto g_z_0_0_xxyzz = cbuffer.data(sh_geom_10_off + 50 * ccomps * dcomps);

            auto g_z_0_0_xxzzz = cbuffer.data(sh_geom_10_off + 51 * ccomps * dcomps);

            auto g_z_0_0_xyyyy = cbuffer.data(sh_geom_10_off + 52 * ccomps * dcomps);

            auto g_z_0_0_xyyyz = cbuffer.data(sh_geom_10_off + 53 * ccomps * dcomps);

            auto g_z_0_0_xyyzz = cbuffer.data(sh_geom_10_off + 54 * ccomps * dcomps);

            auto g_z_0_0_xyzzz = cbuffer.data(sh_geom_10_off + 55 * ccomps * dcomps);

            auto g_z_0_0_xzzzz = cbuffer.data(sh_geom_10_off + 56 * ccomps * dcomps);

            auto g_z_0_0_yyyyy = cbuffer.data(sh_geom_10_off + 57 * ccomps * dcomps);

            auto g_z_0_0_yyyyz = cbuffer.data(sh_geom_10_off + 58 * ccomps * dcomps);

            auto g_z_0_0_yyyzz = cbuffer.data(sh_geom_10_off + 59 * ccomps * dcomps);

            auto g_z_0_0_yyzzz = cbuffer.data(sh_geom_10_off + 60 * ccomps * dcomps);

            auto g_z_0_0_yzzzz = cbuffer.data(sh_geom_10_off + 61 * ccomps * dcomps);

            auto g_z_0_0_zzzzz = cbuffer.data(sh_geom_10_off + 62 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SISS

            const auto si_geom_10_off = idx_geom_10_sixx + i * dcomps + j;

            auto g_x_0_0_xxxxxx = cbuffer.data(si_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_0_xxxxxy = cbuffer.data(si_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_0_xxxxxz = cbuffer.data(si_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_0_xxxxyy = cbuffer.data(si_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_0_xxxxyz = cbuffer.data(si_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_0_xxxxzz = cbuffer.data(si_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_0_xxxyyy = cbuffer.data(si_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_0_xxxyyz = cbuffer.data(si_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_0_xxxyzz = cbuffer.data(si_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_0_xxxzzz = cbuffer.data(si_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_0_xxyyyy = cbuffer.data(si_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_0_xxyyyz = cbuffer.data(si_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_0_xxyyzz = cbuffer.data(si_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_0_xxyzzz = cbuffer.data(si_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_0_xxzzzz = cbuffer.data(si_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_0_xyyyyy = cbuffer.data(si_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_0_xyyyyz = cbuffer.data(si_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_0_xyyyzz = cbuffer.data(si_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_0_xyyzzz = cbuffer.data(si_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_0_xyzzzz = cbuffer.data(si_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_0_xzzzzz = cbuffer.data(si_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_0_yyyyyy = cbuffer.data(si_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_0_yyyyyz = cbuffer.data(si_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_0_yyyyzz = cbuffer.data(si_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_0_yyyzzz = cbuffer.data(si_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_0_yyzzzz = cbuffer.data(si_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_0_yzzzzz = cbuffer.data(si_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_0_zzzzzz = cbuffer.data(si_geom_10_off + 27 * ccomps * dcomps);

            auto g_y_0_0_xxxxxx = cbuffer.data(si_geom_10_off + 28 * ccomps * dcomps);

            auto g_y_0_0_xxxxxy = cbuffer.data(si_geom_10_off + 29 * ccomps * dcomps);

            auto g_y_0_0_xxxxxz = cbuffer.data(si_geom_10_off + 30 * ccomps * dcomps);

            auto g_y_0_0_xxxxyy = cbuffer.data(si_geom_10_off + 31 * ccomps * dcomps);

            auto g_y_0_0_xxxxyz = cbuffer.data(si_geom_10_off + 32 * ccomps * dcomps);

            auto g_y_0_0_xxxxzz = cbuffer.data(si_geom_10_off + 33 * ccomps * dcomps);

            auto g_y_0_0_xxxyyy = cbuffer.data(si_geom_10_off + 34 * ccomps * dcomps);

            auto g_y_0_0_xxxyyz = cbuffer.data(si_geom_10_off + 35 * ccomps * dcomps);

            auto g_y_0_0_xxxyzz = cbuffer.data(si_geom_10_off + 36 * ccomps * dcomps);

            auto g_y_0_0_xxxzzz = cbuffer.data(si_geom_10_off + 37 * ccomps * dcomps);

            auto g_y_0_0_xxyyyy = cbuffer.data(si_geom_10_off + 38 * ccomps * dcomps);

            auto g_y_0_0_xxyyyz = cbuffer.data(si_geom_10_off + 39 * ccomps * dcomps);

            auto g_y_0_0_xxyyzz = cbuffer.data(si_geom_10_off + 40 * ccomps * dcomps);

            auto g_y_0_0_xxyzzz = cbuffer.data(si_geom_10_off + 41 * ccomps * dcomps);

            auto g_y_0_0_xxzzzz = cbuffer.data(si_geom_10_off + 42 * ccomps * dcomps);

            auto g_y_0_0_xyyyyy = cbuffer.data(si_geom_10_off + 43 * ccomps * dcomps);

            auto g_y_0_0_xyyyyz = cbuffer.data(si_geom_10_off + 44 * ccomps * dcomps);

            auto g_y_0_0_xyyyzz = cbuffer.data(si_geom_10_off + 45 * ccomps * dcomps);

            auto g_y_0_0_xyyzzz = cbuffer.data(si_geom_10_off + 46 * ccomps * dcomps);

            auto g_y_0_0_xyzzzz = cbuffer.data(si_geom_10_off + 47 * ccomps * dcomps);

            auto g_y_0_0_xzzzzz = cbuffer.data(si_geom_10_off + 48 * ccomps * dcomps);

            auto g_y_0_0_yyyyyy = cbuffer.data(si_geom_10_off + 49 * ccomps * dcomps);

            auto g_y_0_0_yyyyyz = cbuffer.data(si_geom_10_off + 50 * ccomps * dcomps);

            auto g_y_0_0_yyyyzz = cbuffer.data(si_geom_10_off + 51 * ccomps * dcomps);

            auto g_y_0_0_yyyzzz = cbuffer.data(si_geom_10_off + 52 * ccomps * dcomps);

            auto g_y_0_0_yyzzzz = cbuffer.data(si_geom_10_off + 53 * ccomps * dcomps);

            auto g_y_0_0_yzzzzz = cbuffer.data(si_geom_10_off + 54 * ccomps * dcomps);

            auto g_y_0_0_zzzzzz = cbuffer.data(si_geom_10_off + 55 * ccomps * dcomps);

            auto g_z_0_0_xxxxxx = cbuffer.data(si_geom_10_off + 56 * ccomps * dcomps);

            auto g_z_0_0_xxxxxy = cbuffer.data(si_geom_10_off + 57 * ccomps * dcomps);

            auto g_z_0_0_xxxxxz = cbuffer.data(si_geom_10_off + 58 * ccomps * dcomps);

            auto g_z_0_0_xxxxyy = cbuffer.data(si_geom_10_off + 59 * ccomps * dcomps);

            auto g_z_0_0_xxxxyz = cbuffer.data(si_geom_10_off + 60 * ccomps * dcomps);

            auto g_z_0_0_xxxxzz = cbuffer.data(si_geom_10_off + 61 * ccomps * dcomps);

            auto g_z_0_0_xxxyyy = cbuffer.data(si_geom_10_off + 62 * ccomps * dcomps);

            auto g_z_0_0_xxxyyz = cbuffer.data(si_geom_10_off + 63 * ccomps * dcomps);

            auto g_z_0_0_xxxyzz = cbuffer.data(si_geom_10_off + 64 * ccomps * dcomps);

            auto g_z_0_0_xxxzzz = cbuffer.data(si_geom_10_off + 65 * ccomps * dcomps);

            auto g_z_0_0_xxyyyy = cbuffer.data(si_geom_10_off + 66 * ccomps * dcomps);

            auto g_z_0_0_xxyyyz = cbuffer.data(si_geom_10_off + 67 * ccomps * dcomps);

            auto g_z_0_0_xxyyzz = cbuffer.data(si_geom_10_off + 68 * ccomps * dcomps);

            auto g_z_0_0_xxyzzz = cbuffer.data(si_geom_10_off + 69 * ccomps * dcomps);

            auto g_z_0_0_xxzzzz = cbuffer.data(si_geom_10_off + 70 * ccomps * dcomps);

            auto g_z_0_0_xyyyyy = cbuffer.data(si_geom_10_off + 71 * ccomps * dcomps);

            auto g_z_0_0_xyyyyz = cbuffer.data(si_geom_10_off + 72 * ccomps * dcomps);

            auto g_z_0_0_xyyyzz = cbuffer.data(si_geom_10_off + 73 * ccomps * dcomps);

            auto g_z_0_0_xyyzzz = cbuffer.data(si_geom_10_off + 74 * ccomps * dcomps);

            auto g_z_0_0_xyzzzz = cbuffer.data(si_geom_10_off + 75 * ccomps * dcomps);

            auto g_z_0_0_xzzzzz = cbuffer.data(si_geom_10_off + 76 * ccomps * dcomps);

            auto g_z_0_0_yyyyyy = cbuffer.data(si_geom_10_off + 77 * ccomps * dcomps);

            auto g_z_0_0_yyyyyz = cbuffer.data(si_geom_10_off + 78 * ccomps * dcomps);

            auto g_z_0_0_yyyyzz = cbuffer.data(si_geom_10_off + 79 * ccomps * dcomps);

            auto g_z_0_0_yyyzzz = cbuffer.data(si_geom_10_off + 80 * ccomps * dcomps);

            auto g_z_0_0_yyzzzz = cbuffer.data(si_geom_10_off + 81 * ccomps * dcomps);

            auto g_z_0_0_yzzzzz = cbuffer.data(si_geom_10_off + 82 * ccomps * dcomps);

            auto g_z_0_0_zzzzzz = cbuffer.data(si_geom_10_off + 83 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_phxx

            const auto ph_geom_10_off = idx_geom_10_phxx + i * dcomps + j;

            /// Set up 0-21 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_xxxxx = cbuffer.data(ph_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_x_xxxxy = cbuffer.data(ph_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_x_xxxxz = cbuffer.data(ph_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_x_xxxyy = cbuffer.data(ph_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_x_xxxyz = cbuffer.data(ph_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_x_xxxzz = cbuffer.data(ph_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_x_xxyyy = cbuffer.data(ph_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_x_xxyyz = cbuffer.data(ph_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_x_xxyzz = cbuffer.data(ph_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_x_xxzzz = cbuffer.data(ph_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_x_xyyyy = cbuffer.data(ph_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_x_xyyyz = cbuffer.data(ph_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_x_xyyzz = cbuffer.data(ph_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_x_xyzzz = cbuffer.data(ph_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_x_xzzzz = cbuffer.data(ph_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_x_yyyyy = cbuffer.data(ph_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_x_yyyyz = cbuffer.data(ph_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_x_yyyzz = cbuffer.data(ph_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_x_yyzzz = cbuffer.data(ph_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_x_yzzzz = cbuffer.data(ph_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_x_zzzzz = cbuffer.data(ph_geom_10_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_xxxxx, g_0_xxxxy, g_0_xxxxz, g_0_xxxyy, g_0_xxxyz, g_0_xxxzz, g_0_xxyyy, g_0_xxyyz, g_0_xxyzz, g_0_xxzzz, g_0_xyyyy, g_0_xyyyz, g_0_xyyzz, g_0_xyzzz, g_0_xzzzz, g_0_yyyyy, g_0_yyyyz, g_0_yyyzz, g_0_yyzzz, g_0_yzzzz, g_0_zzzzz, g_x_0_0_xxxxx, g_x_0_0_xxxxxx, g_x_0_0_xxxxxy, g_x_0_0_xxxxxz, g_x_0_0_xxxxy, g_x_0_0_xxxxyy, g_x_0_0_xxxxyz, g_x_0_0_xxxxz, g_x_0_0_xxxxzz, g_x_0_0_xxxyy, g_x_0_0_xxxyyy, g_x_0_0_xxxyyz, g_x_0_0_xxxyz, g_x_0_0_xxxyzz, g_x_0_0_xxxzz, g_x_0_0_xxxzzz, g_x_0_0_xxyyy, g_x_0_0_xxyyyy, g_x_0_0_xxyyyz, g_x_0_0_xxyyz, g_x_0_0_xxyyzz, g_x_0_0_xxyzz, g_x_0_0_xxyzzz, g_x_0_0_xxzzz, g_x_0_0_xxzzzz, g_x_0_0_xyyyy, g_x_0_0_xyyyyy, g_x_0_0_xyyyyz, g_x_0_0_xyyyz, g_x_0_0_xyyyzz, g_x_0_0_xyyzz, g_x_0_0_xyyzzz, g_x_0_0_xyzzz, g_x_0_0_xyzzzz, g_x_0_0_xzzzz, g_x_0_0_xzzzzz, g_x_0_0_yyyyy, g_x_0_0_yyyyz, g_x_0_0_yyyzz, g_x_0_0_yyzzz, g_x_0_0_yzzzz, g_x_0_0_zzzzz, g_x_0_x_xxxxx, g_x_0_x_xxxxy, g_x_0_x_xxxxz, g_x_0_x_xxxyy, g_x_0_x_xxxyz, g_x_0_x_xxxzz, g_x_0_x_xxyyy, g_x_0_x_xxyyz, g_x_0_x_xxyzz, g_x_0_x_xxzzz, g_x_0_x_xyyyy, g_x_0_x_xyyyz, g_x_0_x_xyyzz, g_x_0_x_xyzzz, g_x_0_x_xzzzz, g_x_0_x_yyyyy, g_x_0_x_yyyyz, g_x_0_x_yyyzz, g_x_0_x_yyzzz, g_x_0_x_yzzzz, g_x_0_x_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_xxxxx[k] = -g_0_xxxxx[k] - g_x_0_0_xxxxx[k] * ab_x + g_x_0_0_xxxxxx[k];

                g_x_0_x_xxxxy[k] = -g_0_xxxxy[k] - g_x_0_0_xxxxy[k] * ab_x + g_x_0_0_xxxxxy[k];

                g_x_0_x_xxxxz[k] = -g_0_xxxxz[k] - g_x_0_0_xxxxz[k] * ab_x + g_x_0_0_xxxxxz[k];

                g_x_0_x_xxxyy[k] = -g_0_xxxyy[k] - g_x_0_0_xxxyy[k] * ab_x + g_x_0_0_xxxxyy[k];

                g_x_0_x_xxxyz[k] = -g_0_xxxyz[k] - g_x_0_0_xxxyz[k] * ab_x + g_x_0_0_xxxxyz[k];

                g_x_0_x_xxxzz[k] = -g_0_xxxzz[k] - g_x_0_0_xxxzz[k] * ab_x + g_x_0_0_xxxxzz[k];

                g_x_0_x_xxyyy[k] = -g_0_xxyyy[k] - g_x_0_0_xxyyy[k] * ab_x + g_x_0_0_xxxyyy[k];

                g_x_0_x_xxyyz[k] = -g_0_xxyyz[k] - g_x_0_0_xxyyz[k] * ab_x + g_x_0_0_xxxyyz[k];

                g_x_0_x_xxyzz[k] = -g_0_xxyzz[k] - g_x_0_0_xxyzz[k] * ab_x + g_x_0_0_xxxyzz[k];

                g_x_0_x_xxzzz[k] = -g_0_xxzzz[k] - g_x_0_0_xxzzz[k] * ab_x + g_x_0_0_xxxzzz[k];

                g_x_0_x_xyyyy[k] = -g_0_xyyyy[k] - g_x_0_0_xyyyy[k] * ab_x + g_x_0_0_xxyyyy[k];

                g_x_0_x_xyyyz[k] = -g_0_xyyyz[k] - g_x_0_0_xyyyz[k] * ab_x + g_x_0_0_xxyyyz[k];

                g_x_0_x_xyyzz[k] = -g_0_xyyzz[k] - g_x_0_0_xyyzz[k] * ab_x + g_x_0_0_xxyyzz[k];

                g_x_0_x_xyzzz[k] = -g_0_xyzzz[k] - g_x_0_0_xyzzz[k] * ab_x + g_x_0_0_xxyzzz[k];

                g_x_0_x_xzzzz[k] = -g_0_xzzzz[k] - g_x_0_0_xzzzz[k] * ab_x + g_x_0_0_xxzzzz[k];

                g_x_0_x_yyyyy[k] = -g_0_yyyyy[k] - g_x_0_0_yyyyy[k] * ab_x + g_x_0_0_xyyyyy[k];

                g_x_0_x_yyyyz[k] = -g_0_yyyyz[k] - g_x_0_0_yyyyz[k] * ab_x + g_x_0_0_xyyyyz[k];

                g_x_0_x_yyyzz[k] = -g_0_yyyzz[k] - g_x_0_0_yyyzz[k] * ab_x + g_x_0_0_xyyyzz[k];

                g_x_0_x_yyzzz[k] = -g_0_yyzzz[k] - g_x_0_0_yyzzz[k] * ab_x + g_x_0_0_xyyzzz[k];

                g_x_0_x_yzzzz[k] = -g_0_yzzzz[k] - g_x_0_0_yzzzz[k] * ab_x + g_x_0_0_xyzzzz[k];

                g_x_0_x_zzzzz[k] = -g_0_zzzzz[k] - g_x_0_0_zzzzz[k] * ab_x + g_x_0_0_xzzzzz[k];
            }

            /// Set up 21-42 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_xxxxx = cbuffer.data(ph_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_y_xxxxy = cbuffer.data(ph_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_y_xxxxz = cbuffer.data(ph_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_y_xxxyy = cbuffer.data(ph_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_y_xxxyz = cbuffer.data(ph_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_y_xxxzz = cbuffer.data(ph_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_y_xxyyy = cbuffer.data(ph_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_y_xxyyz = cbuffer.data(ph_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_y_xxyzz = cbuffer.data(ph_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_y_xxzzz = cbuffer.data(ph_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_y_xyyyy = cbuffer.data(ph_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_y_xyyyz = cbuffer.data(ph_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_y_xyyzz = cbuffer.data(ph_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_y_xyzzz = cbuffer.data(ph_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_y_xzzzz = cbuffer.data(ph_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_y_yyyyy = cbuffer.data(ph_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_y_yyyyz = cbuffer.data(ph_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_y_yyyzz = cbuffer.data(ph_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_y_yyzzz = cbuffer.data(ph_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_y_yzzzz = cbuffer.data(ph_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_y_zzzzz = cbuffer.data(ph_geom_10_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_0_xxxxx, g_x_0_0_xxxxxy, g_x_0_0_xxxxy, g_x_0_0_xxxxyy, g_x_0_0_xxxxyz, g_x_0_0_xxxxz, g_x_0_0_xxxyy, g_x_0_0_xxxyyy, g_x_0_0_xxxyyz, g_x_0_0_xxxyz, g_x_0_0_xxxyzz, g_x_0_0_xxxzz, g_x_0_0_xxyyy, g_x_0_0_xxyyyy, g_x_0_0_xxyyyz, g_x_0_0_xxyyz, g_x_0_0_xxyyzz, g_x_0_0_xxyzz, g_x_0_0_xxyzzz, g_x_0_0_xxzzz, g_x_0_0_xyyyy, g_x_0_0_xyyyyy, g_x_0_0_xyyyyz, g_x_0_0_xyyyz, g_x_0_0_xyyyzz, g_x_0_0_xyyzz, g_x_0_0_xyyzzz, g_x_0_0_xyzzz, g_x_0_0_xyzzzz, g_x_0_0_xzzzz, g_x_0_0_yyyyy, g_x_0_0_yyyyyy, g_x_0_0_yyyyyz, g_x_0_0_yyyyz, g_x_0_0_yyyyzz, g_x_0_0_yyyzz, g_x_0_0_yyyzzz, g_x_0_0_yyzzz, g_x_0_0_yyzzzz, g_x_0_0_yzzzz, g_x_0_0_yzzzzz, g_x_0_0_zzzzz, g_x_0_y_xxxxx, g_x_0_y_xxxxy, g_x_0_y_xxxxz, g_x_0_y_xxxyy, g_x_0_y_xxxyz, g_x_0_y_xxxzz, g_x_0_y_xxyyy, g_x_0_y_xxyyz, g_x_0_y_xxyzz, g_x_0_y_xxzzz, g_x_0_y_xyyyy, g_x_0_y_xyyyz, g_x_0_y_xyyzz, g_x_0_y_xyzzz, g_x_0_y_xzzzz, g_x_0_y_yyyyy, g_x_0_y_yyyyz, g_x_0_y_yyyzz, g_x_0_y_yyzzz, g_x_0_y_yzzzz, g_x_0_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_xxxxx[k] = -g_x_0_0_xxxxx[k] * ab_y + g_x_0_0_xxxxxy[k];

                g_x_0_y_xxxxy[k] = -g_x_0_0_xxxxy[k] * ab_y + g_x_0_0_xxxxyy[k];

                g_x_0_y_xxxxz[k] = -g_x_0_0_xxxxz[k] * ab_y + g_x_0_0_xxxxyz[k];

                g_x_0_y_xxxyy[k] = -g_x_0_0_xxxyy[k] * ab_y + g_x_0_0_xxxyyy[k];

                g_x_0_y_xxxyz[k] = -g_x_0_0_xxxyz[k] * ab_y + g_x_0_0_xxxyyz[k];

                g_x_0_y_xxxzz[k] = -g_x_0_0_xxxzz[k] * ab_y + g_x_0_0_xxxyzz[k];

                g_x_0_y_xxyyy[k] = -g_x_0_0_xxyyy[k] * ab_y + g_x_0_0_xxyyyy[k];

                g_x_0_y_xxyyz[k] = -g_x_0_0_xxyyz[k] * ab_y + g_x_0_0_xxyyyz[k];

                g_x_0_y_xxyzz[k] = -g_x_0_0_xxyzz[k] * ab_y + g_x_0_0_xxyyzz[k];

                g_x_0_y_xxzzz[k] = -g_x_0_0_xxzzz[k] * ab_y + g_x_0_0_xxyzzz[k];

                g_x_0_y_xyyyy[k] = -g_x_0_0_xyyyy[k] * ab_y + g_x_0_0_xyyyyy[k];

                g_x_0_y_xyyyz[k] = -g_x_0_0_xyyyz[k] * ab_y + g_x_0_0_xyyyyz[k];

                g_x_0_y_xyyzz[k] = -g_x_0_0_xyyzz[k] * ab_y + g_x_0_0_xyyyzz[k];

                g_x_0_y_xyzzz[k] = -g_x_0_0_xyzzz[k] * ab_y + g_x_0_0_xyyzzz[k];

                g_x_0_y_xzzzz[k] = -g_x_0_0_xzzzz[k] * ab_y + g_x_0_0_xyzzzz[k];

                g_x_0_y_yyyyy[k] = -g_x_0_0_yyyyy[k] * ab_y + g_x_0_0_yyyyyy[k];

                g_x_0_y_yyyyz[k] = -g_x_0_0_yyyyz[k] * ab_y + g_x_0_0_yyyyyz[k];

                g_x_0_y_yyyzz[k] = -g_x_0_0_yyyzz[k] * ab_y + g_x_0_0_yyyyzz[k];

                g_x_0_y_yyzzz[k] = -g_x_0_0_yyzzz[k] * ab_y + g_x_0_0_yyyzzz[k];

                g_x_0_y_yzzzz[k] = -g_x_0_0_yzzzz[k] * ab_y + g_x_0_0_yyzzzz[k];

                g_x_0_y_zzzzz[k] = -g_x_0_0_zzzzz[k] * ab_y + g_x_0_0_yzzzzz[k];
            }

            /// Set up 42-63 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_xxxxx = cbuffer.data(ph_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_z_xxxxy = cbuffer.data(ph_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_z_xxxxz = cbuffer.data(ph_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_z_xxxyy = cbuffer.data(ph_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_z_xxxyz = cbuffer.data(ph_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_z_xxxzz = cbuffer.data(ph_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_z_xxyyy = cbuffer.data(ph_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_z_xxyyz = cbuffer.data(ph_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_z_xxyzz = cbuffer.data(ph_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_z_xxzzz = cbuffer.data(ph_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_z_xyyyy = cbuffer.data(ph_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_z_xyyyz = cbuffer.data(ph_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_z_xyyzz = cbuffer.data(ph_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_z_xyzzz = cbuffer.data(ph_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_z_xzzzz = cbuffer.data(ph_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_z_yyyyy = cbuffer.data(ph_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_z_yyyyz = cbuffer.data(ph_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_z_yyyzz = cbuffer.data(ph_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_z_yyzzz = cbuffer.data(ph_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_z_yzzzz = cbuffer.data(ph_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_z_zzzzz = cbuffer.data(ph_geom_10_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_0_xxxxx, g_x_0_0_xxxxxz, g_x_0_0_xxxxy, g_x_0_0_xxxxyz, g_x_0_0_xxxxz, g_x_0_0_xxxxzz, g_x_0_0_xxxyy, g_x_0_0_xxxyyz, g_x_0_0_xxxyz, g_x_0_0_xxxyzz, g_x_0_0_xxxzz, g_x_0_0_xxxzzz, g_x_0_0_xxyyy, g_x_0_0_xxyyyz, g_x_0_0_xxyyz, g_x_0_0_xxyyzz, g_x_0_0_xxyzz, g_x_0_0_xxyzzz, g_x_0_0_xxzzz, g_x_0_0_xxzzzz, g_x_0_0_xyyyy, g_x_0_0_xyyyyz, g_x_0_0_xyyyz, g_x_0_0_xyyyzz, g_x_0_0_xyyzz, g_x_0_0_xyyzzz, g_x_0_0_xyzzz, g_x_0_0_xyzzzz, g_x_0_0_xzzzz, g_x_0_0_xzzzzz, g_x_0_0_yyyyy, g_x_0_0_yyyyyz, g_x_0_0_yyyyz, g_x_0_0_yyyyzz, g_x_0_0_yyyzz, g_x_0_0_yyyzzz, g_x_0_0_yyzzz, g_x_0_0_yyzzzz, g_x_0_0_yzzzz, g_x_0_0_yzzzzz, g_x_0_0_zzzzz, g_x_0_0_zzzzzz, g_x_0_z_xxxxx, g_x_0_z_xxxxy, g_x_0_z_xxxxz, g_x_0_z_xxxyy, g_x_0_z_xxxyz, g_x_0_z_xxxzz, g_x_0_z_xxyyy, g_x_0_z_xxyyz, g_x_0_z_xxyzz, g_x_0_z_xxzzz, g_x_0_z_xyyyy, g_x_0_z_xyyyz, g_x_0_z_xyyzz, g_x_0_z_xyzzz, g_x_0_z_xzzzz, g_x_0_z_yyyyy, g_x_0_z_yyyyz, g_x_0_z_yyyzz, g_x_0_z_yyzzz, g_x_0_z_yzzzz, g_x_0_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_xxxxx[k] = -g_x_0_0_xxxxx[k] * ab_z + g_x_0_0_xxxxxz[k];

                g_x_0_z_xxxxy[k] = -g_x_0_0_xxxxy[k] * ab_z + g_x_0_0_xxxxyz[k];

                g_x_0_z_xxxxz[k] = -g_x_0_0_xxxxz[k] * ab_z + g_x_0_0_xxxxzz[k];

                g_x_0_z_xxxyy[k] = -g_x_0_0_xxxyy[k] * ab_z + g_x_0_0_xxxyyz[k];

                g_x_0_z_xxxyz[k] = -g_x_0_0_xxxyz[k] * ab_z + g_x_0_0_xxxyzz[k];

                g_x_0_z_xxxzz[k] = -g_x_0_0_xxxzz[k] * ab_z + g_x_0_0_xxxzzz[k];

                g_x_0_z_xxyyy[k] = -g_x_0_0_xxyyy[k] * ab_z + g_x_0_0_xxyyyz[k];

                g_x_0_z_xxyyz[k] = -g_x_0_0_xxyyz[k] * ab_z + g_x_0_0_xxyyzz[k];

                g_x_0_z_xxyzz[k] = -g_x_0_0_xxyzz[k] * ab_z + g_x_0_0_xxyzzz[k];

                g_x_0_z_xxzzz[k] = -g_x_0_0_xxzzz[k] * ab_z + g_x_0_0_xxzzzz[k];

                g_x_0_z_xyyyy[k] = -g_x_0_0_xyyyy[k] * ab_z + g_x_0_0_xyyyyz[k];

                g_x_0_z_xyyyz[k] = -g_x_0_0_xyyyz[k] * ab_z + g_x_0_0_xyyyzz[k];

                g_x_0_z_xyyzz[k] = -g_x_0_0_xyyzz[k] * ab_z + g_x_0_0_xyyzzz[k];

                g_x_0_z_xyzzz[k] = -g_x_0_0_xyzzz[k] * ab_z + g_x_0_0_xyzzzz[k];

                g_x_0_z_xzzzz[k] = -g_x_0_0_xzzzz[k] * ab_z + g_x_0_0_xzzzzz[k];

                g_x_0_z_yyyyy[k] = -g_x_0_0_yyyyy[k] * ab_z + g_x_0_0_yyyyyz[k];

                g_x_0_z_yyyyz[k] = -g_x_0_0_yyyyz[k] * ab_z + g_x_0_0_yyyyzz[k];

                g_x_0_z_yyyzz[k] = -g_x_0_0_yyyzz[k] * ab_z + g_x_0_0_yyyzzz[k];

                g_x_0_z_yyzzz[k] = -g_x_0_0_yyzzz[k] * ab_z + g_x_0_0_yyzzzz[k];

                g_x_0_z_yzzzz[k] = -g_x_0_0_yzzzz[k] * ab_z + g_x_0_0_yzzzzz[k];

                g_x_0_z_zzzzz[k] = -g_x_0_0_zzzzz[k] * ab_z + g_x_0_0_zzzzzz[k];
            }

            /// Set up 63-84 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_xxxxx = cbuffer.data(ph_geom_10_off + 63 * ccomps * dcomps);

            auto g_y_0_x_xxxxy = cbuffer.data(ph_geom_10_off + 64 * ccomps * dcomps);

            auto g_y_0_x_xxxxz = cbuffer.data(ph_geom_10_off + 65 * ccomps * dcomps);

            auto g_y_0_x_xxxyy = cbuffer.data(ph_geom_10_off + 66 * ccomps * dcomps);

            auto g_y_0_x_xxxyz = cbuffer.data(ph_geom_10_off + 67 * ccomps * dcomps);

            auto g_y_0_x_xxxzz = cbuffer.data(ph_geom_10_off + 68 * ccomps * dcomps);

            auto g_y_0_x_xxyyy = cbuffer.data(ph_geom_10_off + 69 * ccomps * dcomps);

            auto g_y_0_x_xxyyz = cbuffer.data(ph_geom_10_off + 70 * ccomps * dcomps);

            auto g_y_0_x_xxyzz = cbuffer.data(ph_geom_10_off + 71 * ccomps * dcomps);

            auto g_y_0_x_xxzzz = cbuffer.data(ph_geom_10_off + 72 * ccomps * dcomps);

            auto g_y_0_x_xyyyy = cbuffer.data(ph_geom_10_off + 73 * ccomps * dcomps);

            auto g_y_0_x_xyyyz = cbuffer.data(ph_geom_10_off + 74 * ccomps * dcomps);

            auto g_y_0_x_xyyzz = cbuffer.data(ph_geom_10_off + 75 * ccomps * dcomps);

            auto g_y_0_x_xyzzz = cbuffer.data(ph_geom_10_off + 76 * ccomps * dcomps);

            auto g_y_0_x_xzzzz = cbuffer.data(ph_geom_10_off + 77 * ccomps * dcomps);

            auto g_y_0_x_yyyyy = cbuffer.data(ph_geom_10_off + 78 * ccomps * dcomps);

            auto g_y_0_x_yyyyz = cbuffer.data(ph_geom_10_off + 79 * ccomps * dcomps);

            auto g_y_0_x_yyyzz = cbuffer.data(ph_geom_10_off + 80 * ccomps * dcomps);

            auto g_y_0_x_yyzzz = cbuffer.data(ph_geom_10_off + 81 * ccomps * dcomps);

            auto g_y_0_x_yzzzz = cbuffer.data(ph_geom_10_off + 82 * ccomps * dcomps);

            auto g_y_0_x_zzzzz = cbuffer.data(ph_geom_10_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_0_xxxxx, g_y_0_0_xxxxxx, g_y_0_0_xxxxxy, g_y_0_0_xxxxxz, g_y_0_0_xxxxy, g_y_0_0_xxxxyy, g_y_0_0_xxxxyz, g_y_0_0_xxxxz, g_y_0_0_xxxxzz, g_y_0_0_xxxyy, g_y_0_0_xxxyyy, g_y_0_0_xxxyyz, g_y_0_0_xxxyz, g_y_0_0_xxxyzz, g_y_0_0_xxxzz, g_y_0_0_xxxzzz, g_y_0_0_xxyyy, g_y_0_0_xxyyyy, g_y_0_0_xxyyyz, g_y_0_0_xxyyz, g_y_0_0_xxyyzz, g_y_0_0_xxyzz, g_y_0_0_xxyzzz, g_y_0_0_xxzzz, g_y_0_0_xxzzzz, g_y_0_0_xyyyy, g_y_0_0_xyyyyy, g_y_0_0_xyyyyz, g_y_0_0_xyyyz, g_y_0_0_xyyyzz, g_y_0_0_xyyzz, g_y_0_0_xyyzzz, g_y_0_0_xyzzz, g_y_0_0_xyzzzz, g_y_0_0_xzzzz, g_y_0_0_xzzzzz, g_y_0_0_yyyyy, g_y_0_0_yyyyz, g_y_0_0_yyyzz, g_y_0_0_yyzzz, g_y_0_0_yzzzz, g_y_0_0_zzzzz, g_y_0_x_xxxxx, g_y_0_x_xxxxy, g_y_0_x_xxxxz, g_y_0_x_xxxyy, g_y_0_x_xxxyz, g_y_0_x_xxxzz, g_y_0_x_xxyyy, g_y_0_x_xxyyz, g_y_0_x_xxyzz, g_y_0_x_xxzzz, g_y_0_x_xyyyy, g_y_0_x_xyyyz, g_y_0_x_xyyzz, g_y_0_x_xyzzz, g_y_0_x_xzzzz, g_y_0_x_yyyyy, g_y_0_x_yyyyz, g_y_0_x_yyyzz, g_y_0_x_yyzzz, g_y_0_x_yzzzz, g_y_0_x_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_xxxxx[k] = -g_y_0_0_xxxxx[k] * ab_x + g_y_0_0_xxxxxx[k];

                g_y_0_x_xxxxy[k] = -g_y_0_0_xxxxy[k] * ab_x + g_y_0_0_xxxxxy[k];

                g_y_0_x_xxxxz[k] = -g_y_0_0_xxxxz[k] * ab_x + g_y_0_0_xxxxxz[k];

                g_y_0_x_xxxyy[k] = -g_y_0_0_xxxyy[k] * ab_x + g_y_0_0_xxxxyy[k];

                g_y_0_x_xxxyz[k] = -g_y_0_0_xxxyz[k] * ab_x + g_y_0_0_xxxxyz[k];

                g_y_0_x_xxxzz[k] = -g_y_0_0_xxxzz[k] * ab_x + g_y_0_0_xxxxzz[k];

                g_y_0_x_xxyyy[k] = -g_y_0_0_xxyyy[k] * ab_x + g_y_0_0_xxxyyy[k];

                g_y_0_x_xxyyz[k] = -g_y_0_0_xxyyz[k] * ab_x + g_y_0_0_xxxyyz[k];

                g_y_0_x_xxyzz[k] = -g_y_0_0_xxyzz[k] * ab_x + g_y_0_0_xxxyzz[k];

                g_y_0_x_xxzzz[k] = -g_y_0_0_xxzzz[k] * ab_x + g_y_0_0_xxxzzz[k];

                g_y_0_x_xyyyy[k] = -g_y_0_0_xyyyy[k] * ab_x + g_y_0_0_xxyyyy[k];

                g_y_0_x_xyyyz[k] = -g_y_0_0_xyyyz[k] * ab_x + g_y_0_0_xxyyyz[k];

                g_y_0_x_xyyzz[k] = -g_y_0_0_xyyzz[k] * ab_x + g_y_0_0_xxyyzz[k];

                g_y_0_x_xyzzz[k] = -g_y_0_0_xyzzz[k] * ab_x + g_y_0_0_xxyzzz[k];

                g_y_0_x_xzzzz[k] = -g_y_0_0_xzzzz[k] * ab_x + g_y_0_0_xxzzzz[k];

                g_y_0_x_yyyyy[k] = -g_y_0_0_yyyyy[k] * ab_x + g_y_0_0_xyyyyy[k];

                g_y_0_x_yyyyz[k] = -g_y_0_0_yyyyz[k] * ab_x + g_y_0_0_xyyyyz[k];

                g_y_0_x_yyyzz[k] = -g_y_0_0_yyyzz[k] * ab_x + g_y_0_0_xyyyzz[k];

                g_y_0_x_yyzzz[k] = -g_y_0_0_yyzzz[k] * ab_x + g_y_0_0_xyyzzz[k];

                g_y_0_x_yzzzz[k] = -g_y_0_0_yzzzz[k] * ab_x + g_y_0_0_xyzzzz[k];

                g_y_0_x_zzzzz[k] = -g_y_0_0_zzzzz[k] * ab_x + g_y_0_0_xzzzzz[k];
            }

            /// Set up 84-105 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_xxxxx = cbuffer.data(ph_geom_10_off + 84 * ccomps * dcomps);

            auto g_y_0_y_xxxxy = cbuffer.data(ph_geom_10_off + 85 * ccomps * dcomps);

            auto g_y_0_y_xxxxz = cbuffer.data(ph_geom_10_off + 86 * ccomps * dcomps);

            auto g_y_0_y_xxxyy = cbuffer.data(ph_geom_10_off + 87 * ccomps * dcomps);

            auto g_y_0_y_xxxyz = cbuffer.data(ph_geom_10_off + 88 * ccomps * dcomps);

            auto g_y_0_y_xxxzz = cbuffer.data(ph_geom_10_off + 89 * ccomps * dcomps);

            auto g_y_0_y_xxyyy = cbuffer.data(ph_geom_10_off + 90 * ccomps * dcomps);

            auto g_y_0_y_xxyyz = cbuffer.data(ph_geom_10_off + 91 * ccomps * dcomps);

            auto g_y_0_y_xxyzz = cbuffer.data(ph_geom_10_off + 92 * ccomps * dcomps);

            auto g_y_0_y_xxzzz = cbuffer.data(ph_geom_10_off + 93 * ccomps * dcomps);

            auto g_y_0_y_xyyyy = cbuffer.data(ph_geom_10_off + 94 * ccomps * dcomps);

            auto g_y_0_y_xyyyz = cbuffer.data(ph_geom_10_off + 95 * ccomps * dcomps);

            auto g_y_0_y_xyyzz = cbuffer.data(ph_geom_10_off + 96 * ccomps * dcomps);

            auto g_y_0_y_xyzzz = cbuffer.data(ph_geom_10_off + 97 * ccomps * dcomps);

            auto g_y_0_y_xzzzz = cbuffer.data(ph_geom_10_off + 98 * ccomps * dcomps);

            auto g_y_0_y_yyyyy = cbuffer.data(ph_geom_10_off + 99 * ccomps * dcomps);

            auto g_y_0_y_yyyyz = cbuffer.data(ph_geom_10_off + 100 * ccomps * dcomps);

            auto g_y_0_y_yyyzz = cbuffer.data(ph_geom_10_off + 101 * ccomps * dcomps);

            auto g_y_0_y_yyzzz = cbuffer.data(ph_geom_10_off + 102 * ccomps * dcomps);

            auto g_y_0_y_yzzzz = cbuffer.data(ph_geom_10_off + 103 * ccomps * dcomps);

            auto g_y_0_y_zzzzz = cbuffer.data(ph_geom_10_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_xxxxx, g_0_xxxxy, g_0_xxxxz, g_0_xxxyy, g_0_xxxyz, g_0_xxxzz, g_0_xxyyy, g_0_xxyyz, g_0_xxyzz, g_0_xxzzz, g_0_xyyyy, g_0_xyyyz, g_0_xyyzz, g_0_xyzzz, g_0_xzzzz, g_0_yyyyy, g_0_yyyyz, g_0_yyyzz, g_0_yyzzz, g_0_yzzzz, g_0_zzzzz, g_y_0_0_xxxxx, g_y_0_0_xxxxxy, g_y_0_0_xxxxy, g_y_0_0_xxxxyy, g_y_0_0_xxxxyz, g_y_0_0_xxxxz, g_y_0_0_xxxyy, g_y_0_0_xxxyyy, g_y_0_0_xxxyyz, g_y_0_0_xxxyz, g_y_0_0_xxxyzz, g_y_0_0_xxxzz, g_y_0_0_xxyyy, g_y_0_0_xxyyyy, g_y_0_0_xxyyyz, g_y_0_0_xxyyz, g_y_0_0_xxyyzz, g_y_0_0_xxyzz, g_y_0_0_xxyzzz, g_y_0_0_xxzzz, g_y_0_0_xyyyy, g_y_0_0_xyyyyy, g_y_0_0_xyyyyz, g_y_0_0_xyyyz, g_y_0_0_xyyyzz, g_y_0_0_xyyzz, g_y_0_0_xyyzzz, g_y_0_0_xyzzz, g_y_0_0_xyzzzz, g_y_0_0_xzzzz, g_y_0_0_yyyyy, g_y_0_0_yyyyyy, g_y_0_0_yyyyyz, g_y_0_0_yyyyz, g_y_0_0_yyyyzz, g_y_0_0_yyyzz, g_y_0_0_yyyzzz, g_y_0_0_yyzzz, g_y_0_0_yyzzzz, g_y_0_0_yzzzz, g_y_0_0_yzzzzz, g_y_0_0_zzzzz, g_y_0_y_xxxxx, g_y_0_y_xxxxy, g_y_0_y_xxxxz, g_y_0_y_xxxyy, g_y_0_y_xxxyz, g_y_0_y_xxxzz, g_y_0_y_xxyyy, g_y_0_y_xxyyz, g_y_0_y_xxyzz, g_y_0_y_xxzzz, g_y_0_y_xyyyy, g_y_0_y_xyyyz, g_y_0_y_xyyzz, g_y_0_y_xyzzz, g_y_0_y_xzzzz, g_y_0_y_yyyyy, g_y_0_y_yyyyz, g_y_0_y_yyyzz, g_y_0_y_yyzzz, g_y_0_y_yzzzz, g_y_0_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_xxxxx[k] = -g_0_xxxxx[k] - g_y_0_0_xxxxx[k] * ab_y + g_y_0_0_xxxxxy[k];

                g_y_0_y_xxxxy[k] = -g_0_xxxxy[k] - g_y_0_0_xxxxy[k] * ab_y + g_y_0_0_xxxxyy[k];

                g_y_0_y_xxxxz[k] = -g_0_xxxxz[k] - g_y_0_0_xxxxz[k] * ab_y + g_y_0_0_xxxxyz[k];

                g_y_0_y_xxxyy[k] = -g_0_xxxyy[k] - g_y_0_0_xxxyy[k] * ab_y + g_y_0_0_xxxyyy[k];

                g_y_0_y_xxxyz[k] = -g_0_xxxyz[k] - g_y_0_0_xxxyz[k] * ab_y + g_y_0_0_xxxyyz[k];

                g_y_0_y_xxxzz[k] = -g_0_xxxzz[k] - g_y_0_0_xxxzz[k] * ab_y + g_y_0_0_xxxyzz[k];

                g_y_0_y_xxyyy[k] = -g_0_xxyyy[k] - g_y_0_0_xxyyy[k] * ab_y + g_y_0_0_xxyyyy[k];

                g_y_0_y_xxyyz[k] = -g_0_xxyyz[k] - g_y_0_0_xxyyz[k] * ab_y + g_y_0_0_xxyyyz[k];

                g_y_0_y_xxyzz[k] = -g_0_xxyzz[k] - g_y_0_0_xxyzz[k] * ab_y + g_y_0_0_xxyyzz[k];

                g_y_0_y_xxzzz[k] = -g_0_xxzzz[k] - g_y_0_0_xxzzz[k] * ab_y + g_y_0_0_xxyzzz[k];

                g_y_0_y_xyyyy[k] = -g_0_xyyyy[k] - g_y_0_0_xyyyy[k] * ab_y + g_y_0_0_xyyyyy[k];

                g_y_0_y_xyyyz[k] = -g_0_xyyyz[k] - g_y_0_0_xyyyz[k] * ab_y + g_y_0_0_xyyyyz[k];

                g_y_0_y_xyyzz[k] = -g_0_xyyzz[k] - g_y_0_0_xyyzz[k] * ab_y + g_y_0_0_xyyyzz[k];

                g_y_0_y_xyzzz[k] = -g_0_xyzzz[k] - g_y_0_0_xyzzz[k] * ab_y + g_y_0_0_xyyzzz[k];

                g_y_0_y_xzzzz[k] = -g_0_xzzzz[k] - g_y_0_0_xzzzz[k] * ab_y + g_y_0_0_xyzzzz[k];

                g_y_0_y_yyyyy[k] = -g_0_yyyyy[k] - g_y_0_0_yyyyy[k] * ab_y + g_y_0_0_yyyyyy[k];

                g_y_0_y_yyyyz[k] = -g_0_yyyyz[k] - g_y_0_0_yyyyz[k] * ab_y + g_y_0_0_yyyyyz[k];

                g_y_0_y_yyyzz[k] = -g_0_yyyzz[k] - g_y_0_0_yyyzz[k] * ab_y + g_y_0_0_yyyyzz[k];

                g_y_0_y_yyzzz[k] = -g_0_yyzzz[k] - g_y_0_0_yyzzz[k] * ab_y + g_y_0_0_yyyzzz[k];

                g_y_0_y_yzzzz[k] = -g_0_yzzzz[k] - g_y_0_0_yzzzz[k] * ab_y + g_y_0_0_yyzzzz[k];

                g_y_0_y_zzzzz[k] = -g_0_zzzzz[k] - g_y_0_0_zzzzz[k] * ab_y + g_y_0_0_yzzzzz[k];
            }

            /// Set up 105-126 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_xxxxx = cbuffer.data(ph_geom_10_off + 105 * ccomps * dcomps);

            auto g_y_0_z_xxxxy = cbuffer.data(ph_geom_10_off + 106 * ccomps * dcomps);

            auto g_y_0_z_xxxxz = cbuffer.data(ph_geom_10_off + 107 * ccomps * dcomps);

            auto g_y_0_z_xxxyy = cbuffer.data(ph_geom_10_off + 108 * ccomps * dcomps);

            auto g_y_0_z_xxxyz = cbuffer.data(ph_geom_10_off + 109 * ccomps * dcomps);

            auto g_y_0_z_xxxzz = cbuffer.data(ph_geom_10_off + 110 * ccomps * dcomps);

            auto g_y_0_z_xxyyy = cbuffer.data(ph_geom_10_off + 111 * ccomps * dcomps);

            auto g_y_0_z_xxyyz = cbuffer.data(ph_geom_10_off + 112 * ccomps * dcomps);

            auto g_y_0_z_xxyzz = cbuffer.data(ph_geom_10_off + 113 * ccomps * dcomps);

            auto g_y_0_z_xxzzz = cbuffer.data(ph_geom_10_off + 114 * ccomps * dcomps);

            auto g_y_0_z_xyyyy = cbuffer.data(ph_geom_10_off + 115 * ccomps * dcomps);

            auto g_y_0_z_xyyyz = cbuffer.data(ph_geom_10_off + 116 * ccomps * dcomps);

            auto g_y_0_z_xyyzz = cbuffer.data(ph_geom_10_off + 117 * ccomps * dcomps);

            auto g_y_0_z_xyzzz = cbuffer.data(ph_geom_10_off + 118 * ccomps * dcomps);

            auto g_y_0_z_xzzzz = cbuffer.data(ph_geom_10_off + 119 * ccomps * dcomps);

            auto g_y_0_z_yyyyy = cbuffer.data(ph_geom_10_off + 120 * ccomps * dcomps);

            auto g_y_0_z_yyyyz = cbuffer.data(ph_geom_10_off + 121 * ccomps * dcomps);

            auto g_y_0_z_yyyzz = cbuffer.data(ph_geom_10_off + 122 * ccomps * dcomps);

            auto g_y_0_z_yyzzz = cbuffer.data(ph_geom_10_off + 123 * ccomps * dcomps);

            auto g_y_0_z_yzzzz = cbuffer.data(ph_geom_10_off + 124 * ccomps * dcomps);

            auto g_y_0_z_zzzzz = cbuffer.data(ph_geom_10_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_0_xxxxx, g_y_0_0_xxxxxz, g_y_0_0_xxxxy, g_y_0_0_xxxxyz, g_y_0_0_xxxxz, g_y_0_0_xxxxzz, g_y_0_0_xxxyy, g_y_0_0_xxxyyz, g_y_0_0_xxxyz, g_y_0_0_xxxyzz, g_y_0_0_xxxzz, g_y_0_0_xxxzzz, g_y_0_0_xxyyy, g_y_0_0_xxyyyz, g_y_0_0_xxyyz, g_y_0_0_xxyyzz, g_y_0_0_xxyzz, g_y_0_0_xxyzzz, g_y_0_0_xxzzz, g_y_0_0_xxzzzz, g_y_0_0_xyyyy, g_y_0_0_xyyyyz, g_y_0_0_xyyyz, g_y_0_0_xyyyzz, g_y_0_0_xyyzz, g_y_0_0_xyyzzz, g_y_0_0_xyzzz, g_y_0_0_xyzzzz, g_y_0_0_xzzzz, g_y_0_0_xzzzzz, g_y_0_0_yyyyy, g_y_0_0_yyyyyz, g_y_0_0_yyyyz, g_y_0_0_yyyyzz, g_y_0_0_yyyzz, g_y_0_0_yyyzzz, g_y_0_0_yyzzz, g_y_0_0_yyzzzz, g_y_0_0_yzzzz, g_y_0_0_yzzzzz, g_y_0_0_zzzzz, g_y_0_0_zzzzzz, g_y_0_z_xxxxx, g_y_0_z_xxxxy, g_y_0_z_xxxxz, g_y_0_z_xxxyy, g_y_0_z_xxxyz, g_y_0_z_xxxzz, g_y_0_z_xxyyy, g_y_0_z_xxyyz, g_y_0_z_xxyzz, g_y_0_z_xxzzz, g_y_0_z_xyyyy, g_y_0_z_xyyyz, g_y_0_z_xyyzz, g_y_0_z_xyzzz, g_y_0_z_xzzzz, g_y_0_z_yyyyy, g_y_0_z_yyyyz, g_y_0_z_yyyzz, g_y_0_z_yyzzz, g_y_0_z_yzzzz, g_y_0_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_xxxxx[k] = -g_y_0_0_xxxxx[k] * ab_z + g_y_0_0_xxxxxz[k];

                g_y_0_z_xxxxy[k] = -g_y_0_0_xxxxy[k] * ab_z + g_y_0_0_xxxxyz[k];

                g_y_0_z_xxxxz[k] = -g_y_0_0_xxxxz[k] * ab_z + g_y_0_0_xxxxzz[k];

                g_y_0_z_xxxyy[k] = -g_y_0_0_xxxyy[k] * ab_z + g_y_0_0_xxxyyz[k];

                g_y_0_z_xxxyz[k] = -g_y_0_0_xxxyz[k] * ab_z + g_y_0_0_xxxyzz[k];

                g_y_0_z_xxxzz[k] = -g_y_0_0_xxxzz[k] * ab_z + g_y_0_0_xxxzzz[k];

                g_y_0_z_xxyyy[k] = -g_y_0_0_xxyyy[k] * ab_z + g_y_0_0_xxyyyz[k];

                g_y_0_z_xxyyz[k] = -g_y_0_0_xxyyz[k] * ab_z + g_y_0_0_xxyyzz[k];

                g_y_0_z_xxyzz[k] = -g_y_0_0_xxyzz[k] * ab_z + g_y_0_0_xxyzzz[k];

                g_y_0_z_xxzzz[k] = -g_y_0_0_xxzzz[k] * ab_z + g_y_0_0_xxzzzz[k];

                g_y_0_z_xyyyy[k] = -g_y_0_0_xyyyy[k] * ab_z + g_y_0_0_xyyyyz[k];

                g_y_0_z_xyyyz[k] = -g_y_0_0_xyyyz[k] * ab_z + g_y_0_0_xyyyzz[k];

                g_y_0_z_xyyzz[k] = -g_y_0_0_xyyzz[k] * ab_z + g_y_0_0_xyyzzz[k];

                g_y_0_z_xyzzz[k] = -g_y_0_0_xyzzz[k] * ab_z + g_y_0_0_xyzzzz[k];

                g_y_0_z_xzzzz[k] = -g_y_0_0_xzzzz[k] * ab_z + g_y_0_0_xzzzzz[k];

                g_y_0_z_yyyyy[k] = -g_y_0_0_yyyyy[k] * ab_z + g_y_0_0_yyyyyz[k];

                g_y_0_z_yyyyz[k] = -g_y_0_0_yyyyz[k] * ab_z + g_y_0_0_yyyyzz[k];

                g_y_0_z_yyyzz[k] = -g_y_0_0_yyyzz[k] * ab_z + g_y_0_0_yyyzzz[k];

                g_y_0_z_yyzzz[k] = -g_y_0_0_yyzzz[k] * ab_z + g_y_0_0_yyzzzz[k];

                g_y_0_z_yzzzz[k] = -g_y_0_0_yzzzz[k] * ab_z + g_y_0_0_yzzzzz[k];

                g_y_0_z_zzzzz[k] = -g_y_0_0_zzzzz[k] * ab_z + g_y_0_0_zzzzzz[k];
            }

            /// Set up 126-147 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_xxxxx = cbuffer.data(ph_geom_10_off + 126 * ccomps * dcomps);

            auto g_z_0_x_xxxxy = cbuffer.data(ph_geom_10_off + 127 * ccomps * dcomps);

            auto g_z_0_x_xxxxz = cbuffer.data(ph_geom_10_off + 128 * ccomps * dcomps);

            auto g_z_0_x_xxxyy = cbuffer.data(ph_geom_10_off + 129 * ccomps * dcomps);

            auto g_z_0_x_xxxyz = cbuffer.data(ph_geom_10_off + 130 * ccomps * dcomps);

            auto g_z_0_x_xxxzz = cbuffer.data(ph_geom_10_off + 131 * ccomps * dcomps);

            auto g_z_0_x_xxyyy = cbuffer.data(ph_geom_10_off + 132 * ccomps * dcomps);

            auto g_z_0_x_xxyyz = cbuffer.data(ph_geom_10_off + 133 * ccomps * dcomps);

            auto g_z_0_x_xxyzz = cbuffer.data(ph_geom_10_off + 134 * ccomps * dcomps);

            auto g_z_0_x_xxzzz = cbuffer.data(ph_geom_10_off + 135 * ccomps * dcomps);

            auto g_z_0_x_xyyyy = cbuffer.data(ph_geom_10_off + 136 * ccomps * dcomps);

            auto g_z_0_x_xyyyz = cbuffer.data(ph_geom_10_off + 137 * ccomps * dcomps);

            auto g_z_0_x_xyyzz = cbuffer.data(ph_geom_10_off + 138 * ccomps * dcomps);

            auto g_z_0_x_xyzzz = cbuffer.data(ph_geom_10_off + 139 * ccomps * dcomps);

            auto g_z_0_x_xzzzz = cbuffer.data(ph_geom_10_off + 140 * ccomps * dcomps);

            auto g_z_0_x_yyyyy = cbuffer.data(ph_geom_10_off + 141 * ccomps * dcomps);

            auto g_z_0_x_yyyyz = cbuffer.data(ph_geom_10_off + 142 * ccomps * dcomps);

            auto g_z_0_x_yyyzz = cbuffer.data(ph_geom_10_off + 143 * ccomps * dcomps);

            auto g_z_0_x_yyzzz = cbuffer.data(ph_geom_10_off + 144 * ccomps * dcomps);

            auto g_z_0_x_yzzzz = cbuffer.data(ph_geom_10_off + 145 * ccomps * dcomps);

            auto g_z_0_x_zzzzz = cbuffer.data(ph_geom_10_off + 146 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_0_xxxxx, g_z_0_0_xxxxxx, g_z_0_0_xxxxxy, g_z_0_0_xxxxxz, g_z_0_0_xxxxy, g_z_0_0_xxxxyy, g_z_0_0_xxxxyz, g_z_0_0_xxxxz, g_z_0_0_xxxxzz, g_z_0_0_xxxyy, g_z_0_0_xxxyyy, g_z_0_0_xxxyyz, g_z_0_0_xxxyz, g_z_0_0_xxxyzz, g_z_0_0_xxxzz, g_z_0_0_xxxzzz, g_z_0_0_xxyyy, g_z_0_0_xxyyyy, g_z_0_0_xxyyyz, g_z_0_0_xxyyz, g_z_0_0_xxyyzz, g_z_0_0_xxyzz, g_z_0_0_xxyzzz, g_z_0_0_xxzzz, g_z_0_0_xxzzzz, g_z_0_0_xyyyy, g_z_0_0_xyyyyy, g_z_0_0_xyyyyz, g_z_0_0_xyyyz, g_z_0_0_xyyyzz, g_z_0_0_xyyzz, g_z_0_0_xyyzzz, g_z_0_0_xyzzz, g_z_0_0_xyzzzz, g_z_0_0_xzzzz, g_z_0_0_xzzzzz, g_z_0_0_yyyyy, g_z_0_0_yyyyz, g_z_0_0_yyyzz, g_z_0_0_yyzzz, g_z_0_0_yzzzz, g_z_0_0_zzzzz, g_z_0_x_xxxxx, g_z_0_x_xxxxy, g_z_0_x_xxxxz, g_z_0_x_xxxyy, g_z_0_x_xxxyz, g_z_0_x_xxxzz, g_z_0_x_xxyyy, g_z_0_x_xxyyz, g_z_0_x_xxyzz, g_z_0_x_xxzzz, g_z_0_x_xyyyy, g_z_0_x_xyyyz, g_z_0_x_xyyzz, g_z_0_x_xyzzz, g_z_0_x_xzzzz, g_z_0_x_yyyyy, g_z_0_x_yyyyz, g_z_0_x_yyyzz, g_z_0_x_yyzzz, g_z_0_x_yzzzz, g_z_0_x_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_xxxxx[k] = -g_z_0_0_xxxxx[k] * ab_x + g_z_0_0_xxxxxx[k];

                g_z_0_x_xxxxy[k] = -g_z_0_0_xxxxy[k] * ab_x + g_z_0_0_xxxxxy[k];

                g_z_0_x_xxxxz[k] = -g_z_0_0_xxxxz[k] * ab_x + g_z_0_0_xxxxxz[k];

                g_z_0_x_xxxyy[k] = -g_z_0_0_xxxyy[k] * ab_x + g_z_0_0_xxxxyy[k];

                g_z_0_x_xxxyz[k] = -g_z_0_0_xxxyz[k] * ab_x + g_z_0_0_xxxxyz[k];

                g_z_0_x_xxxzz[k] = -g_z_0_0_xxxzz[k] * ab_x + g_z_0_0_xxxxzz[k];

                g_z_0_x_xxyyy[k] = -g_z_0_0_xxyyy[k] * ab_x + g_z_0_0_xxxyyy[k];

                g_z_0_x_xxyyz[k] = -g_z_0_0_xxyyz[k] * ab_x + g_z_0_0_xxxyyz[k];

                g_z_0_x_xxyzz[k] = -g_z_0_0_xxyzz[k] * ab_x + g_z_0_0_xxxyzz[k];

                g_z_0_x_xxzzz[k] = -g_z_0_0_xxzzz[k] * ab_x + g_z_0_0_xxxzzz[k];

                g_z_0_x_xyyyy[k] = -g_z_0_0_xyyyy[k] * ab_x + g_z_0_0_xxyyyy[k];

                g_z_0_x_xyyyz[k] = -g_z_0_0_xyyyz[k] * ab_x + g_z_0_0_xxyyyz[k];

                g_z_0_x_xyyzz[k] = -g_z_0_0_xyyzz[k] * ab_x + g_z_0_0_xxyyzz[k];

                g_z_0_x_xyzzz[k] = -g_z_0_0_xyzzz[k] * ab_x + g_z_0_0_xxyzzz[k];

                g_z_0_x_xzzzz[k] = -g_z_0_0_xzzzz[k] * ab_x + g_z_0_0_xxzzzz[k];

                g_z_0_x_yyyyy[k] = -g_z_0_0_yyyyy[k] * ab_x + g_z_0_0_xyyyyy[k];

                g_z_0_x_yyyyz[k] = -g_z_0_0_yyyyz[k] * ab_x + g_z_0_0_xyyyyz[k];

                g_z_0_x_yyyzz[k] = -g_z_0_0_yyyzz[k] * ab_x + g_z_0_0_xyyyzz[k];

                g_z_0_x_yyzzz[k] = -g_z_0_0_yyzzz[k] * ab_x + g_z_0_0_xyyzzz[k];

                g_z_0_x_yzzzz[k] = -g_z_0_0_yzzzz[k] * ab_x + g_z_0_0_xyzzzz[k];

                g_z_0_x_zzzzz[k] = -g_z_0_0_zzzzz[k] * ab_x + g_z_0_0_xzzzzz[k];
            }

            /// Set up 147-168 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_xxxxx = cbuffer.data(ph_geom_10_off + 147 * ccomps * dcomps);

            auto g_z_0_y_xxxxy = cbuffer.data(ph_geom_10_off + 148 * ccomps * dcomps);

            auto g_z_0_y_xxxxz = cbuffer.data(ph_geom_10_off + 149 * ccomps * dcomps);

            auto g_z_0_y_xxxyy = cbuffer.data(ph_geom_10_off + 150 * ccomps * dcomps);

            auto g_z_0_y_xxxyz = cbuffer.data(ph_geom_10_off + 151 * ccomps * dcomps);

            auto g_z_0_y_xxxzz = cbuffer.data(ph_geom_10_off + 152 * ccomps * dcomps);

            auto g_z_0_y_xxyyy = cbuffer.data(ph_geom_10_off + 153 * ccomps * dcomps);

            auto g_z_0_y_xxyyz = cbuffer.data(ph_geom_10_off + 154 * ccomps * dcomps);

            auto g_z_0_y_xxyzz = cbuffer.data(ph_geom_10_off + 155 * ccomps * dcomps);

            auto g_z_0_y_xxzzz = cbuffer.data(ph_geom_10_off + 156 * ccomps * dcomps);

            auto g_z_0_y_xyyyy = cbuffer.data(ph_geom_10_off + 157 * ccomps * dcomps);

            auto g_z_0_y_xyyyz = cbuffer.data(ph_geom_10_off + 158 * ccomps * dcomps);

            auto g_z_0_y_xyyzz = cbuffer.data(ph_geom_10_off + 159 * ccomps * dcomps);

            auto g_z_0_y_xyzzz = cbuffer.data(ph_geom_10_off + 160 * ccomps * dcomps);

            auto g_z_0_y_xzzzz = cbuffer.data(ph_geom_10_off + 161 * ccomps * dcomps);

            auto g_z_0_y_yyyyy = cbuffer.data(ph_geom_10_off + 162 * ccomps * dcomps);

            auto g_z_0_y_yyyyz = cbuffer.data(ph_geom_10_off + 163 * ccomps * dcomps);

            auto g_z_0_y_yyyzz = cbuffer.data(ph_geom_10_off + 164 * ccomps * dcomps);

            auto g_z_0_y_yyzzz = cbuffer.data(ph_geom_10_off + 165 * ccomps * dcomps);

            auto g_z_0_y_yzzzz = cbuffer.data(ph_geom_10_off + 166 * ccomps * dcomps);

            auto g_z_0_y_zzzzz = cbuffer.data(ph_geom_10_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_0_xxxxx, g_z_0_0_xxxxxy, g_z_0_0_xxxxy, g_z_0_0_xxxxyy, g_z_0_0_xxxxyz, g_z_0_0_xxxxz, g_z_0_0_xxxyy, g_z_0_0_xxxyyy, g_z_0_0_xxxyyz, g_z_0_0_xxxyz, g_z_0_0_xxxyzz, g_z_0_0_xxxzz, g_z_0_0_xxyyy, g_z_0_0_xxyyyy, g_z_0_0_xxyyyz, g_z_0_0_xxyyz, g_z_0_0_xxyyzz, g_z_0_0_xxyzz, g_z_0_0_xxyzzz, g_z_0_0_xxzzz, g_z_0_0_xyyyy, g_z_0_0_xyyyyy, g_z_0_0_xyyyyz, g_z_0_0_xyyyz, g_z_0_0_xyyyzz, g_z_0_0_xyyzz, g_z_0_0_xyyzzz, g_z_0_0_xyzzz, g_z_0_0_xyzzzz, g_z_0_0_xzzzz, g_z_0_0_yyyyy, g_z_0_0_yyyyyy, g_z_0_0_yyyyyz, g_z_0_0_yyyyz, g_z_0_0_yyyyzz, g_z_0_0_yyyzz, g_z_0_0_yyyzzz, g_z_0_0_yyzzz, g_z_0_0_yyzzzz, g_z_0_0_yzzzz, g_z_0_0_yzzzzz, g_z_0_0_zzzzz, g_z_0_y_xxxxx, g_z_0_y_xxxxy, g_z_0_y_xxxxz, g_z_0_y_xxxyy, g_z_0_y_xxxyz, g_z_0_y_xxxzz, g_z_0_y_xxyyy, g_z_0_y_xxyyz, g_z_0_y_xxyzz, g_z_0_y_xxzzz, g_z_0_y_xyyyy, g_z_0_y_xyyyz, g_z_0_y_xyyzz, g_z_0_y_xyzzz, g_z_0_y_xzzzz, g_z_0_y_yyyyy, g_z_0_y_yyyyz, g_z_0_y_yyyzz, g_z_0_y_yyzzz, g_z_0_y_yzzzz, g_z_0_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_xxxxx[k] = -g_z_0_0_xxxxx[k] * ab_y + g_z_0_0_xxxxxy[k];

                g_z_0_y_xxxxy[k] = -g_z_0_0_xxxxy[k] * ab_y + g_z_0_0_xxxxyy[k];

                g_z_0_y_xxxxz[k] = -g_z_0_0_xxxxz[k] * ab_y + g_z_0_0_xxxxyz[k];

                g_z_0_y_xxxyy[k] = -g_z_0_0_xxxyy[k] * ab_y + g_z_0_0_xxxyyy[k];

                g_z_0_y_xxxyz[k] = -g_z_0_0_xxxyz[k] * ab_y + g_z_0_0_xxxyyz[k];

                g_z_0_y_xxxzz[k] = -g_z_0_0_xxxzz[k] * ab_y + g_z_0_0_xxxyzz[k];

                g_z_0_y_xxyyy[k] = -g_z_0_0_xxyyy[k] * ab_y + g_z_0_0_xxyyyy[k];

                g_z_0_y_xxyyz[k] = -g_z_0_0_xxyyz[k] * ab_y + g_z_0_0_xxyyyz[k];

                g_z_0_y_xxyzz[k] = -g_z_0_0_xxyzz[k] * ab_y + g_z_0_0_xxyyzz[k];

                g_z_0_y_xxzzz[k] = -g_z_0_0_xxzzz[k] * ab_y + g_z_0_0_xxyzzz[k];

                g_z_0_y_xyyyy[k] = -g_z_0_0_xyyyy[k] * ab_y + g_z_0_0_xyyyyy[k];

                g_z_0_y_xyyyz[k] = -g_z_0_0_xyyyz[k] * ab_y + g_z_0_0_xyyyyz[k];

                g_z_0_y_xyyzz[k] = -g_z_0_0_xyyzz[k] * ab_y + g_z_0_0_xyyyzz[k];

                g_z_0_y_xyzzz[k] = -g_z_0_0_xyzzz[k] * ab_y + g_z_0_0_xyyzzz[k];

                g_z_0_y_xzzzz[k] = -g_z_0_0_xzzzz[k] * ab_y + g_z_0_0_xyzzzz[k];

                g_z_0_y_yyyyy[k] = -g_z_0_0_yyyyy[k] * ab_y + g_z_0_0_yyyyyy[k];

                g_z_0_y_yyyyz[k] = -g_z_0_0_yyyyz[k] * ab_y + g_z_0_0_yyyyyz[k];

                g_z_0_y_yyyzz[k] = -g_z_0_0_yyyzz[k] * ab_y + g_z_0_0_yyyyzz[k];

                g_z_0_y_yyzzz[k] = -g_z_0_0_yyzzz[k] * ab_y + g_z_0_0_yyyzzz[k];

                g_z_0_y_yzzzz[k] = -g_z_0_0_yzzzz[k] * ab_y + g_z_0_0_yyzzzz[k];

                g_z_0_y_zzzzz[k] = -g_z_0_0_zzzzz[k] * ab_y + g_z_0_0_yzzzzz[k];
            }

            /// Set up 168-189 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_xxxxx = cbuffer.data(ph_geom_10_off + 168 * ccomps * dcomps);

            auto g_z_0_z_xxxxy = cbuffer.data(ph_geom_10_off + 169 * ccomps * dcomps);

            auto g_z_0_z_xxxxz = cbuffer.data(ph_geom_10_off + 170 * ccomps * dcomps);

            auto g_z_0_z_xxxyy = cbuffer.data(ph_geom_10_off + 171 * ccomps * dcomps);

            auto g_z_0_z_xxxyz = cbuffer.data(ph_geom_10_off + 172 * ccomps * dcomps);

            auto g_z_0_z_xxxzz = cbuffer.data(ph_geom_10_off + 173 * ccomps * dcomps);

            auto g_z_0_z_xxyyy = cbuffer.data(ph_geom_10_off + 174 * ccomps * dcomps);

            auto g_z_0_z_xxyyz = cbuffer.data(ph_geom_10_off + 175 * ccomps * dcomps);

            auto g_z_0_z_xxyzz = cbuffer.data(ph_geom_10_off + 176 * ccomps * dcomps);

            auto g_z_0_z_xxzzz = cbuffer.data(ph_geom_10_off + 177 * ccomps * dcomps);

            auto g_z_0_z_xyyyy = cbuffer.data(ph_geom_10_off + 178 * ccomps * dcomps);

            auto g_z_0_z_xyyyz = cbuffer.data(ph_geom_10_off + 179 * ccomps * dcomps);

            auto g_z_0_z_xyyzz = cbuffer.data(ph_geom_10_off + 180 * ccomps * dcomps);

            auto g_z_0_z_xyzzz = cbuffer.data(ph_geom_10_off + 181 * ccomps * dcomps);

            auto g_z_0_z_xzzzz = cbuffer.data(ph_geom_10_off + 182 * ccomps * dcomps);

            auto g_z_0_z_yyyyy = cbuffer.data(ph_geom_10_off + 183 * ccomps * dcomps);

            auto g_z_0_z_yyyyz = cbuffer.data(ph_geom_10_off + 184 * ccomps * dcomps);

            auto g_z_0_z_yyyzz = cbuffer.data(ph_geom_10_off + 185 * ccomps * dcomps);

            auto g_z_0_z_yyzzz = cbuffer.data(ph_geom_10_off + 186 * ccomps * dcomps);

            auto g_z_0_z_yzzzz = cbuffer.data(ph_geom_10_off + 187 * ccomps * dcomps);

            auto g_z_0_z_zzzzz = cbuffer.data(ph_geom_10_off + 188 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_xxxxx, g_0_xxxxy, g_0_xxxxz, g_0_xxxyy, g_0_xxxyz, g_0_xxxzz, g_0_xxyyy, g_0_xxyyz, g_0_xxyzz, g_0_xxzzz, g_0_xyyyy, g_0_xyyyz, g_0_xyyzz, g_0_xyzzz, g_0_xzzzz, g_0_yyyyy, g_0_yyyyz, g_0_yyyzz, g_0_yyzzz, g_0_yzzzz, g_0_zzzzz, g_z_0_0_xxxxx, g_z_0_0_xxxxxz, g_z_0_0_xxxxy, g_z_0_0_xxxxyz, g_z_0_0_xxxxz, g_z_0_0_xxxxzz, g_z_0_0_xxxyy, g_z_0_0_xxxyyz, g_z_0_0_xxxyz, g_z_0_0_xxxyzz, g_z_0_0_xxxzz, g_z_0_0_xxxzzz, g_z_0_0_xxyyy, g_z_0_0_xxyyyz, g_z_0_0_xxyyz, g_z_0_0_xxyyzz, g_z_0_0_xxyzz, g_z_0_0_xxyzzz, g_z_0_0_xxzzz, g_z_0_0_xxzzzz, g_z_0_0_xyyyy, g_z_0_0_xyyyyz, g_z_0_0_xyyyz, g_z_0_0_xyyyzz, g_z_0_0_xyyzz, g_z_0_0_xyyzzz, g_z_0_0_xyzzz, g_z_0_0_xyzzzz, g_z_0_0_xzzzz, g_z_0_0_xzzzzz, g_z_0_0_yyyyy, g_z_0_0_yyyyyz, g_z_0_0_yyyyz, g_z_0_0_yyyyzz, g_z_0_0_yyyzz, g_z_0_0_yyyzzz, g_z_0_0_yyzzz, g_z_0_0_yyzzzz, g_z_0_0_yzzzz, g_z_0_0_yzzzzz, g_z_0_0_zzzzz, g_z_0_0_zzzzzz, g_z_0_z_xxxxx, g_z_0_z_xxxxy, g_z_0_z_xxxxz, g_z_0_z_xxxyy, g_z_0_z_xxxyz, g_z_0_z_xxxzz, g_z_0_z_xxyyy, g_z_0_z_xxyyz, g_z_0_z_xxyzz, g_z_0_z_xxzzz, g_z_0_z_xyyyy, g_z_0_z_xyyyz, g_z_0_z_xyyzz, g_z_0_z_xyzzz, g_z_0_z_xzzzz, g_z_0_z_yyyyy, g_z_0_z_yyyyz, g_z_0_z_yyyzz, g_z_0_z_yyzzz, g_z_0_z_yzzzz, g_z_0_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_xxxxx[k] = -g_0_xxxxx[k] - g_z_0_0_xxxxx[k] * ab_z + g_z_0_0_xxxxxz[k];

                g_z_0_z_xxxxy[k] = -g_0_xxxxy[k] - g_z_0_0_xxxxy[k] * ab_z + g_z_0_0_xxxxyz[k];

                g_z_0_z_xxxxz[k] = -g_0_xxxxz[k] - g_z_0_0_xxxxz[k] * ab_z + g_z_0_0_xxxxzz[k];

                g_z_0_z_xxxyy[k] = -g_0_xxxyy[k] - g_z_0_0_xxxyy[k] * ab_z + g_z_0_0_xxxyyz[k];

                g_z_0_z_xxxyz[k] = -g_0_xxxyz[k] - g_z_0_0_xxxyz[k] * ab_z + g_z_0_0_xxxyzz[k];

                g_z_0_z_xxxzz[k] = -g_0_xxxzz[k] - g_z_0_0_xxxzz[k] * ab_z + g_z_0_0_xxxzzz[k];

                g_z_0_z_xxyyy[k] = -g_0_xxyyy[k] - g_z_0_0_xxyyy[k] * ab_z + g_z_0_0_xxyyyz[k];

                g_z_0_z_xxyyz[k] = -g_0_xxyyz[k] - g_z_0_0_xxyyz[k] * ab_z + g_z_0_0_xxyyzz[k];

                g_z_0_z_xxyzz[k] = -g_0_xxyzz[k] - g_z_0_0_xxyzz[k] * ab_z + g_z_0_0_xxyzzz[k];

                g_z_0_z_xxzzz[k] = -g_0_xxzzz[k] - g_z_0_0_xxzzz[k] * ab_z + g_z_0_0_xxzzzz[k];

                g_z_0_z_xyyyy[k] = -g_0_xyyyy[k] - g_z_0_0_xyyyy[k] * ab_z + g_z_0_0_xyyyyz[k];

                g_z_0_z_xyyyz[k] = -g_0_xyyyz[k] - g_z_0_0_xyyyz[k] * ab_z + g_z_0_0_xyyyzz[k];

                g_z_0_z_xyyzz[k] = -g_0_xyyzz[k] - g_z_0_0_xyyzz[k] * ab_z + g_z_0_0_xyyzzz[k];

                g_z_0_z_xyzzz[k] = -g_0_xyzzz[k] - g_z_0_0_xyzzz[k] * ab_z + g_z_0_0_xyzzzz[k];

                g_z_0_z_xzzzz[k] = -g_0_xzzzz[k] - g_z_0_0_xzzzz[k] * ab_z + g_z_0_0_xzzzzz[k];

                g_z_0_z_yyyyy[k] = -g_0_yyyyy[k] - g_z_0_0_yyyyy[k] * ab_z + g_z_0_0_yyyyyz[k];

                g_z_0_z_yyyyz[k] = -g_0_yyyyz[k] - g_z_0_0_yyyyz[k] * ab_z + g_z_0_0_yyyyzz[k];

                g_z_0_z_yyyzz[k] = -g_0_yyyzz[k] - g_z_0_0_yyyzz[k] * ab_z + g_z_0_0_yyyzzz[k];

                g_z_0_z_yyzzz[k] = -g_0_yyzzz[k] - g_z_0_0_yyzzz[k] * ab_z + g_z_0_0_yyzzzz[k];

                g_z_0_z_yzzzz[k] = -g_0_yzzzz[k] - g_z_0_0_yzzzz[k] * ab_z + g_z_0_0_yzzzzz[k];

                g_z_0_z_zzzzz[k] = -g_0_zzzzz[k] - g_z_0_0_zzzzz[k] * ab_z + g_z_0_0_zzzzzz[k];
            }
        }
    }
}

} // erirec namespace

