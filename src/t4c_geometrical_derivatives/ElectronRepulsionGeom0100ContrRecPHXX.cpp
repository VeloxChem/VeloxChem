#include "ElectronRepulsionGeom0100ContrRecPHXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_phxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_phxx,
                                            const size_t idx_shxx,
                                            const size_t idx_geom_01_shxx,
                                            const size_t idx_geom_01_sixx,
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

            /// Set up components of auxilary buffer : SISS

            const auto si_geom_01_off = idx_geom_01_sixx + i * dcomps + j;

            auto g_0_x_0_xxxxxx = cbuffer.data(si_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_0_xxxxxy = cbuffer.data(si_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_0_xxxxxz = cbuffer.data(si_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_0_xxxxyy = cbuffer.data(si_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_0_xxxxyz = cbuffer.data(si_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_0_xxxxzz = cbuffer.data(si_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_0_xxxyyy = cbuffer.data(si_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_0_xxxyyz = cbuffer.data(si_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_0_xxxyzz = cbuffer.data(si_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_0_xxxzzz = cbuffer.data(si_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_0_xxyyyy = cbuffer.data(si_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_0_xxyyyz = cbuffer.data(si_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_0_xxyyzz = cbuffer.data(si_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_0_xxyzzz = cbuffer.data(si_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_0_xxzzzz = cbuffer.data(si_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_0_xyyyyy = cbuffer.data(si_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_0_xyyyyz = cbuffer.data(si_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_0_xyyyzz = cbuffer.data(si_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_0_xyyzzz = cbuffer.data(si_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_0_xyzzzz = cbuffer.data(si_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_0_xzzzzz = cbuffer.data(si_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_0_yyyyyy = cbuffer.data(si_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_0_yyyyyz = cbuffer.data(si_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_0_yyyyzz = cbuffer.data(si_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_0_yyyzzz = cbuffer.data(si_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_0_yyzzzz = cbuffer.data(si_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_0_yzzzzz = cbuffer.data(si_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_0_zzzzzz = cbuffer.data(si_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_y_0_xxxxxx = cbuffer.data(si_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_y_0_xxxxxy = cbuffer.data(si_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_y_0_xxxxxz = cbuffer.data(si_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_y_0_xxxxyy = cbuffer.data(si_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_y_0_xxxxyz = cbuffer.data(si_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_y_0_xxxxzz = cbuffer.data(si_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_y_0_xxxyyy = cbuffer.data(si_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_y_0_xxxyyz = cbuffer.data(si_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_y_0_xxxyzz = cbuffer.data(si_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_y_0_xxxzzz = cbuffer.data(si_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_y_0_xxyyyy = cbuffer.data(si_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_y_0_xxyyyz = cbuffer.data(si_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_y_0_xxyyzz = cbuffer.data(si_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_y_0_xxyzzz = cbuffer.data(si_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_y_0_xxzzzz = cbuffer.data(si_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_y_0_xyyyyy = cbuffer.data(si_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_y_0_xyyyyz = cbuffer.data(si_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_y_0_xyyyzz = cbuffer.data(si_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_y_0_xyyzzz = cbuffer.data(si_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_y_0_xyzzzz = cbuffer.data(si_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_y_0_xzzzzz = cbuffer.data(si_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_y_0_yyyyyy = cbuffer.data(si_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_y_0_yyyyyz = cbuffer.data(si_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_y_0_yyyyzz = cbuffer.data(si_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_y_0_yyyzzz = cbuffer.data(si_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_y_0_yyzzzz = cbuffer.data(si_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_y_0_yzzzzz = cbuffer.data(si_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_y_0_zzzzzz = cbuffer.data(si_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_z_0_xxxxxx = cbuffer.data(si_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_z_0_xxxxxy = cbuffer.data(si_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_z_0_xxxxxz = cbuffer.data(si_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_z_0_xxxxyy = cbuffer.data(si_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_z_0_xxxxyz = cbuffer.data(si_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_z_0_xxxxzz = cbuffer.data(si_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_z_0_xxxyyy = cbuffer.data(si_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_z_0_xxxyyz = cbuffer.data(si_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_z_0_xxxyzz = cbuffer.data(si_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_z_0_xxxzzz = cbuffer.data(si_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_z_0_xxyyyy = cbuffer.data(si_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_z_0_xxyyyz = cbuffer.data(si_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_z_0_xxyyzz = cbuffer.data(si_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_z_0_xxyzzz = cbuffer.data(si_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_z_0_xxzzzz = cbuffer.data(si_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_z_0_xyyyyy = cbuffer.data(si_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_z_0_xyyyyz = cbuffer.data(si_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_z_0_xyyyzz = cbuffer.data(si_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_z_0_xyyzzz = cbuffer.data(si_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_z_0_xyzzzz = cbuffer.data(si_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_z_0_xzzzzz = cbuffer.data(si_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_z_0_yyyyyy = cbuffer.data(si_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_z_0_yyyyyz = cbuffer.data(si_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_z_0_yyyyzz = cbuffer.data(si_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_z_0_yyyzzz = cbuffer.data(si_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_z_0_yyzzzz = cbuffer.data(si_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_z_0_yzzzzz = cbuffer.data(si_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_z_0_zzzzzz = cbuffer.data(si_geom_01_off + 83 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_phxx

            const auto ph_geom_01_off = idx_geom_01_phxx + i * dcomps + j;

            /// Set up 0-21 components of targeted buffer : cbuffer.data(

            auto g_0_x_x_xxxxx = cbuffer.data(ph_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_x_xxxxy = cbuffer.data(ph_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_x_xxxxz = cbuffer.data(ph_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_x_xxxyy = cbuffer.data(ph_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_x_xxxyz = cbuffer.data(ph_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_x_xxxzz = cbuffer.data(ph_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_x_xxyyy = cbuffer.data(ph_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_x_xxyyz = cbuffer.data(ph_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_x_xxyzz = cbuffer.data(ph_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_x_xxzzz = cbuffer.data(ph_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_x_xyyyy = cbuffer.data(ph_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_x_xyyyz = cbuffer.data(ph_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_x_xyyzz = cbuffer.data(ph_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_x_xyzzz = cbuffer.data(ph_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_x_xzzzz = cbuffer.data(ph_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_x_yyyyy = cbuffer.data(ph_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_x_yyyyz = cbuffer.data(ph_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_x_yyyzz = cbuffer.data(ph_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_x_yyzzz = cbuffer.data(ph_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_x_yzzzz = cbuffer.data(ph_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_x_zzzzz = cbuffer.data(ph_geom_01_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xxxxx, g_0_x_0_xxxxxx, g_0_x_0_xxxxxy, g_0_x_0_xxxxxz, g_0_x_0_xxxxy, g_0_x_0_xxxxyy, g_0_x_0_xxxxyz, g_0_x_0_xxxxz, g_0_x_0_xxxxzz, g_0_x_0_xxxyy, g_0_x_0_xxxyyy, g_0_x_0_xxxyyz, g_0_x_0_xxxyz, g_0_x_0_xxxyzz, g_0_x_0_xxxzz, g_0_x_0_xxxzzz, g_0_x_0_xxyyy, g_0_x_0_xxyyyy, g_0_x_0_xxyyyz, g_0_x_0_xxyyz, g_0_x_0_xxyyzz, g_0_x_0_xxyzz, g_0_x_0_xxyzzz, g_0_x_0_xxzzz, g_0_x_0_xxzzzz, g_0_x_0_xyyyy, g_0_x_0_xyyyyy, g_0_x_0_xyyyyz, g_0_x_0_xyyyz, g_0_x_0_xyyyzz, g_0_x_0_xyyzz, g_0_x_0_xyyzzz, g_0_x_0_xyzzz, g_0_x_0_xyzzzz, g_0_x_0_xzzzz, g_0_x_0_xzzzzz, g_0_x_0_yyyyy, g_0_x_0_yyyyz, g_0_x_0_yyyzz, g_0_x_0_yyzzz, g_0_x_0_yzzzz, g_0_x_0_zzzzz, g_0_x_x_xxxxx, g_0_x_x_xxxxy, g_0_x_x_xxxxz, g_0_x_x_xxxyy, g_0_x_x_xxxyz, g_0_x_x_xxxzz, g_0_x_x_xxyyy, g_0_x_x_xxyyz, g_0_x_x_xxyzz, g_0_x_x_xxzzz, g_0_x_x_xyyyy, g_0_x_x_xyyyz, g_0_x_x_xyyzz, g_0_x_x_xyzzz, g_0_x_x_xzzzz, g_0_x_x_yyyyy, g_0_x_x_yyyyz, g_0_x_x_yyyzz, g_0_x_x_yyzzz, g_0_x_x_yzzzz, g_0_x_x_zzzzz, g_0_xxxxx, g_0_xxxxy, g_0_xxxxz, g_0_xxxyy, g_0_xxxyz, g_0_xxxzz, g_0_xxyyy, g_0_xxyyz, g_0_xxyzz, g_0_xxzzz, g_0_xyyyy, g_0_xyyyz, g_0_xyyzz, g_0_xyzzz, g_0_xzzzz, g_0_yyyyy, g_0_yyyyz, g_0_yyyzz, g_0_yyzzz, g_0_yzzzz, g_0_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_x_xxxxx[k] = g_0_xxxxx[k] - g_0_x_0_xxxxx[k] * ab_x + g_0_x_0_xxxxxx[k];

                g_0_x_x_xxxxy[k] = g_0_xxxxy[k] - g_0_x_0_xxxxy[k] * ab_x + g_0_x_0_xxxxxy[k];

                g_0_x_x_xxxxz[k] = g_0_xxxxz[k] - g_0_x_0_xxxxz[k] * ab_x + g_0_x_0_xxxxxz[k];

                g_0_x_x_xxxyy[k] = g_0_xxxyy[k] - g_0_x_0_xxxyy[k] * ab_x + g_0_x_0_xxxxyy[k];

                g_0_x_x_xxxyz[k] = g_0_xxxyz[k] - g_0_x_0_xxxyz[k] * ab_x + g_0_x_0_xxxxyz[k];

                g_0_x_x_xxxzz[k] = g_0_xxxzz[k] - g_0_x_0_xxxzz[k] * ab_x + g_0_x_0_xxxxzz[k];

                g_0_x_x_xxyyy[k] = g_0_xxyyy[k] - g_0_x_0_xxyyy[k] * ab_x + g_0_x_0_xxxyyy[k];

                g_0_x_x_xxyyz[k] = g_0_xxyyz[k] - g_0_x_0_xxyyz[k] * ab_x + g_0_x_0_xxxyyz[k];

                g_0_x_x_xxyzz[k] = g_0_xxyzz[k] - g_0_x_0_xxyzz[k] * ab_x + g_0_x_0_xxxyzz[k];

                g_0_x_x_xxzzz[k] = g_0_xxzzz[k] - g_0_x_0_xxzzz[k] * ab_x + g_0_x_0_xxxzzz[k];

                g_0_x_x_xyyyy[k] = g_0_xyyyy[k] - g_0_x_0_xyyyy[k] * ab_x + g_0_x_0_xxyyyy[k];

                g_0_x_x_xyyyz[k] = g_0_xyyyz[k] - g_0_x_0_xyyyz[k] * ab_x + g_0_x_0_xxyyyz[k];

                g_0_x_x_xyyzz[k] = g_0_xyyzz[k] - g_0_x_0_xyyzz[k] * ab_x + g_0_x_0_xxyyzz[k];

                g_0_x_x_xyzzz[k] = g_0_xyzzz[k] - g_0_x_0_xyzzz[k] * ab_x + g_0_x_0_xxyzzz[k];

                g_0_x_x_xzzzz[k] = g_0_xzzzz[k] - g_0_x_0_xzzzz[k] * ab_x + g_0_x_0_xxzzzz[k];

                g_0_x_x_yyyyy[k] = g_0_yyyyy[k] - g_0_x_0_yyyyy[k] * ab_x + g_0_x_0_xyyyyy[k];

                g_0_x_x_yyyyz[k] = g_0_yyyyz[k] - g_0_x_0_yyyyz[k] * ab_x + g_0_x_0_xyyyyz[k];

                g_0_x_x_yyyzz[k] = g_0_yyyzz[k] - g_0_x_0_yyyzz[k] * ab_x + g_0_x_0_xyyyzz[k];

                g_0_x_x_yyzzz[k] = g_0_yyzzz[k] - g_0_x_0_yyzzz[k] * ab_x + g_0_x_0_xyyzzz[k];

                g_0_x_x_yzzzz[k] = g_0_yzzzz[k] - g_0_x_0_yzzzz[k] * ab_x + g_0_x_0_xyzzzz[k];

                g_0_x_x_zzzzz[k] = g_0_zzzzz[k] - g_0_x_0_zzzzz[k] * ab_x + g_0_x_0_xzzzzz[k];
            }

            /// Set up 21-42 components of targeted buffer : cbuffer.data(

            auto g_0_x_y_xxxxx = cbuffer.data(ph_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_y_xxxxy = cbuffer.data(ph_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_y_xxxxz = cbuffer.data(ph_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_y_xxxyy = cbuffer.data(ph_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_y_xxxyz = cbuffer.data(ph_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_y_xxxzz = cbuffer.data(ph_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_y_xxyyy = cbuffer.data(ph_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_y_xxyyz = cbuffer.data(ph_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_y_xxyzz = cbuffer.data(ph_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_y_xxzzz = cbuffer.data(ph_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_y_xyyyy = cbuffer.data(ph_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_y_xyyyz = cbuffer.data(ph_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_y_xyyzz = cbuffer.data(ph_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_y_xyzzz = cbuffer.data(ph_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_y_xzzzz = cbuffer.data(ph_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_y_yyyyy = cbuffer.data(ph_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_y_yyyyz = cbuffer.data(ph_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_y_yyyzz = cbuffer.data(ph_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_y_yyzzz = cbuffer.data(ph_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_y_yzzzz = cbuffer.data(ph_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_y_zzzzz = cbuffer.data(ph_geom_01_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xxxxx, g_0_x_0_xxxxxy, g_0_x_0_xxxxy, g_0_x_0_xxxxyy, g_0_x_0_xxxxyz, g_0_x_0_xxxxz, g_0_x_0_xxxyy, g_0_x_0_xxxyyy, g_0_x_0_xxxyyz, g_0_x_0_xxxyz, g_0_x_0_xxxyzz, g_0_x_0_xxxzz, g_0_x_0_xxyyy, g_0_x_0_xxyyyy, g_0_x_0_xxyyyz, g_0_x_0_xxyyz, g_0_x_0_xxyyzz, g_0_x_0_xxyzz, g_0_x_0_xxyzzz, g_0_x_0_xxzzz, g_0_x_0_xyyyy, g_0_x_0_xyyyyy, g_0_x_0_xyyyyz, g_0_x_0_xyyyz, g_0_x_0_xyyyzz, g_0_x_0_xyyzz, g_0_x_0_xyyzzz, g_0_x_0_xyzzz, g_0_x_0_xyzzzz, g_0_x_0_xzzzz, g_0_x_0_yyyyy, g_0_x_0_yyyyyy, g_0_x_0_yyyyyz, g_0_x_0_yyyyz, g_0_x_0_yyyyzz, g_0_x_0_yyyzz, g_0_x_0_yyyzzz, g_0_x_0_yyzzz, g_0_x_0_yyzzzz, g_0_x_0_yzzzz, g_0_x_0_yzzzzz, g_0_x_0_zzzzz, g_0_x_y_xxxxx, g_0_x_y_xxxxy, g_0_x_y_xxxxz, g_0_x_y_xxxyy, g_0_x_y_xxxyz, g_0_x_y_xxxzz, g_0_x_y_xxyyy, g_0_x_y_xxyyz, g_0_x_y_xxyzz, g_0_x_y_xxzzz, g_0_x_y_xyyyy, g_0_x_y_xyyyz, g_0_x_y_xyyzz, g_0_x_y_xyzzz, g_0_x_y_xzzzz, g_0_x_y_yyyyy, g_0_x_y_yyyyz, g_0_x_y_yyyzz, g_0_x_y_yyzzz, g_0_x_y_yzzzz, g_0_x_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_y_xxxxx[k] = -g_0_x_0_xxxxx[k] * ab_y + g_0_x_0_xxxxxy[k];

                g_0_x_y_xxxxy[k] = -g_0_x_0_xxxxy[k] * ab_y + g_0_x_0_xxxxyy[k];

                g_0_x_y_xxxxz[k] = -g_0_x_0_xxxxz[k] * ab_y + g_0_x_0_xxxxyz[k];

                g_0_x_y_xxxyy[k] = -g_0_x_0_xxxyy[k] * ab_y + g_0_x_0_xxxyyy[k];

                g_0_x_y_xxxyz[k] = -g_0_x_0_xxxyz[k] * ab_y + g_0_x_0_xxxyyz[k];

                g_0_x_y_xxxzz[k] = -g_0_x_0_xxxzz[k] * ab_y + g_0_x_0_xxxyzz[k];

                g_0_x_y_xxyyy[k] = -g_0_x_0_xxyyy[k] * ab_y + g_0_x_0_xxyyyy[k];

                g_0_x_y_xxyyz[k] = -g_0_x_0_xxyyz[k] * ab_y + g_0_x_0_xxyyyz[k];

                g_0_x_y_xxyzz[k] = -g_0_x_0_xxyzz[k] * ab_y + g_0_x_0_xxyyzz[k];

                g_0_x_y_xxzzz[k] = -g_0_x_0_xxzzz[k] * ab_y + g_0_x_0_xxyzzz[k];

                g_0_x_y_xyyyy[k] = -g_0_x_0_xyyyy[k] * ab_y + g_0_x_0_xyyyyy[k];

                g_0_x_y_xyyyz[k] = -g_0_x_0_xyyyz[k] * ab_y + g_0_x_0_xyyyyz[k];

                g_0_x_y_xyyzz[k] = -g_0_x_0_xyyzz[k] * ab_y + g_0_x_0_xyyyzz[k];

                g_0_x_y_xyzzz[k] = -g_0_x_0_xyzzz[k] * ab_y + g_0_x_0_xyyzzz[k];

                g_0_x_y_xzzzz[k] = -g_0_x_0_xzzzz[k] * ab_y + g_0_x_0_xyzzzz[k];

                g_0_x_y_yyyyy[k] = -g_0_x_0_yyyyy[k] * ab_y + g_0_x_0_yyyyyy[k];

                g_0_x_y_yyyyz[k] = -g_0_x_0_yyyyz[k] * ab_y + g_0_x_0_yyyyyz[k];

                g_0_x_y_yyyzz[k] = -g_0_x_0_yyyzz[k] * ab_y + g_0_x_0_yyyyzz[k];

                g_0_x_y_yyzzz[k] = -g_0_x_0_yyzzz[k] * ab_y + g_0_x_0_yyyzzz[k];

                g_0_x_y_yzzzz[k] = -g_0_x_0_yzzzz[k] * ab_y + g_0_x_0_yyzzzz[k];

                g_0_x_y_zzzzz[k] = -g_0_x_0_zzzzz[k] * ab_y + g_0_x_0_yzzzzz[k];
            }

            /// Set up 42-63 components of targeted buffer : cbuffer.data(

            auto g_0_x_z_xxxxx = cbuffer.data(ph_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_z_xxxxy = cbuffer.data(ph_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_z_xxxxz = cbuffer.data(ph_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_z_xxxyy = cbuffer.data(ph_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_z_xxxyz = cbuffer.data(ph_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_z_xxxzz = cbuffer.data(ph_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_z_xxyyy = cbuffer.data(ph_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_z_xxyyz = cbuffer.data(ph_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_z_xxyzz = cbuffer.data(ph_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_z_xxzzz = cbuffer.data(ph_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_z_xyyyy = cbuffer.data(ph_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_z_xyyyz = cbuffer.data(ph_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_z_xyyzz = cbuffer.data(ph_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_z_xyzzz = cbuffer.data(ph_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_z_xzzzz = cbuffer.data(ph_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_z_yyyyy = cbuffer.data(ph_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_z_yyyyz = cbuffer.data(ph_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_z_yyyzz = cbuffer.data(ph_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_z_yyzzz = cbuffer.data(ph_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_z_yzzzz = cbuffer.data(ph_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_z_zzzzz = cbuffer.data(ph_geom_01_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xxxxx, g_0_x_0_xxxxxz, g_0_x_0_xxxxy, g_0_x_0_xxxxyz, g_0_x_0_xxxxz, g_0_x_0_xxxxzz, g_0_x_0_xxxyy, g_0_x_0_xxxyyz, g_0_x_0_xxxyz, g_0_x_0_xxxyzz, g_0_x_0_xxxzz, g_0_x_0_xxxzzz, g_0_x_0_xxyyy, g_0_x_0_xxyyyz, g_0_x_0_xxyyz, g_0_x_0_xxyyzz, g_0_x_0_xxyzz, g_0_x_0_xxyzzz, g_0_x_0_xxzzz, g_0_x_0_xxzzzz, g_0_x_0_xyyyy, g_0_x_0_xyyyyz, g_0_x_0_xyyyz, g_0_x_0_xyyyzz, g_0_x_0_xyyzz, g_0_x_0_xyyzzz, g_0_x_0_xyzzz, g_0_x_0_xyzzzz, g_0_x_0_xzzzz, g_0_x_0_xzzzzz, g_0_x_0_yyyyy, g_0_x_0_yyyyyz, g_0_x_0_yyyyz, g_0_x_0_yyyyzz, g_0_x_0_yyyzz, g_0_x_0_yyyzzz, g_0_x_0_yyzzz, g_0_x_0_yyzzzz, g_0_x_0_yzzzz, g_0_x_0_yzzzzz, g_0_x_0_zzzzz, g_0_x_0_zzzzzz, g_0_x_z_xxxxx, g_0_x_z_xxxxy, g_0_x_z_xxxxz, g_0_x_z_xxxyy, g_0_x_z_xxxyz, g_0_x_z_xxxzz, g_0_x_z_xxyyy, g_0_x_z_xxyyz, g_0_x_z_xxyzz, g_0_x_z_xxzzz, g_0_x_z_xyyyy, g_0_x_z_xyyyz, g_0_x_z_xyyzz, g_0_x_z_xyzzz, g_0_x_z_xzzzz, g_0_x_z_yyyyy, g_0_x_z_yyyyz, g_0_x_z_yyyzz, g_0_x_z_yyzzz, g_0_x_z_yzzzz, g_0_x_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_z_xxxxx[k] = -g_0_x_0_xxxxx[k] * ab_z + g_0_x_0_xxxxxz[k];

                g_0_x_z_xxxxy[k] = -g_0_x_0_xxxxy[k] * ab_z + g_0_x_0_xxxxyz[k];

                g_0_x_z_xxxxz[k] = -g_0_x_0_xxxxz[k] * ab_z + g_0_x_0_xxxxzz[k];

                g_0_x_z_xxxyy[k] = -g_0_x_0_xxxyy[k] * ab_z + g_0_x_0_xxxyyz[k];

                g_0_x_z_xxxyz[k] = -g_0_x_0_xxxyz[k] * ab_z + g_0_x_0_xxxyzz[k];

                g_0_x_z_xxxzz[k] = -g_0_x_0_xxxzz[k] * ab_z + g_0_x_0_xxxzzz[k];

                g_0_x_z_xxyyy[k] = -g_0_x_0_xxyyy[k] * ab_z + g_0_x_0_xxyyyz[k];

                g_0_x_z_xxyyz[k] = -g_0_x_0_xxyyz[k] * ab_z + g_0_x_0_xxyyzz[k];

                g_0_x_z_xxyzz[k] = -g_0_x_0_xxyzz[k] * ab_z + g_0_x_0_xxyzzz[k];

                g_0_x_z_xxzzz[k] = -g_0_x_0_xxzzz[k] * ab_z + g_0_x_0_xxzzzz[k];

                g_0_x_z_xyyyy[k] = -g_0_x_0_xyyyy[k] * ab_z + g_0_x_0_xyyyyz[k];

                g_0_x_z_xyyyz[k] = -g_0_x_0_xyyyz[k] * ab_z + g_0_x_0_xyyyzz[k];

                g_0_x_z_xyyzz[k] = -g_0_x_0_xyyzz[k] * ab_z + g_0_x_0_xyyzzz[k];

                g_0_x_z_xyzzz[k] = -g_0_x_0_xyzzz[k] * ab_z + g_0_x_0_xyzzzz[k];

                g_0_x_z_xzzzz[k] = -g_0_x_0_xzzzz[k] * ab_z + g_0_x_0_xzzzzz[k];

                g_0_x_z_yyyyy[k] = -g_0_x_0_yyyyy[k] * ab_z + g_0_x_0_yyyyyz[k];

                g_0_x_z_yyyyz[k] = -g_0_x_0_yyyyz[k] * ab_z + g_0_x_0_yyyyzz[k];

                g_0_x_z_yyyzz[k] = -g_0_x_0_yyyzz[k] * ab_z + g_0_x_0_yyyzzz[k];

                g_0_x_z_yyzzz[k] = -g_0_x_0_yyzzz[k] * ab_z + g_0_x_0_yyzzzz[k];

                g_0_x_z_yzzzz[k] = -g_0_x_0_yzzzz[k] * ab_z + g_0_x_0_yzzzzz[k];

                g_0_x_z_zzzzz[k] = -g_0_x_0_zzzzz[k] * ab_z + g_0_x_0_zzzzzz[k];
            }

            /// Set up 63-84 components of targeted buffer : cbuffer.data(

            auto g_0_y_x_xxxxx = cbuffer.data(ph_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_y_x_xxxxy = cbuffer.data(ph_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_y_x_xxxxz = cbuffer.data(ph_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_y_x_xxxyy = cbuffer.data(ph_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_y_x_xxxyz = cbuffer.data(ph_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_y_x_xxxzz = cbuffer.data(ph_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_y_x_xxyyy = cbuffer.data(ph_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_y_x_xxyyz = cbuffer.data(ph_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_y_x_xxyzz = cbuffer.data(ph_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_y_x_xxzzz = cbuffer.data(ph_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_y_x_xyyyy = cbuffer.data(ph_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_y_x_xyyyz = cbuffer.data(ph_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_y_x_xyyzz = cbuffer.data(ph_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_y_x_xyzzz = cbuffer.data(ph_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_y_x_xzzzz = cbuffer.data(ph_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_y_x_yyyyy = cbuffer.data(ph_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_y_x_yyyyz = cbuffer.data(ph_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_y_x_yyyzz = cbuffer.data(ph_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_y_x_yyzzz = cbuffer.data(ph_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_y_x_yzzzz = cbuffer.data(ph_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_y_x_zzzzz = cbuffer.data(ph_geom_01_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_0_xxxxx, g_0_y_0_xxxxxx, g_0_y_0_xxxxxy, g_0_y_0_xxxxxz, g_0_y_0_xxxxy, g_0_y_0_xxxxyy, g_0_y_0_xxxxyz, g_0_y_0_xxxxz, g_0_y_0_xxxxzz, g_0_y_0_xxxyy, g_0_y_0_xxxyyy, g_0_y_0_xxxyyz, g_0_y_0_xxxyz, g_0_y_0_xxxyzz, g_0_y_0_xxxzz, g_0_y_0_xxxzzz, g_0_y_0_xxyyy, g_0_y_0_xxyyyy, g_0_y_0_xxyyyz, g_0_y_0_xxyyz, g_0_y_0_xxyyzz, g_0_y_0_xxyzz, g_0_y_0_xxyzzz, g_0_y_0_xxzzz, g_0_y_0_xxzzzz, g_0_y_0_xyyyy, g_0_y_0_xyyyyy, g_0_y_0_xyyyyz, g_0_y_0_xyyyz, g_0_y_0_xyyyzz, g_0_y_0_xyyzz, g_0_y_0_xyyzzz, g_0_y_0_xyzzz, g_0_y_0_xyzzzz, g_0_y_0_xzzzz, g_0_y_0_xzzzzz, g_0_y_0_yyyyy, g_0_y_0_yyyyz, g_0_y_0_yyyzz, g_0_y_0_yyzzz, g_0_y_0_yzzzz, g_0_y_0_zzzzz, g_0_y_x_xxxxx, g_0_y_x_xxxxy, g_0_y_x_xxxxz, g_0_y_x_xxxyy, g_0_y_x_xxxyz, g_0_y_x_xxxzz, g_0_y_x_xxyyy, g_0_y_x_xxyyz, g_0_y_x_xxyzz, g_0_y_x_xxzzz, g_0_y_x_xyyyy, g_0_y_x_xyyyz, g_0_y_x_xyyzz, g_0_y_x_xyzzz, g_0_y_x_xzzzz, g_0_y_x_yyyyy, g_0_y_x_yyyyz, g_0_y_x_yyyzz, g_0_y_x_yyzzz, g_0_y_x_yzzzz, g_0_y_x_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_x_xxxxx[k] = -g_0_y_0_xxxxx[k] * ab_x + g_0_y_0_xxxxxx[k];

                g_0_y_x_xxxxy[k] = -g_0_y_0_xxxxy[k] * ab_x + g_0_y_0_xxxxxy[k];

                g_0_y_x_xxxxz[k] = -g_0_y_0_xxxxz[k] * ab_x + g_0_y_0_xxxxxz[k];

                g_0_y_x_xxxyy[k] = -g_0_y_0_xxxyy[k] * ab_x + g_0_y_0_xxxxyy[k];

                g_0_y_x_xxxyz[k] = -g_0_y_0_xxxyz[k] * ab_x + g_0_y_0_xxxxyz[k];

                g_0_y_x_xxxzz[k] = -g_0_y_0_xxxzz[k] * ab_x + g_0_y_0_xxxxzz[k];

                g_0_y_x_xxyyy[k] = -g_0_y_0_xxyyy[k] * ab_x + g_0_y_0_xxxyyy[k];

                g_0_y_x_xxyyz[k] = -g_0_y_0_xxyyz[k] * ab_x + g_0_y_0_xxxyyz[k];

                g_0_y_x_xxyzz[k] = -g_0_y_0_xxyzz[k] * ab_x + g_0_y_0_xxxyzz[k];

                g_0_y_x_xxzzz[k] = -g_0_y_0_xxzzz[k] * ab_x + g_0_y_0_xxxzzz[k];

                g_0_y_x_xyyyy[k] = -g_0_y_0_xyyyy[k] * ab_x + g_0_y_0_xxyyyy[k];

                g_0_y_x_xyyyz[k] = -g_0_y_0_xyyyz[k] * ab_x + g_0_y_0_xxyyyz[k];

                g_0_y_x_xyyzz[k] = -g_0_y_0_xyyzz[k] * ab_x + g_0_y_0_xxyyzz[k];

                g_0_y_x_xyzzz[k] = -g_0_y_0_xyzzz[k] * ab_x + g_0_y_0_xxyzzz[k];

                g_0_y_x_xzzzz[k] = -g_0_y_0_xzzzz[k] * ab_x + g_0_y_0_xxzzzz[k];

                g_0_y_x_yyyyy[k] = -g_0_y_0_yyyyy[k] * ab_x + g_0_y_0_xyyyyy[k];

                g_0_y_x_yyyyz[k] = -g_0_y_0_yyyyz[k] * ab_x + g_0_y_0_xyyyyz[k];

                g_0_y_x_yyyzz[k] = -g_0_y_0_yyyzz[k] * ab_x + g_0_y_0_xyyyzz[k];

                g_0_y_x_yyzzz[k] = -g_0_y_0_yyzzz[k] * ab_x + g_0_y_0_xyyzzz[k];

                g_0_y_x_yzzzz[k] = -g_0_y_0_yzzzz[k] * ab_x + g_0_y_0_xyzzzz[k];

                g_0_y_x_zzzzz[k] = -g_0_y_0_zzzzz[k] * ab_x + g_0_y_0_xzzzzz[k];
            }

            /// Set up 84-105 components of targeted buffer : cbuffer.data(

            auto g_0_y_y_xxxxx = cbuffer.data(ph_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_y_y_xxxxy = cbuffer.data(ph_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_y_y_xxxxz = cbuffer.data(ph_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_y_y_xxxyy = cbuffer.data(ph_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_y_y_xxxyz = cbuffer.data(ph_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_y_y_xxxzz = cbuffer.data(ph_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_y_y_xxyyy = cbuffer.data(ph_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_y_y_xxyyz = cbuffer.data(ph_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_y_y_xxyzz = cbuffer.data(ph_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_y_y_xxzzz = cbuffer.data(ph_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_y_y_xyyyy = cbuffer.data(ph_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_y_y_xyyyz = cbuffer.data(ph_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_y_y_xyyzz = cbuffer.data(ph_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_y_y_xyzzz = cbuffer.data(ph_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_y_y_xzzzz = cbuffer.data(ph_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_y_y_yyyyy = cbuffer.data(ph_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_y_y_yyyyz = cbuffer.data(ph_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_y_y_yyyzz = cbuffer.data(ph_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_y_y_yyzzz = cbuffer.data(ph_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_y_y_yzzzz = cbuffer.data(ph_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_y_y_zzzzz = cbuffer.data(ph_geom_01_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_xxxxx, g_0_xxxxy, g_0_xxxxz, g_0_xxxyy, g_0_xxxyz, g_0_xxxzz, g_0_xxyyy, g_0_xxyyz, g_0_xxyzz, g_0_xxzzz, g_0_xyyyy, g_0_xyyyz, g_0_xyyzz, g_0_xyzzz, g_0_xzzzz, g_0_y_0_xxxxx, g_0_y_0_xxxxxy, g_0_y_0_xxxxy, g_0_y_0_xxxxyy, g_0_y_0_xxxxyz, g_0_y_0_xxxxz, g_0_y_0_xxxyy, g_0_y_0_xxxyyy, g_0_y_0_xxxyyz, g_0_y_0_xxxyz, g_0_y_0_xxxyzz, g_0_y_0_xxxzz, g_0_y_0_xxyyy, g_0_y_0_xxyyyy, g_0_y_0_xxyyyz, g_0_y_0_xxyyz, g_0_y_0_xxyyzz, g_0_y_0_xxyzz, g_0_y_0_xxyzzz, g_0_y_0_xxzzz, g_0_y_0_xyyyy, g_0_y_0_xyyyyy, g_0_y_0_xyyyyz, g_0_y_0_xyyyz, g_0_y_0_xyyyzz, g_0_y_0_xyyzz, g_0_y_0_xyyzzz, g_0_y_0_xyzzz, g_0_y_0_xyzzzz, g_0_y_0_xzzzz, g_0_y_0_yyyyy, g_0_y_0_yyyyyy, g_0_y_0_yyyyyz, g_0_y_0_yyyyz, g_0_y_0_yyyyzz, g_0_y_0_yyyzz, g_0_y_0_yyyzzz, g_0_y_0_yyzzz, g_0_y_0_yyzzzz, g_0_y_0_yzzzz, g_0_y_0_yzzzzz, g_0_y_0_zzzzz, g_0_y_y_xxxxx, g_0_y_y_xxxxy, g_0_y_y_xxxxz, g_0_y_y_xxxyy, g_0_y_y_xxxyz, g_0_y_y_xxxzz, g_0_y_y_xxyyy, g_0_y_y_xxyyz, g_0_y_y_xxyzz, g_0_y_y_xxzzz, g_0_y_y_xyyyy, g_0_y_y_xyyyz, g_0_y_y_xyyzz, g_0_y_y_xyzzz, g_0_y_y_xzzzz, g_0_y_y_yyyyy, g_0_y_y_yyyyz, g_0_y_y_yyyzz, g_0_y_y_yyzzz, g_0_y_y_yzzzz, g_0_y_y_zzzzz, g_0_yyyyy, g_0_yyyyz, g_0_yyyzz, g_0_yyzzz, g_0_yzzzz, g_0_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_y_xxxxx[k] = g_0_xxxxx[k] - g_0_y_0_xxxxx[k] * ab_y + g_0_y_0_xxxxxy[k];

                g_0_y_y_xxxxy[k] = g_0_xxxxy[k] - g_0_y_0_xxxxy[k] * ab_y + g_0_y_0_xxxxyy[k];

                g_0_y_y_xxxxz[k] = g_0_xxxxz[k] - g_0_y_0_xxxxz[k] * ab_y + g_0_y_0_xxxxyz[k];

                g_0_y_y_xxxyy[k] = g_0_xxxyy[k] - g_0_y_0_xxxyy[k] * ab_y + g_0_y_0_xxxyyy[k];

                g_0_y_y_xxxyz[k] = g_0_xxxyz[k] - g_0_y_0_xxxyz[k] * ab_y + g_0_y_0_xxxyyz[k];

                g_0_y_y_xxxzz[k] = g_0_xxxzz[k] - g_0_y_0_xxxzz[k] * ab_y + g_0_y_0_xxxyzz[k];

                g_0_y_y_xxyyy[k] = g_0_xxyyy[k] - g_0_y_0_xxyyy[k] * ab_y + g_0_y_0_xxyyyy[k];

                g_0_y_y_xxyyz[k] = g_0_xxyyz[k] - g_0_y_0_xxyyz[k] * ab_y + g_0_y_0_xxyyyz[k];

                g_0_y_y_xxyzz[k] = g_0_xxyzz[k] - g_0_y_0_xxyzz[k] * ab_y + g_0_y_0_xxyyzz[k];

                g_0_y_y_xxzzz[k] = g_0_xxzzz[k] - g_0_y_0_xxzzz[k] * ab_y + g_0_y_0_xxyzzz[k];

                g_0_y_y_xyyyy[k] = g_0_xyyyy[k] - g_0_y_0_xyyyy[k] * ab_y + g_0_y_0_xyyyyy[k];

                g_0_y_y_xyyyz[k] = g_0_xyyyz[k] - g_0_y_0_xyyyz[k] * ab_y + g_0_y_0_xyyyyz[k];

                g_0_y_y_xyyzz[k] = g_0_xyyzz[k] - g_0_y_0_xyyzz[k] * ab_y + g_0_y_0_xyyyzz[k];

                g_0_y_y_xyzzz[k] = g_0_xyzzz[k] - g_0_y_0_xyzzz[k] * ab_y + g_0_y_0_xyyzzz[k];

                g_0_y_y_xzzzz[k] = g_0_xzzzz[k] - g_0_y_0_xzzzz[k] * ab_y + g_0_y_0_xyzzzz[k];

                g_0_y_y_yyyyy[k] = g_0_yyyyy[k] - g_0_y_0_yyyyy[k] * ab_y + g_0_y_0_yyyyyy[k];

                g_0_y_y_yyyyz[k] = g_0_yyyyz[k] - g_0_y_0_yyyyz[k] * ab_y + g_0_y_0_yyyyyz[k];

                g_0_y_y_yyyzz[k] = g_0_yyyzz[k] - g_0_y_0_yyyzz[k] * ab_y + g_0_y_0_yyyyzz[k];

                g_0_y_y_yyzzz[k] = g_0_yyzzz[k] - g_0_y_0_yyzzz[k] * ab_y + g_0_y_0_yyyzzz[k];

                g_0_y_y_yzzzz[k] = g_0_yzzzz[k] - g_0_y_0_yzzzz[k] * ab_y + g_0_y_0_yyzzzz[k];

                g_0_y_y_zzzzz[k] = g_0_zzzzz[k] - g_0_y_0_zzzzz[k] * ab_y + g_0_y_0_yzzzzz[k];
            }

            /// Set up 105-126 components of targeted buffer : cbuffer.data(

            auto g_0_y_z_xxxxx = cbuffer.data(ph_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_y_z_xxxxy = cbuffer.data(ph_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_y_z_xxxxz = cbuffer.data(ph_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_y_z_xxxyy = cbuffer.data(ph_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_y_z_xxxyz = cbuffer.data(ph_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_y_z_xxxzz = cbuffer.data(ph_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_y_z_xxyyy = cbuffer.data(ph_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_y_z_xxyyz = cbuffer.data(ph_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_y_z_xxyzz = cbuffer.data(ph_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_y_z_xxzzz = cbuffer.data(ph_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_y_z_xyyyy = cbuffer.data(ph_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_y_z_xyyyz = cbuffer.data(ph_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_y_z_xyyzz = cbuffer.data(ph_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_y_z_xyzzz = cbuffer.data(ph_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_y_z_xzzzz = cbuffer.data(ph_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_y_z_yyyyy = cbuffer.data(ph_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_y_z_yyyyz = cbuffer.data(ph_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_y_z_yyyzz = cbuffer.data(ph_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_y_z_yyzzz = cbuffer.data(ph_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_y_z_yzzzz = cbuffer.data(ph_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_y_z_zzzzz = cbuffer.data(ph_geom_01_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_0_xxxxx, g_0_y_0_xxxxxz, g_0_y_0_xxxxy, g_0_y_0_xxxxyz, g_0_y_0_xxxxz, g_0_y_0_xxxxzz, g_0_y_0_xxxyy, g_0_y_0_xxxyyz, g_0_y_0_xxxyz, g_0_y_0_xxxyzz, g_0_y_0_xxxzz, g_0_y_0_xxxzzz, g_0_y_0_xxyyy, g_0_y_0_xxyyyz, g_0_y_0_xxyyz, g_0_y_0_xxyyzz, g_0_y_0_xxyzz, g_0_y_0_xxyzzz, g_0_y_0_xxzzz, g_0_y_0_xxzzzz, g_0_y_0_xyyyy, g_0_y_0_xyyyyz, g_0_y_0_xyyyz, g_0_y_0_xyyyzz, g_0_y_0_xyyzz, g_0_y_0_xyyzzz, g_0_y_0_xyzzz, g_0_y_0_xyzzzz, g_0_y_0_xzzzz, g_0_y_0_xzzzzz, g_0_y_0_yyyyy, g_0_y_0_yyyyyz, g_0_y_0_yyyyz, g_0_y_0_yyyyzz, g_0_y_0_yyyzz, g_0_y_0_yyyzzz, g_0_y_0_yyzzz, g_0_y_0_yyzzzz, g_0_y_0_yzzzz, g_0_y_0_yzzzzz, g_0_y_0_zzzzz, g_0_y_0_zzzzzz, g_0_y_z_xxxxx, g_0_y_z_xxxxy, g_0_y_z_xxxxz, g_0_y_z_xxxyy, g_0_y_z_xxxyz, g_0_y_z_xxxzz, g_0_y_z_xxyyy, g_0_y_z_xxyyz, g_0_y_z_xxyzz, g_0_y_z_xxzzz, g_0_y_z_xyyyy, g_0_y_z_xyyyz, g_0_y_z_xyyzz, g_0_y_z_xyzzz, g_0_y_z_xzzzz, g_0_y_z_yyyyy, g_0_y_z_yyyyz, g_0_y_z_yyyzz, g_0_y_z_yyzzz, g_0_y_z_yzzzz, g_0_y_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_z_xxxxx[k] = -g_0_y_0_xxxxx[k] * ab_z + g_0_y_0_xxxxxz[k];

                g_0_y_z_xxxxy[k] = -g_0_y_0_xxxxy[k] * ab_z + g_0_y_0_xxxxyz[k];

                g_0_y_z_xxxxz[k] = -g_0_y_0_xxxxz[k] * ab_z + g_0_y_0_xxxxzz[k];

                g_0_y_z_xxxyy[k] = -g_0_y_0_xxxyy[k] * ab_z + g_0_y_0_xxxyyz[k];

                g_0_y_z_xxxyz[k] = -g_0_y_0_xxxyz[k] * ab_z + g_0_y_0_xxxyzz[k];

                g_0_y_z_xxxzz[k] = -g_0_y_0_xxxzz[k] * ab_z + g_0_y_0_xxxzzz[k];

                g_0_y_z_xxyyy[k] = -g_0_y_0_xxyyy[k] * ab_z + g_0_y_0_xxyyyz[k];

                g_0_y_z_xxyyz[k] = -g_0_y_0_xxyyz[k] * ab_z + g_0_y_0_xxyyzz[k];

                g_0_y_z_xxyzz[k] = -g_0_y_0_xxyzz[k] * ab_z + g_0_y_0_xxyzzz[k];

                g_0_y_z_xxzzz[k] = -g_0_y_0_xxzzz[k] * ab_z + g_0_y_0_xxzzzz[k];

                g_0_y_z_xyyyy[k] = -g_0_y_0_xyyyy[k] * ab_z + g_0_y_0_xyyyyz[k];

                g_0_y_z_xyyyz[k] = -g_0_y_0_xyyyz[k] * ab_z + g_0_y_0_xyyyzz[k];

                g_0_y_z_xyyzz[k] = -g_0_y_0_xyyzz[k] * ab_z + g_0_y_0_xyyzzz[k];

                g_0_y_z_xyzzz[k] = -g_0_y_0_xyzzz[k] * ab_z + g_0_y_0_xyzzzz[k];

                g_0_y_z_xzzzz[k] = -g_0_y_0_xzzzz[k] * ab_z + g_0_y_0_xzzzzz[k];

                g_0_y_z_yyyyy[k] = -g_0_y_0_yyyyy[k] * ab_z + g_0_y_0_yyyyyz[k];

                g_0_y_z_yyyyz[k] = -g_0_y_0_yyyyz[k] * ab_z + g_0_y_0_yyyyzz[k];

                g_0_y_z_yyyzz[k] = -g_0_y_0_yyyzz[k] * ab_z + g_0_y_0_yyyzzz[k];

                g_0_y_z_yyzzz[k] = -g_0_y_0_yyzzz[k] * ab_z + g_0_y_0_yyzzzz[k];

                g_0_y_z_yzzzz[k] = -g_0_y_0_yzzzz[k] * ab_z + g_0_y_0_yzzzzz[k];

                g_0_y_z_zzzzz[k] = -g_0_y_0_zzzzz[k] * ab_z + g_0_y_0_zzzzzz[k];
            }

            /// Set up 126-147 components of targeted buffer : cbuffer.data(

            auto g_0_z_x_xxxxx = cbuffer.data(ph_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_z_x_xxxxy = cbuffer.data(ph_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_z_x_xxxxz = cbuffer.data(ph_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_z_x_xxxyy = cbuffer.data(ph_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_z_x_xxxyz = cbuffer.data(ph_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_z_x_xxxzz = cbuffer.data(ph_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_z_x_xxyyy = cbuffer.data(ph_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_z_x_xxyyz = cbuffer.data(ph_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_z_x_xxyzz = cbuffer.data(ph_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_z_x_xxzzz = cbuffer.data(ph_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_z_x_xyyyy = cbuffer.data(ph_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_z_x_xyyyz = cbuffer.data(ph_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_z_x_xyyzz = cbuffer.data(ph_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_z_x_xyzzz = cbuffer.data(ph_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_z_x_xzzzz = cbuffer.data(ph_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_z_x_yyyyy = cbuffer.data(ph_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_z_x_yyyyz = cbuffer.data(ph_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_z_x_yyyzz = cbuffer.data(ph_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_z_x_yyzzz = cbuffer.data(ph_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_z_x_yzzzz = cbuffer.data(ph_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_z_x_zzzzz = cbuffer.data(ph_geom_01_off + 146 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_0_xxxxx, g_0_z_0_xxxxxx, g_0_z_0_xxxxxy, g_0_z_0_xxxxxz, g_0_z_0_xxxxy, g_0_z_0_xxxxyy, g_0_z_0_xxxxyz, g_0_z_0_xxxxz, g_0_z_0_xxxxzz, g_0_z_0_xxxyy, g_0_z_0_xxxyyy, g_0_z_0_xxxyyz, g_0_z_0_xxxyz, g_0_z_0_xxxyzz, g_0_z_0_xxxzz, g_0_z_0_xxxzzz, g_0_z_0_xxyyy, g_0_z_0_xxyyyy, g_0_z_0_xxyyyz, g_0_z_0_xxyyz, g_0_z_0_xxyyzz, g_0_z_0_xxyzz, g_0_z_0_xxyzzz, g_0_z_0_xxzzz, g_0_z_0_xxzzzz, g_0_z_0_xyyyy, g_0_z_0_xyyyyy, g_0_z_0_xyyyyz, g_0_z_0_xyyyz, g_0_z_0_xyyyzz, g_0_z_0_xyyzz, g_0_z_0_xyyzzz, g_0_z_0_xyzzz, g_0_z_0_xyzzzz, g_0_z_0_xzzzz, g_0_z_0_xzzzzz, g_0_z_0_yyyyy, g_0_z_0_yyyyz, g_0_z_0_yyyzz, g_0_z_0_yyzzz, g_0_z_0_yzzzz, g_0_z_0_zzzzz, g_0_z_x_xxxxx, g_0_z_x_xxxxy, g_0_z_x_xxxxz, g_0_z_x_xxxyy, g_0_z_x_xxxyz, g_0_z_x_xxxzz, g_0_z_x_xxyyy, g_0_z_x_xxyyz, g_0_z_x_xxyzz, g_0_z_x_xxzzz, g_0_z_x_xyyyy, g_0_z_x_xyyyz, g_0_z_x_xyyzz, g_0_z_x_xyzzz, g_0_z_x_xzzzz, g_0_z_x_yyyyy, g_0_z_x_yyyyz, g_0_z_x_yyyzz, g_0_z_x_yyzzz, g_0_z_x_yzzzz, g_0_z_x_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_x_xxxxx[k] = -g_0_z_0_xxxxx[k] * ab_x + g_0_z_0_xxxxxx[k];

                g_0_z_x_xxxxy[k] = -g_0_z_0_xxxxy[k] * ab_x + g_0_z_0_xxxxxy[k];

                g_0_z_x_xxxxz[k] = -g_0_z_0_xxxxz[k] * ab_x + g_0_z_0_xxxxxz[k];

                g_0_z_x_xxxyy[k] = -g_0_z_0_xxxyy[k] * ab_x + g_0_z_0_xxxxyy[k];

                g_0_z_x_xxxyz[k] = -g_0_z_0_xxxyz[k] * ab_x + g_0_z_0_xxxxyz[k];

                g_0_z_x_xxxzz[k] = -g_0_z_0_xxxzz[k] * ab_x + g_0_z_0_xxxxzz[k];

                g_0_z_x_xxyyy[k] = -g_0_z_0_xxyyy[k] * ab_x + g_0_z_0_xxxyyy[k];

                g_0_z_x_xxyyz[k] = -g_0_z_0_xxyyz[k] * ab_x + g_0_z_0_xxxyyz[k];

                g_0_z_x_xxyzz[k] = -g_0_z_0_xxyzz[k] * ab_x + g_0_z_0_xxxyzz[k];

                g_0_z_x_xxzzz[k] = -g_0_z_0_xxzzz[k] * ab_x + g_0_z_0_xxxzzz[k];

                g_0_z_x_xyyyy[k] = -g_0_z_0_xyyyy[k] * ab_x + g_0_z_0_xxyyyy[k];

                g_0_z_x_xyyyz[k] = -g_0_z_0_xyyyz[k] * ab_x + g_0_z_0_xxyyyz[k];

                g_0_z_x_xyyzz[k] = -g_0_z_0_xyyzz[k] * ab_x + g_0_z_0_xxyyzz[k];

                g_0_z_x_xyzzz[k] = -g_0_z_0_xyzzz[k] * ab_x + g_0_z_0_xxyzzz[k];

                g_0_z_x_xzzzz[k] = -g_0_z_0_xzzzz[k] * ab_x + g_0_z_0_xxzzzz[k];

                g_0_z_x_yyyyy[k] = -g_0_z_0_yyyyy[k] * ab_x + g_0_z_0_xyyyyy[k];

                g_0_z_x_yyyyz[k] = -g_0_z_0_yyyyz[k] * ab_x + g_0_z_0_xyyyyz[k];

                g_0_z_x_yyyzz[k] = -g_0_z_0_yyyzz[k] * ab_x + g_0_z_0_xyyyzz[k];

                g_0_z_x_yyzzz[k] = -g_0_z_0_yyzzz[k] * ab_x + g_0_z_0_xyyzzz[k];

                g_0_z_x_yzzzz[k] = -g_0_z_0_yzzzz[k] * ab_x + g_0_z_0_xyzzzz[k];

                g_0_z_x_zzzzz[k] = -g_0_z_0_zzzzz[k] * ab_x + g_0_z_0_xzzzzz[k];
            }

            /// Set up 147-168 components of targeted buffer : cbuffer.data(

            auto g_0_z_y_xxxxx = cbuffer.data(ph_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_z_y_xxxxy = cbuffer.data(ph_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_z_y_xxxxz = cbuffer.data(ph_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_z_y_xxxyy = cbuffer.data(ph_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_z_y_xxxyz = cbuffer.data(ph_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_z_y_xxxzz = cbuffer.data(ph_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_z_y_xxyyy = cbuffer.data(ph_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_z_y_xxyyz = cbuffer.data(ph_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_z_y_xxyzz = cbuffer.data(ph_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_z_y_xxzzz = cbuffer.data(ph_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_z_y_xyyyy = cbuffer.data(ph_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_z_y_xyyyz = cbuffer.data(ph_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_z_y_xyyzz = cbuffer.data(ph_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_z_y_xyzzz = cbuffer.data(ph_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_z_y_xzzzz = cbuffer.data(ph_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_z_y_yyyyy = cbuffer.data(ph_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_z_y_yyyyz = cbuffer.data(ph_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_z_y_yyyzz = cbuffer.data(ph_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_z_y_yyzzz = cbuffer.data(ph_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_z_y_yzzzz = cbuffer.data(ph_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_z_y_zzzzz = cbuffer.data(ph_geom_01_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_0_xxxxx, g_0_z_0_xxxxxy, g_0_z_0_xxxxy, g_0_z_0_xxxxyy, g_0_z_0_xxxxyz, g_0_z_0_xxxxz, g_0_z_0_xxxyy, g_0_z_0_xxxyyy, g_0_z_0_xxxyyz, g_0_z_0_xxxyz, g_0_z_0_xxxyzz, g_0_z_0_xxxzz, g_0_z_0_xxyyy, g_0_z_0_xxyyyy, g_0_z_0_xxyyyz, g_0_z_0_xxyyz, g_0_z_0_xxyyzz, g_0_z_0_xxyzz, g_0_z_0_xxyzzz, g_0_z_0_xxzzz, g_0_z_0_xyyyy, g_0_z_0_xyyyyy, g_0_z_0_xyyyyz, g_0_z_0_xyyyz, g_0_z_0_xyyyzz, g_0_z_0_xyyzz, g_0_z_0_xyyzzz, g_0_z_0_xyzzz, g_0_z_0_xyzzzz, g_0_z_0_xzzzz, g_0_z_0_yyyyy, g_0_z_0_yyyyyy, g_0_z_0_yyyyyz, g_0_z_0_yyyyz, g_0_z_0_yyyyzz, g_0_z_0_yyyzz, g_0_z_0_yyyzzz, g_0_z_0_yyzzz, g_0_z_0_yyzzzz, g_0_z_0_yzzzz, g_0_z_0_yzzzzz, g_0_z_0_zzzzz, g_0_z_y_xxxxx, g_0_z_y_xxxxy, g_0_z_y_xxxxz, g_0_z_y_xxxyy, g_0_z_y_xxxyz, g_0_z_y_xxxzz, g_0_z_y_xxyyy, g_0_z_y_xxyyz, g_0_z_y_xxyzz, g_0_z_y_xxzzz, g_0_z_y_xyyyy, g_0_z_y_xyyyz, g_0_z_y_xyyzz, g_0_z_y_xyzzz, g_0_z_y_xzzzz, g_0_z_y_yyyyy, g_0_z_y_yyyyz, g_0_z_y_yyyzz, g_0_z_y_yyzzz, g_0_z_y_yzzzz, g_0_z_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_y_xxxxx[k] = -g_0_z_0_xxxxx[k] * ab_y + g_0_z_0_xxxxxy[k];

                g_0_z_y_xxxxy[k] = -g_0_z_0_xxxxy[k] * ab_y + g_0_z_0_xxxxyy[k];

                g_0_z_y_xxxxz[k] = -g_0_z_0_xxxxz[k] * ab_y + g_0_z_0_xxxxyz[k];

                g_0_z_y_xxxyy[k] = -g_0_z_0_xxxyy[k] * ab_y + g_0_z_0_xxxyyy[k];

                g_0_z_y_xxxyz[k] = -g_0_z_0_xxxyz[k] * ab_y + g_0_z_0_xxxyyz[k];

                g_0_z_y_xxxzz[k] = -g_0_z_0_xxxzz[k] * ab_y + g_0_z_0_xxxyzz[k];

                g_0_z_y_xxyyy[k] = -g_0_z_0_xxyyy[k] * ab_y + g_0_z_0_xxyyyy[k];

                g_0_z_y_xxyyz[k] = -g_0_z_0_xxyyz[k] * ab_y + g_0_z_0_xxyyyz[k];

                g_0_z_y_xxyzz[k] = -g_0_z_0_xxyzz[k] * ab_y + g_0_z_0_xxyyzz[k];

                g_0_z_y_xxzzz[k] = -g_0_z_0_xxzzz[k] * ab_y + g_0_z_0_xxyzzz[k];

                g_0_z_y_xyyyy[k] = -g_0_z_0_xyyyy[k] * ab_y + g_0_z_0_xyyyyy[k];

                g_0_z_y_xyyyz[k] = -g_0_z_0_xyyyz[k] * ab_y + g_0_z_0_xyyyyz[k];

                g_0_z_y_xyyzz[k] = -g_0_z_0_xyyzz[k] * ab_y + g_0_z_0_xyyyzz[k];

                g_0_z_y_xyzzz[k] = -g_0_z_0_xyzzz[k] * ab_y + g_0_z_0_xyyzzz[k];

                g_0_z_y_xzzzz[k] = -g_0_z_0_xzzzz[k] * ab_y + g_0_z_0_xyzzzz[k];

                g_0_z_y_yyyyy[k] = -g_0_z_0_yyyyy[k] * ab_y + g_0_z_0_yyyyyy[k];

                g_0_z_y_yyyyz[k] = -g_0_z_0_yyyyz[k] * ab_y + g_0_z_0_yyyyyz[k];

                g_0_z_y_yyyzz[k] = -g_0_z_0_yyyzz[k] * ab_y + g_0_z_0_yyyyzz[k];

                g_0_z_y_yyzzz[k] = -g_0_z_0_yyzzz[k] * ab_y + g_0_z_0_yyyzzz[k];

                g_0_z_y_yzzzz[k] = -g_0_z_0_yzzzz[k] * ab_y + g_0_z_0_yyzzzz[k];

                g_0_z_y_zzzzz[k] = -g_0_z_0_zzzzz[k] * ab_y + g_0_z_0_yzzzzz[k];
            }

            /// Set up 168-189 components of targeted buffer : cbuffer.data(

            auto g_0_z_z_xxxxx = cbuffer.data(ph_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_z_z_xxxxy = cbuffer.data(ph_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_z_z_xxxxz = cbuffer.data(ph_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_z_z_xxxyy = cbuffer.data(ph_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_z_z_xxxyz = cbuffer.data(ph_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_z_z_xxxzz = cbuffer.data(ph_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_z_z_xxyyy = cbuffer.data(ph_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_z_z_xxyyz = cbuffer.data(ph_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_z_z_xxyzz = cbuffer.data(ph_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_z_z_xxzzz = cbuffer.data(ph_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_z_z_xyyyy = cbuffer.data(ph_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_z_z_xyyyz = cbuffer.data(ph_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_z_z_xyyzz = cbuffer.data(ph_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_z_z_xyzzz = cbuffer.data(ph_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_z_z_xzzzz = cbuffer.data(ph_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_z_z_yyyyy = cbuffer.data(ph_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_z_z_yyyyz = cbuffer.data(ph_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_z_z_yyyzz = cbuffer.data(ph_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_z_z_yyzzz = cbuffer.data(ph_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_z_z_yzzzz = cbuffer.data(ph_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_z_z_zzzzz = cbuffer.data(ph_geom_01_off + 188 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_xxxxx, g_0_xxxxy, g_0_xxxxz, g_0_xxxyy, g_0_xxxyz, g_0_xxxzz, g_0_xxyyy, g_0_xxyyz, g_0_xxyzz, g_0_xxzzz, g_0_xyyyy, g_0_xyyyz, g_0_xyyzz, g_0_xyzzz, g_0_xzzzz, g_0_yyyyy, g_0_yyyyz, g_0_yyyzz, g_0_yyzzz, g_0_yzzzz, g_0_z_0_xxxxx, g_0_z_0_xxxxxz, g_0_z_0_xxxxy, g_0_z_0_xxxxyz, g_0_z_0_xxxxz, g_0_z_0_xxxxzz, g_0_z_0_xxxyy, g_0_z_0_xxxyyz, g_0_z_0_xxxyz, g_0_z_0_xxxyzz, g_0_z_0_xxxzz, g_0_z_0_xxxzzz, g_0_z_0_xxyyy, g_0_z_0_xxyyyz, g_0_z_0_xxyyz, g_0_z_0_xxyyzz, g_0_z_0_xxyzz, g_0_z_0_xxyzzz, g_0_z_0_xxzzz, g_0_z_0_xxzzzz, g_0_z_0_xyyyy, g_0_z_0_xyyyyz, g_0_z_0_xyyyz, g_0_z_0_xyyyzz, g_0_z_0_xyyzz, g_0_z_0_xyyzzz, g_0_z_0_xyzzz, g_0_z_0_xyzzzz, g_0_z_0_xzzzz, g_0_z_0_xzzzzz, g_0_z_0_yyyyy, g_0_z_0_yyyyyz, g_0_z_0_yyyyz, g_0_z_0_yyyyzz, g_0_z_0_yyyzz, g_0_z_0_yyyzzz, g_0_z_0_yyzzz, g_0_z_0_yyzzzz, g_0_z_0_yzzzz, g_0_z_0_yzzzzz, g_0_z_0_zzzzz, g_0_z_0_zzzzzz, g_0_z_z_xxxxx, g_0_z_z_xxxxy, g_0_z_z_xxxxz, g_0_z_z_xxxyy, g_0_z_z_xxxyz, g_0_z_z_xxxzz, g_0_z_z_xxyyy, g_0_z_z_xxyyz, g_0_z_z_xxyzz, g_0_z_z_xxzzz, g_0_z_z_xyyyy, g_0_z_z_xyyyz, g_0_z_z_xyyzz, g_0_z_z_xyzzz, g_0_z_z_xzzzz, g_0_z_z_yyyyy, g_0_z_z_yyyyz, g_0_z_z_yyyzz, g_0_z_z_yyzzz, g_0_z_z_yzzzz, g_0_z_z_zzzzz, g_0_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_z_xxxxx[k] = g_0_xxxxx[k] - g_0_z_0_xxxxx[k] * ab_z + g_0_z_0_xxxxxz[k];

                g_0_z_z_xxxxy[k] = g_0_xxxxy[k] - g_0_z_0_xxxxy[k] * ab_z + g_0_z_0_xxxxyz[k];

                g_0_z_z_xxxxz[k] = g_0_xxxxz[k] - g_0_z_0_xxxxz[k] * ab_z + g_0_z_0_xxxxzz[k];

                g_0_z_z_xxxyy[k] = g_0_xxxyy[k] - g_0_z_0_xxxyy[k] * ab_z + g_0_z_0_xxxyyz[k];

                g_0_z_z_xxxyz[k] = g_0_xxxyz[k] - g_0_z_0_xxxyz[k] * ab_z + g_0_z_0_xxxyzz[k];

                g_0_z_z_xxxzz[k] = g_0_xxxzz[k] - g_0_z_0_xxxzz[k] * ab_z + g_0_z_0_xxxzzz[k];

                g_0_z_z_xxyyy[k] = g_0_xxyyy[k] - g_0_z_0_xxyyy[k] * ab_z + g_0_z_0_xxyyyz[k];

                g_0_z_z_xxyyz[k] = g_0_xxyyz[k] - g_0_z_0_xxyyz[k] * ab_z + g_0_z_0_xxyyzz[k];

                g_0_z_z_xxyzz[k] = g_0_xxyzz[k] - g_0_z_0_xxyzz[k] * ab_z + g_0_z_0_xxyzzz[k];

                g_0_z_z_xxzzz[k] = g_0_xxzzz[k] - g_0_z_0_xxzzz[k] * ab_z + g_0_z_0_xxzzzz[k];

                g_0_z_z_xyyyy[k] = g_0_xyyyy[k] - g_0_z_0_xyyyy[k] * ab_z + g_0_z_0_xyyyyz[k];

                g_0_z_z_xyyyz[k] = g_0_xyyyz[k] - g_0_z_0_xyyyz[k] * ab_z + g_0_z_0_xyyyzz[k];

                g_0_z_z_xyyzz[k] = g_0_xyyzz[k] - g_0_z_0_xyyzz[k] * ab_z + g_0_z_0_xyyzzz[k];

                g_0_z_z_xyzzz[k] = g_0_xyzzz[k] - g_0_z_0_xyzzz[k] * ab_z + g_0_z_0_xyzzzz[k];

                g_0_z_z_xzzzz[k] = g_0_xzzzz[k] - g_0_z_0_xzzzz[k] * ab_z + g_0_z_0_xzzzzz[k];

                g_0_z_z_yyyyy[k] = g_0_yyyyy[k] - g_0_z_0_yyyyy[k] * ab_z + g_0_z_0_yyyyyz[k];

                g_0_z_z_yyyyz[k] = g_0_yyyyz[k] - g_0_z_0_yyyyz[k] * ab_z + g_0_z_0_yyyyzz[k];

                g_0_z_z_yyyzz[k] = g_0_yyyzz[k] - g_0_z_0_yyyzz[k] * ab_z + g_0_z_0_yyyzzz[k];

                g_0_z_z_yyzzz[k] = g_0_yyzzz[k] - g_0_z_0_yyzzz[k] * ab_z + g_0_z_0_yyzzzz[k];

                g_0_z_z_yzzzz[k] = g_0_yzzzz[k] - g_0_z_0_yzzzz[k] * ab_z + g_0_z_0_yzzzzz[k];

                g_0_z_z_zzzzz[k] = g_0_zzzzz[k] - g_0_z_0_zzzzz[k] * ab_z + g_0_z_0_zzzzzz[k];
            }
        }
    }
}

} // erirec namespace

