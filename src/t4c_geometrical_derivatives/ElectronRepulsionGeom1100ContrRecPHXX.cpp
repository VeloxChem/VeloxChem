#include "ElectronRepulsionGeom1100ContrRecPHXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom11_hrr_electron_repulsion_phxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_11_phxx,
                                            const size_t idx_geom_01_shxx,
                                            const size_t idx_geom_10_shxx,
                                            const size_t idx_geom_11_shxx,
                                            const size_t idx_geom_11_sixx,
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

            /// Set up components of auxilary buffer : SHSS

            const auto sh_geom_11_off = idx_geom_11_shxx + i * dcomps + j;

            auto g_x_x_0_xxxxx = cbuffer.data(sh_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_0_xxxxy = cbuffer.data(sh_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_0_xxxxz = cbuffer.data(sh_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_0_xxxyy = cbuffer.data(sh_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_0_xxxyz = cbuffer.data(sh_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_0_xxxzz = cbuffer.data(sh_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_0_xxyyy = cbuffer.data(sh_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_0_xxyyz = cbuffer.data(sh_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_0_xxyzz = cbuffer.data(sh_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_0_xxzzz = cbuffer.data(sh_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_0_xyyyy = cbuffer.data(sh_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_0_xyyyz = cbuffer.data(sh_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_0_xyyzz = cbuffer.data(sh_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_0_xyzzz = cbuffer.data(sh_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_0_xzzzz = cbuffer.data(sh_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_0_yyyyy = cbuffer.data(sh_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_0_yyyyz = cbuffer.data(sh_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_0_yyyzz = cbuffer.data(sh_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_x_0_yyzzz = cbuffer.data(sh_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_0_yzzzz = cbuffer.data(sh_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_x_0_zzzzz = cbuffer.data(sh_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_y_0_xxxxx = cbuffer.data(sh_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_y_0_xxxxy = cbuffer.data(sh_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_y_0_xxxxz = cbuffer.data(sh_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_y_0_xxxyy = cbuffer.data(sh_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_y_0_xxxyz = cbuffer.data(sh_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_y_0_xxxzz = cbuffer.data(sh_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_y_0_xxyyy = cbuffer.data(sh_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_y_0_xxyyz = cbuffer.data(sh_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_y_0_xxyzz = cbuffer.data(sh_geom_11_off + 29 * ccomps * dcomps);

            auto g_x_y_0_xxzzz = cbuffer.data(sh_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_y_0_xyyyy = cbuffer.data(sh_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_y_0_xyyyz = cbuffer.data(sh_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_y_0_xyyzz = cbuffer.data(sh_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_y_0_xyzzz = cbuffer.data(sh_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_y_0_xzzzz = cbuffer.data(sh_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_y_0_yyyyy = cbuffer.data(sh_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_y_0_yyyyz = cbuffer.data(sh_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_y_0_yyyzz = cbuffer.data(sh_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_y_0_yyzzz = cbuffer.data(sh_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_y_0_yzzzz = cbuffer.data(sh_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_y_0_zzzzz = cbuffer.data(sh_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_z_0_xxxxx = cbuffer.data(sh_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_z_0_xxxxy = cbuffer.data(sh_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_z_0_xxxxz = cbuffer.data(sh_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_z_0_xxxyy = cbuffer.data(sh_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_z_0_xxxyz = cbuffer.data(sh_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_z_0_xxxzz = cbuffer.data(sh_geom_11_off + 47 * ccomps * dcomps);

            auto g_x_z_0_xxyyy = cbuffer.data(sh_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_z_0_xxyyz = cbuffer.data(sh_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_z_0_xxyzz = cbuffer.data(sh_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_z_0_xxzzz = cbuffer.data(sh_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_z_0_xyyyy = cbuffer.data(sh_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_z_0_xyyyz = cbuffer.data(sh_geom_11_off + 53 * ccomps * dcomps);

            auto g_x_z_0_xyyzz = cbuffer.data(sh_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_z_0_xyzzz = cbuffer.data(sh_geom_11_off + 55 * ccomps * dcomps);

            auto g_x_z_0_xzzzz = cbuffer.data(sh_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_z_0_yyyyy = cbuffer.data(sh_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_z_0_yyyyz = cbuffer.data(sh_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_z_0_yyyzz = cbuffer.data(sh_geom_11_off + 59 * ccomps * dcomps);

            auto g_x_z_0_yyzzz = cbuffer.data(sh_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_z_0_yzzzz = cbuffer.data(sh_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_z_0_zzzzz = cbuffer.data(sh_geom_11_off + 62 * ccomps * dcomps);

            auto g_y_x_0_xxxxx = cbuffer.data(sh_geom_11_off + 63 * ccomps * dcomps);

            auto g_y_x_0_xxxxy = cbuffer.data(sh_geom_11_off + 64 * ccomps * dcomps);

            auto g_y_x_0_xxxxz = cbuffer.data(sh_geom_11_off + 65 * ccomps * dcomps);

            auto g_y_x_0_xxxyy = cbuffer.data(sh_geom_11_off + 66 * ccomps * dcomps);

            auto g_y_x_0_xxxyz = cbuffer.data(sh_geom_11_off + 67 * ccomps * dcomps);

            auto g_y_x_0_xxxzz = cbuffer.data(sh_geom_11_off + 68 * ccomps * dcomps);

            auto g_y_x_0_xxyyy = cbuffer.data(sh_geom_11_off + 69 * ccomps * dcomps);

            auto g_y_x_0_xxyyz = cbuffer.data(sh_geom_11_off + 70 * ccomps * dcomps);

            auto g_y_x_0_xxyzz = cbuffer.data(sh_geom_11_off + 71 * ccomps * dcomps);

            auto g_y_x_0_xxzzz = cbuffer.data(sh_geom_11_off + 72 * ccomps * dcomps);

            auto g_y_x_0_xyyyy = cbuffer.data(sh_geom_11_off + 73 * ccomps * dcomps);

            auto g_y_x_0_xyyyz = cbuffer.data(sh_geom_11_off + 74 * ccomps * dcomps);

            auto g_y_x_0_xyyzz = cbuffer.data(sh_geom_11_off + 75 * ccomps * dcomps);

            auto g_y_x_0_xyzzz = cbuffer.data(sh_geom_11_off + 76 * ccomps * dcomps);

            auto g_y_x_0_xzzzz = cbuffer.data(sh_geom_11_off + 77 * ccomps * dcomps);

            auto g_y_x_0_yyyyy = cbuffer.data(sh_geom_11_off + 78 * ccomps * dcomps);

            auto g_y_x_0_yyyyz = cbuffer.data(sh_geom_11_off + 79 * ccomps * dcomps);

            auto g_y_x_0_yyyzz = cbuffer.data(sh_geom_11_off + 80 * ccomps * dcomps);

            auto g_y_x_0_yyzzz = cbuffer.data(sh_geom_11_off + 81 * ccomps * dcomps);

            auto g_y_x_0_yzzzz = cbuffer.data(sh_geom_11_off + 82 * ccomps * dcomps);

            auto g_y_x_0_zzzzz = cbuffer.data(sh_geom_11_off + 83 * ccomps * dcomps);

            auto g_y_y_0_xxxxx = cbuffer.data(sh_geom_11_off + 84 * ccomps * dcomps);

            auto g_y_y_0_xxxxy = cbuffer.data(sh_geom_11_off + 85 * ccomps * dcomps);

            auto g_y_y_0_xxxxz = cbuffer.data(sh_geom_11_off + 86 * ccomps * dcomps);

            auto g_y_y_0_xxxyy = cbuffer.data(sh_geom_11_off + 87 * ccomps * dcomps);

            auto g_y_y_0_xxxyz = cbuffer.data(sh_geom_11_off + 88 * ccomps * dcomps);

            auto g_y_y_0_xxxzz = cbuffer.data(sh_geom_11_off + 89 * ccomps * dcomps);

            auto g_y_y_0_xxyyy = cbuffer.data(sh_geom_11_off + 90 * ccomps * dcomps);

            auto g_y_y_0_xxyyz = cbuffer.data(sh_geom_11_off + 91 * ccomps * dcomps);

            auto g_y_y_0_xxyzz = cbuffer.data(sh_geom_11_off + 92 * ccomps * dcomps);

            auto g_y_y_0_xxzzz = cbuffer.data(sh_geom_11_off + 93 * ccomps * dcomps);

            auto g_y_y_0_xyyyy = cbuffer.data(sh_geom_11_off + 94 * ccomps * dcomps);

            auto g_y_y_0_xyyyz = cbuffer.data(sh_geom_11_off + 95 * ccomps * dcomps);

            auto g_y_y_0_xyyzz = cbuffer.data(sh_geom_11_off + 96 * ccomps * dcomps);

            auto g_y_y_0_xyzzz = cbuffer.data(sh_geom_11_off + 97 * ccomps * dcomps);

            auto g_y_y_0_xzzzz = cbuffer.data(sh_geom_11_off + 98 * ccomps * dcomps);

            auto g_y_y_0_yyyyy = cbuffer.data(sh_geom_11_off + 99 * ccomps * dcomps);

            auto g_y_y_0_yyyyz = cbuffer.data(sh_geom_11_off + 100 * ccomps * dcomps);

            auto g_y_y_0_yyyzz = cbuffer.data(sh_geom_11_off + 101 * ccomps * dcomps);

            auto g_y_y_0_yyzzz = cbuffer.data(sh_geom_11_off + 102 * ccomps * dcomps);

            auto g_y_y_0_yzzzz = cbuffer.data(sh_geom_11_off + 103 * ccomps * dcomps);

            auto g_y_y_0_zzzzz = cbuffer.data(sh_geom_11_off + 104 * ccomps * dcomps);

            auto g_y_z_0_xxxxx = cbuffer.data(sh_geom_11_off + 105 * ccomps * dcomps);

            auto g_y_z_0_xxxxy = cbuffer.data(sh_geom_11_off + 106 * ccomps * dcomps);

            auto g_y_z_0_xxxxz = cbuffer.data(sh_geom_11_off + 107 * ccomps * dcomps);

            auto g_y_z_0_xxxyy = cbuffer.data(sh_geom_11_off + 108 * ccomps * dcomps);

            auto g_y_z_0_xxxyz = cbuffer.data(sh_geom_11_off + 109 * ccomps * dcomps);

            auto g_y_z_0_xxxzz = cbuffer.data(sh_geom_11_off + 110 * ccomps * dcomps);

            auto g_y_z_0_xxyyy = cbuffer.data(sh_geom_11_off + 111 * ccomps * dcomps);

            auto g_y_z_0_xxyyz = cbuffer.data(sh_geom_11_off + 112 * ccomps * dcomps);

            auto g_y_z_0_xxyzz = cbuffer.data(sh_geom_11_off + 113 * ccomps * dcomps);

            auto g_y_z_0_xxzzz = cbuffer.data(sh_geom_11_off + 114 * ccomps * dcomps);

            auto g_y_z_0_xyyyy = cbuffer.data(sh_geom_11_off + 115 * ccomps * dcomps);

            auto g_y_z_0_xyyyz = cbuffer.data(sh_geom_11_off + 116 * ccomps * dcomps);

            auto g_y_z_0_xyyzz = cbuffer.data(sh_geom_11_off + 117 * ccomps * dcomps);

            auto g_y_z_0_xyzzz = cbuffer.data(sh_geom_11_off + 118 * ccomps * dcomps);

            auto g_y_z_0_xzzzz = cbuffer.data(sh_geom_11_off + 119 * ccomps * dcomps);

            auto g_y_z_0_yyyyy = cbuffer.data(sh_geom_11_off + 120 * ccomps * dcomps);

            auto g_y_z_0_yyyyz = cbuffer.data(sh_geom_11_off + 121 * ccomps * dcomps);

            auto g_y_z_0_yyyzz = cbuffer.data(sh_geom_11_off + 122 * ccomps * dcomps);

            auto g_y_z_0_yyzzz = cbuffer.data(sh_geom_11_off + 123 * ccomps * dcomps);

            auto g_y_z_0_yzzzz = cbuffer.data(sh_geom_11_off + 124 * ccomps * dcomps);

            auto g_y_z_0_zzzzz = cbuffer.data(sh_geom_11_off + 125 * ccomps * dcomps);

            auto g_z_x_0_xxxxx = cbuffer.data(sh_geom_11_off + 126 * ccomps * dcomps);

            auto g_z_x_0_xxxxy = cbuffer.data(sh_geom_11_off + 127 * ccomps * dcomps);

            auto g_z_x_0_xxxxz = cbuffer.data(sh_geom_11_off + 128 * ccomps * dcomps);

            auto g_z_x_0_xxxyy = cbuffer.data(sh_geom_11_off + 129 * ccomps * dcomps);

            auto g_z_x_0_xxxyz = cbuffer.data(sh_geom_11_off + 130 * ccomps * dcomps);

            auto g_z_x_0_xxxzz = cbuffer.data(sh_geom_11_off + 131 * ccomps * dcomps);

            auto g_z_x_0_xxyyy = cbuffer.data(sh_geom_11_off + 132 * ccomps * dcomps);

            auto g_z_x_0_xxyyz = cbuffer.data(sh_geom_11_off + 133 * ccomps * dcomps);

            auto g_z_x_0_xxyzz = cbuffer.data(sh_geom_11_off + 134 * ccomps * dcomps);

            auto g_z_x_0_xxzzz = cbuffer.data(sh_geom_11_off + 135 * ccomps * dcomps);

            auto g_z_x_0_xyyyy = cbuffer.data(sh_geom_11_off + 136 * ccomps * dcomps);

            auto g_z_x_0_xyyyz = cbuffer.data(sh_geom_11_off + 137 * ccomps * dcomps);

            auto g_z_x_0_xyyzz = cbuffer.data(sh_geom_11_off + 138 * ccomps * dcomps);

            auto g_z_x_0_xyzzz = cbuffer.data(sh_geom_11_off + 139 * ccomps * dcomps);

            auto g_z_x_0_xzzzz = cbuffer.data(sh_geom_11_off + 140 * ccomps * dcomps);

            auto g_z_x_0_yyyyy = cbuffer.data(sh_geom_11_off + 141 * ccomps * dcomps);

            auto g_z_x_0_yyyyz = cbuffer.data(sh_geom_11_off + 142 * ccomps * dcomps);

            auto g_z_x_0_yyyzz = cbuffer.data(sh_geom_11_off + 143 * ccomps * dcomps);

            auto g_z_x_0_yyzzz = cbuffer.data(sh_geom_11_off + 144 * ccomps * dcomps);

            auto g_z_x_0_yzzzz = cbuffer.data(sh_geom_11_off + 145 * ccomps * dcomps);

            auto g_z_x_0_zzzzz = cbuffer.data(sh_geom_11_off + 146 * ccomps * dcomps);

            auto g_z_y_0_xxxxx = cbuffer.data(sh_geom_11_off + 147 * ccomps * dcomps);

            auto g_z_y_0_xxxxy = cbuffer.data(sh_geom_11_off + 148 * ccomps * dcomps);

            auto g_z_y_0_xxxxz = cbuffer.data(sh_geom_11_off + 149 * ccomps * dcomps);

            auto g_z_y_0_xxxyy = cbuffer.data(sh_geom_11_off + 150 * ccomps * dcomps);

            auto g_z_y_0_xxxyz = cbuffer.data(sh_geom_11_off + 151 * ccomps * dcomps);

            auto g_z_y_0_xxxzz = cbuffer.data(sh_geom_11_off + 152 * ccomps * dcomps);

            auto g_z_y_0_xxyyy = cbuffer.data(sh_geom_11_off + 153 * ccomps * dcomps);

            auto g_z_y_0_xxyyz = cbuffer.data(sh_geom_11_off + 154 * ccomps * dcomps);

            auto g_z_y_0_xxyzz = cbuffer.data(sh_geom_11_off + 155 * ccomps * dcomps);

            auto g_z_y_0_xxzzz = cbuffer.data(sh_geom_11_off + 156 * ccomps * dcomps);

            auto g_z_y_0_xyyyy = cbuffer.data(sh_geom_11_off + 157 * ccomps * dcomps);

            auto g_z_y_0_xyyyz = cbuffer.data(sh_geom_11_off + 158 * ccomps * dcomps);

            auto g_z_y_0_xyyzz = cbuffer.data(sh_geom_11_off + 159 * ccomps * dcomps);

            auto g_z_y_0_xyzzz = cbuffer.data(sh_geom_11_off + 160 * ccomps * dcomps);

            auto g_z_y_0_xzzzz = cbuffer.data(sh_geom_11_off + 161 * ccomps * dcomps);

            auto g_z_y_0_yyyyy = cbuffer.data(sh_geom_11_off + 162 * ccomps * dcomps);

            auto g_z_y_0_yyyyz = cbuffer.data(sh_geom_11_off + 163 * ccomps * dcomps);

            auto g_z_y_0_yyyzz = cbuffer.data(sh_geom_11_off + 164 * ccomps * dcomps);

            auto g_z_y_0_yyzzz = cbuffer.data(sh_geom_11_off + 165 * ccomps * dcomps);

            auto g_z_y_0_yzzzz = cbuffer.data(sh_geom_11_off + 166 * ccomps * dcomps);

            auto g_z_y_0_zzzzz = cbuffer.data(sh_geom_11_off + 167 * ccomps * dcomps);

            auto g_z_z_0_xxxxx = cbuffer.data(sh_geom_11_off + 168 * ccomps * dcomps);

            auto g_z_z_0_xxxxy = cbuffer.data(sh_geom_11_off + 169 * ccomps * dcomps);

            auto g_z_z_0_xxxxz = cbuffer.data(sh_geom_11_off + 170 * ccomps * dcomps);

            auto g_z_z_0_xxxyy = cbuffer.data(sh_geom_11_off + 171 * ccomps * dcomps);

            auto g_z_z_0_xxxyz = cbuffer.data(sh_geom_11_off + 172 * ccomps * dcomps);

            auto g_z_z_0_xxxzz = cbuffer.data(sh_geom_11_off + 173 * ccomps * dcomps);

            auto g_z_z_0_xxyyy = cbuffer.data(sh_geom_11_off + 174 * ccomps * dcomps);

            auto g_z_z_0_xxyyz = cbuffer.data(sh_geom_11_off + 175 * ccomps * dcomps);

            auto g_z_z_0_xxyzz = cbuffer.data(sh_geom_11_off + 176 * ccomps * dcomps);

            auto g_z_z_0_xxzzz = cbuffer.data(sh_geom_11_off + 177 * ccomps * dcomps);

            auto g_z_z_0_xyyyy = cbuffer.data(sh_geom_11_off + 178 * ccomps * dcomps);

            auto g_z_z_0_xyyyz = cbuffer.data(sh_geom_11_off + 179 * ccomps * dcomps);

            auto g_z_z_0_xyyzz = cbuffer.data(sh_geom_11_off + 180 * ccomps * dcomps);

            auto g_z_z_0_xyzzz = cbuffer.data(sh_geom_11_off + 181 * ccomps * dcomps);

            auto g_z_z_0_xzzzz = cbuffer.data(sh_geom_11_off + 182 * ccomps * dcomps);

            auto g_z_z_0_yyyyy = cbuffer.data(sh_geom_11_off + 183 * ccomps * dcomps);

            auto g_z_z_0_yyyyz = cbuffer.data(sh_geom_11_off + 184 * ccomps * dcomps);

            auto g_z_z_0_yyyzz = cbuffer.data(sh_geom_11_off + 185 * ccomps * dcomps);

            auto g_z_z_0_yyzzz = cbuffer.data(sh_geom_11_off + 186 * ccomps * dcomps);

            auto g_z_z_0_yzzzz = cbuffer.data(sh_geom_11_off + 187 * ccomps * dcomps);

            auto g_z_z_0_zzzzz = cbuffer.data(sh_geom_11_off + 188 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SISS

            const auto si_geom_11_off = idx_geom_11_sixx + i * dcomps + j;

            auto g_x_x_0_xxxxxx = cbuffer.data(si_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_0_xxxxxy = cbuffer.data(si_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_0_xxxxxz = cbuffer.data(si_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_0_xxxxyy = cbuffer.data(si_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_0_xxxxyz = cbuffer.data(si_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_0_xxxxzz = cbuffer.data(si_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_0_xxxyyy = cbuffer.data(si_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_0_xxxyyz = cbuffer.data(si_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_0_xxxyzz = cbuffer.data(si_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_0_xxxzzz = cbuffer.data(si_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_0_xxyyyy = cbuffer.data(si_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_0_xxyyyz = cbuffer.data(si_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_0_xxyyzz = cbuffer.data(si_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_0_xxyzzz = cbuffer.data(si_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_0_xxzzzz = cbuffer.data(si_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_0_xyyyyy = cbuffer.data(si_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_0_xyyyyz = cbuffer.data(si_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_0_xyyyzz = cbuffer.data(si_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_x_0_xyyzzz = cbuffer.data(si_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_0_xyzzzz = cbuffer.data(si_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_x_0_xzzzzz = cbuffer.data(si_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_x_0_yyyyyy = cbuffer.data(si_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_0_yyyyyz = cbuffer.data(si_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_0_yyyyzz = cbuffer.data(si_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_x_0_yyyzzz = cbuffer.data(si_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_0_yyzzzz = cbuffer.data(si_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_0_yzzzzz = cbuffer.data(si_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_x_0_zzzzzz = cbuffer.data(si_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_y_0_xxxxxx = cbuffer.data(si_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_y_0_xxxxxy = cbuffer.data(si_geom_11_off + 29 * ccomps * dcomps);

            auto g_x_y_0_xxxxxz = cbuffer.data(si_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_y_0_xxxxyy = cbuffer.data(si_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_y_0_xxxxyz = cbuffer.data(si_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_y_0_xxxxzz = cbuffer.data(si_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_y_0_xxxyyy = cbuffer.data(si_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_y_0_xxxyyz = cbuffer.data(si_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_y_0_xxxyzz = cbuffer.data(si_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_y_0_xxxzzz = cbuffer.data(si_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_y_0_xxyyyy = cbuffer.data(si_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_y_0_xxyyyz = cbuffer.data(si_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_y_0_xxyyzz = cbuffer.data(si_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_y_0_xxyzzz = cbuffer.data(si_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_y_0_xxzzzz = cbuffer.data(si_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_y_0_xyyyyy = cbuffer.data(si_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_y_0_xyyyyz = cbuffer.data(si_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_y_0_xyyyzz = cbuffer.data(si_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_y_0_xyyzzz = cbuffer.data(si_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_y_0_xyzzzz = cbuffer.data(si_geom_11_off + 47 * ccomps * dcomps);

            auto g_x_y_0_xzzzzz = cbuffer.data(si_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_y_0_yyyyyy = cbuffer.data(si_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_y_0_yyyyyz = cbuffer.data(si_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_y_0_yyyyzz = cbuffer.data(si_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_y_0_yyyzzz = cbuffer.data(si_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_y_0_yyzzzz = cbuffer.data(si_geom_11_off + 53 * ccomps * dcomps);

            auto g_x_y_0_yzzzzz = cbuffer.data(si_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_y_0_zzzzzz = cbuffer.data(si_geom_11_off + 55 * ccomps * dcomps);

            auto g_x_z_0_xxxxxx = cbuffer.data(si_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_z_0_xxxxxy = cbuffer.data(si_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_z_0_xxxxxz = cbuffer.data(si_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_z_0_xxxxyy = cbuffer.data(si_geom_11_off + 59 * ccomps * dcomps);

            auto g_x_z_0_xxxxyz = cbuffer.data(si_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_z_0_xxxxzz = cbuffer.data(si_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_z_0_xxxyyy = cbuffer.data(si_geom_11_off + 62 * ccomps * dcomps);

            auto g_x_z_0_xxxyyz = cbuffer.data(si_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_z_0_xxxyzz = cbuffer.data(si_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_z_0_xxxzzz = cbuffer.data(si_geom_11_off + 65 * ccomps * dcomps);

            auto g_x_z_0_xxyyyy = cbuffer.data(si_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_z_0_xxyyyz = cbuffer.data(si_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_z_0_xxyyzz = cbuffer.data(si_geom_11_off + 68 * ccomps * dcomps);

            auto g_x_z_0_xxyzzz = cbuffer.data(si_geom_11_off + 69 * ccomps * dcomps);

            auto g_x_z_0_xxzzzz = cbuffer.data(si_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_z_0_xyyyyy = cbuffer.data(si_geom_11_off + 71 * ccomps * dcomps);

            auto g_x_z_0_xyyyyz = cbuffer.data(si_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_z_0_xyyyzz = cbuffer.data(si_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_z_0_xyyzzz = cbuffer.data(si_geom_11_off + 74 * ccomps * dcomps);

            auto g_x_z_0_xyzzzz = cbuffer.data(si_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_z_0_xzzzzz = cbuffer.data(si_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_z_0_yyyyyy = cbuffer.data(si_geom_11_off + 77 * ccomps * dcomps);

            auto g_x_z_0_yyyyyz = cbuffer.data(si_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_z_0_yyyyzz = cbuffer.data(si_geom_11_off + 79 * ccomps * dcomps);

            auto g_x_z_0_yyyzzz = cbuffer.data(si_geom_11_off + 80 * ccomps * dcomps);

            auto g_x_z_0_yyzzzz = cbuffer.data(si_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_z_0_yzzzzz = cbuffer.data(si_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_z_0_zzzzzz = cbuffer.data(si_geom_11_off + 83 * ccomps * dcomps);

            auto g_y_x_0_xxxxxx = cbuffer.data(si_geom_11_off + 84 * ccomps * dcomps);

            auto g_y_x_0_xxxxxy = cbuffer.data(si_geom_11_off + 85 * ccomps * dcomps);

            auto g_y_x_0_xxxxxz = cbuffer.data(si_geom_11_off + 86 * ccomps * dcomps);

            auto g_y_x_0_xxxxyy = cbuffer.data(si_geom_11_off + 87 * ccomps * dcomps);

            auto g_y_x_0_xxxxyz = cbuffer.data(si_geom_11_off + 88 * ccomps * dcomps);

            auto g_y_x_0_xxxxzz = cbuffer.data(si_geom_11_off + 89 * ccomps * dcomps);

            auto g_y_x_0_xxxyyy = cbuffer.data(si_geom_11_off + 90 * ccomps * dcomps);

            auto g_y_x_0_xxxyyz = cbuffer.data(si_geom_11_off + 91 * ccomps * dcomps);

            auto g_y_x_0_xxxyzz = cbuffer.data(si_geom_11_off + 92 * ccomps * dcomps);

            auto g_y_x_0_xxxzzz = cbuffer.data(si_geom_11_off + 93 * ccomps * dcomps);

            auto g_y_x_0_xxyyyy = cbuffer.data(si_geom_11_off + 94 * ccomps * dcomps);

            auto g_y_x_0_xxyyyz = cbuffer.data(si_geom_11_off + 95 * ccomps * dcomps);

            auto g_y_x_0_xxyyzz = cbuffer.data(si_geom_11_off + 96 * ccomps * dcomps);

            auto g_y_x_0_xxyzzz = cbuffer.data(si_geom_11_off + 97 * ccomps * dcomps);

            auto g_y_x_0_xxzzzz = cbuffer.data(si_geom_11_off + 98 * ccomps * dcomps);

            auto g_y_x_0_xyyyyy = cbuffer.data(si_geom_11_off + 99 * ccomps * dcomps);

            auto g_y_x_0_xyyyyz = cbuffer.data(si_geom_11_off + 100 * ccomps * dcomps);

            auto g_y_x_0_xyyyzz = cbuffer.data(si_geom_11_off + 101 * ccomps * dcomps);

            auto g_y_x_0_xyyzzz = cbuffer.data(si_geom_11_off + 102 * ccomps * dcomps);

            auto g_y_x_0_xyzzzz = cbuffer.data(si_geom_11_off + 103 * ccomps * dcomps);

            auto g_y_x_0_xzzzzz = cbuffer.data(si_geom_11_off + 104 * ccomps * dcomps);

            auto g_y_x_0_yyyyyy = cbuffer.data(si_geom_11_off + 105 * ccomps * dcomps);

            auto g_y_x_0_yyyyyz = cbuffer.data(si_geom_11_off + 106 * ccomps * dcomps);

            auto g_y_x_0_yyyyzz = cbuffer.data(si_geom_11_off + 107 * ccomps * dcomps);

            auto g_y_x_0_yyyzzz = cbuffer.data(si_geom_11_off + 108 * ccomps * dcomps);

            auto g_y_x_0_yyzzzz = cbuffer.data(si_geom_11_off + 109 * ccomps * dcomps);

            auto g_y_x_0_yzzzzz = cbuffer.data(si_geom_11_off + 110 * ccomps * dcomps);

            auto g_y_x_0_zzzzzz = cbuffer.data(si_geom_11_off + 111 * ccomps * dcomps);

            auto g_y_y_0_xxxxxx = cbuffer.data(si_geom_11_off + 112 * ccomps * dcomps);

            auto g_y_y_0_xxxxxy = cbuffer.data(si_geom_11_off + 113 * ccomps * dcomps);

            auto g_y_y_0_xxxxxz = cbuffer.data(si_geom_11_off + 114 * ccomps * dcomps);

            auto g_y_y_0_xxxxyy = cbuffer.data(si_geom_11_off + 115 * ccomps * dcomps);

            auto g_y_y_0_xxxxyz = cbuffer.data(si_geom_11_off + 116 * ccomps * dcomps);

            auto g_y_y_0_xxxxzz = cbuffer.data(si_geom_11_off + 117 * ccomps * dcomps);

            auto g_y_y_0_xxxyyy = cbuffer.data(si_geom_11_off + 118 * ccomps * dcomps);

            auto g_y_y_0_xxxyyz = cbuffer.data(si_geom_11_off + 119 * ccomps * dcomps);

            auto g_y_y_0_xxxyzz = cbuffer.data(si_geom_11_off + 120 * ccomps * dcomps);

            auto g_y_y_0_xxxzzz = cbuffer.data(si_geom_11_off + 121 * ccomps * dcomps);

            auto g_y_y_0_xxyyyy = cbuffer.data(si_geom_11_off + 122 * ccomps * dcomps);

            auto g_y_y_0_xxyyyz = cbuffer.data(si_geom_11_off + 123 * ccomps * dcomps);

            auto g_y_y_0_xxyyzz = cbuffer.data(si_geom_11_off + 124 * ccomps * dcomps);

            auto g_y_y_0_xxyzzz = cbuffer.data(si_geom_11_off + 125 * ccomps * dcomps);

            auto g_y_y_0_xxzzzz = cbuffer.data(si_geom_11_off + 126 * ccomps * dcomps);

            auto g_y_y_0_xyyyyy = cbuffer.data(si_geom_11_off + 127 * ccomps * dcomps);

            auto g_y_y_0_xyyyyz = cbuffer.data(si_geom_11_off + 128 * ccomps * dcomps);

            auto g_y_y_0_xyyyzz = cbuffer.data(si_geom_11_off + 129 * ccomps * dcomps);

            auto g_y_y_0_xyyzzz = cbuffer.data(si_geom_11_off + 130 * ccomps * dcomps);

            auto g_y_y_0_xyzzzz = cbuffer.data(si_geom_11_off + 131 * ccomps * dcomps);

            auto g_y_y_0_xzzzzz = cbuffer.data(si_geom_11_off + 132 * ccomps * dcomps);

            auto g_y_y_0_yyyyyy = cbuffer.data(si_geom_11_off + 133 * ccomps * dcomps);

            auto g_y_y_0_yyyyyz = cbuffer.data(si_geom_11_off + 134 * ccomps * dcomps);

            auto g_y_y_0_yyyyzz = cbuffer.data(si_geom_11_off + 135 * ccomps * dcomps);

            auto g_y_y_0_yyyzzz = cbuffer.data(si_geom_11_off + 136 * ccomps * dcomps);

            auto g_y_y_0_yyzzzz = cbuffer.data(si_geom_11_off + 137 * ccomps * dcomps);

            auto g_y_y_0_yzzzzz = cbuffer.data(si_geom_11_off + 138 * ccomps * dcomps);

            auto g_y_y_0_zzzzzz = cbuffer.data(si_geom_11_off + 139 * ccomps * dcomps);

            auto g_y_z_0_xxxxxx = cbuffer.data(si_geom_11_off + 140 * ccomps * dcomps);

            auto g_y_z_0_xxxxxy = cbuffer.data(si_geom_11_off + 141 * ccomps * dcomps);

            auto g_y_z_0_xxxxxz = cbuffer.data(si_geom_11_off + 142 * ccomps * dcomps);

            auto g_y_z_0_xxxxyy = cbuffer.data(si_geom_11_off + 143 * ccomps * dcomps);

            auto g_y_z_0_xxxxyz = cbuffer.data(si_geom_11_off + 144 * ccomps * dcomps);

            auto g_y_z_0_xxxxzz = cbuffer.data(si_geom_11_off + 145 * ccomps * dcomps);

            auto g_y_z_0_xxxyyy = cbuffer.data(si_geom_11_off + 146 * ccomps * dcomps);

            auto g_y_z_0_xxxyyz = cbuffer.data(si_geom_11_off + 147 * ccomps * dcomps);

            auto g_y_z_0_xxxyzz = cbuffer.data(si_geom_11_off + 148 * ccomps * dcomps);

            auto g_y_z_0_xxxzzz = cbuffer.data(si_geom_11_off + 149 * ccomps * dcomps);

            auto g_y_z_0_xxyyyy = cbuffer.data(si_geom_11_off + 150 * ccomps * dcomps);

            auto g_y_z_0_xxyyyz = cbuffer.data(si_geom_11_off + 151 * ccomps * dcomps);

            auto g_y_z_0_xxyyzz = cbuffer.data(si_geom_11_off + 152 * ccomps * dcomps);

            auto g_y_z_0_xxyzzz = cbuffer.data(si_geom_11_off + 153 * ccomps * dcomps);

            auto g_y_z_0_xxzzzz = cbuffer.data(si_geom_11_off + 154 * ccomps * dcomps);

            auto g_y_z_0_xyyyyy = cbuffer.data(si_geom_11_off + 155 * ccomps * dcomps);

            auto g_y_z_0_xyyyyz = cbuffer.data(si_geom_11_off + 156 * ccomps * dcomps);

            auto g_y_z_0_xyyyzz = cbuffer.data(si_geom_11_off + 157 * ccomps * dcomps);

            auto g_y_z_0_xyyzzz = cbuffer.data(si_geom_11_off + 158 * ccomps * dcomps);

            auto g_y_z_0_xyzzzz = cbuffer.data(si_geom_11_off + 159 * ccomps * dcomps);

            auto g_y_z_0_xzzzzz = cbuffer.data(si_geom_11_off + 160 * ccomps * dcomps);

            auto g_y_z_0_yyyyyy = cbuffer.data(si_geom_11_off + 161 * ccomps * dcomps);

            auto g_y_z_0_yyyyyz = cbuffer.data(si_geom_11_off + 162 * ccomps * dcomps);

            auto g_y_z_0_yyyyzz = cbuffer.data(si_geom_11_off + 163 * ccomps * dcomps);

            auto g_y_z_0_yyyzzz = cbuffer.data(si_geom_11_off + 164 * ccomps * dcomps);

            auto g_y_z_0_yyzzzz = cbuffer.data(si_geom_11_off + 165 * ccomps * dcomps);

            auto g_y_z_0_yzzzzz = cbuffer.data(si_geom_11_off + 166 * ccomps * dcomps);

            auto g_y_z_0_zzzzzz = cbuffer.data(si_geom_11_off + 167 * ccomps * dcomps);

            auto g_z_x_0_xxxxxx = cbuffer.data(si_geom_11_off + 168 * ccomps * dcomps);

            auto g_z_x_0_xxxxxy = cbuffer.data(si_geom_11_off + 169 * ccomps * dcomps);

            auto g_z_x_0_xxxxxz = cbuffer.data(si_geom_11_off + 170 * ccomps * dcomps);

            auto g_z_x_0_xxxxyy = cbuffer.data(si_geom_11_off + 171 * ccomps * dcomps);

            auto g_z_x_0_xxxxyz = cbuffer.data(si_geom_11_off + 172 * ccomps * dcomps);

            auto g_z_x_0_xxxxzz = cbuffer.data(si_geom_11_off + 173 * ccomps * dcomps);

            auto g_z_x_0_xxxyyy = cbuffer.data(si_geom_11_off + 174 * ccomps * dcomps);

            auto g_z_x_0_xxxyyz = cbuffer.data(si_geom_11_off + 175 * ccomps * dcomps);

            auto g_z_x_0_xxxyzz = cbuffer.data(si_geom_11_off + 176 * ccomps * dcomps);

            auto g_z_x_0_xxxzzz = cbuffer.data(si_geom_11_off + 177 * ccomps * dcomps);

            auto g_z_x_0_xxyyyy = cbuffer.data(si_geom_11_off + 178 * ccomps * dcomps);

            auto g_z_x_0_xxyyyz = cbuffer.data(si_geom_11_off + 179 * ccomps * dcomps);

            auto g_z_x_0_xxyyzz = cbuffer.data(si_geom_11_off + 180 * ccomps * dcomps);

            auto g_z_x_0_xxyzzz = cbuffer.data(si_geom_11_off + 181 * ccomps * dcomps);

            auto g_z_x_0_xxzzzz = cbuffer.data(si_geom_11_off + 182 * ccomps * dcomps);

            auto g_z_x_0_xyyyyy = cbuffer.data(si_geom_11_off + 183 * ccomps * dcomps);

            auto g_z_x_0_xyyyyz = cbuffer.data(si_geom_11_off + 184 * ccomps * dcomps);

            auto g_z_x_0_xyyyzz = cbuffer.data(si_geom_11_off + 185 * ccomps * dcomps);

            auto g_z_x_0_xyyzzz = cbuffer.data(si_geom_11_off + 186 * ccomps * dcomps);

            auto g_z_x_0_xyzzzz = cbuffer.data(si_geom_11_off + 187 * ccomps * dcomps);

            auto g_z_x_0_xzzzzz = cbuffer.data(si_geom_11_off + 188 * ccomps * dcomps);

            auto g_z_x_0_yyyyyy = cbuffer.data(si_geom_11_off + 189 * ccomps * dcomps);

            auto g_z_x_0_yyyyyz = cbuffer.data(si_geom_11_off + 190 * ccomps * dcomps);

            auto g_z_x_0_yyyyzz = cbuffer.data(si_geom_11_off + 191 * ccomps * dcomps);

            auto g_z_x_0_yyyzzz = cbuffer.data(si_geom_11_off + 192 * ccomps * dcomps);

            auto g_z_x_0_yyzzzz = cbuffer.data(si_geom_11_off + 193 * ccomps * dcomps);

            auto g_z_x_0_yzzzzz = cbuffer.data(si_geom_11_off + 194 * ccomps * dcomps);

            auto g_z_x_0_zzzzzz = cbuffer.data(si_geom_11_off + 195 * ccomps * dcomps);

            auto g_z_y_0_xxxxxx = cbuffer.data(si_geom_11_off + 196 * ccomps * dcomps);

            auto g_z_y_0_xxxxxy = cbuffer.data(si_geom_11_off + 197 * ccomps * dcomps);

            auto g_z_y_0_xxxxxz = cbuffer.data(si_geom_11_off + 198 * ccomps * dcomps);

            auto g_z_y_0_xxxxyy = cbuffer.data(si_geom_11_off + 199 * ccomps * dcomps);

            auto g_z_y_0_xxxxyz = cbuffer.data(si_geom_11_off + 200 * ccomps * dcomps);

            auto g_z_y_0_xxxxzz = cbuffer.data(si_geom_11_off + 201 * ccomps * dcomps);

            auto g_z_y_0_xxxyyy = cbuffer.data(si_geom_11_off + 202 * ccomps * dcomps);

            auto g_z_y_0_xxxyyz = cbuffer.data(si_geom_11_off + 203 * ccomps * dcomps);

            auto g_z_y_0_xxxyzz = cbuffer.data(si_geom_11_off + 204 * ccomps * dcomps);

            auto g_z_y_0_xxxzzz = cbuffer.data(si_geom_11_off + 205 * ccomps * dcomps);

            auto g_z_y_0_xxyyyy = cbuffer.data(si_geom_11_off + 206 * ccomps * dcomps);

            auto g_z_y_0_xxyyyz = cbuffer.data(si_geom_11_off + 207 * ccomps * dcomps);

            auto g_z_y_0_xxyyzz = cbuffer.data(si_geom_11_off + 208 * ccomps * dcomps);

            auto g_z_y_0_xxyzzz = cbuffer.data(si_geom_11_off + 209 * ccomps * dcomps);

            auto g_z_y_0_xxzzzz = cbuffer.data(si_geom_11_off + 210 * ccomps * dcomps);

            auto g_z_y_0_xyyyyy = cbuffer.data(si_geom_11_off + 211 * ccomps * dcomps);

            auto g_z_y_0_xyyyyz = cbuffer.data(si_geom_11_off + 212 * ccomps * dcomps);

            auto g_z_y_0_xyyyzz = cbuffer.data(si_geom_11_off + 213 * ccomps * dcomps);

            auto g_z_y_0_xyyzzz = cbuffer.data(si_geom_11_off + 214 * ccomps * dcomps);

            auto g_z_y_0_xyzzzz = cbuffer.data(si_geom_11_off + 215 * ccomps * dcomps);

            auto g_z_y_0_xzzzzz = cbuffer.data(si_geom_11_off + 216 * ccomps * dcomps);

            auto g_z_y_0_yyyyyy = cbuffer.data(si_geom_11_off + 217 * ccomps * dcomps);

            auto g_z_y_0_yyyyyz = cbuffer.data(si_geom_11_off + 218 * ccomps * dcomps);

            auto g_z_y_0_yyyyzz = cbuffer.data(si_geom_11_off + 219 * ccomps * dcomps);

            auto g_z_y_0_yyyzzz = cbuffer.data(si_geom_11_off + 220 * ccomps * dcomps);

            auto g_z_y_0_yyzzzz = cbuffer.data(si_geom_11_off + 221 * ccomps * dcomps);

            auto g_z_y_0_yzzzzz = cbuffer.data(si_geom_11_off + 222 * ccomps * dcomps);

            auto g_z_y_0_zzzzzz = cbuffer.data(si_geom_11_off + 223 * ccomps * dcomps);

            auto g_z_z_0_xxxxxx = cbuffer.data(si_geom_11_off + 224 * ccomps * dcomps);

            auto g_z_z_0_xxxxxy = cbuffer.data(si_geom_11_off + 225 * ccomps * dcomps);

            auto g_z_z_0_xxxxxz = cbuffer.data(si_geom_11_off + 226 * ccomps * dcomps);

            auto g_z_z_0_xxxxyy = cbuffer.data(si_geom_11_off + 227 * ccomps * dcomps);

            auto g_z_z_0_xxxxyz = cbuffer.data(si_geom_11_off + 228 * ccomps * dcomps);

            auto g_z_z_0_xxxxzz = cbuffer.data(si_geom_11_off + 229 * ccomps * dcomps);

            auto g_z_z_0_xxxyyy = cbuffer.data(si_geom_11_off + 230 * ccomps * dcomps);

            auto g_z_z_0_xxxyyz = cbuffer.data(si_geom_11_off + 231 * ccomps * dcomps);

            auto g_z_z_0_xxxyzz = cbuffer.data(si_geom_11_off + 232 * ccomps * dcomps);

            auto g_z_z_0_xxxzzz = cbuffer.data(si_geom_11_off + 233 * ccomps * dcomps);

            auto g_z_z_0_xxyyyy = cbuffer.data(si_geom_11_off + 234 * ccomps * dcomps);

            auto g_z_z_0_xxyyyz = cbuffer.data(si_geom_11_off + 235 * ccomps * dcomps);

            auto g_z_z_0_xxyyzz = cbuffer.data(si_geom_11_off + 236 * ccomps * dcomps);

            auto g_z_z_0_xxyzzz = cbuffer.data(si_geom_11_off + 237 * ccomps * dcomps);

            auto g_z_z_0_xxzzzz = cbuffer.data(si_geom_11_off + 238 * ccomps * dcomps);

            auto g_z_z_0_xyyyyy = cbuffer.data(si_geom_11_off + 239 * ccomps * dcomps);

            auto g_z_z_0_xyyyyz = cbuffer.data(si_geom_11_off + 240 * ccomps * dcomps);

            auto g_z_z_0_xyyyzz = cbuffer.data(si_geom_11_off + 241 * ccomps * dcomps);

            auto g_z_z_0_xyyzzz = cbuffer.data(si_geom_11_off + 242 * ccomps * dcomps);

            auto g_z_z_0_xyzzzz = cbuffer.data(si_geom_11_off + 243 * ccomps * dcomps);

            auto g_z_z_0_xzzzzz = cbuffer.data(si_geom_11_off + 244 * ccomps * dcomps);

            auto g_z_z_0_yyyyyy = cbuffer.data(si_geom_11_off + 245 * ccomps * dcomps);

            auto g_z_z_0_yyyyyz = cbuffer.data(si_geom_11_off + 246 * ccomps * dcomps);

            auto g_z_z_0_yyyyzz = cbuffer.data(si_geom_11_off + 247 * ccomps * dcomps);

            auto g_z_z_0_yyyzzz = cbuffer.data(si_geom_11_off + 248 * ccomps * dcomps);

            auto g_z_z_0_yyzzzz = cbuffer.data(si_geom_11_off + 249 * ccomps * dcomps);

            auto g_z_z_0_yzzzzz = cbuffer.data(si_geom_11_off + 250 * ccomps * dcomps);

            auto g_z_z_0_zzzzzz = cbuffer.data(si_geom_11_off + 251 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_phxx

            const auto ph_geom_11_off = idx_geom_11_phxx + i * dcomps + j;

            /// Set up 0-21 components of targeted buffer : cbuffer.data(

            auto g_x_x_x_xxxxx = cbuffer.data(ph_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_x_xxxxy = cbuffer.data(ph_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_x_xxxxz = cbuffer.data(ph_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_x_xxxyy = cbuffer.data(ph_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_x_xxxyz = cbuffer.data(ph_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_x_xxxzz = cbuffer.data(ph_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_x_xxyyy = cbuffer.data(ph_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_x_xxyyz = cbuffer.data(ph_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_x_xxyzz = cbuffer.data(ph_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_x_xxzzz = cbuffer.data(ph_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_x_xyyyy = cbuffer.data(ph_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_x_xyyyz = cbuffer.data(ph_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_x_xyyzz = cbuffer.data(ph_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_x_xyzzz = cbuffer.data(ph_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_x_xzzzz = cbuffer.data(ph_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_x_yyyyy = cbuffer.data(ph_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_x_yyyyz = cbuffer.data(ph_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_x_yyyzz = cbuffer.data(ph_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_x_x_yyzzz = cbuffer.data(ph_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_x_yzzzz = cbuffer.data(ph_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_x_x_zzzzz = cbuffer.data(ph_geom_11_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xxxxx, g_0_x_0_xxxxy, g_0_x_0_xxxxz, g_0_x_0_xxxyy, g_0_x_0_xxxyz, g_0_x_0_xxxzz, g_0_x_0_xxyyy, g_0_x_0_xxyyz, g_0_x_0_xxyzz, g_0_x_0_xxzzz, g_0_x_0_xyyyy, g_0_x_0_xyyyz, g_0_x_0_xyyzz, g_0_x_0_xyzzz, g_0_x_0_xzzzz, g_0_x_0_yyyyy, g_0_x_0_yyyyz, g_0_x_0_yyyzz, g_0_x_0_yyzzz, g_0_x_0_yzzzz, g_0_x_0_zzzzz, g_x_0_0_xxxxx, g_x_0_0_xxxxy, g_x_0_0_xxxxz, g_x_0_0_xxxyy, g_x_0_0_xxxyz, g_x_0_0_xxxzz, g_x_0_0_xxyyy, g_x_0_0_xxyyz, g_x_0_0_xxyzz, g_x_0_0_xxzzz, g_x_0_0_xyyyy, g_x_0_0_xyyyz, g_x_0_0_xyyzz, g_x_0_0_xyzzz, g_x_0_0_xzzzz, g_x_0_0_yyyyy, g_x_0_0_yyyyz, g_x_0_0_yyyzz, g_x_0_0_yyzzz, g_x_0_0_yzzzz, g_x_0_0_zzzzz, g_x_x_0_xxxxx, g_x_x_0_xxxxxx, g_x_x_0_xxxxxy, g_x_x_0_xxxxxz, g_x_x_0_xxxxy, g_x_x_0_xxxxyy, g_x_x_0_xxxxyz, g_x_x_0_xxxxz, g_x_x_0_xxxxzz, g_x_x_0_xxxyy, g_x_x_0_xxxyyy, g_x_x_0_xxxyyz, g_x_x_0_xxxyz, g_x_x_0_xxxyzz, g_x_x_0_xxxzz, g_x_x_0_xxxzzz, g_x_x_0_xxyyy, g_x_x_0_xxyyyy, g_x_x_0_xxyyyz, g_x_x_0_xxyyz, g_x_x_0_xxyyzz, g_x_x_0_xxyzz, g_x_x_0_xxyzzz, g_x_x_0_xxzzz, g_x_x_0_xxzzzz, g_x_x_0_xyyyy, g_x_x_0_xyyyyy, g_x_x_0_xyyyyz, g_x_x_0_xyyyz, g_x_x_0_xyyyzz, g_x_x_0_xyyzz, g_x_x_0_xyyzzz, g_x_x_0_xyzzz, g_x_x_0_xyzzzz, g_x_x_0_xzzzz, g_x_x_0_xzzzzz, g_x_x_0_yyyyy, g_x_x_0_yyyyz, g_x_x_0_yyyzz, g_x_x_0_yyzzz, g_x_x_0_yzzzz, g_x_x_0_zzzzz, g_x_x_x_xxxxx, g_x_x_x_xxxxy, g_x_x_x_xxxxz, g_x_x_x_xxxyy, g_x_x_x_xxxyz, g_x_x_x_xxxzz, g_x_x_x_xxyyy, g_x_x_x_xxyyz, g_x_x_x_xxyzz, g_x_x_x_xxzzz, g_x_x_x_xyyyy, g_x_x_x_xyyyz, g_x_x_x_xyyzz, g_x_x_x_xyzzz, g_x_x_x_xzzzz, g_x_x_x_yyyyy, g_x_x_x_yyyyz, g_x_x_x_yyyzz, g_x_x_x_yyzzz, g_x_x_x_yzzzz, g_x_x_x_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_x_xxxxx[k] = -g_0_x_0_xxxxx[k] + g_x_0_0_xxxxx[k] - g_x_x_0_xxxxx[k] * ab_x + g_x_x_0_xxxxxx[k];

                g_x_x_x_xxxxy[k] = -g_0_x_0_xxxxy[k] + g_x_0_0_xxxxy[k] - g_x_x_0_xxxxy[k] * ab_x + g_x_x_0_xxxxxy[k];

                g_x_x_x_xxxxz[k] = -g_0_x_0_xxxxz[k] + g_x_0_0_xxxxz[k] - g_x_x_0_xxxxz[k] * ab_x + g_x_x_0_xxxxxz[k];

                g_x_x_x_xxxyy[k] = -g_0_x_0_xxxyy[k] + g_x_0_0_xxxyy[k] - g_x_x_0_xxxyy[k] * ab_x + g_x_x_0_xxxxyy[k];

                g_x_x_x_xxxyz[k] = -g_0_x_0_xxxyz[k] + g_x_0_0_xxxyz[k] - g_x_x_0_xxxyz[k] * ab_x + g_x_x_0_xxxxyz[k];

                g_x_x_x_xxxzz[k] = -g_0_x_0_xxxzz[k] + g_x_0_0_xxxzz[k] - g_x_x_0_xxxzz[k] * ab_x + g_x_x_0_xxxxzz[k];

                g_x_x_x_xxyyy[k] = -g_0_x_0_xxyyy[k] + g_x_0_0_xxyyy[k] - g_x_x_0_xxyyy[k] * ab_x + g_x_x_0_xxxyyy[k];

                g_x_x_x_xxyyz[k] = -g_0_x_0_xxyyz[k] + g_x_0_0_xxyyz[k] - g_x_x_0_xxyyz[k] * ab_x + g_x_x_0_xxxyyz[k];

                g_x_x_x_xxyzz[k] = -g_0_x_0_xxyzz[k] + g_x_0_0_xxyzz[k] - g_x_x_0_xxyzz[k] * ab_x + g_x_x_0_xxxyzz[k];

                g_x_x_x_xxzzz[k] = -g_0_x_0_xxzzz[k] + g_x_0_0_xxzzz[k] - g_x_x_0_xxzzz[k] * ab_x + g_x_x_0_xxxzzz[k];

                g_x_x_x_xyyyy[k] = -g_0_x_0_xyyyy[k] + g_x_0_0_xyyyy[k] - g_x_x_0_xyyyy[k] * ab_x + g_x_x_0_xxyyyy[k];

                g_x_x_x_xyyyz[k] = -g_0_x_0_xyyyz[k] + g_x_0_0_xyyyz[k] - g_x_x_0_xyyyz[k] * ab_x + g_x_x_0_xxyyyz[k];

                g_x_x_x_xyyzz[k] = -g_0_x_0_xyyzz[k] + g_x_0_0_xyyzz[k] - g_x_x_0_xyyzz[k] * ab_x + g_x_x_0_xxyyzz[k];

                g_x_x_x_xyzzz[k] = -g_0_x_0_xyzzz[k] + g_x_0_0_xyzzz[k] - g_x_x_0_xyzzz[k] * ab_x + g_x_x_0_xxyzzz[k];

                g_x_x_x_xzzzz[k] = -g_0_x_0_xzzzz[k] + g_x_0_0_xzzzz[k] - g_x_x_0_xzzzz[k] * ab_x + g_x_x_0_xxzzzz[k];

                g_x_x_x_yyyyy[k] = -g_0_x_0_yyyyy[k] + g_x_0_0_yyyyy[k] - g_x_x_0_yyyyy[k] * ab_x + g_x_x_0_xyyyyy[k];

                g_x_x_x_yyyyz[k] = -g_0_x_0_yyyyz[k] + g_x_0_0_yyyyz[k] - g_x_x_0_yyyyz[k] * ab_x + g_x_x_0_xyyyyz[k];

                g_x_x_x_yyyzz[k] = -g_0_x_0_yyyzz[k] + g_x_0_0_yyyzz[k] - g_x_x_0_yyyzz[k] * ab_x + g_x_x_0_xyyyzz[k];

                g_x_x_x_yyzzz[k] = -g_0_x_0_yyzzz[k] + g_x_0_0_yyzzz[k] - g_x_x_0_yyzzz[k] * ab_x + g_x_x_0_xyyzzz[k];

                g_x_x_x_yzzzz[k] = -g_0_x_0_yzzzz[k] + g_x_0_0_yzzzz[k] - g_x_x_0_yzzzz[k] * ab_x + g_x_x_0_xyzzzz[k];

                g_x_x_x_zzzzz[k] = -g_0_x_0_zzzzz[k] + g_x_0_0_zzzzz[k] - g_x_x_0_zzzzz[k] * ab_x + g_x_x_0_xzzzzz[k];
            }

            /// Set up 21-42 components of targeted buffer : cbuffer.data(

            auto g_x_x_y_xxxxx = cbuffer.data(ph_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_y_xxxxy = cbuffer.data(ph_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_y_xxxxz = cbuffer.data(ph_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_x_y_xxxyy = cbuffer.data(ph_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_y_xxxyz = cbuffer.data(ph_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_y_xxxzz = cbuffer.data(ph_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_x_y_xxyyy = cbuffer.data(ph_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_x_y_xxyyz = cbuffer.data(ph_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_x_y_xxyzz = cbuffer.data(ph_geom_11_off + 29 * ccomps * dcomps);

            auto g_x_x_y_xxzzz = cbuffer.data(ph_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_x_y_xyyyy = cbuffer.data(ph_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_x_y_xyyyz = cbuffer.data(ph_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_x_y_xyyzz = cbuffer.data(ph_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_x_y_xyzzz = cbuffer.data(ph_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_x_y_xzzzz = cbuffer.data(ph_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_x_y_yyyyy = cbuffer.data(ph_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_x_y_yyyyz = cbuffer.data(ph_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_x_y_yyyzz = cbuffer.data(ph_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_x_y_yyzzz = cbuffer.data(ph_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_x_y_yzzzz = cbuffer.data(ph_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_x_y_zzzzz = cbuffer.data(ph_geom_11_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_0_xxxxx, g_x_x_0_xxxxxy, g_x_x_0_xxxxy, g_x_x_0_xxxxyy, g_x_x_0_xxxxyz, g_x_x_0_xxxxz, g_x_x_0_xxxyy, g_x_x_0_xxxyyy, g_x_x_0_xxxyyz, g_x_x_0_xxxyz, g_x_x_0_xxxyzz, g_x_x_0_xxxzz, g_x_x_0_xxyyy, g_x_x_0_xxyyyy, g_x_x_0_xxyyyz, g_x_x_0_xxyyz, g_x_x_0_xxyyzz, g_x_x_0_xxyzz, g_x_x_0_xxyzzz, g_x_x_0_xxzzz, g_x_x_0_xyyyy, g_x_x_0_xyyyyy, g_x_x_0_xyyyyz, g_x_x_0_xyyyz, g_x_x_0_xyyyzz, g_x_x_0_xyyzz, g_x_x_0_xyyzzz, g_x_x_0_xyzzz, g_x_x_0_xyzzzz, g_x_x_0_xzzzz, g_x_x_0_yyyyy, g_x_x_0_yyyyyy, g_x_x_0_yyyyyz, g_x_x_0_yyyyz, g_x_x_0_yyyyzz, g_x_x_0_yyyzz, g_x_x_0_yyyzzz, g_x_x_0_yyzzz, g_x_x_0_yyzzzz, g_x_x_0_yzzzz, g_x_x_0_yzzzzz, g_x_x_0_zzzzz, g_x_x_y_xxxxx, g_x_x_y_xxxxy, g_x_x_y_xxxxz, g_x_x_y_xxxyy, g_x_x_y_xxxyz, g_x_x_y_xxxzz, g_x_x_y_xxyyy, g_x_x_y_xxyyz, g_x_x_y_xxyzz, g_x_x_y_xxzzz, g_x_x_y_xyyyy, g_x_x_y_xyyyz, g_x_x_y_xyyzz, g_x_x_y_xyzzz, g_x_x_y_xzzzz, g_x_x_y_yyyyy, g_x_x_y_yyyyz, g_x_x_y_yyyzz, g_x_x_y_yyzzz, g_x_x_y_yzzzz, g_x_x_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_y_xxxxx[k] = -g_x_x_0_xxxxx[k] * ab_y + g_x_x_0_xxxxxy[k];

                g_x_x_y_xxxxy[k] = -g_x_x_0_xxxxy[k] * ab_y + g_x_x_0_xxxxyy[k];

                g_x_x_y_xxxxz[k] = -g_x_x_0_xxxxz[k] * ab_y + g_x_x_0_xxxxyz[k];

                g_x_x_y_xxxyy[k] = -g_x_x_0_xxxyy[k] * ab_y + g_x_x_0_xxxyyy[k];

                g_x_x_y_xxxyz[k] = -g_x_x_0_xxxyz[k] * ab_y + g_x_x_0_xxxyyz[k];

                g_x_x_y_xxxzz[k] = -g_x_x_0_xxxzz[k] * ab_y + g_x_x_0_xxxyzz[k];

                g_x_x_y_xxyyy[k] = -g_x_x_0_xxyyy[k] * ab_y + g_x_x_0_xxyyyy[k];

                g_x_x_y_xxyyz[k] = -g_x_x_0_xxyyz[k] * ab_y + g_x_x_0_xxyyyz[k];

                g_x_x_y_xxyzz[k] = -g_x_x_0_xxyzz[k] * ab_y + g_x_x_0_xxyyzz[k];

                g_x_x_y_xxzzz[k] = -g_x_x_0_xxzzz[k] * ab_y + g_x_x_0_xxyzzz[k];

                g_x_x_y_xyyyy[k] = -g_x_x_0_xyyyy[k] * ab_y + g_x_x_0_xyyyyy[k];

                g_x_x_y_xyyyz[k] = -g_x_x_0_xyyyz[k] * ab_y + g_x_x_0_xyyyyz[k];

                g_x_x_y_xyyzz[k] = -g_x_x_0_xyyzz[k] * ab_y + g_x_x_0_xyyyzz[k];

                g_x_x_y_xyzzz[k] = -g_x_x_0_xyzzz[k] * ab_y + g_x_x_0_xyyzzz[k];

                g_x_x_y_xzzzz[k] = -g_x_x_0_xzzzz[k] * ab_y + g_x_x_0_xyzzzz[k];

                g_x_x_y_yyyyy[k] = -g_x_x_0_yyyyy[k] * ab_y + g_x_x_0_yyyyyy[k];

                g_x_x_y_yyyyz[k] = -g_x_x_0_yyyyz[k] * ab_y + g_x_x_0_yyyyyz[k];

                g_x_x_y_yyyzz[k] = -g_x_x_0_yyyzz[k] * ab_y + g_x_x_0_yyyyzz[k];

                g_x_x_y_yyzzz[k] = -g_x_x_0_yyzzz[k] * ab_y + g_x_x_0_yyyzzz[k];

                g_x_x_y_yzzzz[k] = -g_x_x_0_yzzzz[k] * ab_y + g_x_x_0_yyzzzz[k];

                g_x_x_y_zzzzz[k] = -g_x_x_0_zzzzz[k] * ab_y + g_x_x_0_yzzzzz[k];
            }

            /// Set up 42-63 components of targeted buffer : cbuffer.data(

            auto g_x_x_z_xxxxx = cbuffer.data(ph_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_x_z_xxxxy = cbuffer.data(ph_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_x_z_xxxxz = cbuffer.data(ph_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_x_z_xxxyy = cbuffer.data(ph_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_x_z_xxxyz = cbuffer.data(ph_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_x_z_xxxzz = cbuffer.data(ph_geom_11_off + 47 * ccomps * dcomps);

            auto g_x_x_z_xxyyy = cbuffer.data(ph_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_x_z_xxyyz = cbuffer.data(ph_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_x_z_xxyzz = cbuffer.data(ph_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_x_z_xxzzz = cbuffer.data(ph_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_x_z_xyyyy = cbuffer.data(ph_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_x_z_xyyyz = cbuffer.data(ph_geom_11_off + 53 * ccomps * dcomps);

            auto g_x_x_z_xyyzz = cbuffer.data(ph_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_x_z_xyzzz = cbuffer.data(ph_geom_11_off + 55 * ccomps * dcomps);

            auto g_x_x_z_xzzzz = cbuffer.data(ph_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_x_z_yyyyy = cbuffer.data(ph_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_x_z_yyyyz = cbuffer.data(ph_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_x_z_yyyzz = cbuffer.data(ph_geom_11_off + 59 * ccomps * dcomps);

            auto g_x_x_z_yyzzz = cbuffer.data(ph_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_x_z_yzzzz = cbuffer.data(ph_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_x_z_zzzzz = cbuffer.data(ph_geom_11_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_0_xxxxx, g_x_x_0_xxxxxz, g_x_x_0_xxxxy, g_x_x_0_xxxxyz, g_x_x_0_xxxxz, g_x_x_0_xxxxzz, g_x_x_0_xxxyy, g_x_x_0_xxxyyz, g_x_x_0_xxxyz, g_x_x_0_xxxyzz, g_x_x_0_xxxzz, g_x_x_0_xxxzzz, g_x_x_0_xxyyy, g_x_x_0_xxyyyz, g_x_x_0_xxyyz, g_x_x_0_xxyyzz, g_x_x_0_xxyzz, g_x_x_0_xxyzzz, g_x_x_0_xxzzz, g_x_x_0_xxzzzz, g_x_x_0_xyyyy, g_x_x_0_xyyyyz, g_x_x_0_xyyyz, g_x_x_0_xyyyzz, g_x_x_0_xyyzz, g_x_x_0_xyyzzz, g_x_x_0_xyzzz, g_x_x_0_xyzzzz, g_x_x_0_xzzzz, g_x_x_0_xzzzzz, g_x_x_0_yyyyy, g_x_x_0_yyyyyz, g_x_x_0_yyyyz, g_x_x_0_yyyyzz, g_x_x_0_yyyzz, g_x_x_0_yyyzzz, g_x_x_0_yyzzz, g_x_x_0_yyzzzz, g_x_x_0_yzzzz, g_x_x_0_yzzzzz, g_x_x_0_zzzzz, g_x_x_0_zzzzzz, g_x_x_z_xxxxx, g_x_x_z_xxxxy, g_x_x_z_xxxxz, g_x_x_z_xxxyy, g_x_x_z_xxxyz, g_x_x_z_xxxzz, g_x_x_z_xxyyy, g_x_x_z_xxyyz, g_x_x_z_xxyzz, g_x_x_z_xxzzz, g_x_x_z_xyyyy, g_x_x_z_xyyyz, g_x_x_z_xyyzz, g_x_x_z_xyzzz, g_x_x_z_xzzzz, g_x_x_z_yyyyy, g_x_x_z_yyyyz, g_x_x_z_yyyzz, g_x_x_z_yyzzz, g_x_x_z_yzzzz, g_x_x_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_z_xxxxx[k] = -g_x_x_0_xxxxx[k] * ab_z + g_x_x_0_xxxxxz[k];

                g_x_x_z_xxxxy[k] = -g_x_x_0_xxxxy[k] * ab_z + g_x_x_0_xxxxyz[k];

                g_x_x_z_xxxxz[k] = -g_x_x_0_xxxxz[k] * ab_z + g_x_x_0_xxxxzz[k];

                g_x_x_z_xxxyy[k] = -g_x_x_0_xxxyy[k] * ab_z + g_x_x_0_xxxyyz[k];

                g_x_x_z_xxxyz[k] = -g_x_x_0_xxxyz[k] * ab_z + g_x_x_0_xxxyzz[k];

                g_x_x_z_xxxzz[k] = -g_x_x_0_xxxzz[k] * ab_z + g_x_x_0_xxxzzz[k];

                g_x_x_z_xxyyy[k] = -g_x_x_0_xxyyy[k] * ab_z + g_x_x_0_xxyyyz[k];

                g_x_x_z_xxyyz[k] = -g_x_x_0_xxyyz[k] * ab_z + g_x_x_0_xxyyzz[k];

                g_x_x_z_xxyzz[k] = -g_x_x_0_xxyzz[k] * ab_z + g_x_x_0_xxyzzz[k];

                g_x_x_z_xxzzz[k] = -g_x_x_0_xxzzz[k] * ab_z + g_x_x_0_xxzzzz[k];

                g_x_x_z_xyyyy[k] = -g_x_x_0_xyyyy[k] * ab_z + g_x_x_0_xyyyyz[k];

                g_x_x_z_xyyyz[k] = -g_x_x_0_xyyyz[k] * ab_z + g_x_x_0_xyyyzz[k];

                g_x_x_z_xyyzz[k] = -g_x_x_0_xyyzz[k] * ab_z + g_x_x_0_xyyzzz[k];

                g_x_x_z_xyzzz[k] = -g_x_x_0_xyzzz[k] * ab_z + g_x_x_0_xyzzzz[k];

                g_x_x_z_xzzzz[k] = -g_x_x_0_xzzzz[k] * ab_z + g_x_x_0_xzzzzz[k];

                g_x_x_z_yyyyy[k] = -g_x_x_0_yyyyy[k] * ab_z + g_x_x_0_yyyyyz[k];

                g_x_x_z_yyyyz[k] = -g_x_x_0_yyyyz[k] * ab_z + g_x_x_0_yyyyzz[k];

                g_x_x_z_yyyzz[k] = -g_x_x_0_yyyzz[k] * ab_z + g_x_x_0_yyyzzz[k];

                g_x_x_z_yyzzz[k] = -g_x_x_0_yyzzz[k] * ab_z + g_x_x_0_yyzzzz[k];

                g_x_x_z_yzzzz[k] = -g_x_x_0_yzzzz[k] * ab_z + g_x_x_0_yzzzzz[k];

                g_x_x_z_zzzzz[k] = -g_x_x_0_zzzzz[k] * ab_z + g_x_x_0_zzzzzz[k];
            }

            /// Set up 63-84 components of targeted buffer : cbuffer.data(

            auto g_x_y_x_xxxxx = cbuffer.data(ph_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_y_x_xxxxy = cbuffer.data(ph_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_y_x_xxxxz = cbuffer.data(ph_geom_11_off + 65 * ccomps * dcomps);

            auto g_x_y_x_xxxyy = cbuffer.data(ph_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_y_x_xxxyz = cbuffer.data(ph_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_y_x_xxxzz = cbuffer.data(ph_geom_11_off + 68 * ccomps * dcomps);

            auto g_x_y_x_xxyyy = cbuffer.data(ph_geom_11_off + 69 * ccomps * dcomps);

            auto g_x_y_x_xxyyz = cbuffer.data(ph_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_y_x_xxyzz = cbuffer.data(ph_geom_11_off + 71 * ccomps * dcomps);

            auto g_x_y_x_xxzzz = cbuffer.data(ph_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_y_x_xyyyy = cbuffer.data(ph_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_y_x_xyyyz = cbuffer.data(ph_geom_11_off + 74 * ccomps * dcomps);

            auto g_x_y_x_xyyzz = cbuffer.data(ph_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_y_x_xyzzz = cbuffer.data(ph_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_y_x_xzzzz = cbuffer.data(ph_geom_11_off + 77 * ccomps * dcomps);

            auto g_x_y_x_yyyyy = cbuffer.data(ph_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_y_x_yyyyz = cbuffer.data(ph_geom_11_off + 79 * ccomps * dcomps);

            auto g_x_y_x_yyyzz = cbuffer.data(ph_geom_11_off + 80 * ccomps * dcomps);

            auto g_x_y_x_yyzzz = cbuffer.data(ph_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_y_x_yzzzz = cbuffer.data(ph_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_y_x_zzzzz = cbuffer.data(ph_geom_11_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_0_xxxxx, g_0_y_0_xxxxy, g_0_y_0_xxxxz, g_0_y_0_xxxyy, g_0_y_0_xxxyz, g_0_y_0_xxxzz, g_0_y_0_xxyyy, g_0_y_0_xxyyz, g_0_y_0_xxyzz, g_0_y_0_xxzzz, g_0_y_0_xyyyy, g_0_y_0_xyyyz, g_0_y_0_xyyzz, g_0_y_0_xyzzz, g_0_y_0_xzzzz, g_0_y_0_yyyyy, g_0_y_0_yyyyz, g_0_y_0_yyyzz, g_0_y_0_yyzzz, g_0_y_0_yzzzz, g_0_y_0_zzzzz, g_x_y_0_xxxxx, g_x_y_0_xxxxxx, g_x_y_0_xxxxxy, g_x_y_0_xxxxxz, g_x_y_0_xxxxy, g_x_y_0_xxxxyy, g_x_y_0_xxxxyz, g_x_y_0_xxxxz, g_x_y_0_xxxxzz, g_x_y_0_xxxyy, g_x_y_0_xxxyyy, g_x_y_0_xxxyyz, g_x_y_0_xxxyz, g_x_y_0_xxxyzz, g_x_y_0_xxxzz, g_x_y_0_xxxzzz, g_x_y_0_xxyyy, g_x_y_0_xxyyyy, g_x_y_0_xxyyyz, g_x_y_0_xxyyz, g_x_y_0_xxyyzz, g_x_y_0_xxyzz, g_x_y_0_xxyzzz, g_x_y_0_xxzzz, g_x_y_0_xxzzzz, g_x_y_0_xyyyy, g_x_y_0_xyyyyy, g_x_y_0_xyyyyz, g_x_y_0_xyyyz, g_x_y_0_xyyyzz, g_x_y_0_xyyzz, g_x_y_0_xyyzzz, g_x_y_0_xyzzz, g_x_y_0_xyzzzz, g_x_y_0_xzzzz, g_x_y_0_xzzzzz, g_x_y_0_yyyyy, g_x_y_0_yyyyz, g_x_y_0_yyyzz, g_x_y_0_yyzzz, g_x_y_0_yzzzz, g_x_y_0_zzzzz, g_x_y_x_xxxxx, g_x_y_x_xxxxy, g_x_y_x_xxxxz, g_x_y_x_xxxyy, g_x_y_x_xxxyz, g_x_y_x_xxxzz, g_x_y_x_xxyyy, g_x_y_x_xxyyz, g_x_y_x_xxyzz, g_x_y_x_xxzzz, g_x_y_x_xyyyy, g_x_y_x_xyyyz, g_x_y_x_xyyzz, g_x_y_x_xyzzz, g_x_y_x_xzzzz, g_x_y_x_yyyyy, g_x_y_x_yyyyz, g_x_y_x_yyyzz, g_x_y_x_yyzzz, g_x_y_x_yzzzz, g_x_y_x_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_x_xxxxx[k] = -g_0_y_0_xxxxx[k] - g_x_y_0_xxxxx[k] * ab_x + g_x_y_0_xxxxxx[k];

                g_x_y_x_xxxxy[k] = -g_0_y_0_xxxxy[k] - g_x_y_0_xxxxy[k] * ab_x + g_x_y_0_xxxxxy[k];

                g_x_y_x_xxxxz[k] = -g_0_y_0_xxxxz[k] - g_x_y_0_xxxxz[k] * ab_x + g_x_y_0_xxxxxz[k];

                g_x_y_x_xxxyy[k] = -g_0_y_0_xxxyy[k] - g_x_y_0_xxxyy[k] * ab_x + g_x_y_0_xxxxyy[k];

                g_x_y_x_xxxyz[k] = -g_0_y_0_xxxyz[k] - g_x_y_0_xxxyz[k] * ab_x + g_x_y_0_xxxxyz[k];

                g_x_y_x_xxxzz[k] = -g_0_y_0_xxxzz[k] - g_x_y_0_xxxzz[k] * ab_x + g_x_y_0_xxxxzz[k];

                g_x_y_x_xxyyy[k] = -g_0_y_0_xxyyy[k] - g_x_y_0_xxyyy[k] * ab_x + g_x_y_0_xxxyyy[k];

                g_x_y_x_xxyyz[k] = -g_0_y_0_xxyyz[k] - g_x_y_0_xxyyz[k] * ab_x + g_x_y_0_xxxyyz[k];

                g_x_y_x_xxyzz[k] = -g_0_y_0_xxyzz[k] - g_x_y_0_xxyzz[k] * ab_x + g_x_y_0_xxxyzz[k];

                g_x_y_x_xxzzz[k] = -g_0_y_0_xxzzz[k] - g_x_y_0_xxzzz[k] * ab_x + g_x_y_0_xxxzzz[k];

                g_x_y_x_xyyyy[k] = -g_0_y_0_xyyyy[k] - g_x_y_0_xyyyy[k] * ab_x + g_x_y_0_xxyyyy[k];

                g_x_y_x_xyyyz[k] = -g_0_y_0_xyyyz[k] - g_x_y_0_xyyyz[k] * ab_x + g_x_y_0_xxyyyz[k];

                g_x_y_x_xyyzz[k] = -g_0_y_0_xyyzz[k] - g_x_y_0_xyyzz[k] * ab_x + g_x_y_0_xxyyzz[k];

                g_x_y_x_xyzzz[k] = -g_0_y_0_xyzzz[k] - g_x_y_0_xyzzz[k] * ab_x + g_x_y_0_xxyzzz[k];

                g_x_y_x_xzzzz[k] = -g_0_y_0_xzzzz[k] - g_x_y_0_xzzzz[k] * ab_x + g_x_y_0_xxzzzz[k];

                g_x_y_x_yyyyy[k] = -g_0_y_0_yyyyy[k] - g_x_y_0_yyyyy[k] * ab_x + g_x_y_0_xyyyyy[k];

                g_x_y_x_yyyyz[k] = -g_0_y_0_yyyyz[k] - g_x_y_0_yyyyz[k] * ab_x + g_x_y_0_xyyyyz[k];

                g_x_y_x_yyyzz[k] = -g_0_y_0_yyyzz[k] - g_x_y_0_yyyzz[k] * ab_x + g_x_y_0_xyyyzz[k];

                g_x_y_x_yyzzz[k] = -g_0_y_0_yyzzz[k] - g_x_y_0_yyzzz[k] * ab_x + g_x_y_0_xyyzzz[k];

                g_x_y_x_yzzzz[k] = -g_0_y_0_yzzzz[k] - g_x_y_0_yzzzz[k] * ab_x + g_x_y_0_xyzzzz[k];

                g_x_y_x_zzzzz[k] = -g_0_y_0_zzzzz[k] - g_x_y_0_zzzzz[k] * ab_x + g_x_y_0_xzzzzz[k];
            }

            /// Set up 84-105 components of targeted buffer : cbuffer.data(

            auto g_x_y_y_xxxxx = cbuffer.data(ph_geom_11_off + 84 * ccomps * dcomps);

            auto g_x_y_y_xxxxy = cbuffer.data(ph_geom_11_off + 85 * ccomps * dcomps);

            auto g_x_y_y_xxxxz = cbuffer.data(ph_geom_11_off + 86 * ccomps * dcomps);

            auto g_x_y_y_xxxyy = cbuffer.data(ph_geom_11_off + 87 * ccomps * dcomps);

            auto g_x_y_y_xxxyz = cbuffer.data(ph_geom_11_off + 88 * ccomps * dcomps);

            auto g_x_y_y_xxxzz = cbuffer.data(ph_geom_11_off + 89 * ccomps * dcomps);

            auto g_x_y_y_xxyyy = cbuffer.data(ph_geom_11_off + 90 * ccomps * dcomps);

            auto g_x_y_y_xxyyz = cbuffer.data(ph_geom_11_off + 91 * ccomps * dcomps);

            auto g_x_y_y_xxyzz = cbuffer.data(ph_geom_11_off + 92 * ccomps * dcomps);

            auto g_x_y_y_xxzzz = cbuffer.data(ph_geom_11_off + 93 * ccomps * dcomps);

            auto g_x_y_y_xyyyy = cbuffer.data(ph_geom_11_off + 94 * ccomps * dcomps);

            auto g_x_y_y_xyyyz = cbuffer.data(ph_geom_11_off + 95 * ccomps * dcomps);

            auto g_x_y_y_xyyzz = cbuffer.data(ph_geom_11_off + 96 * ccomps * dcomps);

            auto g_x_y_y_xyzzz = cbuffer.data(ph_geom_11_off + 97 * ccomps * dcomps);

            auto g_x_y_y_xzzzz = cbuffer.data(ph_geom_11_off + 98 * ccomps * dcomps);

            auto g_x_y_y_yyyyy = cbuffer.data(ph_geom_11_off + 99 * ccomps * dcomps);

            auto g_x_y_y_yyyyz = cbuffer.data(ph_geom_11_off + 100 * ccomps * dcomps);

            auto g_x_y_y_yyyzz = cbuffer.data(ph_geom_11_off + 101 * ccomps * dcomps);

            auto g_x_y_y_yyzzz = cbuffer.data(ph_geom_11_off + 102 * ccomps * dcomps);

            auto g_x_y_y_yzzzz = cbuffer.data(ph_geom_11_off + 103 * ccomps * dcomps);

            auto g_x_y_y_zzzzz = cbuffer.data(ph_geom_11_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_0_xxxxx, g_x_0_0_xxxxy, g_x_0_0_xxxxz, g_x_0_0_xxxyy, g_x_0_0_xxxyz, g_x_0_0_xxxzz, g_x_0_0_xxyyy, g_x_0_0_xxyyz, g_x_0_0_xxyzz, g_x_0_0_xxzzz, g_x_0_0_xyyyy, g_x_0_0_xyyyz, g_x_0_0_xyyzz, g_x_0_0_xyzzz, g_x_0_0_xzzzz, g_x_0_0_yyyyy, g_x_0_0_yyyyz, g_x_0_0_yyyzz, g_x_0_0_yyzzz, g_x_0_0_yzzzz, g_x_0_0_zzzzz, g_x_y_0_xxxxx, g_x_y_0_xxxxxy, g_x_y_0_xxxxy, g_x_y_0_xxxxyy, g_x_y_0_xxxxyz, g_x_y_0_xxxxz, g_x_y_0_xxxyy, g_x_y_0_xxxyyy, g_x_y_0_xxxyyz, g_x_y_0_xxxyz, g_x_y_0_xxxyzz, g_x_y_0_xxxzz, g_x_y_0_xxyyy, g_x_y_0_xxyyyy, g_x_y_0_xxyyyz, g_x_y_0_xxyyz, g_x_y_0_xxyyzz, g_x_y_0_xxyzz, g_x_y_0_xxyzzz, g_x_y_0_xxzzz, g_x_y_0_xyyyy, g_x_y_0_xyyyyy, g_x_y_0_xyyyyz, g_x_y_0_xyyyz, g_x_y_0_xyyyzz, g_x_y_0_xyyzz, g_x_y_0_xyyzzz, g_x_y_0_xyzzz, g_x_y_0_xyzzzz, g_x_y_0_xzzzz, g_x_y_0_yyyyy, g_x_y_0_yyyyyy, g_x_y_0_yyyyyz, g_x_y_0_yyyyz, g_x_y_0_yyyyzz, g_x_y_0_yyyzz, g_x_y_0_yyyzzz, g_x_y_0_yyzzz, g_x_y_0_yyzzzz, g_x_y_0_yzzzz, g_x_y_0_yzzzzz, g_x_y_0_zzzzz, g_x_y_y_xxxxx, g_x_y_y_xxxxy, g_x_y_y_xxxxz, g_x_y_y_xxxyy, g_x_y_y_xxxyz, g_x_y_y_xxxzz, g_x_y_y_xxyyy, g_x_y_y_xxyyz, g_x_y_y_xxyzz, g_x_y_y_xxzzz, g_x_y_y_xyyyy, g_x_y_y_xyyyz, g_x_y_y_xyyzz, g_x_y_y_xyzzz, g_x_y_y_xzzzz, g_x_y_y_yyyyy, g_x_y_y_yyyyz, g_x_y_y_yyyzz, g_x_y_y_yyzzz, g_x_y_y_yzzzz, g_x_y_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_y_xxxxx[k] = g_x_0_0_xxxxx[k] - g_x_y_0_xxxxx[k] * ab_y + g_x_y_0_xxxxxy[k];

                g_x_y_y_xxxxy[k] = g_x_0_0_xxxxy[k] - g_x_y_0_xxxxy[k] * ab_y + g_x_y_0_xxxxyy[k];

                g_x_y_y_xxxxz[k] = g_x_0_0_xxxxz[k] - g_x_y_0_xxxxz[k] * ab_y + g_x_y_0_xxxxyz[k];

                g_x_y_y_xxxyy[k] = g_x_0_0_xxxyy[k] - g_x_y_0_xxxyy[k] * ab_y + g_x_y_0_xxxyyy[k];

                g_x_y_y_xxxyz[k] = g_x_0_0_xxxyz[k] - g_x_y_0_xxxyz[k] * ab_y + g_x_y_0_xxxyyz[k];

                g_x_y_y_xxxzz[k] = g_x_0_0_xxxzz[k] - g_x_y_0_xxxzz[k] * ab_y + g_x_y_0_xxxyzz[k];

                g_x_y_y_xxyyy[k] = g_x_0_0_xxyyy[k] - g_x_y_0_xxyyy[k] * ab_y + g_x_y_0_xxyyyy[k];

                g_x_y_y_xxyyz[k] = g_x_0_0_xxyyz[k] - g_x_y_0_xxyyz[k] * ab_y + g_x_y_0_xxyyyz[k];

                g_x_y_y_xxyzz[k] = g_x_0_0_xxyzz[k] - g_x_y_0_xxyzz[k] * ab_y + g_x_y_0_xxyyzz[k];

                g_x_y_y_xxzzz[k] = g_x_0_0_xxzzz[k] - g_x_y_0_xxzzz[k] * ab_y + g_x_y_0_xxyzzz[k];

                g_x_y_y_xyyyy[k] = g_x_0_0_xyyyy[k] - g_x_y_0_xyyyy[k] * ab_y + g_x_y_0_xyyyyy[k];

                g_x_y_y_xyyyz[k] = g_x_0_0_xyyyz[k] - g_x_y_0_xyyyz[k] * ab_y + g_x_y_0_xyyyyz[k];

                g_x_y_y_xyyzz[k] = g_x_0_0_xyyzz[k] - g_x_y_0_xyyzz[k] * ab_y + g_x_y_0_xyyyzz[k];

                g_x_y_y_xyzzz[k] = g_x_0_0_xyzzz[k] - g_x_y_0_xyzzz[k] * ab_y + g_x_y_0_xyyzzz[k];

                g_x_y_y_xzzzz[k] = g_x_0_0_xzzzz[k] - g_x_y_0_xzzzz[k] * ab_y + g_x_y_0_xyzzzz[k];

                g_x_y_y_yyyyy[k] = g_x_0_0_yyyyy[k] - g_x_y_0_yyyyy[k] * ab_y + g_x_y_0_yyyyyy[k];

                g_x_y_y_yyyyz[k] = g_x_0_0_yyyyz[k] - g_x_y_0_yyyyz[k] * ab_y + g_x_y_0_yyyyyz[k];

                g_x_y_y_yyyzz[k] = g_x_0_0_yyyzz[k] - g_x_y_0_yyyzz[k] * ab_y + g_x_y_0_yyyyzz[k];

                g_x_y_y_yyzzz[k] = g_x_0_0_yyzzz[k] - g_x_y_0_yyzzz[k] * ab_y + g_x_y_0_yyyzzz[k];

                g_x_y_y_yzzzz[k] = g_x_0_0_yzzzz[k] - g_x_y_0_yzzzz[k] * ab_y + g_x_y_0_yyzzzz[k];

                g_x_y_y_zzzzz[k] = g_x_0_0_zzzzz[k] - g_x_y_0_zzzzz[k] * ab_y + g_x_y_0_yzzzzz[k];
            }

            /// Set up 105-126 components of targeted buffer : cbuffer.data(

            auto g_x_y_z_xxxxx = cbuffer.data(ph_geom_11_off + 105 * ccomps * dcomps);

            auto g_x_y_z_xxxxy = cbuffer.data(ph_geom_11_off + 106 * ccomps * dcomps);

            auto g_x_y_z_xxxxz = cbuffer.data(ph_geom_11_off + 107 * ccomps * dcomps);

            auto g_x_y_z_xxxyy = cbuffer.data(ph_geom_11_off + 108 * ccomps * dcomps);

            auto g_x_y_z_xxxyz = cbuffer.data(ph_geom_11_off + 109 * ccomps * dcomps);

            auto g_x_y_z_xxxzz = cbuffer.data(ph_geom_11_off + 110 * ccomps * dcomps);

            auto g_x_y_z_xxyyy = cbuffer.data(ph_geom_11_off + 111 * ccomps * dcomps);

            auto g_x_y_z_xxyyz = cbuffer.data(ph_geom_11_off + 112 * ccomps * dcomps);

            auto g_x_y_z_xxyzz = cbuffer.data(ph_geom_11_off + 113 * ccomps * dcomps);

            auto g_x_y_z_xxzzz = cbuffer.data(ph_geom_11_off + 114 * ccomps * dcomps);

            auto g_x_y_z_xyyyy = cbuffer.data(ph_geom_11_off + 115 * ccomps * dcomps);

            auto g_x_y_z_xyyyz = cbuffer.data(ph_geom_11_off + 116 * ccomps * dcomps);

            auto g_x_y_z_xyyzz = cbuffer.data(ph_geom_11_off + 117 * ccomps * dcomps);

            auto g_x_y_z_xyzzz = cbuffer.data(ph_geom_11_off + 118 * ccomps * dcomps);

            auto g_x_y_z_xzzzz = cbuffer.data(ph_geom_11_off + 119 * ccomps * dcomps);

            auto g_x_y_z_yyyyy = cbuffer.data(ph_geom_11_off + 120 * ccomps * dcomps);

            auto g_x_y_z_yyyyz = cbuffer.data(ph_geom_11_off + 121 * ccomps * dcomps);

            auto g_x_y_z_yyyzz = cbuffer.data(ph_geom_11_off + 122 * ccomps * dcomps);

            auto g_x_y_z_yyzzz = cbuffer.data(ph_geom_11_off + 123 * ccomps * dcomps);

            auto g_x_y_z_yzzzz = cbuffer.data(ph_geom_11_off + 124 * ccomps * dcomps);

            auto g_x_y_z_zzzzz = cbuffer.data(ph_geom_11_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_0_xxxxx, g_x_y_0_xxxxxz, g_x_y_0_xxxxy, g_x_y_0_xxxxyz, g_x_y_0_xxxxz, g_x_y_0_xxxxzz, g_x_y_0_xxxyy, g_x_y_0_xxxyyz, g_x_y_0_xxxyz, g_x_y_0_xxxyzz, g_x_y_0_xxxzz, g_x_y_0_xxxzzz, g_x_y_0_xxyyy, g_x_y_0_xxyyyz, g_x_y_0_xxyyz, g_x_y_0_xxyyzz, g_x_y_0_xxyzz, g_x_y_0_xxyzzz, g_x_y_0_xxzzz, g_x_y_0_xxzzzz, g_x_y_0_xyyyy, g_x_y_0_xyyyyz, g_x_y_0_xyyyz, g_x_y_0_xyyyzz, g_x_y_0_xyyzz, g_x_y_0_xyyzzz, g_x_y_0_xyzzz, g_x_y_0_xyzzzz, g_x_y_0_xzzzz, g_x_y_0_xzzzzz, g_x_y_0_yyyyy, g_x_y_0_yyyyyz, g_x_y_0_yyyyz, g_x_y_0_yyyyzz, g_x_y_0_yyyzz, g_x_y_0_yyyzzz, g_x_y_0_yyzzz, g_x_y_0_yyzzzz, g_x_y_0_yzzzz, g_x_y_0_yzzzzz, g_x_y_0_zzzzz, g_x_y_0_zzzzzz, g_x_y_z_xxxxx, g_x_y_z_xxxxy, g_x_y_z_xxxxz, g_x_y_z_xxxyy, g_x_y_z_xxxyz, g_x_y_z_xxxzz, g_x_y_z_xxyyy, g_x_y_z_xxyyz, g_x_y_z_xxyzz, g_x_y_z_xxzzz, g_x_y_z_xyyyy, g_x_y_z_xyyyz, g_x_y_z_xyyzz, g_x_y_z_xyzzz, g_x_y_z_xzzzz, g_x_y_z_yyyyy, g_x_y_z_yyyyz, g_x_y_z_yyyzz, g_x_y_z_yyzzz, g_x_y_z_yzzzz, g_x_y_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_z_xxxxx[k] = -g_x_y_0_xxxxx[k] * ab_z + g_x_y_0_xxxxxz[k];

                g_x_y_z_xxxxy[k] = -g_x_y_0_xxxxy[k] * ab_z + g_x_y_0_xxxxyz[k];

                g_x_y_z_xxxxz[k] = -g_x_y_0_xxxxz[k] * ab_z + g_x_y_0_xxxxzz[k];

                g_x_y_z_xxxyy[k] = -g_x_y_0_xxxyy[k] * ab_z + g_x_y_0_xxxyyz[k];

                g_x_y_z_xxxyz[k] = -g_x_y_0_xxxyz[k] * ab_z + g_x_y_0_xxxyzz[k];

                g_x_y_z_xxxzz[k] = -g_x_y_0_xxxzz[k] * ab_z + g_x_y_0_xxxzzz[k];

                g_x_y_z_xxyyy[k] = -g_x_y_0_xxyyy[k] * ab_z + g_x_y_0_xxyyyz[k];

                g_x_y_z_xxyyz[k] = -g_x_y_0_xxyyz[k] * ab_z + g_x_y_0_xxyyzz[k];

                g_x_y_z_xxyzz[k] = -g_x_y_0_xxyzz[k] * ab_z + g_x_y_0_xxyzzz[k];

                g_x_y_z_xxzzz[k] = -g_x_y_0_xxzzz[k] * ab_z + g_x_y_0_xxzzzz[k];

                g_x_y_z_xyyyy[k] = -g_x_y_0_xyyyy[k] * ab_z + g_x_y_0_xyyyyz[k];

                g_x_y_z_xyyyz[k] = -g_x_y_0_xyyyz[k] * ab_z + g_x_y_0_xyyyzz[k];

                g_x_y_z_xyyzz[k] = -g_x_y_0_xyyzz[k] * ab_z + g_x_y_0_xyyzzz[k];

                g_x_y_z_xyzzz[k] = -g_x_y_0_xyzzz[k] * ab_z + g_x_y_0_xyzzzz[k];

                g_x_y_z_xzzzz[k] = -g_x_y_0_xzzzz[k] * ab_z + g_x_y_0_xzzzzz[k];

                g_x_y_z_yyyyy[k] = -g_x_y_0_yyyyy[k] * ab_z + g_x_y_0_yyyyyz[k];

                g_x_y_z_yyyyz[k] = -g_x_y_0_yyyyz[k] * ab_z + g_x_y_0_yyyyzz[k];

                g_x_y_z_yyyzz[k] = -g_x_y_0_yyyzz[k] * ab_z + g_x_y_0_yyyzzz[k];

                g_x_y_z_yyzzz[k] = -g_x_y_0_yyzzz[k] * ab_z + g_x_y_0_yyzzzz[k];

                g_x_y_z_yzzzz[k] = -g_x_y_0_yzzzz[k] * ab_z + g_x_y_0_yzzzzz[k];

                g_x_y_z_zzzzz[k] = -g_x_y_0_zzzzz[k] * ab_z + g_x_y_0_zzzzzz[k];
            }

            /// Set up 126-147 components of targeted buffer : cbuffer.data(

            auto g_x_z_x_xxxxx = cbuffer.data(ph_geom_11_off + 126 * ccomps * dcomps);

            auto g_x_z_x_xxxxy = cbuffer.data(ph_geom_11_off + 127 * ccomps * dcomps);

            auto g_x_z_x_xxxxz = cbuffer.data(ph_geom_11_off + 128 * ccomps * dcomps);

            auto g_x_z_x_xxxyy = cbuffer.data(ph_geom_11_off + 129 * ccomps * dcomps);

            auto g_x_z_x_xxxyz = cbuffer.data(ph_geom_11_off + 130 * ccomps * dcomps);

            auto g_x_z_x_xxxzz = cbuffer.data(ph_geom_11_off + 131 * ccomps * dcomps);

            auto g_x_z_x_xxyyy = cbuffer.data(ph_geom_11_off + 132 * ccomps * dcomps);

            auto g_x_z_x_xxyyz = cbuffer.data(ph_geom_11_off + 133 * ccomps * dcomps);

            auto g_x_z_x_xxyzz = cbuffer.data(ph_geom_11_off + 134 * ccomps * dcomps);

            auto g_x_z_x_xxzzz = cbuffer.data(ph_geom_11_off + 135 * ccomps * dcomps);

            auto g_x_z_x_xyyyy = cbuffer.data(ph_geom_11_off + 136 * ccomps * dcomps);

            auto g_x_z_x_xyyyz = cbuffer.data(ph_geom_11_off + 137 * ccomps * dcomps);

            auto g_x_z_x_xyyzz = cbuffer.data(ph_geom_11_off + 138 * ccomps * dcomps);

            auto g_x_z_x_xyzzz = cbuffer.data(ph_geom_11_off + 139 * ccomps * dcomps);

            auto g_x_z_x_xzzzz = cbuffer.data(ph_geom_11_off + 140 * ccomps * dcomps);

            auto g_x_z_x_yyyyy = cbuffer.data(ph_geom_11_off + 141 * ccomps * dcomps);

            auto g_x_z_x_yyyyz = cbuffer.data(ph_geom_11_off + 142 * ccomps * dcomps);

            auto g_x_z_x_yyyzz = cbuffer.data(ph_geom_11_off + 143 * ccomps * dcomps);

            auto g_x_z_x_yyzzz = cbuffer.data(ph_geom_11_off + 144 * ccomps * dcomps);

            auto g_x_z_x_yzzzz = cbuffer.data(ph_geom_11_off + 145 * ccomps * dcomps);

            auto g_x_z_x_zzzzz = cbuffer.data(ph_geom_11_off + 146 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_0_xxxxx, g_0_z_0_xxxxy, g_0_z_0_xxxxz, g_0_z_0_xxxyy, g_0_z_0_xxxyz, g_0_z_0_xxxzz, g_0_z_0_xxyyy, g_0_z_0_xxyyz, g_0_z_0_xxyzz, g_0_z_0_xxzzz, g_0_z_0_xyyyy, g_0_z_0_xyyyz, g_0_z_0_xyyzz, g_0_z_0_xyzzz, g_0_z_0_xzzzz, g_0_z_0_yyyyy, g_0_z_0_yyyyz, g_0_z_0_yyyzz, g_0_z_0_yyzzz, g_0_z_0_yzzzz, g_0_z_0_zzzzz, g_x_z_0_xxxxx, g_x_z_0_xxxxxx, g_x_z_0_xxxxxy, g_x_z_0_xxxxxz, g_x_z_0_xxxxy, g_x_z_0_xxxxyy, g_x_z_0_xxxxyz, g_x_z_0_xxxxz, g_x_z_0_xxxxzz, g_x_z_0_xxxyy, g_x_z_0_xxxyyy, g_x_z_0_xxxyyz, g_x_z_0_xxxyz, g_x_z_0_xxxyzz, g_x_z_0_xxxzz, g_x_z_0_xxxzzz, g_x_z_0_xxyyy, g_x_z_0_xxyyyy, g_x_z_0_xxyyyz, g_x_z_0_xxyyz, g_x_z_0_xxyyzz, g_x_z_0_xxyzz, g_x_z_0_xxyzzz, g_x_z_0_xxzzz, g_x_z_0_xxzzzz, g_x_z_0_xyyyy, g_x_z_0_xyyyyy, g_x_z_0_xyyyyz, g_x_z_0_xyyyz, g_x_z_0_xyyyzz, g_x_z_0_xyyzz, g_x_z_0_xyyzzz, g_x_z_0_xyzzz, g_x_z_0_xyzzzz, g_x_z_0_xzzzz, g_x_z_0_xzzzzz, g_x_z_0_yyyyy, g_x_z_0_yyyyz, g_x_z_0_yyyzz, g_x_z_0_yyzzz, g_x_z_0_yzzzz, g_x_z_0_zzzzz, g_x_z_x_xxxxx, g_x_z_x_xxxxy, g_x_z_x_xxxxz, g_x_z_x_xxxyy, g_x_z_x_xxxyz, g_x_z_x_xxxzz, g_x_z_x_xxyyy, g_x_z_x_xxyyz, g_x_z_x_xxyzz, g_x_z_x_xxzzz, g_x_z_x_xyyyy, g_x_z_x_xyyyz, g_x_z_x_xyyzz, g_x_z_x_xyzzz, g_x_z_x_xzzzz, g_x_z_x_yyyyy, g_x_z_x_yyyyz, g_x_z_x_yyyzz, g_x_z_x_yyzzz, g_x_z_x_yzzzz, g_x_z_x_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_x_xxxxx[k] = -g_0_z_0_xxxxx[k] - g_x_z_0_xxxxx[k] * ab_x + g_x_z_0_xxxxxx[k];

                g_x_z_x_xxxxy[k] = -g_0_z_0_xxxxy[k] - g_x_z_0_xxxxy[k] * ab_x + g_x_z_0_xxxxxy[k];

                g_x_z_x_xxxxz[k] = -g_0_z_0_xxxxz[k] - g_x_z_0_xxxxz[k] * ab_x + g_x_z_0_xxxxxz[k];

                g_x_z_x_xxxyy[k] = -g_0_z_0_xxxyy[k] - g_x_z_0_xxxyy[k] * ab_x + g_x_z_0_xxxxyy[k];

                g_x_z_x_xxxyz[k] = -g_0_z_0_xxxyz[k] - g_x_z_0_xxxyz[k] * ab_x + g_x_z_0_xxxxyz[k];

                g_x_z_x_xxxzz[k] = -g_0_z_0_xxxzz[k] - g_x_z_0_xxxzz[k] * ab_x + g_x_z_0_xxxxzz[k];

                g_x_z_x_xxyyy[k] = -g_0_z_0_xxyyy[k] - g_x_z_0_xxyyy[k] * ab_x + g_x_z_0_xxxyyy[k];

                g_x_z_x_xxyyz[k] = -g_0_z_0_xxyyz[k] - g_x_z_0_xxyyz[k] * ab_x + g_x_z_0_xxxyyz[k];

                g_x_z_x_xxyzz[k] = -g_0_z_0_xxyzz[k] - g_x_z_0_xxyzz[k] * ab_x + g_x_z_0_xxxyzz[k];

                g_x_z_x_xxzzz[k] = -g_0_z_0_xxzzz[k] - g_x_z_0_xxzzz[k] * ab_x + g_x_z_0_xxxzzz[k];

                g_x_z_x_xyyyy[k] = -g_0_z_0_xyyyy[k] - g_x_z_0_xyyyy[k] * ab_x + g_x_z_0_xxyyyy[k];

                g_x_z_x_xyyyz[k] = -g_0_z_0_xyyyz[k] - g_x_z_0_xyyyz[k] * ab_x + g_x_z_0_xxyyyz[k];

                g_x_z_x_xyyzz[k] = -g_0_z_0_xyyzz[k] - g_x_z_0_xyyzz[k] * ab_x + g_x_z_0_xxyyzz[k];

                g_x_z_x_xyzzz[k] = -g_0_z_0_xyzzz[k] - g_x_z_0_xyzzz[k] * ab_x + g_x_z_0_xxyzzz[k];

                g_x_z_x_xzzzz[k] = -g_0_z_0_xzzzz[k] - g_x_z_0_xzzzz[k] * ab_x + g_x_z_0_xxzzzz[k];

                g_x_z_x_yyyyy[k] = -g_0_z_0_yyyyy[k] - g_x_z_0_yyyyy[k] * ab_x + g_x_z_0_xyyyyy[k];

                g_x_z_x_yyyyz[k] = -g_0_z_0_yyyyz[k] - g_x_z_0_yyyyz[k] * ab_x + g_x_z_0_xyyyyz[k];

                g_x_z_x_yyyzz[k] = -g_0_z_0_yyyzz[k] - g_x_z_0_yyyzz[k] * ab_x + g_x_z_0_xyyyzz[k];

                g_x_z_x_yyzzz[k] = -g_0_z_0_yyzzz[k] - g_x_z_0_yyzzz[k] * ab_x + g_x_z_0_xyyzzz[k];

                g_x_z_x_yzzzz[k] = -g_0_z_0_yzzzz[k] - g_x_z_0_yzzzz[k] * ab_x + g_x_z_0_xyzzzz[k];

                g_x_z_x_zzzzz[k] = -g_0_z_0_zzzzz[k] - g_x_z_0_zzzzz[k] * ab_x + g_x_z_0_xzzzzz[k];
            }

            /// Set up 147-168 components of targeted buffer : cbuffer.data(

            auto g_x_z_y_xxxxx = cbuffer.data(ph_geom_11_off + 147 * ccomps * dcomps);

            auto g_x_z_y_xxxxy = cbuffer.data(ph_geom_11_off + 148 * ccomps * dcomps);

            auto g_x_z_y_xxxxz = cbuffer.data(ph_geom_11_off + 149 * ccomps * dcomps);

            auto g_x_z_y_xxxyy = cbuffer.data(ph_geom_11_off + 150 * ccomps * dcomps);

            auto g_x_z_y_xxxyz = cbuffer.data(ph_geom_11_off + 151 * ccomps * dcomps);

            auto g_x_z_y_xxxzz = cbuffer.data(ph_geom_11_off + 152 * ccomps * dcomps);

            auto g_x_z_y_xxyyy = cbuffer.data(ph_geom_11_off + 153 * ccomps * dcomps);

            auto g_x_z_y_xxyyz = cbuffer.data(ph_geom_11_off + 154 * ccomps * dcomps);

            auto g_x_z_y_xxyzz = cbuffer.data(ph_geom_11_off + 155 * ccomps * dcomps);

            auto g_x_z_y_xxzzz = cbuffer.data(ph_geom_11_off + 156 * ccomps * dcomps);

            auto g_x_z_y_xyyyy = cbuffer.data(ph_geom_11_off + 157 * ccomps * dcomps);

            auto g_x_z_y_xyyyz = cbuffer.data(ph_geom_11_off + 158 * ccomps * dcomps);

            auto g_x_z_y_xyyzz = cbuffer.data(ph_geom_11_off + 159 * ccomps * dcomps);

            auto g_x_z_y_xyzzz = cbuffer.data(ph_geom_11_off + 160 * ccomps * dcomps);

            auto g_x_z_y_xzzzz = cbuffer.data(ph_geom_11_off + 161 * ccomps * dcomps);

            auto g_x_z_y_yyyyy = cbuffer.data(ph_geom_11_off + 162 * ccomps * dcomps);

            auto g_x_z_y_yyyyz = cbuffer.data(ph_geom_11_off + 163 * ccomps * dcomps);

            auto g_x_z_y_yyyzz = cbuffer.data(ph_geom_11_off + 164 * ccomps * dcomps);

            auto g_x_z_y_yyzzz = cbuffer.data(ph_geom_11_off + 165 * ccomps * dcomps);

            auto g_x_z_y_yzzzz = cbuffer.data(ph_geom_11_off + 166 * ccomps * dcomps);

            auto g_x_z_y_zzzzz = cbuffer.data(ph_geom_11_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_0_xxxxx, g_x_z_0_xxxxxy, g_x_z_0_xxxxy, g_x_z_0_xxxxyy, g_x_z_0_xxxxyz, g_x_z_0_xxxxz, g_x_z_0_xxxyy, g_x_z_0_xxxyyy, g_x_z_0_xxxyyz, g_x_z_0_xxxyz, g_x_z_0_xxxyzz, g_x_z_0_xxxzz, g_x_z_0_xxyyy, g_x_z_0_xxyyyy, g_x_z_0_xxyyyz, g_x_z_0_xxyyz, g_x_z_0_xxyyzz, g_x_z_0_xxyzz, g_x_z_0_xxyzzz, g_x_z_0_xxzzz, g_x_z_0_xyyyy, g_x_z_0_xyyyyy, g_x_z_0_xyyyyz, g_x_z_0_xyyyz, g_x_z_0_xyyyzz, g_x_z_0_xyyzz, g_x_z_0_xyyzzz, g_x_z_0_xyzzz, g_x_z_0_xyzzzz, g_x_z_0_xzzzz, g_x_z_0_yyyyy, g_x_z_0_yyyyyy, g_x_z_0_yyyyyz, g_x_z_0_yyyyz, g_x_z_0_yyyyzz, g_x_z_0_yyyzz, g_x_z_0_yyyzzz, g_x_z_0_yyzzz, g_x_z_0_yyzzzz, g_x_z_0_yzzzz, g_x_z_0_yzzzzz, g_x_z_0_zzzzz, g_x_z_y_xxxxx, g_x_z_y_xxxxy, g_x_z_y_xxxxz, g_x_z_y_xxxyy, g_x_z_y_xxxyz, g_x_z_y_xxxzz, g_x_z_y_xxyyy, g_x_z_y_xxyyz, g_x_z_y_xxyzz, g_x_z_y_xxzzz, g_x_z_y_xyyyy, g_x_z_y_xyyyz, g_x_z_y_xyyzz, g_x_z_y_xyzzz, g_x_z_y_xzzzz, g_x_z_y_yyyyy, g_x_z_y_yyyyz, g_x_z_y_yyyzz, g_x_z_y_yyzzz, g_x_z_y_yzzzz, g_x_z_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_y_xxxxx[k] = -g_x_z_0_xxxxx[k] * ab_y + g_x_z_0_xxxxxy[k];

                g_x_z_y_xxxxy[k] = -g_x_z_0_xxxxy[k] * ab_y + g_x_z_0_xxxxyy[k];

                g_x_z_y_xxxxz[k] = -g_x_z_0_xxxxz[k] * ab_y + g_x_z_0_xxxxyz[k];

                g_x_z_y_xxxyy[k] = -g_x_z_0_xxxyy[k] * ab_y + g_x_z_0_xxxyyy[k];

                g_x_z_y_xxxyz[k] = -g_x_z_0_xxxyz[k] * ab_y + g_x_z_0_xxxyyz[k];

                g_x_z_y_xxxzz[k] = -g_x_z_0_xxxzz[k] * ab_y + g_x_z_0_xxxyzz[k];

                g_x_z_y_xxyyy[k] = -g_x_z_0_xxyyy[k] * ab_y + g_x_z_0_xxyyyy[k];

                g_x_z_y_xxyyz[k] = -g_x_z_0_xxyyz[k] * ab_y + g_x_z_0_xxyyyz[k];

                g_x_z_y_xxyzz[k] = -g_x_z_0_xxyzz[k] * ab_y + g_x_z_0_xxyyzz[k];

                g_x_z_y_xxzzz[k] = -g_x_z_0_xxzzz[k] * ab_y + g_x_z_0_xxyzzz[k];

                g_x_z_y_xyyyy[k] = -g_x_z_0_xyyyy[k] * ab_y + g_x_z_0_xyyyyy[k];

                g_x_z_y_xyyyz[k] = -g_x_z_0_xyyyz[k] * ab_y + g_x_z_0_xyyyyz[k];

                g_x_z_y_xyyzz[k] = -g_x_z_0_xyyzz[k] * ab_y + g_x_z_0_xyyyzz[k];

                g_x_z_y_xyzzz[k] = -g_x_z_0_xyzzz[k] * ab_y + g_x_z_0_xyyzzz[k];

                g_x_z_y_xzzzz[k] = -g_x_z_0_xzzzz[k] * ab_y + g_x_z_0_xyzzzz[k];

                g_x_z_y_yyyyy[k] = -g_x_z_0_yyyyy[k] * ab_y + g_x_z_0_yyyyyy[k];

                g_x_z_y_yyyyz[k] = -g_x_z_0_yyyyz[k] * ab_y + g_x_z_0_yyyyyz[k];

                g_x_z_y_yyyzz[k] = -g_x_z_0_yyyzz[k] * ab_y + g_x_z_0_yyyyzz[k];

                g_x_z_y_yyzzz[k] = -g_x_z_0_yyzzz[k] * ab_y + g_x_z_0_yyyzzz[k];

                g_x_z_y_yzzzz[k] = -g_x_z_0_yzzzz[k] * ab_y + g_x_z_0_yyzzzz[k];

                g_x_z_y_zzzzz[k] = -g_x_z_0_zzzzz[k] * ab_y + g_x_z_0_yzzzzz[k];
            }

            /// Set up 168-189 components of targeted buffer : cbuffer.data(

            auto g_x_z_z_xxxxx = cbuffer.data(ph_geom_11_off + 168 * ccomps * dcomps);

            auto g_x_z_z_xxxxy = cbuffer.data(ph_geom_11_off + 169 * ccomps * dcomps);

            auto g_x_z_z_xxxxz = cbuffer.data(ph_geom_11_off + 170 * ccomps * dcomps);

            auto g_x_z_z_xxxyy = cbuffer.data(ph_geom_11_off + 171 * ccomps * dcomps);

            auto g_x_z_z_xxxyz = cbuffer.data(ph_geom_11_off + 172 * ccomps * dcomps);

            auto g_x_z_z_xxxzz = cbuffer.data(ph_geom_11_off + 173 * ccomps * dcomps);

            auto g_x_z_z_xxyyy = cbuffer.data(ph_geom_11_off + 174 * ccomps * dcomps);

            auto g_x_z_z_xxyyz = cbuffer.data(ph_geom_11_off + 175 * ccomps * dcomps);

            auto g_x_z_z_xxyzz = cbuffer.data(ph_geom_11_off + 176 * ccomps * dcomps);

            auto g_x_z_z_xxzzz = cbuffer.data(ph_geom_11_off + 177 * ccomps * dcomps);

            auto g_x_z_z_xyyyy = cbuffer.data(ph_geom_11_off + 178 * ccomps * dcomps);

            auto g_x_z_z_xyyyz = cbuffer.data(ph_geom_11_off + 179 * ccomps * dcomps);

            auto g_x_z_z_xyyzz = cbuffer.data(ph_geom_11_off + 180 * ccomps * dcomps);

            auto g_x_z_z_xyzzz = cbuffer.data(ph_geom_11_off + 181 * ccomps * dcomps);

            auto g_x_z_z_xzzzz = cbuffer.data(ph_geom_11_off + 182 * ccomps * dcomps);

            auto g_x_z_z_yyyyy = cbuffer.data(ph_geom_11_off + 183 * ccomps * dcomps);

            auto g_x_z_z_yyyyz = cbuffer.data(ph_geom_11_off + 184 * ccomps * dcomps);

            auto g_x_z_z_yyyzz = cbuffer.data(ph_geom_11_off + 185 * ccomps * dcomps);

            auto g_x_z_z_yyzzz = cbuffer.data(ph_geom_11_off + 186 * ccomps * dcomps);

            auto g_x_z_z_yzzzz = cbuffer.data(ph_geom_11_off + 187 * ccomps * dcomps);

            auto g_x_z_z_zzzzz = cbuffer.data(ph_geom_11_off + 188 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_0_xxxxx, g_x_0_0_xxxxy, g_x_0_0_xxxxz, g_x_0_0_xxxyy, g_x_0_0_xxxyz, g_x_0_0_xxxzz, g_x_0_0_xxyyy, g_x_0_0_xxyyz, g_x_0_0_xxyzz, g_x_0_0_xxzzz, g_x_0_0_xyyyy, g_x_0_0_xyyyz, g_x_0_0_xyyzz, g_x_0_0_xyzzz, g_x_0_0_xzzzz, g_x_0_0_yyyyy, g_x_0_0_yyyyz, g_x_0_0_yyyzz, g_x_0_0_yyzzz, g_x_0_0_yzzzz, g_x_0_0_zzzzz, g_x_z_0_xxxxx, g_x_z_0_xxxxxz, g_x_z_0_xxxxy, g_x_z_0_xxxxyz, g_x_z_0_xxxxz, g_x_z_0_xxxxzz, g_x_z_0_xxxyy, g_x_z_0_xxxyyz, g_x_z_0_xxxyz, g_x_z_0_xxxyzz, g_x_z_0_xxxzz, g_x_z_0_xxxzzz, g_x_z_0_xxyyy, g_x_z_0_xxyyyz, g_x_z_0_xxyyz, g_x_z_0_xxyyzz, g_x_z_0_xxyzz, g_x_z_0_xxyzzz, g_x_z_0_xxzzz, g_x_z_0_xxzzzz, g_x_z_0_xyyyy, g_x_z_0_xyyyyz, g_x_z_0_xyyyz, g_x_z_0_xyyyzz, g_x_z_0_xyyzz, g_x_z_0_xyyzzz, g_x_z_0_xyzzz, g_x_z_0_xyzzzz, g_x_z_0_xzzzz, g_x_z_0_xzzzzz, g_x_z_0_yyyyy, g_x_z_0_yyyyyz, g_x_z_0_yyyyz, g_x_z_0_yyyyzz, g_x_z_0_yyyzz, g_x_z_0_yyyzzz, g_x_z_0_yyzzz, g_x_z_0_yyzzzz, g_x_z_0_yzzzz, g_x_z_0_yzzzzz, g_x_z_0_zzzzz, g_x_z_0_zzzzzz, g_x_z_z_xxxxx, g_x_z_z_xxxxy, g_x_z_z_xxxxz, g_x_z_z_xxxyy, g_x_z_z_xxxyz, g_x_z_z_xxxzz, g_x_z_z_xxyyy, g_x_z_z_xxyyz, g_x_z_z_xxyzz, g_x_z_z_xxzzz, g_x_z_z_xyyyy, g_x_z_z_xyyyz, g_x_z_z_xyyzz, g_x_z_z_xyzzz, g_x_z_z_xzzzz, g_x_z_z_yyyyy, g_x_z_z_yyyyz, g_x_z_z_yyyzz, g_x_z_z_yyzzz, g_x_z_z_yzzzz, g_x_z_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_z_xxxxx[k] = g_x_0_0_xxxxx[k] - g_x_z_0_xxxxx[k] * ab_z + g_x_z_0_xxxxxz[k];

                g_x_z_z_xxxxy[k] = g_x_0_0_xxxxy[k] - g_x_z_0_xxxxy[k] * ab_z + g_x_z_0_xxxxyz[k];

                g_x_z_z_xxxxz[k] = g_x_0_0_xxxxz[k] - g_x_z_0_xxxxz[k] * ab_z + g_x_z_0_xxxxzz[k];

                g_x_z_z_xxxyy[k] = g_x_0_0_xxxyy[k] - g_x_z_0_xxxyy[k] * ab_z + g_x_z_0_xxxyyz[k];

                g_x_z_z_xxxyz[k] = g_x_0_0_xxxyz[k] - g_x_z_0_xxxyz[k] * ab_z + g_x_z_0_xxxyzz[k];

                g_x_z_z_xxxzz[k] = g_x_0_0_xxxzz[k] - g_x_z_0_xxxzz[k] * ab_z + g_x_z_0_xxxzzz[k];

                g_x_z_z_xxyyy[k] = g_x_0_0_xxyyy[k] - g_x_z_0_xxyyy[k] * ab_z + g_x_z_0_xxyyyz[k];

                g_x_z_z_xxyyz[k] = g_x_0_0_xxyyz[k] - g_x_z_0_xxyyz[k] * ab_z + g_x_z_0_xxyyzz[k];

                g_x_z_z_xxyzz[k] = g_x_0_0_xxyzz[k] - g_x_z_0_xxyzz[k] * ab_z + g_x_z_0_xxyzzz[k];

                g_x_z_z_xxzzz[k] = g_x_0_0_xxzzz[k] - g_x_z_0_xxzzz[k] * ab_z + g_x_z_0_xxzzzz[k];

                g_x_z_z_xyyyy[k] = g_x_0_0_xyyyy[k] - g_x_z_0_xyyyy[k] * ab_z + g_x_z_0_xyyyyz[k];

                g_x_z_z_xyyyz[k] = g_x_0_0_xyyyz[k] - g_x_z_0_xyyyz[k] * ab_z + g_x_z_0_xyyyzz[k];

                g_x_z_z_xyyzz[k] = g_x_0_0_xyyzz[k] - g_x_z_0_xyyzz[k] * ab_z + g_x_z_0_xyyzzz[k];

                g_x_z_z_xyzzz[k] = g_x_0_0_xyzzz[k] - g_x_z_0_xyzzz[k] * ab_z + g_x_z_0_xyzzzz[k];

                g_x_z_z_xzzzz[k] = g_x_0_0_xzzzz[k] - g_x_z_0_xzzzz[k] * ab_z + g_x_z_0_xzzzzz[k];

                g_x_z_z_yyyyy[k] = g_x_0_0_yyyyy[k] - g_x_z_0_yyyyy[k] * ab_z + g_x_z_0_yyyyyz[k];

                g_x_z_z_yyyyz[k] = g_x_0_0_yyyyz[k] - g_x_z_0_yyyyz[k] * ab_z + g_x_z_0_yyyyzz[k];

                g_x_z_z_yyyzz[k] = g_x_0_0_yyyzz[k] - g_x_z_0_yyyzz[k] * ab_z + g_x_z_0_yyyzzz[k];

                g_x_z_z_yyzzz[k] = g_x_0_0_yyzzz[k] - g_x_z_0_yyzzz[k] * ab_z + g_x_z_0_yyzzzz[k];

                g_x_z_z_yzzzz[k] = g_x_0_0_yzzzz[k] - g_x_z_0_yzzzz[k] * ab_z + g_x_z_0_yzzzzz[k];

                g_x_z_z_zzzzz[k] = g_x_0_0_zzzzz[k] - g_x_z_0_zzzzz[k] * ab_z + g_x_z_0_zzzzzz[k];
            }

            /// Set up 189-210 components of targeted buffer : cbuffer.data(

            auto g_y_x_x_xxxxx = cbuffer.data(ph_geom_11_off + 189 * ccomps * dcomps);

            auto g_y_x_x_xxxxy = cbuffer.data(ph_geom_11_off + 190 * ccomps * dcomps);

            auto g_y_x_x_xxxxz = cbuffer.data(ph_geom_11_off + 191 * ccomps * dcomps);

            auto g_y_x_x_xxxyy = cbuffer.data(ph_geom_11_off + 192 * ccomps * dcomps);

            auto g_y_x_x_xxxyz = cbuffer.data(ph_geom_11_off + 193 * ccomps * dcomps);

            auto g_y_x_x_xxxzz = cbuffer.data(ph_geom_11_off + 194 * ccomps * dcomps);

            auto g_y_x_x_xxyyy = cbuffer.data(ph_geom_11_off + 195 * ccomps * dcomps);

            auto g_y_x_x_xxyyz = cbuffer.data(ph_geom_11_off + 196 * ccomps * dcomps);

            auto g_y_x_x_xxyzz = cbuffer.data(ph_geom_11_off + 197 * ccomps * dcomps);

            auto g_y_x_x_xxzzz = cbuffer.data(ph_geom_11_off + 198 * ccomps * dcomps);

            auto g_y_x_x_xyyyy = cbuffer.data(ph_geom_11_off + 199 * ccomps * dcomps);

            auto g_y_x_x_xyyyz = cbuffer.data(ph_geom_11_off + 200 * ccomps * dcomps);

            auto g_y_x_x_xyyzz = cbuffer.data(ph_geom_11_off + 201 * ccomps * dcomps);

            auto g_y_x_x_xyzzz = cbuffer.data(ph_geom_11_off + 202 * ccomps * dcomps);

            auto g_y_x_x_xzzzz = cbuffer.data(ph_geom_11_off + 203 * ccomps * dcomps);

            auto g_y_x_x_yyyyy = cbuffer.data(ph_geom_11_off + 204 * ccomps * dcomps);

            auto g_y_x_x_yyyyz = cbuffer.data(ph_geom_11_off + 205 * ccomps * dcomps);

            auto g_y_x_x_yyyzz = cbuffer.data(ph_geom_11_off + 206 * ccomps * dcomps);

            auto g_y_x_x_yyzzz = cbuffer.data(ph_geom_11_off + 207 * ccomps * dcomps);

            auto g_y_x_x_yzzzz = cbuffer.data(ph_geom_11_off + 208 * ccomps * dcomps);

            auto g_y_x_x_zzzzz = cbuffer.data(ph_geom_11_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_0_xxxxx, g_y_0_0_xxxxy, g_y_0_0_xxxxz, g_y_0_0_xxxyy, g_y_0_0_xxxyz, g_y_0_0_xxxzz, g_y_0_0_xxyyy, g_y_0_0_xxyyz, g_y_0_0_xxyzz, g_y_0_0_xxzzz, g_y_0_0_xyyyy, g_y_0_0_xyyyz, g_y_0_0_xyyzz, g_y_0_0_xyzzz, g_y_0_0_xzzzz, g_y_0_0_yyyyy, g_y_0_0_yyyyz, g_y_0_0_yyyzz, g_y_0_0_yyzzz, g_y_0_0_yzzzz, g_y_0_0_zzzzz, g_y_x_0_xxxxx, g_y_x_0_xxxxxx, g_y_x_0_xxxxxy, g_y_x_0_xxxxxz, g_y_x_0_xxxxy, g_y_x_0_xxxxyy, g_y_x_0_xxxxyz, g_y_x_0_xxxxz, g_y_x_0_xxxxzz, g_y_x_0_xxxyy, g_y_x_0_xxxyyy, g_y_x_0_xxxyyz, g_y_x_0_xxxyz, g_y_x_0_xxxyzz, g_y_x_0_xxxzz, g_y_x_0_xxxzzz, g_y_x_0_xxyyy, g_y_x_0_xxyyyy, g_y_x_0_xxyyyz, g_y_x_0_xxyyz, g_y_x_0_xxyyzz, g_y_x_0_xxyzz, g_y_x_0_xxyzzz, g_y_x_0_xxzzz, g_y_x_0_xxzzzz, g_y_x_0_xyyyy, g_y_x_0_xyyyyy, g_y_x_0_xyyyyz, g_y_x_0_xyyyz, g_y_x_0_xyyyzz, g_y_x_0_xyyzz, g_y_x_0_xyyzzz, g_y_x_0_xyzzz, g_y_x_0_xyzzzz, g_y_x_0_xzzzz, g_y_x_0_xzzzzz, g_y_x_0_yyyyy, g_y_x_0_yyyyz, g_y_x_0_yyyzz, g_y_x_0_yyzzz, g_y_x_0_yzzzz, g_y_x_0_zzzzz, g_y_x_x_xxxxx, g_y_x_x_xxxxy, g_y_x_x_xxxxz, g_y_x_x_xxxyy, g_y_x_x_xxxyz, g_y_x_x_xxxzz, g_y_x_x_xxyyy, g_y_x_x_xxyyz, g_y_x_x_xxyzz, g_y_x_x_xxzzz, g_y_x_x_xyyyy, g_y_x_x_xyyyz, g_y_x_x_xyyzz, g_y_x_x_xyzzz, g_y_x_x_xzzzz, g_y_x_x_yyyyy, g_y_x_x_yyyyz, g_y_x_x_yyyzz, g_y_x_x_yyzzz, g_y_x_x_yzzzz, g_y_x_x_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_x_xxxxx[k] = g_y_0_0_xxxxx[k] - g_y_x_0_xxxxx[k] * ab_x + g_y_x_0_xxxxxx[k];

                g_y_x_x_xxxxy[k] = g_y_0_0_xxxxy[k] - g_y_x_0_xxxxy[k] * ab_x + g_y_x_0_xxxxxy[k];

                g_y_x_x_xxxxz[k] = g_y_0_0_xxxxz[k] - g_y_x_0_xxxxz[k] * ab_x + g_y_x_0_xxxxxz[k];

                g_y_x_x_xxxyy[k] = g_y_0_0_xxxyy[k] - g_y_x_0_xxxyy[k] * ab_x + g_y_x_0_xxxxyy[k];

                g_y_x_x_xxxyz[k] = g_y_0_0_xxxyz[k] - g_y_x_0_xxxyz[k] * ab_x + g_y_x_0_xxxxyz[k];

                g_y_x_x_xxxzz[k] = g_y_0_0_xxxzz[k] - g_y_x_0_xxxzz[k] * ab_x + g_y_x_0_xxxxzz[k];

                g_y_x_x_xxyyy[k] = g_y_0_0_xxyyy[k] - g_y_x_0_xxyyy[k] * ab_x + g_y_x_0_xxxyyy[k];

                g_y_x_x_xxyyz[k] = g_y_0_0_xxyyz[k] - g_y_x_0_xxyyz[k] * ab_x + g_y_x_0_xxxyyz[k];

                g_y_x_x_xxyzz[k] = g_y_0_0_xxyzz[k] - g_y_x_0_xxyzz[k] * ab_x + g_y_x_0_xxxyzz[k];

                g_y_x_x_xxzzz[k] = g_y_0_0_xxzzz[k] - g_y_x_0_xxzzz[k] * ab_x + g_y_x_0_xxxzzz[k];

                g_y_x_x_xyyyy[k] = g_y_0_0_xyyyy[k] - g_y_x_0_xyyyy[k] * ab_x + g_y_x_0_xxyyyy[k];

                g_y_x_x_xyyyz[k] = g_y_0_0_xyyyz[k] - g_y_x_0_xyyyz[k] * ab_x + g_y_x_0_xxyyyz[k];

                g_y_x_x_xyyzz[k] = g_y_0_0_xyyzz[k] - g_y_x_0_xyyzz[k] * ab_x + g_y_x_0_xxyyzz[k];

                g_y_x_x_xyzzz[k] = g_y_0_0_xyzzz[k] - g_y_x_0_xyzzz[k] * ab_x + g_y_x_0_xxyzzz[k];

                g_y_x_x_xzzzz[k] = g_y_0_0_xzzzz[k] - g_y_x_0_xzzzz[k] * ab_x + g_y_x_0_xxzzzz[k];

                g_y_x_x_yyyyy[k] = g_y_0_0_yyyyy[k] - g_y_x_0_yyyyy[k] * ab_x + g_y_x_0_xyyyyy[k];

                g_y_x_x_yyyyz[k] = g_y_0_0_yyyyz[k] - g_y_x_0_yyyyz[k] * ab_x + g_y_x_0_xyyyyz[k];

                g_y_x_x_yyyzz[k] = g_y_0_0_yyyzz[k] - g_y_x_0_yyyzz[k] * ab_x + g_y_x_0_xyyyzz[k];

                g_y_x_x_yyzzz[k] = g_y_0_0_yyzzz[k] - g_y_x_0_yyzzz[k] * ab_x + g_y_x_0_xyyzzz[k];

                g_y_x_x_yzzzz[k] = g_y_0_0_yzzzz[k] - g_y_x_0_yzzzz[k] * ab_x + g_y_x_0_xyzzzz[k];

                g_y_x_x_zzzzz[k] = g_y_0_0_zzzzz[k] - g_y_x_0_zzzzz[k] * ab_x + g_y_x_0_xzzzzz[k];
            }

            /// Set up 210-231 components of targeted buffer : cbuffer.data(

            auto g_y_x_y_xxxxx = cbuffer.data(ph_geom_11_off + 210 * ccomps * dcomps);

            auto g_y_x_y_xxxxy = cbuffer.data(ph_geom_11_off + 211 * ccomps * dcomps);

            auto g_y_x_y_xxxxz = cbuffer.data(ph_geom_11_off + 212 * ccomps * dcomps);

            auto g_y_x_y_xxxyy = cbuffer.data(ph_geom_11_off + 213 * ccomps * dcomps);

            auto g_y_x_y_xxxyz = cbuffer.data(ph_geom_11_off + 214 * ccomps * dcomps);

            auto g_y_x_y_xxxzz = cbuffer.data(ph_geom_11_off + 215 * ccomps * dcomps);

            auto g_y_x_y_xxyyy = cbuffer.data(ph_geom_11_off + 216 * ccomps * dcomps);

            auto g_y_x_y_xxyyz = cbuffer.data(ph_geom_11_off + 217 * ccomps * dcomps);

            auto g_y_x_y_xxyzz = cbuffer.data(ph_geom_11_off + 218 * ccomps * dcomps);

            auto g_y_x_y_xxzzz = cbuffer.data(ph_geom_11_off + 219 * ccomps * dcomps);

            auto g_y_x_y_xyyyy = cbuffer.data(ph_geom_11_off + 220 * ccomps * dcomps);

            auto g_y_x_y_xyyyz = cbuffer.data(ph_geom_11_off + 221 * ccomps * dcomps);

            auto g_y_x_y_xyyzz = cbuffer.data(ph_geom_11_off + 222 * ccomps * dcomps);

            auto g_y_x_y_xyzzz = cbuffer.data(ph_geom_11_off + 223 * ccomps * dcomps);

            auto g_y_x_y_xzzzz = cbuffer.data(ph_geom_11_off + 224 * ccomps * dcomps);

            auto g_y_x_y_yyyyy = cbuffer.data(ph_geom_11_off + 225 * ccomps * dcomps);

            auto g_y_x_y_yyyyz = cbuffer.data(ph_geom_11_off + 226 * ccomps * dcomps);

            auto g_y_x_y_yyyzz = cbuffer.data(ph_geom_11_off + 227 * ccomps * dcomps);

            auto g_y_x_y_yyzzz = cbuffer.data(ph_geom_11_off + 228 * ccomps * dcomps);

            auto g_y_x_y_yzzzz = cbuffer.data(ph_geom_11_off + 229 * ccomps * dcomps);

            auto g_y_x_y_zzzzz = cbuffer.data(ph_geom_11_off + 230 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xxxxx, g_0_x_0_xxxxy, g_0_x_0_xxxxz, g_0_x_0_xxxyy, g_0_x_0_xxxyz, g_0_x_0_xxxzz, g_0_x_0_xxyyy, g_0_x_0_xxyyz, g_0_x_0_xxyzz, g_0_x_0_xxzzz, g_0_x_0_xyyyy, g_0_x_0_xyyyz, g_0_x_0_xyyzz, g_0_x_0_xyzzz, g_0_x_0_xzzzz, g_0_x_0_yyyyy, g_0_x_0_yyyyz, g_0_x_0_yyyzz, g_0_x_0_yyzzz, g_0_x_0_yzzzz, g_0_x_0_zzzzz, g_y_x_0_xxxxx, g_y_x_0_xxxxxy, g_y_x_0_xxxxy, g_y_x_0_xxxxyy, g_y_x_0_xxxxyz, g_y_x_0_xxxxz, g_y_x_0_xxxyy, g_y_x_0_xxxyyy, g_y_x_0_xxxyyz, g_y_x_0_xxxyz, g_y_x_0_xxxyzz, g_y_x_0_xxxzz, g_y_x_0_xxyyy, g_y_x_0_xxyyyy, g_y_x_0_xxyyyz, g_y_x_0_xxyyz, g_y_x_0_xxyyzz, g_y_x_0_xxyzz, g_y_x_0_xxyzzz, g_y_x_0_xxzzz, g_y_x_0_xyyyy, g_y_x_0_xyyyyy, g_y_x_0_xyyyyz, g_y_x_0_xyyyz, g_y_x_0_xyyyzz, g_y_x_0_xyyzz, g_y_x_0_xyyzzz, g_y_x_0_xyzzz, g_y_x_0_xyzzzz, g_y_x_0_xzzzz, g_y_x_0_yyyyy, g_y_x_0_yyyyyy, g_y_x_0_yyyyyz, g_y_x_0_yyyyz, g_y_x_0_yyyyzz, g_y_x_0_yyyzz, g_y_x_0_yyyzzz, g_y_x_0_yyzzz, g_y_x_0_yyzzzz, g_y_x_0_yzzzz, g_y_x_0_yzzzzz, g_y_x_0_zzzzz, g_y_x_y_xxxxx, g_y_x_y_xxxxy, g_y_x_y_xxxxz, g_y_x_y_xxxyy, g_y_x_y_xxxyz, g_y_x_y_xxxzz, g_y_x_y_xxyyy, g_y_x_y_xxyyz, g_y_x_y_xxyzz, g_y_x_y_xxzzz, g_y_x_y_xyyyy, g_y_x_y_xyyyz, g_y_x_y_xyyzz, g_y_x_y_xyzzz, g_y_x_y_xzzzz, g_y_x_y_yyyyy, g_y_x_y_yyyyz, g_y_x_y_yyyzz, g_y_x_y_yyzzz, g_y_x_y_yzzzz, g_y_x_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_y_xxxxx[k] = -g_0_x_0_xxxxx[k] - g_y_x_0_xxxxx[k] * ab_y + g_y_x_0_xxxxxy[k];

                g_y_x_y_xxxxy[k] = -g_0_x_0_xxxxy[k] - g_y_x_0_xxxxy[k] * ab_y + g_y_x_0_xxxxyy[k];

                g_y_x_y_xxxxz[k] = -g_0_x_0_xxxxz[k] - g_y_x_0_xxxxz[k] * ab_y + g_y_x_0_xxxxyz[k];

                g_y_x_y_xxxyy[k] = -g_0_x_0_xxxyy[k] - g_y_x_0_xxxyy[k] * ab_y + g_y_x_0_xxxyyy[k];

                g_y_x_y_xxxyz[k] = -g_0_x_0_xxxyz[k] - g_y_x_0_xxxyz[k] * ab_y + g_y_x_0_xxxyyz[k];

                g_y_x_y_xxxzz[k] = -g_0_x_0_xxxzz[k] - g_y_x_0_xxxzz[k] * ab_y + g_y_x_0_xxxyzz[k];

                g_y_x_y_xxyyy[k] = -g_0_x_0_xxyyy[k] - g_y_x_0_xxyyy[k] * ab_y + g_y_x_0_xxyyyy[k];

                g_y_x_y_xxyyz[k] = -g_0_x_0_xxyyz[k] - g_y_x_0_xxyyz[k] * ab_y + g_y_x_0_xxyyyz[k];

                g_y_x_y_xxyzz[k] = -g_0_x_0_xxyzz[k] - g_y_x_0_xxyzz[k] * ab_y + g_y_x_0_xxyyzz[k];

                g_y_x_y_xxzzz[k] = -g_0_x_0_xxzzz[k] - g_y_x_0_xxzzz[k] * ab_y + g_y_x_0_xxyzzz[k];

                g_y_x_y_xyyyy[k] = -g_0_x_0_xyyyy[k] - g_y_x_0_xyyyy[k] * ab_y + g_y_x_0_xyyyyy[k];

                g_y_x_y_xyyyz[k] = -g_0_x_0_xyyyz[k] - g_y_x_0_xyyyz[k] * ab_y + g_y_x_0_xyyyyz[k];

                g_y_x_y_xyyzz[k] = -g_0_x_0_xyyzz[k] - g_y_x_0_xyyzz[k] * ab_y + g_y_x_0_xyyyzz[k];

                g_y_x_y_xyzzz[k] = -g_0_x_0_xyzzz[k] - g_y_x_0_xyzzz[k] * ab_y + g_y_x_0_xyyzzz[k];

                g_y_x_y_xzzzz[k] = -g_0_x_0_xzzzz[k] - g_y_x_0_xzzzz[k] * ab_y + g_y_x_0_xyzzzz[k];

                g_y_x_y_yyyyy[k] = -g_0_x_0_yyyyy[k] - g_y_x_0_yyyyy[k] * ab_y + g_y_x_0_yyyyyy[k];

                g_y_x_y_yyyyz[k] = -g_0_x_0_yyyyz[k] - g_y_x_0_yyyyz[k] * ab_y + g_y_x_0_yyyyyz[k];

                g_y_x_y_yyyzz[k] = -g_0_x_0_yyyzz[k] - g_y_x_0_yyyzz[k] * ab_y + g_y_x_0_yyyyzz[k];

                g_y_x_y_yyzzz[k] = -g_0_x_0_yyzzz[k] - g_y_x_0_yyzzz[k] * ab_y + g_y_x_0_yyyzzz[k];

                g_y_x_y_yzzzz[k] = -g_0_x_0_yzzzz[k] - g_y_x_0_yzzzz[k] * ab_y + g_y_x_0_yyzzzz[k];

                g_y_x_y_zzzzz[k] = -g_0_x_0_zzzzz[k] - g_y_x_0_zzzzz[k] * ab_y + g_y_x_0_yzzzzz[k];
            }

            /// Set up 231-252 components of targeted buffer : cbuffer.data(

            auto g_y_x_z_xxxxx = cbuffer.data(ph_geom_11_off + 231 * ccomps * dcomps);

            auto g_y_x_z_xxxxy = cbuffer.data(ph_geom_11_off + 232 * ccomps * dcomps);

            auto g_y_x_z_xxxxz = cbuffer.data(ph_geom_11_off + 233 * ccomps * dcomps);

            auto g_y_x_z_xxxyy = cbuffer.data(ph_geom_11_off + 234 * ccomps * dcomps);

            auto g_y_x_z_xxxyz = cbuffer.data(ph_geom_11_off + 235 * ccomps * dcomps);

            auto g_y_x_z_xxxzz = cbuffer.data(ph_geom_11_off + 236 * ccomps * dcomps);

            auto g_y_x_z_xxyyy = cbuffer.data(ph_geom_11_off + 237 * ccomps * dcomps);

            auto g_y_x_z_xxyyz = cbuffer.data(ph_geom_11_off + 238 * ccomps * dcomps);

            auto g_y_x_z_xxyzz = cbuffer.data(ph_geom_11_off + 239 * ccomps * dcomps);

            auto g_y_x_z_xxzzz = cbuffer.data(ph_geom_11_off + 240 * ccomps * dcomps);

            auto g_y_x_z_xyyyy = cbuffer.data(ph_geom_11_off + 241 * ccomps * dcomps);

            auto g_y_x_z_xyyyz = cbuffer.data(ph_geom_11_off + 242 * ccomps * dcomps);

            auto g_y_x_z_xyyzz = cbuffer.data(ph_geom_11_off + 243 * ccomps * dcomps);

            auto g_y_x_z_xyzzz = cbuffer.data(ph_geom_11_off + 244 * ccomps * dcomps);

            auto g_y_x_z_xzzzz = cbuffer.data(ph_geom_11_off + 245 * ccomps * dcomps);

            auto g_y_x_z_yyyyy = cbuffer.data(ph_geom_11_off + 246 * ccomps * dcomps);

            auto g_y_x_z_yyyyz = cbuffer.data(ph_geom_11_off + 247 * ccomps * dcomps);

            auto g_y_x_z_yyyzz = cbuffer.data(ph_geom_11_off + 248 * ccomps * dcomps);

            auto g_y_x_z_yyzzz = cbuffer.data(ph_geom_11_off + 249 * ccomps * dcomps);

            auto g_y_x_z_yzzzz = cbuffer.data(ph_geom_11_off + 250 * ccomps * dcomps);

            auto g_y_x_z_zzzzz = cbuffer.data(ph_geom_11_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_0_xxxxx, g_y_x_0_xxxxxz, g_y_x_0_xxxxy, g_y_x_0_xxxxyz, g_y_x_0_xxxxz, g_y_x_0_xxxxzz, g_y_x_0_xxxyy, g_y_x_0_xxxyyz, g_y_x_0_xxxyz, g_y_x_0_xxxyzz, g_y_x_0_xxxzz, g_y_x_0_xxxzzz, g_y_x_0_xxyyy, g_y_x_0_xxyyyz, g_y_x_0_xxyyz, g_y_x_0_xxyyzz, g_y_x_0_xxyzz, g_y_x_0_xxyzzz, g_y_x_0_xxzzz, g_y_x_0_xxzzzz, g_y_x_0_xyyyy, g_y_x_0_xyyyyz, g_y_x_0_xyyyz, g_y_x_0_xyyyzz, g_y_x_0_xyyzz, g_y_x_0_xyyzzz, g_y_x_0_xyzzz, g_y_x_0_xyzzzz, g_y_x_0_xzzzz, g_y_x_0_xzzzzz, g_y_x_0_yyyyy, g_y_x_0_yyyyyz, g_y_x_0_yyyyz, g_y_x_0_yyyyzz, g_y_x_0_yyyzz, g_y_x_0_yyyzzz, g_y_x_0_yyzzz, g_y_x_0_yyzzzz, g_y_x_0_yzzzz, g_y_x_0_yzzzzz, g_y_x_0_zzzzz, g_y_x_0_zzzzzz, g_y_x_z_xxxxx, g_y_x_z_xxxxy, g_y_x_z_xxxxz, g_y_x_z_xxxyy, g_y_x_z_xxxyz, g_y_x_z_xxxzz, g_y_x_z_xxyyy, g_y_x_z_xxyyz, g_y_x_z_xxyzz, g_y_x_z_xxzzz, g_y_x_z_xyyyy, g_y_x_z_xyyyz, g_y_x_z_xyyzz, g_y_x_z_xyzzz, g_y_x_z_xzzzz, g_y_x_z_yyyyy, g_y_x_z_yyyyz, g_y_x_z_yyyzz, g_y_x_z_yyzzz, g_y_x_z_yzzzz, g_y_x_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_z_xxxxx[k] = -g_y_x_0_xxxxx[k] * ab_z + g_y_x_0_xxxxxz[k];

                g_y_x_z_xxxxy[k] = -g_y_x_0_xxxxy[k] * ab_z + g_y_x_0_xxxxyz[k];

                g_y_x_z_xxxxz[k] = -g_y_x_0_xxxxz[k] * ab_z + g_y_x_0_xxxxzz[k];

                g_y_x_z_xxxyy[k] = -g_y_x_0_xxxyy[k] * ab_z + g_y_x_0_xxxyyz[k];

                g_y_x_z_xxxyz[k] = -g_y_x_0_xxxyz[k] * ab_z + g_y_x_0_xxxyzz[k];

                g_y_x_z_xxxzz[k] = -g_y_x_0_xxxzz[k] * ab_z + g_y_x_0_xxxzzz[k];

                g_y_x_z_xxyyy[k] = -g_y_x_0_xxyyy[k] * ab_z + g_y_x_0_xxyyyz[k];

                g_y_x_z_xxyyz[k] = -g_y_x_0_xxyyz[k] * ab_z + g_y_x_0_xxyyzz[k];

                g_y_x_z_xxyzz[k] = -g_y_x_0_xxyzz[k] * ab_z + g_y_x_0_xxyzzz[k];

                g_y_x_z_xxzzz[k] = -g_y_x_0_xxzzz[k] * ab_z + g_y_x_0_xxzzzz[k];

                g_y_x_z_xyyyy[k] = -g_y_x_0_xyyyy[k] * ab_z + g_y_x_0_xyyyyz[k];

                g_y_x_z_xyyyz[k] = -g_y_x_0_xyyyz[k] * ab_z + g_y_x_0_xyyyzz[k];

                g_y_x_z_xyyzz[k] = -g_y_x_0_xyyzz[k] * ab_z + g_y_x_0_xyyzzz[k];

                g_y_x_z_xyzzz[k] = -g_y_x_0_xyzzz[k] * ab_z + g_y_x_0_xyzzzz[k];

                g_y_x_z_xzzzz[k] = -g_y_x_0_xzzzz[k] * ab_z + g_y_x_0_xzzzzz[k];

                g_y_x_z_yyyyy[k] = -g_y_x_0_yyyyy[k] * ab_z + g_y_x_0_yyyyyz[k];

                g_y_x_z_yyyyz[k] = -g_y_x_0_yyyyz[k] * ab_z + g_y_x_0_yyyyzz[k];

                g_y_x_z_yyyzz[k] = -g_y_x_0_yyyzz[k] * ab_z + g_y_x_0_yyyzzz[k];

                g_y_x_z_yyzzz[k] = -g_y_x_0_yyzzz[k] * ab_z + g_y_x_0_yyzzzz[k];

                g_y_x_z_yzzzz[k] = -g_y_x_0_yzzzz[k] * ab_z + g_y_x_0_yzzzzz[k];

                g_y_x_z_zzzzz[k] = -g_y_x_0_zzzzz[k] * ab_z + g_y_x_0_zzzzzz[k];
            }

            /// Set up 252-273 components of targeted buffer : cbuffer.data(

            auto g_y_y_x_xxxxx = cbuffer.data(ph_geom_11_off + 252 * ccomps * dcomps);

            auto g_y_y_x_xxxxy = cbuffer.data(ph_geom_11_off + 253 * ccomps * dcomps);

            auto g_y_y_x_xxxxz = cbuffer.data(ph_geom_11_off + 254 * ccomps * dcomps);

            auto g_y_y_x_xxxyy = cbuffer.data(ph_geom_11_off + 255 * ccomps * dcomps);

            auto g_y_y_x_xxxyz = cbuffer.data(ph_geom_11_off + 256 * ccomps * dcomps);

            auto g_y_y_x_xxxzz = cbuffer.data(ph_geom_11_off + 257 * ccomps * dcomps);

            auto g_y_y_x_xxyyy = cbuffer.data(ph_geom_11_off + 258 * ccomps * dcomps);

            auto g_y_y_x_xxyyz = cbuffer.data(ph_geom_11_off + 259 * ccomps * dcomps);

            auto g_y_y_x_xxyzz = cbuffer.data(ph_geom_11_off + 260 * ccomps * dcomps);

            auto g_y_y_x_xxzzz = cbuffer.data(ph_geom_11_off + 261 * ccomps * dcomps);

            auto g_y_y_x_xyyyy = cbuffer.data(ph_geom_11_off + 262 * ccomps * dcomps);

            auto g_y_y_x_xyyyz = cbuffer.data(ph_geom_11_off + 263 * ccomps * dcomps);

            auto g_y_y_x_xyyzz = cbuffer.data(ph_geom_11_off + 264 * ccomps * dcomps);

            auto g_y_y_x_xyzzz = cbuffer.data(ph_geom_11_off + 265 * ccomps * dcomps);

            auto g_y_y_x_xzzzz = cbuffer.data(ph_geom_11_off + 266 * ccomps * dcomps);

            auto g_y_y_x_yyyyy = cbuffer.data(ph_geom_11_off + 267 * ccomps * dcomps);

            auto g_y_y_x_yyyyz = cbuffer.data(ph_geom_11_off + 268 * ccomps * dcomps);

            auto g_y_y_x_yyyzz = cbuffer.data(ph_geom_11_off + 269 * ccomps * dcomps);

            auto g_y_y_x_yyzzz = cbuffer.data(ph_geom_11_off + 270 * ccomps * dcomps);

            auto g_y_y_x_yzzzz = cbuffer.data(ph_geom_11_off + 271 * ccomps * dcomps);

            auto g_y_y_x_zzzzz = cbuffer.data(ph_geom_11_off + 272 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_0_xxxxx, g_y_y_0_xxxxxx, g_y_y_0_xxxxxy, g_y_y_0_xxxxxz, g_y_y_0_xxxxy, g_y_y_0_xxxxyy, g_y_y_0_xxxxyz, g_y_y_0_xxxxz, g_y_y_0_xxxxzz, g_y_y_0_xxxyy, g_y_y_0_xxxyyy, g_y_y_0_xxxyyz, g_y_y_0_xxxyz, g_y_y_0_xxxyzz, g_y_y_0_xxxzz, g_y_y_0_xxxzzz, g_y_y_0_xxyyy, g_y_y_0_xxyyyy, g_y_y_0_xxyyyz, g_y_y_0_xxyyz, g_y_y_0_xxyyzz, g_y_y_0_xxyzz, g_y_y_0_xxyzzz, g_y_y_0_xxzzz, g_y_y_0_xxzzzz, g_y_y_0_xyyyy, g_y_y_0_xyyyyy, g_y_y_0_xyyyyz, g_y_y_0_xyyyz, g_y_y_0_xyyyzz, g_y_y_0_xyyzz, g_y_y_0_xyyzzz, g_y_y_0_xyzzz, g_y_y_0_xyzzzz, g_y_y_0_xzzzz, g_y_y_0_xzzzzz, g_y_y_0_yyyyy, g_y_y_0_yyyyz, g_y_y_0_yyyzz, g_y_y_0_yyzzz, g_y_y_0_yzzzz, g_y_y_0_zzzzz, g_y_y_x_xxxxx, g_y_y_x_xxxxy, g_y_y_x_xxxxz, g_y_y_x_xxxyy, g_y_y_x_xxxyz, g_y_y_x_xxxzz, g_y_y_x_xxyyy, g_y_y_x_xxyyz, g_y_y_x_xxyzz, g_y_y_x_xxzzz, g_y_y_x_xyyyy, g_y_y_x_xyyyz, g_y_y_x_xyyzz, g_y_y_x_xyzzz, g_y_y_x_xzzzz, g_y_y_x_yyyyy, g_y_y_x_yyyyz, g_y_y_x_yyyzz, g_y_y_x_yyzzz, g_y_y_x_yzzzz, g_y_y_x_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_x_xxxxx[k] = -g_y_y_0_xxxxx[k] * ab_x + g_y_y_0_xxxxxx[k];

                g_y_y_x_xxxxy[k] = -g_y_y_0_xxxxy[k] * ab_x + g_y_y_0_xxxxxy[k];

                g_y_y_x_xxxxz[k] = -g_y_y_0_xxxxz[k] * ab_x + g_y_y_0_xxxxxz[k];

                g_y_y_x_xxxyy[k] = -g_y_y_0_xxxyy[k] * ab_x + g_y_y_0_xxxxyy[k];

                g_y_y_x_xxxyz[k] = -g_y_y_0_xxxyz[k] * ab_x + g_y_y_0_xxxxyz[k];

                g_y_y_x_xxxzz[k] = -g_y_y_0_xxxzz[k] * ab_x + g_y_y_0_xxxxzz[k];

                g_y_y_x_xxyyy[k] = -g_y_y_0_xxyyy[k] * ab_x + g_y_y_0_xxxyyy[k];

                g_y_y_x_xxyyz[k] = -g_y_y_0_xxyyz[k] * ab_x + g_y_y_0_xxxyyz[k];

                g_y_y_x_xxyzz[k] = -g_y_y_0_xxyzz[k] * ab_x + g_y_y_0_xxxyzz[k];

                g_y_y_x_xxzzz[k] = -g_y_y_0_xxzzz[k] * ab_x + g_y_y_0_xxxzzz[k];

                g_y_y_x_xyyyy[k] = -g_y_y_0_xyyyy[k] * ab_x + g_y_y_0_xxyyyy[k];

                g_y_y_x_xyyyz[k] = -g_y_y_0_xyyyz[k] * ab_x + g_y_y_0_xxyyyz[k];

                g_y_y_x_xyyzz[k] = -g_y_y_0_xyyzz[k] * ab_x + g_y_y_0_xxyyzz[k];

                g_y_y_x_xyzzz[k] = -g_y_y_0_xyzzz[k] * ab_x + g_y_y_0_xxyzzz[k];

                g_y_y_x_xzzzz[k] = -g_y_y_0_xzzzz[k] * ab_x + g_y_y_0_xxzzzz[k];

                g_y_y_x_yyyyy[k] = -g_y_y_0_yyyyy[k] * ab_x + g_y_y_0_xyyyyy[k];

                g_y_y_x_yyyyz[k] = -g_y_y_0_yyyyz[k] * ab_x + g_y_y_0_xyyyyz[k];

                g_y_y_x_yyyzz[k] = -g_y_y_0_yyyzz[k] * ab_x + g_y_y_0_xyyyzz[k];

                g_y_y_x_yyzzz[k] = -g_y_y_0_yyzzz[k] * ab_x + g_y_y_0_xyyzzz[k];

                g_y_y_x_yzzzz[k] = -g_y_y_0_yzzzz[k] * ab_x + g_y_y_0_xyzzzz[k];

                g_y_y_x_zzzzz[k] = -g_y_y_0_zzzzz[k] * ab_x + g_y_y_0_xzzzzz[k];
            }

            /// Set up 273-294 components of targeted buffer : cbuffer.data(

            auto g_y_y_y_xxxxx = cbuffer.data(ph_geom_11_off + 273 * ccomps * dcomps);

            auto g_y_y_y_xxxxy = cbuffer.data(ph_geom_11_off + 274 * ccomps * dcomps);

            auto g_y_y_y_xxxxz = cbuffer.data(ph_geom_11_off + 275 * ccomps * dcomps);

            auto g_y_y_y_xxxyy = cbuffer.data(ph_geom_11_off + 276 * ccomps * dcomps);

            auto g_y_y_y_xxxyz = cbuffer.data(ph_geom_11_off + 277 * ccomps * dcomps);

            auto g_y_y_y_xxxzz = cbuffer.data(ph_geom_11_off + 278 * ccomps * dcomps);

            auto g_y_y_y_xxyyy = cbuffer.data(ph_geom_11_off + 279 * ccomps * dcomps);

            auto g_y_y_y_xxyyz = cbuffer.data(ph_geom_11_off + 280 * ccomps * dcomps);

            auto g_y_y_y_xxyzz = cbuffer.data(ph_geom_11_off + 281 * ccomps * dcomps);

            auto g_y_y_y_xxzzz = cbuffer.data(ph_geom_11_off + 282 * ccomps * dcomps);

            auto g_y_y_y_xyyyy = cbuffer.data(ph_geom_11_off + 283 * ccomps * dcomps);

            auto g_y_y_y_xyyyz = cbuffer.data(ph_geom_11_off + 284 * ccomps * dcomps);

            auto g_y_y_y_xyyzz = cbuffer.data(ph_geom_11_off + 285 * ccomps * dcomps);

            auto g_y_y_y_xyzzz = cbuffer.data(ph_geom_11_off + 286 * ccomps * dcomps);

            auto g_y_y_y_xzzzz = cbuffer.data(ph_geom_11_off + 287 * ccomps * dcomps);

            auto g_y_y_y_yyyyy = cbuffer.data(ph_geom_11_off + 288 * ccomps * dcomps);

            auto g_y_y_y_yyyyz = cbuffer.data(ph_geom_11_off + 289 * ccomps * dcomps);

            auto g_y_y_y_yyyzz = cbuffer.data(ph_geom_11_off + 290 * ccomps * dcomps);

            auto g_y_y_y_yyzzz = cbuffer.data(ph_geom_11_off + 291 * ccomps * dcomps);

            auto g_y_y_y_yzzzz = cbuffer.data(ph_geom_11_off + 292 * ccomps * dcomps);

            auto g_y_y_y_zzzzz = cbuffer.data(ph_geom_11_off + 293 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_0_xxxxx, g_0_y_0_xxxxy, g_0_y_0_xxxxz, g_0_y_0_xxxyy, g_0_y_0_xxxyz, g_0_y_0_xxxzz, g_0_y_0_xxyyy, g_0_y_0_xxyyz, g_0_y_0_xxyzz, g_0_y_0_xxzzz, g_0_y_0_xyyyy, g_0_y_0_xyyyz, g_0_y_0_xyyzz, g_0_y_0_xyzzz, g_0_y_0_xzzzz, g_0_y_0_yyyyy, g_0_y_0_yyyyz, g_0_y_0_yyyzz, g_0_y_0_yyzzz, g_0_y_0_yzzzz, g_0_y_0_zzzzz, g_y_0_0_xxxxx, g_y_0_0_xxxxy, g_y_0_0_xxxxz, g_y_0_0_xxxyy, g_y_0_0_xxxyz, g_y_0_0_xxxzz, g_y_0_0_xxyyy, g_y_0_0_xxyyz, g_y_0_0_xxyzz, g_y_0_0_xxzzz, g_y_0_0_xyyyy, g_y_0_0_xyyyz, g_y_0_0_xyyzz, g_y_0_0_xyzzz, g_y_0_0_xzzzz, g_y_0_0_yyyyy, g_y_0_0_yyyyz, g_y_0_0_yyyzz, g_y_0_0_yyzzz, g_y_0_0_yzzzz, g_y_0_0_zzzzz, g_y_y_0_xxxxx, g_y_y_0_xxxxxy, g_y_y_0_xxxxy, g_y_y_0_xxxxyy, g_y_y_0_xxxxyz, g_y_y_0_xxxxz, g_y_y_0_xxxyy, g_y_y_0_xxxyyy, g_y_y_0_xxxyyz, g_y_y_0_xxxyz, g_y_y_0_xxxyzz, g_y_y_0_xxxzz, g_y_y_0_xxyyy, g_y_y_0_xxyyyy, g_y_y_0_xxyyyz, g_y_y_0_xxyyz, g_y_y_0_xxyyzz, g_y_y_0_xxyzz, g_y_y_0_xxyzzz, g_y_y_0_xxzzz, g_y_y_0_xyyyy, g_y_y_0_xyyyyy, g_y_y_0_xyyyyz, g_y_y_0_xyyyz, g_y_y_0_xyyyzz, g_y_y_0_xyyzz, g_y_y_0_xyyzzz, g_y_y_0_xyzzz, g_y_y_0_xyzzzz, g_y_y_0_xzzzz, g_y_y_0_yyyyy, g_y_y_0_yyyyyy, g_y_y_0_yyyyyz, g_y_y_0_yyyyz, g_y_y_0_yyyyzz, g_y_y_0_yyyzz, g_y_y_0_yyyzzz, g_y_y_0_yyzzz, g_y_y_0_yyzzzz, g_y_y_0_yzzzz, g_y_y_0_yzzzzz, g_y_y_0_zzzzz, g_y_y_y_xxxxx, g_y_y_y_xxxxy, g_y_y_y_xxxxz, g_y_y_y_xxxyy, g_y_y_y_xxxyz, g_y_y_y_xxxzz, g_y_y_y_xxyyy, g_y_y_y_xxyyz, g_y_y_y_xxyzz, g_y_y_y_xxzzz, g_y_y_y_xyyyy, g_y_y_y_xyyyz, g_y_y_y_xyyzz, g_y_y_y_xyzzz, g_y_y_y_xzzzz, g_y_y_y_yyyyy, g_y_y_y_yyyyz, g_y_y_y_yyyzz, g_y_y_y_yyzzz, g_y_y_y_yzzzz, g_y_y_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_y_xxxxx[k] = -g_0_y_0_xxxxx[k] + g_y_0_0_xxxxx[k] - g_y_y_0_xxxxx[k] * ab_y + g_y_y_0_xxxxxy[k];

                g_y_y_y_xxxxy[k] = -g_0_y_0_xxxxy[k] + g_y_0_0_xxxxy[k] - g_y_y_0_xxxxy[k] * ab_y + g_y_y_0_xxxxyy[k];

                g_y_y_y_xxxxz[k] = -g_0_y_0_xxxxz[k] + g_y_0_0_xxxxz[k] - g_y_y_0_xxxxz[k] * ab_y + g_y_y_0_xxxxyz[k];

                g_y_y_y_xxxyy[k] = -g_0_y_0_xxxyy[k] + g_y_0_0_xxxyy[k] - g_y_y_0_xxxyy[k] * ab_y + g_y_y_0_xxxyyy[k];

                g_y_y_y_xxxyz[k] = -g_0_y_0_xxxyz[k] + g_y_0_0_xxxyz[k] - g_y_y_0_xxxyz[k] * ab_y + g_y_y_0_xxxyyz[k];

                g_y_y_y_xxxzz[k] = -g_0_y_0_xxxzz[k] + g_y_0_0_xxxzz[k] - g_y_y_0_xxxzz[k] * ab_y + g_y_y_0_xxxyzz[k];

                g_y_y_y_xxyyy[k] = -g_0_y_0_xxyyy[k] + g_y_0_0_xxyyy[k] - g_y_y_0_xxyyy[k] * ab_y + g_y_y_0_xxyyyy[k];

                g_y_y_y_xxyyz[k] = -g_0_y_0_xxyyz[k] + g_y_0_0_xxyyz[k] - g_y_y_0_xxyyz[k] * ab_y + g_y_y_0_xxyyyz[k];

                g_y_y_y_xxyzz[k] = -g_0_y_0_xxyzz[k] + g_y_0_0_xxyzz[k] - g_y_y_0_xxyzz[k] * ab_y + g_y_y_0_xxyyzz[k];

                g_y_y_y_xxzzz[k] = -g_0_y_0_xxzzz[k] + g_y_0_0_xxzzz[k] - g_y_y_0_xxzzz[k] * ab_y + g_y_y_0_xxyzzz[k];

                g_y_y_y_xyyyy[k] = -g_0_y_0_xyyyy[k] + g_y_0_0_xyyyy[k] - g_y_y_0_xyyyy[k] * ab_y + g_y_y_0_xyyyyy[k];

                g_y_y_y_xyyyz[k] = -g_0_y_0_xyyyz[k] + g_y_0_0_xyyyz[k] - g_y_y_0_xyyyz[k] * ab_y + g_y_y_0_xyyyyz[k];

                g_y_y_y_xyyzz[k] = -g_0_y_0_xyyzz[k] + g_y_0_0_xyyzz[k] - g_y_y_0_xyyzz[k] * ab_y + g_y_y_0_xyyyzz[k];

                g_y_y_y_xyzzz[k] = -g_0_y_0_xyzzz[k] + g_y_0_0_xyzzz[k] - g_y_y_0_xyzzz[k] * ab_y + g_y_y_0_xyyzzz[k];

                g_y_y_y_xzzzz[k] = -g_0_y_0_xzzzz[k] + g_y_0_0_xzzzz[k] - g_y_y_0_xzzzz[k] * ab_y + g_y_y_0_xyzzzz[k];

                g_y_y_y_yyyyy[k] = -g_0_y_0_yyyyy[k] + g_y_0_0_yyyyy[k] - g_y_y_0_yyyyy[k] * ab_y + g_y_y_0_yyyyyy[k];

                g_y_y_y_yyyyz[k] = -g_0_y_0_yyyyz[k] + g_y_0_0_yyyyz[k] - g_y_y_0_yyyyz[k] * ab_y + g_y_y_0_yyyyyz[k];

                g_y_y_y_yyyzz[k] = -g_0_y_0_yyyzz[k] + g_y_0_0_yyyzz[k] - g_y_y_0_yyyzz[k] * ab_y + g_y_y_0_yyyyzz[k];

                g_y_y_y_yyzzz[k] = -g_0_y_0_yyzzz[k] + g_y_0_0_yyzzz[k] - g_y_y_0_yyzzz[k] * ab_y + g_y_y_0_yyyzzz[k];

                g_y_y_y_yzzzz[k] = -g_0_y_0_yzzzz[k] + g_y_0_0_yzzzz[k] - g_y_y_0_yzzzz[k] * ab_y + g_y_y_0_yyzzzz[k];

                g_y_y_y_zzzzz[k] = -g_0_y_0_zzzzz[k] + g_y_0_0_zzzzz[k] - g_y_y_0_zzzzz[k] * ab_y + g_y_y_0_yzzzzz[k];
            }

            /// Set up 294-315 components of targeted buffer : cbuffer.data(

            auto g_y_y_z_xxxxx = cbuffer.data(ph_geom_11_off + 294 * ccomps * dcomps);

            auto g_y_y_z_xxxxy = cbuffer.data(ph_geom_11_off + 295 * ccomps * dcomps);

            auto g_y_y_z_xxxxz = cbuffer.data(ph_geom_11_off + 296 * ccomps * dcomps);

            auto g_y_y_z_xxxyy = cbuffer.data(ph_geom_11_off + 297 * ccomps * dcomps);

            auto g_y_y_z_xxxyz = cbuffer.data(ph_geom_11_off + 298 * ccomps * dcomps);

            auto g_y_y_z_xxxzz = cbuffer.data(ph_geom_11_off + 299 * ccomps * dcomps);

            auto g_y_y_z_xxyyy = cbuffer.data(ph_geom_11_off + 300 * ccomps * dcomps);

            auto g_y_y_z_xxyyz = cbuffer.data(ph_geom_11_off + 301 * ccomps * dcomps);

            auto g_y_y_z_xxyzz = cbuffer.data(ph_geom_11_off + 302 * ccomps * dcomps);

            auto g_y_y_z_xxzzz = cbuffer.data(ph_geom_11_off + 303 * ccomps * dcomps);

            auto g_y_y_z_xyyyy = cbuffer.data(ph_geom_11_off + 304 * ccomps * dcomps);

            auto g_y_y_z_xyyyz = cbuffer.data(ph_geom_11_off + 305 * ccomps * dcomps);

            auto g_y_y_z_xyyzz = cbuffer.data(ph_geom_11_off + 306 * ccomps * dcomps);

            auto g_y_y_z_xyzzz = cbuffer.data(ph_geom_11_off + 307 * ccomps * dcomps);

            auto g_y_y_z_xzzzz = cbuffer.data(ph_geom_11_off + 308 * ccomps * dcomps);

            auto g_y_y_z_yyyyy = cbuffer.data(ph_geom_11_off + 309 * ccomps * dcomps);

            auto g_y_y_z_yyyyz = cbuffer.data(ph_geom_11_off + 310 * ccomps * dcomps);

            auto g_y_y_z_yyyzz = cbuffer.data(ph_geom_11_off + 311 * ccomps * dcomps);

            auto g_y_y_z_yyzzz = cbuffer.data(ph_geom_11_off + 312 * ccomps * dcomps);

            auto g_y_y_z_yzzzz = cbuffer.data(ph_geom_11_off + 313 * ccomps * dcomps);

            auto g_y_y_z_zzzzz = cbuffer.data(ph_geom_11_off + 314 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_0_xxxxx, g_y_y_0_xxxxxz, g_y_y_0_xxxxy, g_y_y_0_xxxxyz, g_y_y_0_xxxxz, g_y_y_0_xxxxzz, g_y_y_0_xxxyy, g_y_y_0_xxxyyz, g_y_y_0_xxxyz, g_y_y_0_xxxyzz, g_y_y_0_xxxzz, g_y_y_0_xxxzzz, g_y_y_0_xxyyy, g_y_y_0_xxyyyz, g_y_y_0_xxyyz, g_y_y_0_xxyyzz, g_y_y_0_xxyzz, g_y_y_0_xxyzzz, g_y_y_0_xxzzz, g_y_y_0_xxzzzz, g_y_y_0_xyyyy, g_y_y_0_xyyyyz, g_y_y_0_xyyyz, g_y_y_0_xyyyzz, g_y_y_0_xyyzz, g_y_y_0_xyyzzz, g_y_y_0_xyzzz, g_y_y_0_xyzzzz, g_y_y_0_xzzzz, g_y_y_0_xzzzzz, g_y_y_0_yyyyy, g_y_y_0_yyyyyz, g_y_y_0_yyyyz, g_y_y_0_yyyyzz, g_y_y_0_yyyzz, g_y_y_0_yyyzzz, g_y_y_0_yyzzz, g_y_y_0_yyzzzz, g_y_y_0_yzzzz, g_y_y_0_yzzzzz, g_y_y_0_zzzzz, g_y_y_0_zzzzzz, g_y_y_z_xxxxx, g_y_y_z_xxxxy, g_y_y_z_xxxxz, g_y_y_z_xxxyy, g_y_y_z_xxxyz, g_y_y_z_xxxzz, g_y_y_z_xxyyy, g_y_y_z_xxyyz, g_y_y_z_xxyzz, g_y_y_z_xxzzz, g_y_y_z_xyyyy, g_y_y_z_xyyyz, g_y_y_z_xyyzz, g_y_y_z_xyzzz, g_y_y_z_xzzzz, g_y_y_z_yyyyy, g_y_y_z_yyyyz, g_y_y_z_yyyzz, g_y_y_z_yyzzz, g_y_y_z_yzzzz, g_y_y_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_z_xxxxx[k] = -g_y_y_0_xxxxx[k] * ab_z + g_y_y_0_xxxxxz[k];

                g_y_y_z_xxxxy[k] = -g_y_y_0_xxxxy[k] * ab_z + g_y_y_0_xxxxyz[k];

                g_y_y_z_xxxxz[k] = -g_y_y_0_xxxxz[k] * ab_z + g_y_y_0_xxxxzz[k];

                g_y_y_z_xxxyy[k] = -g_y_y_0_xxxyy[k] * ab_z + g_y_y_0_xxxyyz[k];

                g_y_y_z_xxxyz[k] = -g_y_y_0_xxxyz[k] * ab_z + g_y_y_0_xxxyzz[k];

                g_y_y_z_xxxzz[k] = -g_y_y_0_xxxzz[k] * ab_z + g_y_y_0_xxxzzz[k];

                g_y_y_z_xxyyy[k] = -g_y_y_0_xxyyy[k] * ab_z + g_y_y_0_xxyyyz[k];

                g_y_y_z_xxyyz[k] = -g_y_y_0_xxyyz[k] * ab_z + g_y_y_0_xxyyzz[k];

                g_y_y_z_xxyzz[k] = -g_y_y_0_xxyzz[k] * ab_z + g_y_y_0_xxyzzz[k];

                g_y_y_z_xxzzz[k] = -g_y_y_0_xxzzz[k] * ab_z + g_y_y_0_xxzzzz[k];

                g_y_y_z_xyyyy[k] = -g_y_y_0_xyyyy[k] * ab_z + g_y_y_0_xyyyyz[k];

                g_y_y_z_xyyyz[k] = -g_y_y_0_xyyyz[k] * ab_z + g_y_y_0_xyyyzz[k];

                g_y_y_z_xyyzz[k] = -g_y_y_0_xyyzz[k] * ab_z + g_y_y_0_xyyzzz[k];

                g_y_y_z_xyzzz[k] = -g_y_y_0_xyzzz[k] * ab_z + g_y_y_0_xyzzzz[k];

                g_y_y_z_xzzzz[k] = -g_y_y_0_xzzzz[k] * ab_z + g_y_y_0_xzzzzz[k];

                g_y_y_z_yyyyy[k] = -g_y_y_0_yyyyy[k] * ab_z + g_y_y_0_yyyyyz[k];

                g_y_y_z_yyyyz[k] = -g_y_y_0_yyyyz[k] * ab_z + g_y_y_0_yyyyzz[k];

                g_y_y_z_yyyzz[k] = -g_y_y_0_yyyzz[k] * ab_z + g_y_y_0_yyyzzz[k];

                g_y_y_z_yyzzz[k] = -g_y_y_0_yyzzz[k] * ab_z + g_y_y_0_yyzzzz[k];

                g_y_y_z_yzzzz[k] = -g_y_y_0_yzzzz[k] * ab_z + g_y_y_0_yzzzzz[k];

                g_y_y_z_zzzzz[k] = -g_y_y_0_zzzzz[k] * ab_z + g_y_y_0_zzzzzz[k];
            }

            /// Set up 315-336 components of targeted buffer : cbuffer.data(

            auto g_y_z_x_xxxxx = cbuffer.data(ph_geom_11_off + 315 * ccomps * dcomps);

            auto g_y_z_x_xxxxy = cbuffer.data(ph_geom_11_off + 316 * ccomps * dcomps);

            auto g_y_z_x_xxxxz = cbuffer.data(ph_geom_11_off + 317 * ccomps * dcomps);

            auto g_y_z_x_xxxyy = cbuffer.data(ph_geom_11_off + 318 * ccomps * dcomps);

            auto g_y_z_x_xxxyz = cbuffer.data(ph_geom_11_off + 319 * ccomps * dcomps);

            auto g_y_z_x_xxxzz = cbuffer.data(ph_geom_11_off + 320 * ccomps * dcomps);

            auto g_y_z_x_xxyyy = cbuffer.data(ph_geom_11_off + 321 * ccomps * dcomps);

            auto g_y_z_x_xxyyz = cbuffer.data(ph_geom_11_off + 322 * ccomps * dcomps);

            auto g_y_z_x_xxyzz = cbuffer.data(ph_geom_11_off + 323 * ccomps * dcomps);

            auto g_y_z_x_xxzzz = cbuffer.data(ph_geom_11_off + 324 * ccomps * dcomps);

            auto g_y_z_x_xyyyy = cbuffer.data(ph_geom_11_off + 325 * ccomps * dcomps);

            auto g_y_z_x_xyyyz = cbuffer.data(ph_geom_11_off + 326 * ccomps * dcomps);

            auto g_y_z_x_xyyzz = cbuffer.data(ph_geom_11_off + 327 * ccomps * dcomps);

            auto g_y_z_x_xyzzz = cbuffer.data(ph_geom_11_off + 328 * ccomps * dcomps);

            auto g_y_z_x_xzzzz = cbuffer.data(ph_geom_11_off + 329 * ccomps * dcomps);

            auto g_y_z_x_yyyyy = cbuffer.data(ph_geom_11_off + 330 * ccomps * dcomps);

            auto g_y_z_x_yyyyz = cbuffer.data(ph_geom_11_off + 331 * ccomps * dcomps);

            auto g_y_z_x_yyyzz = cbuffer.data(ph_geom_11_off + 332 * ccomps * dcomps);

            auto g_y_z_x_yyzzz = cbuffer.data(ph_geom_11_off + 333 * ccomps * dcomps);

            auto g_y_z_x_yzzzz = cbuffer.data(ph_geom_11_off + 334 * ccomps * dcomps);

            auto g_y_z_x_zzzzz = cbuffer.data(ph_geom_11_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_0_xxxxx, g_y_z_0_xxxxxx, g_y_z_0_xxxxxy, g_y_z_0_xxxxxz, g_y_z_0_xxxxy, g_y_z_0_xxxxyy, g_y_z_0_xxxxyz, g_y_z_0_xxxxz, g_y_z_0_xxxxzz, g_y_z_0_xxxyy, g_y_z_0_xxxyyy, g_y_z_0_xxxyyz, g_y_z_0_xxxyz, g_y_z_0_xxxyzz, g_y_z_0_xxxzz, g_y_z_0_xxxzzz, g_y_z_0_xxyyy, g_y_z_0_xxyyyy, g_y_z_0_xxyyyz, g_y_z_0_xxyyz, g_y_z_0_xxyyzz, g_y_z_0_xxyzz, g_y_z_0_xxyzzz, g_y_z_0_xxzzz, g_y_z_0_xxzzzz, g_y_z_0_xyyyy, g_y_z_0_xyyyyy, g_y_z_0_xyyyyz, g_y_z_0_xyyyz, g_y_z_0_xyyyzz, g_y_z_0_xyyzz, g_y_z_0_xyyzzz, g_y_z_0_xyzzz, g_y_z_0_xyzzzz, g_y_z_0_xzzzz, g_y_z_0_xzzzzz, g_y_z_0_yyyyy, g_y_z_0_yyyyz, g_y_z_0_yyyzz, g_y_z_0_yyzzz, g_y_z_0_yzzzz, g_y_z_0_zzzzz, g_y_z_x_xxxxx, g_y_z_x_xxxxy, g_y_z_x_xxxxz, g_y_z_x_xxxyy, g_y_z_x_xxxyz, g_y_z_x_xxxzz, g_y_z_x_xxyyy, g_y_z_x_xxyyz, g_y_z_x_xxyzz, g_y_z_x_xxzzz, g_y_z_x_xyyyy, g_y_z_x_xyyyz, g_y_z_x_xyyzz, g_y_z_x_xyzzz, g_y_z_x_xzzzz, g_y_z_x_yyyyy, g_y_z_x_yyyyz, g_y_z_x_yyyzz, g_y_z_x_yyzzz, g_y_z_x_yzzzz, g_y_z_x_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_x_xxxxx[k] = -g_y_z_0_xxxxx[k] * ab_x + g_y_z_0_xxxxxx[k];

                g_y_z_x_xxxxy[k] = -g_y_z_0_xxxxy[k] * ab_x + g_y_z_0_xxxxxy[k];

                g_y_z_x_xxxxz[k] = -g_y_z_0_xxxxz[k] * ab_x + g_y_z_0_xxxxxz[k];

                g_y_z_x_xxxyy[k] = -g_y_z_0_xxxyy[k] * ab_x + g_y_z_0_xxxxyy[k];

                g_y_z_x_xxxyz[k] = -g_y_z_0_xxxyz[k] * ab_x + g_y_z_0_xxxxyz[k];

                g_y_z_x_xxxzz[k] = -g_y_z_0_xxxzz[k] * ab_x + g_y_z_0_xxxxzz[k];

                g_y_z_x_xxyyy[k] = -g_y_z_0_xxyyy[k] * ab_x + g_y_z_0_xxxyyy[k];

                g_y_z_x_xxyyz[k] = -g_y_z_0_xxyyz[k] * ab_x + g_y_z_0_xxxyyz[k];

                g_y_z_x_xxyzz[k] = -g_y_z_0_xxyzz[k] * ab_x + g_y_z_0_xxxyzz[k];

                g_y_z_x_xxzzz[k] = -g_y_z_0_xxzzz[k] * ab_x + g_y_z_0_xxxzzz[k];

                g_y_z_x_xyyyy[k] = -g_y_z_0_xyyyy[k] * ab_x + g_y_z_0_xxyyyy[k];

                g_y_z_x_xyyyz[k] = -g_y_z_0_xyyyz[k] * ab_x + g_y_z_0_xxyyyz[k];

                g_y_z_x_xyyzz[k] = -g_y_z_0_xyyzz[k] * ab_x + g_y_z_0_xxyyzz[k];

                g_y_z_x_xyzzz[k] = -g_y_z_0_xyzzz[k] * ab_x + g_y_z_0_xxyzzz[k];

                g_y_z_x_xzzzz[k] = -g_y_z_0_xzzzz[k] * ab_x + g_y_z_0_xxzzzz[k];

                g_y_z_x_yyyyy[k] = -g_y_z_0_yyyyy[k] * ab_x + g_y_z_0_xyyyyy[k];

                g_y_z_x_yyyyz[k] = -g_y_z_0_yyyyz[k] * ab_x + g_y_z_0_xyyyyz[k];

                g_y_z_x_yyyzz[k] = -g_y_z_0_yyyzz[k] * ab_x + g_y_z_0_xyyyzz[k];

                g_y_z_x_yyzzz[k] = -g_y_z_0_yyzzz[k] * ab_x + g_y_z_0_xyyzzz[k];

                g_y_z_x_yzzzz[k] = -g_y_z_0_yzzzz[k] * ab_x + g_y_z_0_xyzzzz[k];

                g_y_z_x_zzzzz[k] = -g_y_z_0_zzzzz[k] * ab_x + g_y_z_0_xzzzzz[k];
            }

            /// Set up 336-357 components of targeted buffer : cbuffer.data(

            auto g_y_z_y_xxxxx = cbuffer.data(ph_geom_11_off + 336 * ccomps * dcomps);

            auto g_y_z_y_xxxxy = cbuffer.data(ph_geom_11_off + 337 * ccomps * dcomps);

            auto g_y_z_y_xxxxz = cbuffer.data(ph_geom_11_off + 338 * ccomps * dcomps);

            auto g_y_z_y_xxxyy = cbuffer.data(ph_geom_11_off + 339 * ccomps * dcomps);

            auto g_y_z_y_xxxyz = cbuffer.data(ph_geom_11_off + 340 * ccomps * dcomps);

            auto g_y_z_y_xxxzz = cbuffer.data(ph_geom_11_off + 341 * ccomps * dcomps);

            auto g_y_z_y_xxyyy = cbuffer.data(ph_geom_11_off + 342 * ccomps * dcomps);

            auto g_y_z_y_xxyyz = cbuffer.data(ph_geom_11_off + 343 * ccomps * dcomps);

            auto g_y_z_y_xxyzz = cbuffer.data(ph_geom_11_off + 344 * ccomps * dcomps);

            auto g_y_z_y_xxzzz = cbuffer.data(ph_geom_11_off + 345 * ccomps * dcomps);

            auto g_y_z_y_xyyyy = cbuffer.data(ph_geom_11_off + 346 * ccomps * dcomps);

            auto g_y_z_y_xyyyz = cbuffer.data(ph_geom_11_off + 347 * ccomps * dcomps);

            auto g_y_z_y_xyyzz = cbuffer.data(ph_geom_11_off + 348 * ccomps * dcomps);

            auto g_y_z_y_xyzzz = cbuffer.data(ph_geom_11_off + 349 * ccomps * dcomps);

            auto g_y_z_y_xzzzz = cbuffer.data(ph_geom_11_off + 350 * ccomps * dcomps);

            auto g_y_z_y_yyyyy = cbuffer.data(ph_geom_11_off + 351 * ccomps * dcomps);

            auto g_y_z_y_yyyyz = cbuffer.data(ph_geom_11_off + 352 * ccomps * dcomps);

            auto g_y_z_y_yyyzz = cbuffer.data(ph_geom_11_off + 353 * ccomps * dcomps);

            auto g_y_z_y_yyzzz = cbuffer.data(ph_geom_11_off + 354 * ccomps * dcomps);

            auto g_y_z_y_yzzzz = cbuffer.data(ph_geom_11_off + 355 * ccomps * dcomps);

            auto g_y_z_y_zzzzz = cbuffer.data(ph_geom_11_off + 356 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_0_xxxxx, g_0_z_0_xxxxy, g_0_z_0_xxxxz, g_0_z_0_xxxyy, g_0_z_0_xxxyz, g_0_z_0_xxxzz, g_0_z_0_xxyyy, g_0_z_0_xxyyz, g_0_z_0_xxyzz, g_0_z_0_xxzzz, g_0_z_0_xyyyy, g_0_z_0_xyyyz, g_0_z_0_xyyzz, g_0_z_0_xyzzz, g_0_z_0_xzzzz, g_0_z_0_yyyyy, g_0_z_0_yyyyz, g_0_z_0_yyyzz, g_0_z_0_yyzzz, g_0_z_0_yzzzz, g_0_z_0_zzzzz, g_y_z_0_xxxxx, g_y_z_0_xxxxxy, g_y_z_0_xxxxy, g_y_z_0_xxxxyy, g_y_z_0_xxxxyz, g_y_z_0_xxxxz, g_y_z_0_xxxyy, g_y_z_0_xxxyyy, g_y_z_0_xxxyyz, g_y_z_0_xxxyz, g_y_z_0_xxxyzz, g_y_z_0_xxxzz, g_y_z_0_xxyyy, g_y_z_0_xxyyyy, g_y_z_0_xxyyyz, g_y_z_0_xxyyz, g_y_z_0_xxyyzz, g_y_z_0_xxyzz, g_y_z_0_xxyzzz, g_y_z_0_xxzzz, g_y_z_0_xyyyy, g_y_z_0_xyyyyy, g_y_z_0_xyyyyz, g_y_z_0_xyyyz, g_y_z_0_xyyyzz, g_y_z_0_xyyzz, g_y_z_0_xyyzzz, g_y_z_0_xyzzz, g_y_z_0_xyzzzz, g_y_z_0_xzzzz, g_y_z_0_yyyyy, g_y_z_0_yyyyyy, g_y_z_0_yyyyyz, g_y_z_0_yyyyz, g_y_z_0_yyyyzz, g_y_z_0_yyyzz, g_y_z_0_yyyzzz, g_y_z_0_yyzzz, g_y_z_0_yyzzzz, g_y_z_0_yzzzz, g_y_z_0_yzzzzz, g_y_z_0_zzzzz, g_y_z_y_xxxxx, g_y_z_y_xxxxy, g_y_z_y_xxxxz, g_y_z_y_xxxyy, g_y_z_y_xxxyz, g_y_z_y_xxxzz, g_y_z_y_xxyyy, g_y_z_y_xxyyz, g_y_z_y_xxyzz, g_y_z_y_xxzzz, g_y_z_y_xyyyy, g_y_z_y_xyyyz, g_y_z_y_xyyzz, g_y_z_y_xyzzz, g_y_z_y_xzzzz, g_y_z_y_yyyyy, g_y_z_y_yyyyz, g_y_z_y_yyyzz, g_y_z_y_yyzzz, g_y_z_y_yzzzz, g_y_z_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_y_xxxxx[k] = -g_0_z_0_xxxxx[k] - g_y_z_0_xxxxx[k] * ab_y + g_y_z_0_xxxxxy[k];

                g_y_z_y_xxxxy[k] = -g_0_z_0_xxxxy[k] - g_y_z_0_xxxxy[k] * ab_y + g_y_z_0_xxxxyy[k];

                g_y_z_y_xxxxz[k] = -g_0_z_0_xxxxz[k] - g_y_z_0_xxxxz[k] * ab_y + g_y_z_0_xxxxyz[k];

                g_y_z_y_xxxyy[k] = -g_0_z_0_xxxyy[k] - g_y_z_0_xxxyy[k] * ab_y + g_y_z_0_xxxyyy[k];

                g_y_z_y_xxxyz[k] = -g_0_z_0_xxxyz[k] - g_y_z_0_xxxyz[k] * ab_y + g_y_z_0_xxxyyz[k];

                g_y_z_y_xxxzz[k] = -g_0_z_0_xxxzz[k] - g_y_z_0_xxxzz[k] * ab_y + g_y_z_0_xxxyzz[k];

                g_y_z_y_xxyyy[k] = -g_0_z_0_xxyyy[k] - g_y_z_0_xxyyy[k] * ab_y + g_y_z_0_xxyyyy[k];

                g_y_z_y_xxyyz[k] = -g_0_z_0_xxyyz[k] - g_y_z_0_xxyyz[k] * ab_y + g_y_z_0_xxyyyz[k];

                g_y_z_y_xxyzz[k] = -g_0_z_0_xxyzz[k] - g_y_z_0_xxyzz[k] * ab_y + g_y_z_0_xxyyzz[k];

                g_y_z_y_xxzzz[k] = -g_0_z_0_xxzzz[k] - g_y_z_0_xxzzz[k] * ab_y + g_y_z_0_xxyzzz[k];

                g_y_z_y_xyyyy[k] = -g_0_z_0_xyyyy[k] - g_y_z_0_xyyyy[k] * ab_y + g_y_z_0_xyyyyy[k];

                g_y_z_y_xyyyz[k] = -g_0_z_0_xyyyz[k] - g_y_z_0_xyyyz[k] * ab_y + g_y_z_0_xyyyyz[k];

                g_y_z_y_xyyzz[k] = -g_0_z_0_xyyzz[k] - g_y_z_0_xyyzz[k] * ab_y + g_y_z_0_xyyyzz[k];

                g_y_z_y_xyzzz[k] = -g_0_z_0_xyzzz[k] - g_y_z_0_xyzzz[k] * ab_y + g_y_z_0_xyyzzz[k];

                g_y_z_y_xzzzz[k] = -g_0_z_0_xzzzz[k] - g_y_z_0_xzzzz[k] * ab_y + g_y_z_0_xyzzzz[k];

                g_y_z_y_yyyyy[k] = -g_0_z_0_yyyyy[k] - g_y_z_0_yyyyy[k] * ab_y + g_y_z_0_yyyyyy[k];

                g_y_z_y_yyyyz[k] = -g_0_z_0_yyyyz[k] - g_y_z_0_yyyyz[k] * ab_y + g_y_z_0_yyyyyz[k];

                g_y_z_y_yyyzz[k] = -g_0_z_0_yyyzz[k] - g_y_z_0_yyyzz[k] * ab_y + g_y_z_0_yyyyzz[k];

                g_y_z_y_yyzzz[k] = -g_0_z_0_yyzzz[k] - g_y_z_0_yyzzz[k] * ab_y + g_y_z_0_yyyzzz[k];

                g_y_z_y_yzzzz[k] = -g_0_z_0_yzzzz[k] - g_y_z_0_yzzzz[k] * ab_y + g_y_z_0_yyzzzz[k];

                g_y_z_y_zzzzz[k] = -g_0_z_0_zzzzz[k] - g_y_z_0_zzzzz[k] * ab_y + g_y_z_0_yzzzzz[k];
            }

            /// Set up 357-378 components of targeted buffer : cbuffer.data(

            auto g_y_z_z_xxxxx = cbuffer.data(ph_geom_11_off + 357 * ccomps * dcomps);

            auto g_y_z_z_xxxxy = cbuffer.data(ph_geom_11_off + 358 * ccomps * dcomps);

            auto g_y_z_z_xxxxz = cbuffer.data(ph_geom_11_off + 359 * ccomps * dcomps);

            auto g_y_z_z_xxxyy = cbuffer.data(ph_geom_11_off + 360 * ccomps * dcomps);

            auto g_y_z_z_xxxyz = cbuffer.data(ph_geom_11_off + 361 * ccomps * dcomps);

            auto g_y_z_z_xxxzz = cbuffer.data(ph_geom_11_off + 362 * ccomps * dcomps);

            auto g_y_z_z_xxyyy = cbuffer.data(ph_geom_11_off + 363 * ccomps * dcomps);

            auto g_y_z_z_xxyyz = cbuffer.data(ph_geom_11_off + 364 * ccomps * dcomps);

            auto g_y_z_z_xxyzz = cbuffer.data(ph_geom_11_off + 365 * ccomps * dcomps);

            auto g_y_z_z_xxzzz = cbuffer.data(ph_geom_11_off + 366 * ccomps * dcomps);

            auto g_y_z_z_xyyyy = cbuffer.data(ph_geom_11_off + 367 * ccomps * dcomps);

            auto g_y_z_z_xyyyz = cbuffer.data(ph_geom_11_off + 368 * ccomps * dcomps);

            auto g_y_z_z_xyyzz = cbuffer.data(ph_geom_11_off + 369 * ccomps * dcomps);

            auto g_y_z_z_xyzzz = cbuffer.data(ph_geom_11_off + 370 * ccomps * dcomps);

            auto g_y_z_z_xzzzz = cbuffer.data(ph_geom_11_off + 371 * ccomps * dcomps);

            auto g_y_z_z_yyyyy = cbuffer.data(ph_geom_11_off + 372 * ccomps * dcomps);

            auto g_y_z_z_yyyyz = cbuffer.data(ph_geom_11_off + 373 * ccomps * dcomps);

            auto g_y_z_z_yyyzz = cbuffer.data(ph_geom_11_off + 374 * ccomps * dcomps);

            auto g_y_z_z_yyzzz = cbuffer.data(ph_geom_11_off + 375 * ccomps * dcomps);

            auto g_y_z_z_yzzzz = cbuffer.data(ph_geom_11_off + 376 * ccomps * dcomps);

            auto g_y_z_z_zzzzz = cbuffer.data(ph_geom_11_off + 377 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_0_xxxxx, g_y_0_0_xxxxy, g_y_0_0_xxxxz, g_y_0_0_xxxyy, g_y_0_0_xxxyz, g_y_0_0_xxxzz, g_y_0_0_xxyyy, g_y_0_0_xxyyz, g_y_0_0_xxyzz, g_y_0_0_xxzzz, g_y_0_0_xyyyy, g_y_0_0_xyyyz, g_y_0_0_xyyzz, g_y_0_0_xyzzz, g_y_0_0_xzzzz, g_y_0_0_yyyyy, g_y_0_0_yyyyz, g_y_0_0_yyyzz, g_y_0_0_yyzzz, g_y_0_0_yzzzz, g_y_0_0_zzzzz, g_y_z_0_xxxxx, g_y_z_0_xxxxxz, g_y_z_0_xxxxy, g_y_z_0_xxxxyz, g_y_z_0_xxxxz, g_y_z_0_xxxxzz, g_y_z_0_xxxyy, g_y_z_0_xxxyyz, g_y_z_0_xxxyz, g_y_z_0_xxxyzz, g_y_z_0_xxxzz, g_y_z_0_xxxzzz, g_y_z_0_xxyyy, g_y_z_0_xxyyyz, g_y_z_0_xxyyz, g_y_z_0_xxyyzz, g_y_z_0_xxyzz, g_y_z_0_xxyzzz, g_y_z_0_xxzzz, g_y_z_0_xxzzzz, g_y_z_0_xyyyy, g_y_z_0_xyyyyz, g_y_z_0_xyyyz, g_y_z_0_xyyyzz, g_y_z_0_xyyzz, g_y_z_0_xyyzzz, g_y_z_0_xyzzz, g_y_z_0_xyzzzz, g_y_z_0_xzzzz, g_y_z_0_xzzzzz, g_y_z_0_yyyyy, g_y_z_0_yyyyyz, g_y_z_0_yyyyz, g_y_z_0_yyyyzz, g_y_z_0_yyyzz, g_y_z_0_yyyzzz, g_y_z_0_yyzzz, g_y_z_0_yyzzzz, g_y_z_0_yzzzz, g_y_z_0_yzzzzz, g_y_z_0_zzzzz, g_y_z_0_zzzzzz, g_y_z_z_xxxxx, g_y_z_z_xxxxy, g_y_z_z_xxxxz, g_y_z_z_xxxyy, g_y_z_z_xxxyz, g_y_z_z_xxxzz, g_y_z_z_xxyyy, g_y_z_z_xxyyz, g_y_z_z_xxyzz, g_y_z_z_xxzzz, g_y_z_z_xyyyy, g_y_z_z_xyyyz, g_y_z_z_xyyzz, g_y_z_z_xyzzz, g_y_z_z_xzzzz, g_y_z_z_yyyyy, g_y_z_z_yyyyz, g_y_z_z_yyyzz, g_y_z_z_yyzzz, g_y_z_z_yzzzz, g_y_z_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_z_xxxxx[k] = g_y_0_0_xxxxx[k] - g_y_z_0_xxxxx[k] * ab_z + g_y_z_0_xxxxxz[k];

                g_y_z_z_xxxxy[k] = g_y_0_0_xxxxy[k] - g_y_z_0_xxxxy[k] * ab_z + g_y_z_0_xxxxyz[k];

                g_y_z_z_xxxxz[k] = g_y_0_0_xxxxz[k] - g_y_z_0_xxxxz[k] * ab_z + g_y_z_0_xxxxzz[k];

                g_y_z_z_xxxyy[k] = g_y_0_0_xxxyy[k] - g_y_z_0_xxxyy[k] * ab_z + g_y_z_0_xxxyyz[k];

                g_y_z_z_xxxyz[k] = g_y_0_0_xxxyz[k] - g_y_z_0_xxxyz[k] * ab_z + g_y_z_0_xxxyzz[k];

                g_y_z_z_xxxzz[k] = g_y_0_0_xxxzz[k] - g_y_z_0_xxxzz[k] * ab_z + g_y_z_0_xxxzzz[k];

                g_y_z_z_xxyyy[k] = g_y_0_0_xxyyy[k] - g_y_z_0_xxyyy[k] * ab_z + g_y_z_0_xxyyyz[k];

                g_y_z_z_xxyyz[k] = g_y_0_0_xxyyz[k] - g_y_z_0_xxyyz[k] * ab_z + g_y_z_0_xxyyzz[k];

                g_y_z_z_xxyzz[k] = g_y_0_0_xxyzz[k] - g_y_z_0_xxyzz[k] * ab_z + g_y_z_0_xxyzzz[k];

                g_y_z_z_xxzzz[k] = g_y_0_0_xxzzz[k] - g_y_z_0_xxzzz[k] * ab_z + g_y_z_0_xxzzzz[k];

                g_y_z_z_xyyyy[k] = g_y_0_0_xyyyy[k] - g_y_z_0_xyyyy[k] * ab_z + g_y_z_0_xyyyyz[k];

                g_y_z_z_xyyyz[k] = g_y_0_0_xyyyz[k] - g_y_z_0_xyyyz[k] * ab_z + g_y_z_0_xyyyzz[k];

                g_y_z_z_xyyzz[k] = g_y_0_0_xyyzz[k] - g_y_z_0_xyyzz[k] * ab_z + g_y_z_0_xyyzzz[k];

                g_y_z_z_xyzzz[k] = g_y_0_0_xyzzz[k] - g_y_z_0_xyzzz[k] * ab_z + g_y_z_0_xyzzzz[k];

                g_y_z_z_xzzzz[k] = g_y_0_0_xzzzz[k] - g_y_z_0_xzzzz[k] * ab_z + g_y_z_0_xzzzzz[k];

                g_y_z_z_yyyyy[k] = g_y_0_0_yyyyy[k] - g_y_z_0_yyyyy[k] * ab_z + g_y_z_0_yyyyyz[k];

                g_y_z_z_yyyyz[k] = g_y_0_0_yyyyz[k] - g_y_z_0_yyyyz[k] * ab_z + g_y_z_0_yyyyzz[k];

                g_y_z_z_yyyzz[k] = g_y_0_0_yyyzz[k] - g_y_z_0_yyyzz[k] * ab_z + g_y_z_0_yyyzzz[k];

                g_y_z_z_yyzzz[k] = g_y_0_0_yyzzz[k] - g_y_z_0_yyzzz[k] * ab_z + g_y_z_0_yyzzzz[k];

                g_y_z_z_yzzzz[k] = g_y_0_0_yzzzz[k] - g_y_z_0_yzzzz[k] * ab_z + g_y_z_0_yzzzzz[k];

                g_y_z_z_zzzzz[k] = g_y_0_0_zzzzz[k] - g_y_z_0_zzzzz[k] * ab_z + g_y_z_0_zzzzzz[k];
            }

            /// Set up 378-399 components of targeted buffer : cbuffer.data(

            auto g_z_x_x_xxxxx = cbuffer.data(ph_geom_11_off + 378 * ccomps * dcomps);

            auto g_z_x_x_xxxxy = cbuffer.data(ph_geom_11_off + 379 * ccomps * dcomps);

            auto g_z_x_x_xxxxz = cbuffer.data(ph_geom_11_off + 380 * ccomps * dcomps);

            auto g_z_x_x_xxxyy = cbuffer.data(ph_geom_11_off + 381 * ccomps * dcomps);

            auto g_z_x_x_xxxyz = cbuffer.data(ph_geom_11_off + 382 * ccomps * dcomps);

            auto g_z_x_x_xxxzz = cbuffer.data(ph_geom_11_off + 383 * ccomps * dcomps);

            auto g_z_x_x_xxyyy = cbuffer.data(ph_geom_11_off + 384 * ccomps * dcomps);

            auto g_z_x_x_xxyyz = cbuffer.data(ph_geom_11_off + 385 * ccomps * dcomps);

            auto g_z_x_x_xxyzz = cbuffer.data(ph_geom_11_off + 386 * ccomps * dcomps);

            auto g_z_x_x_xxzzz = cbuffer.data(ph_geom_11_off + 387 * ccomps * dcomps);

            auto g_z_x_x_xyyyy = cbuffer.data(ph_geom_11_off + 388 * ccomps * dcomps);

            auto g_z_x_x_xyyyz = cbuffer.data(ph_geom_11_off + 389 * ccomps * dcomps);

            auto g_z_x_x_xyyzz = cbuffer.data(ph_geom_11_off + 390 * ccomps * dcomps);

            auto g_z_x_x_xyzzz = cbuffer.data(ph_geom_11_off + 391 * ccomps * dcomps);

            auto g_z_x_x_xzzzz = cbuffer.data(ph_geom_11_off + 392 * ccomps * dcomps);

            auto g_z_x_x_yyyyy = cbuffer.data(ph_geom_11_off + 393 * ccomps * dcomps);

            auto g_z_x_x_yyyyz = cbuffer.data(ph_geom_11_off + 394 * ccomps * dcomps);

            auto g_z_x_x_yyyzz = cbuffer.data(ph_geom_11_off + 395 * ccomps * dcomps);

            auto g_z_x_x_yyzzz = cbuffer.data(ph_geom_11_off + 396 * ccomps * dcomps);

            auto g_z_x_x_yzzzz = cbuffer.data(ph_geom_11_off + 397 * ccomps * dcomps);

            auto g_z_x_x_zzzzz = cbuffer.data(ph_geom_11_off + 398 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_0_xxxxx, g_z_0_0_xxxxy, g_z_0_0_xxxxz, g_z_0_0_xxxyy, g_z_0_0_xxxyz, g_z_0_0_xxxzz, g_z_0_0_xxyyy, g_z_0_0_xxyyz, g_z_0_0_xxyzz, g_z_0_0_xxzzz, g_z_0_0_xyyyy, g_z_0_0_xyyyz, g_z_0_0_xyyzz, g_z_0_0_xyzzz, g_z_0_0_xzzzz, g_z_0_0_yyyyy, g_z_0_0_yyyyz, g_z_0_0_yyyzz, g_z_0_0_yyzzz, g_z_0_0_yzzzz, g_z_0_0_zzzzz, g_z_x_0_xxxxx, g_z_x_0_xxxxxx, g_z_x_0_xxxxxy, g_z_x_0_xxxxxz, g_z_x_0_xxxxy, g_z_x_0_xxxxyy, g_z_x_0_xxxxyz, g_z_x_0_xxxxz, g_z_x_0_xxxxzz, g_z_x_0_xxxyy, g_z_x_0_xxxyyy, g_z_x_0_xxxyyz, g_z_x_0_xxxyz, g_z_x_0_xxxyzz, g_z_x_0_xxxzz, g_z_x_0_xxxzzz, g_z_x_0_xxyyy, g_z_x_0_xxyyyy, g_z_x_0_xxyyyz, g_z_x_0_xxyyz, g_z_x_0_xxyyzz, g_z_x_0_xxyzz, g_z_x_0_xxyzzz, g_z_x_0_xxzzz, g_z_x_0_xxzzzz, g_z_x_0_xyyyy, g_z_x_0_xyyyyy, g_z_x_0_xyyyyz, g_z_x_0_xyyyz, g_z_x_0_xyyyzz, g_z_x_0_xyyzz, g_z_x_0_xyyzzz, g_z_x_0_xyzzz, g_z_x_0_xyzzzz, g_z_x_0_xzzzz, g_z_x_0_xzzzzz, g_z_x_0_yyyyy, g_z_x_0_yyyyz, g_z_x_0_yyyzz, g_z_x_0_yyzzz, g_z_x_0_yzzzz, g_z_x_0_zzzzz, g_z_x_x_xxxxx, g_z_x_x_xxxxy, g_z_x_x_xxxxz, g_z_x_x_xxxyy, g_z_x_x_xxxyz, g_z_x_x_xxxzz, g_z_x_x_xxyyy, g_z_x_x_xxyyz, g_z_x_x_xxyzz, g_z_x_x_xxzzz, g_z_x_x_xyyyy, g_z_x_x_xyyyz, g_z_x_x_xyyzz, g_z_x_x_xyzzz, g_z_x_x_xzzzz, g_z_x_x_yyyyy, g_z_x_x_yyyyz, g_z_x_x_yyyzz, g_z_x_x_yyzzz, g_z_x_x_yzzzz, g_z_x_x_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_x_xxxxx[k] = g_z_0_0_xxxxx[k] - g_z_x_0_xxxxx[k] * ab_x + g_z_x_0_xxxxxx[k];

                g_z_x_x_xxxxy[k] = g_z_0_0_xxxxy[k] - g_z_x_0_xxxxy[k] * ab_x + g_z_x_0_xxxxxy[k];

                g_z_x_x_xxxxz[k] = g_z_0_0_xxxxz[k] - g_z_x_0_xxxxz[k] * ab_x + g_z_x_0_xxxxxz[k];

                g_z_x_x_xxxyy[k] = g_z_0_0_xxxyy[k] - g_z_x_0_xxxyy[k] * ab_x + g_z_x_0_xxxxyy[k];

                g_z_x_x_xxxyz[k] = g_z_0_0_xxxyz[k] - g_z_x_0_xxxyz[k] * ab_x + g_z_x_0_xxxxyz[k];

                g_z_x_x_xxxzz[k] = g_z_0_0_xxxzz[k] - g_z_x_0_xxxzz[k] * ab_x + g_z_x_0_xxxxzz[k];

                g_z_x_x_xxyyy[k] = g_z_0_0_xxyyy[k] - g_z_x_0_xxyyy[k] * ab_x + g_z_x_0_xxxyyy[k];

                g_z_x_x_xxyyz[k] = g_z_0_0_xxyyz[k] - g_z_x_0_xxyyz[k] * ab_x + g_z_x_0_xxxyyz[k];

                g_z_x_x_xxyzz[k] = g_z_0_0_xxyzz[k] - g_z_x_0_xxyzz[k] * ab_x + g_z_x_0_xxxyzz[k];

                g_z_x_x_xxzzz[k] = g_z_0_0_xxzzz[k] - g_z_x_0_xxzzz[k] * ab_x + g_z_x_0_xxxzzz[k];

                g_z_x_x_xyyyy[k] = g_z_0_0_xyyyy[k] - g_z_x_0_xyyyy[k] * ab_x + g_z_x_0_xxyyyy[k];

                g_z_x_x_xyyyz[k] = g_z_0_0_xyyyz[k] - g_z_x_0_xyyyz[k] * ab_x + g_z_x_0_xxyyyz[k];

                g_z_x_x_xyyzz[k] = g_z_0_0_xyyzz[k] - g_z_x_0_xyyzz[k] * ab_x + g_z_x_0_xxyyzz[k];

                g_z_x_x_xyzzz[k] = g_z_0_0_xyzzz[k] - g_z_x_0_xyzzz[k] * ab_x + g_z_x_0_xxyzzz[k];

                g_z_x_x_xzzzz[k] = g_z_0_0_xzzzz[k] - g_z_x_0_xzzzz[k] * ab_x + g_z_x_0_xxzzzz[k];

                g_z_x_x_yyyyy[k] = g_z_0_0_yyyyy[k] - g_z_x_0_yyyyy[k] * ab_x + g_z_x_0_xyyyyy[k];

                g_z_x_x_yyyyz[k] = g_z_0_0_yyyyz[k] - g_z_x_0_yyyyz[k] * ab_x + g_z_x_0_xyyyyz[k];

                g_z_x_x_yyyzz[k] = g_z_0_0_yyyzz[k] - g_z_x_0_yyyzz[k] * ab_x + g_z_x_0_xyyyzz[k];

                g_z_x_x_yyzzz[k] = g_z_0_0_yyzzz[k] - g_z_x_0_yyzzz[k] * ab_x + g_z_x_0_xyyzzz[k];

                g_z_x_x_yzzzz[k] = g_z_0_0_yzzzz[k] - g_z_x_0_yzzzz[k] * ab_x + g_z_x_0_xyzzzz[k];

                g_z_x_x_zzzzz[k] = g_z_0_0_zzzzz[k] - g_z_x_0_zzzzz[k] * ab_x + g_z_x_0_xzzzzz[k];
            }

            /// Set up 399-420 components of targeted buffer : cbuffer.data(

            auto g_z_x_y_xxxxx = cbuffer.data(ph_geom_11_off + 399 * ccomps * dcomps);

            auto g_z_x_y_xxxxy = cbuffer.data(ph_geom_11_off + 400 * ccomps * dcomps);

            auto g_z_x_y_xxxxz = cbuffer.data(ph_geom_11_off + 401 * ccomps * dcomps);

            auto g_z_x_y_xxxyy = cbuffer.data(ph_geom_11_off + 402 * ccomps * dcomps);

            auto g_z_x_y_xxxyz = cbuffer.data(ph_geom_11_off + 403 * ccomps * dcomps);

            auto g_z_x_y_xxxzz = cbuffer.data(ph_geom_11_off + 404 * ccomps * dcomps);

            auto g_z_x_y_xxyyy = cbuffer.data(ph_geom_11_off + 405 * ccomps * dcomps);

            auto g_z_x_y_xxyyz = cbuffer.data(ph_geom_11_off + 406 * ccomps * dcomps);

            auto g_z_x_y_xxyzz = cbuffer.data(ph_geom_11_off + 407 * ccomps * dcomps);

            auto g_z_x_y_xxzzz = cbuffer.data(ph_geom_11_off + 408 * ccomps * dcomps);

            auto g_z_x_y_xyyyy = cbuffer.data(ph_geom_11_off + 409 * ccomps * dcomps);

            auto g_z_x_y_xyyyz = cbuffer.data(ph_geom_11_off + 410 * ccomps * dcomps);

            auto g_z_x_y_xyyzz = cbuffer.data(ph_geom_11_off + 411 * ccomps * dcomps);

            auto g_z_x_y_xyzzz = cbuffer.data(ph_geom_11_off + 412 * ccomps * dcomps);

            auto g_z_x_y_xzzzz = cbuffer.data(ph_geom_11_off + 413 * ccomps * dcomps);

            auto g_z_x_y_yyyyy = cbuffer.data(ph_geom_11_off + 414 * ccomps * dcomps);

            auto g_z_x_y_yyyyz = cbuffer.data(ph_geom_11_off + 415 * ccomps * dcomps);

            auto g_z_x_y_yyyzz = cbuffer.data(ph_geom_11_off + 416 * ccomps * dcomps);

            auto g_z_x_y_yyzzz = cbuffer.data(ph_geom_11_off + 417 * ccomps * dcomps);

            auto g_z_x_y_yzzzz = cbuffer.data(ph_geom_11_off + 418 * ccomps * dcomps);

            auto g_z_x_y_zzzzz = cbuffer.data(ph_geom_11_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_0_xxxxx, g_z_x_0_xxxxxy, g_z_x_0_xxxxy, g_z_x_0_xxxxyy, g_z_x_0_xxxxyz, g_z_x_0_xxxxz, g_z_x_0_xxxyy, g_z_x_0_xxxyyy, g_z_x_0_xxxyyz, g_z_x_0_xxxyz, g_z_x_0_xxxyzz, g_z_x_0_xxxzz, g_z_x_0_xxyyy, g_z_x_0_xxyyyy, g_z_x_0_xxyyyz, g_z_x_0_xxyyz, g_z_x_0_xxyyzz, g_z_x_0_xxyzz, g_z_x_0_xxyzzz, g_z_x_0_xxzzz, g_z_x_0_xyyyy, g_z_x_0_xyyyyy, g_z_x_0_xyyyyz, g_z_x_0_xyyyz, g_z_x_0_xyyyzz, g_z_x_0_xyyzz, g_z_x_0_xyyzzz, g_z_x_0_xyzzz, g_z_x_0_xyzzzz, g_z_x_0_xzzzz, g_z_x_0_yyyyy, g_z_x_0_yyyyyy, g_z_x_0_yyyyyz, g_z_x_0_yyyyz, g_z_x_0_yyyyzz, g_z_x_0_yyyzz, g_z_x_0_yyyzzz, g_z_x_0_yyzzz, g_z_x_0_yyzzzz, g_z_x_0_yzzzz, g_z_x_0_yzzzzz, g_z_x_0_zzzzz, g_z_x_y_xxxxx, g_z_x_y_xxxxy, g_z_x_y_xxxxz, g_z_x_y_xxxyy, g_z_x_y_xxxyz, g_z_x_y_xxxzz, g_z_x_y_xxyyy, g_z_x_y_xxyyz, g_z_x_y_xxyzz, g_z_x_y_xxzzz, g_z_x_y_xyyyy, g_z_x_y_xyyyz, g_z_x_y_xyyzz, g_z_x_y_xyzzz, g_z_x_y_xzzzz, g_z_x_y_yyyyy, g_z_x_y_yyyyz, g_z_x_y_yyyzz, g_z_x_y_yyzzz, g_z_x_y_yzzzz, g_z_x_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_y_xxxxx[k] = -g_z_x_0_xxxxx[k] * ab_y + g_z_x_0_xxxxxy[k];

                g_z_x_y_xxxxy[k] = -g_z_x_0_xxxxy[k] * ab_y + g_z_x_0_xxxxyy[k];

                g_z_x_y_xxxxz[k] = -g_z_x_0_xxxxz[k] * ab_y + g_z_x_0_xxxxyz[k];

                g_z_x_y_xxxyy[k] = -g_z_x_0_xxxyy[k] * ab_y + g_z_x_0_xxxyyy[k];

                g_z_x_y_xxxyz[k] = -g_z_x_0_xxxyz[k] * ab_y + g_z_x_0_xxxyyz[k];

                g_z_x_y_xxxzz[k] = -g_z_x_0_xxxzz[k] * ab_y + g_z_x_0_xxxyzz[k];

                g_z_x_y_xxyyy[k] = -g_z_x_0_xxyyy[k] * ab_y + g_z_x_0_xxyyyy[k];

                g_z_x_y_xxyyz[k] = -g_z_x_0_xxyyz[k] * ab_y + g_z_x_0_xxyyyz[k];

                g_z_x_y_xxyzz[k] = -g_z_x_0_xxyzz[k] * ab_y + g_z_x_0_xxyyzz[k];

                g_z_x_y_xxzzz[k] = -g_z_x_0_xxzzz[k] * ab_y + g_z_x_0_xxyzzz[k];

                g_z_x_y_xyyyy[k] = -g_z_x_0_xyyyy[k] * ab_y + g_z_x_0_xyyyyy[k];

                g_z_x_y_xyyyz[k] = -g_z_x_0_xyyyz[k] * ab_y + g_z_x_0_xyyyyz[k];

                g_z_x_y_xyyzz[k] = -g_z_x_0_xyyzz[k] * ab_y + g_z_x_0_xyyyzz[k];

                g_z_x_y_xyzzz[k] = -g_z_x_0_xyzzz[k] * ab_y + g_z_x_0_xyyzzz[k];

                g_z_x_y_xzzzz[k] = -g_z_x_0_xzzzz[k] * ab_y + g_z_x_0_xyzzzz[k];

                g_z_x_y_yyyyy[k] = -g_z_x_0_yyyyy[k] * ab_y + g_z_x_0_yyyyyy[k];

                g_z_x_y_yyyyz[k] = -g_z_x_0_yyyyz[k] * ab_y + g_z_x_0_yyyyyz[k];

                g_z_x_y_yyyzz[k] = -g_z_x_0_yyyzz[k] * ab_y + g_z_x_0_yyyyzz[k];

                g_z_x_y_yyzzz[k] = -g_z_x_0_yyzzz[k] * ab_y + g_z_x_0_yyyzzz[k];

                g_z_x_y_yzzzz[k] = -g_z_x_0_yzzzz[k] * ab_y + g_z_x_0_yyzzzz[k];

                g_z_x_y_zzzzz[k] = -g_z_x_0_zzzzz[k] * ab_y + g_z_x_0_yzzzzz[k];
            }

            /// Set up 420-441 components of targeted buffer : cbuffer.data(

            auto g_z_x_z_xxxxx = cbuffer.data(ph_geom_11_off + 420 * ccomps * dcomps);

            auto g_z_x_z_xxxxy = cbuffer.data(ph_geom_11_off + 421 * ccomps * dcomps);

            auto g_z_x_z_xxxxz = cbuffer.data(ph_geom_11_off + 422 * ccomps * dcomps);

            auto g_z_x_z_xxxyy = cbuffer.data(ph_geom_11_off + 423 * ccomps * dcomps);

            auto g_z_x_z_xxxyz = cbuffer.data(ph_geom_11_off + 424 * ccomps * dcomps);

            auto g_z_x_z_xxxzz = cbuffer.data(ph_geom_11_off + 425 * ccomps * dcomps);

            auto g_z_x_z_xxyyy = cbuffer.data(ph_geom_11_off + 426 * ccomps * dcomps);

            auto g_z_x_z_xxyyz = cbuffer.data(ph_geom_11_off + 427 * ccomps * dcomps);

            auto g_z_x_z_xxyzz = cbuffer.data(ph_geom_11_off + 428 * ccomps * dcomps);

            auto g_z_x_z_xxzzz = cbuffer.data(ph_geom_11_off + 429 * ccomps * dcomps);

            auto g_z_x_z_xyyyy = cbuffer.data(ph_geom_11_off + 430 * ccomps * dcomps);

            auto g_z_x_z_xyyyz = cbuffer.data(ph_geom_11_off + 431 * ccomps * dcomps);

            auto g_z_x_z_xyyzz = cbuffer.data(ph_geom_11_off + 432 * ccomps * dcomps);

            auto g_z_x_z_xyzzz = cbuffer.data(ph_geom_11_off + 433 * ccomps * dcomps);

            auto g_z_x_z_xzzzz = cbuffer.data(ph_geom_11_off + 434 * ccomps * dcomps);

            auto g_z_x_z_yyyyy = cbuffer.data(ph_geom_11_off + 435 * ccomps * dcomps);

            auto g_z_x_z_yyyyz = cbuffer.data(ph_geom_11_off + 436 * ccomps * dcomps);

            auto g_z_x_z_yyyzz = cbuffer.data(ph_geom_11_off + 437 * ccomps * dcomps);

            auto g_z_x_z_yyzzz = cbuffer.data(ph_geom_11_off + 438 * ccomps * dcomps);

            auto g_z_x_z_yzzzz = cbuffer.data(ph_geom_11_off + 439 * ccomps * dcomps);

            auto g_z_x_z_zzzzz = cbuffer.data(ph_geom_11_off + 440 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xxxxx, g_0_x_0_xxxxy, g_0_x_0_xxxxz, g_0_x_0_xxxyy, g_0_x_0_xxxyz, g_0_x_0_xxxzz, g_0_x_0_xxyyy, g_0_x_0_xxyyz, g_0_x_0_xxyzz, g_0_x_0_xxzzz, g_0_x_0_xyyyy, g_0_x_0_xyyyz, g_0_x_0_xyyzz, g_0_x_0_xyzzz, g_0_x_0_xzzzz, g_0_x_0_yyyyy, g_0_x_0_yyyyz, g_0_x_0_yyyzz, g_0_x_0_yyzzz, g_0_x_0_yzzzz, g_0_x_0_zzzzz, g_z_x_0_xxxxx, g_z_x_0_xxxxxz, g_z_x_0_xxxxy, g_z_x_0_xxxxyz, g_z_x_0_xxxxz, g_z_x_0_xxxxzz, g_z_x_0_xxxyy, g_z_x_0_xxxyyz, g_z_x_0_xxxyz, g_z_x_0_xxxyzz, g_z_x_0_xxxzz, g_z_x_0_xxxzzz, g_z_x_0_xxyyy, g_z_x_0_xxyyyz, g_z_x_0_xxyyz, g_z_x_0_xxyyzz, g_z_x_0_xxyzz, g_z_x_0_xxyzzz, g_z_x_0_xxzzz, g_z_x_0_xxzzzz, g_z_x_0_xyyyy, g_z_x_0_xyyyyz, g_z_x_0_xyyyz, g_z_x_0_xyyyzz, g_z_x_0_xyyzz, g_z_x_0_xyyzzz, g_z_x_0_xyzzz, g_z_x_0_xyzzzz, g_z_x_0_xzzzz, g_z_x_0_xzzzzz, g_z_x_0_yyyyy, g_z_x_0_yyyyyz, g_z_x_0_yyyyz, g_z_x_0_yyyyzz, g_z_x_0_yyyzz, g_z_x_0_yyyzzz, g_z_x_0_yyzzz, g_z_x_0_yyzzzz, g_z_x_0_yzzzz, g_z_x_0_yzzzzz, g_z_x_0_zzzzz, g_z_x_0_zzzzzz, g_z_x_z_xxxxx, g_z_x_z_xxxxy, g_z_x_z_xxxxz, g_z_x_z_xxxyy, g_z_x_z_xxxyz, g_z_x_z_xxxzz, g_z_x_z_xxyyy, g_z_x_z_xxyyz, g_z_x_z_xxyzz, g_z_x_z_xxzzz, g_z_x_z_xyyyy, g_z_x_z_xyyyz, g_z_x_z_xyyzz, g_z_x_z_xyzzz, g_z_x_z_xzzzz, g_z_x_z_yyyyy, g_z_x_z_yyyyz, g_z_x_z_yyyzz, g_z_x_z_yyzzz, g_z_x_z_yzzzz, g_z_x_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_z_xxxxx[k] = -g_0_x_0_xxxxx[k] - g_z_x_0_xxxxx[k] * ab_z + g_z_x_0_xxxxxz[k];

                g_z_x_z_xxxxy[k] = -g_0_x_0_xxxxy[k] - g_z_x_0_xxxxy[k] * ab_z + g_z_x_0_xxxxyz[k];

                g_z_x_z_xxxxz[k] = -g_0_x_0_xxxxz[k] - g_z_x_0_xxxxz[k] * ab_z + g_z_x_0_xxxxzz[k];

                g_z_x_z_xxxyy[k] = -g_0_x_0_xxxyy[k] - g_z_x_0_xxxyy[k] * ab_z + g_z_x_0_xxxyyz[k];

                g_z_x_z_xxxyz[k] = -g_0_x_0_xxxyz[k] - g_z_x_0_xxxyz[k] * ab_z + g_z_x_0_xxxyzz[k];

                g_z_x_z_xxxzz[k] = -g_0_x_0_xxxzz[k] - g_z_x_0_xxxzz[k] * ab_z + g_z_x_0_xxxzzz[k];

                g_z_x_z_xxyyy[k] = -g_0_x_0_xxyyy[k] - g_z_x_0_xxyyy[k] * ab_z + g_z_x_0_xxyyyz[k];

                g_z_x_z_xxyyz[k] = -g_0_x_0_xxyyz[k] - g_z_x_0_xxyyz[k] * ab_z + g_z_x_0_xxyyzz[k];

                g_z_x_z_xxyzz[k] = -g_0_x_0_xxyzz[k] - g_z_x_0_xxyzz[k] * ab_z + g_z_x_0_xxyzzz[k];

                g_z_x_z_xxzzz[k] = -g_0_x_0_xxzzz[k] - g_z_x_0_xxzzz[k] * ab_z + g_z_x_0_xxzzzz[k];

                g_z_x_z_xyyyy[k] = -g_0_x_0_xyyyy[k] - g_z_x_0_xyyyy[k] * ab_z + g_z_x_0_xyyyyz[k];

                g_z_x_z_xyyyz[k] = -g_0_x_0_xyyyz[k] - g_z_x_0_xyyyz[k] * ab_z + g_z_x_0_xyyyzz[k];

                g_z_x_z_xyyzz[k] = -g_0_x_0_xyyzz[k] - g_z_x_0_xyyzz[k] * ab_z + g_z_x_0_xyyzzz[k];

                g_z_x_z_xyzzz[k] = -g_0_x_0_xyzzz[k] - g_z_x_0_xyzzz[k] * ab_z + g_z_x_0_xyzzzz[k];

                g_z_x_z_xzzzz[k] = -g_0_x_0_xzzzz[k] - g_z_x_0_xzzzz[k] * ab_z + g_z_x_0_xzzzzz[k];

                g_z_x_z_yyyyy[k] = -g_0_x_0_yyyyy[k] - g_z_x_0_yyyyy[k] * ab_z + g_z_x_0_yyyyyz[k];

                g_z_x_z_yyyyz[k] = -g_0_x_0_yyyyz[k] - g_z_x_0_yyyyz[k] * ab_z + g_z_x_0_yyyyzz[k];

                g_z_x_z_yyyzz[k] = -g_0_x_0_yyyzz[k] - g_z_x_0_yyyzz[k] * ab_z + g_z_x_0_yyyzzz[k];

                g_z_x_z_yyzzz[k] = -g_0_x_0_yyzzz[k] - g_z_x_0_yyzzz[k] * ab_z + g_z_x_0_yyzzzz[k];

                g_z_x_z_yzzzz[k] = -g_0_x_0_yzzzz[k] - g_z_x_0_yzzzz[k] * ab_z + g_z_x_0_yzzzzz[k];

                g_z_x_z_zzzzz[k] = -g_0_x_0_zzzzz[k] - g_z_x_0_zzzzz[k] * ab_z + g_z_x_0_zzzzzz[k];
            }

            /// Set up 441-462 components of targeted buffer : cbuffer.data(

            auto g_z_y_x_xxxxx = cbuffer.data(ph_geom_11_off + 441 * ccomps * dcomps);

            auto g_z_y_x_xxxxy = cbuffer.data(ph_geom_11_off + 442 * ccomps * dcomps);

            auto g_z_y_x_xxxxz = cbuffer.data(ph_geom_11_off + 443 * ccomps * dcomps);

            auto g_z_y_x_xxxyy = cbuffer.data(ph_geom_11_off + 444 * ccomps * dcomps);

            auto g_z_y_x_xxxyz = cbuffer.data(ph_geom_11_off + 445 * ccomps * dcomps);

            auto g_z_y_x_xxxzz = cbuffer.data(ph_geom_11_off + 446 * ccomps * dcomps);

            auto g_z_y_x_xxyyy = cbuffer.data(ph_geom_11_off + 447 * ccomps * dcomps);

            auto g_z_y_x_xxyyz = cbuffer.data(ph_geom_11_off + 448 * ccomps * dcomps);

            auto g_z_y_x_xxyzz = cbuffer.data(ph_geom_11_off + 449 * ccomps * dcomps);

            auto g_z_y_x_xxzzz = cbuffer.data(ph_geom_11_off + 450 * ccomps * dcomps);

            auto g_z_y_x_xyyyy = cbuffer.data(ph_geom_11_off + 451 * ccomps * dcomps);

            auto g_z_y_x_xyyyz = cbuffer.data(ph_geom_11_off + 452 * ccomps * dcomps);

            auto g_z_y_x_xyyzz = cbuffer.data(ph_geom_11_off + 453 * ccomps * dcomps);

            auto g_z_y_x_xyzzz = cbuffer.data(ph_geom_11_off + 454 * ccomps * dcomps);

            auto g_z_y_x_xzzzz = cbuffer.data(ph_geom_11_off + 455 * ccomps * dcomps);

            auto g_z_y_x_yyyyy = cbuffer.data(ph_geom_11_off + 456 * ccomps * dcomps);

            auto g_z_y_x_yyyyz = cbuffer.data(ph_geom_11_off + 457 * ccomps * dcomps);

            auto g_z_y_x_yyyzz = cbuffer.data(ph_geom_11_off + 458 * ccomps * dcomps);

            auto g_z_y_x_yyzzz = cbuffer.data(ph_geom_11_off + 459 * ccomps * dcomps);

            auto g_z_y_x_yzzzz = cbuffer.data(ph_geom_11_off + 460 * ccomps * dcomps);

            auto g_z_y_x_zzzzz = cbuffer.data(ph_geom_11_off + 461 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_0_xxxxx, g_z_y_0_xxxxxx, g_z_y_0_xxxxxy, g_z_y_0_xxxxxz, g_z_y_0_xxxxy, g_z_y_0_xxxxyy, g_z_y_0_xxxxyz, g_z_y_0_xxxxz, g_z_y_0_xxxxzz, g_z_y_0_xxxyy, g_z_y_0_xxxyyy, g_z_y_0_xxxyyz, g_z_y_0_xxxyz, g_z_y_0_xxxyzz, g_z_y_0_xxxzz, g_z_y_0_xxxzzz, g_z_y_0_xxyyy, g_z_y_0_xxyyyy, g_z_y_0_xxyyyz, g_z_y_0_xxyyz, g_z_y_0_xxyyzz, g_z_y_0_xxyzz, g_z_y_0_xxyzzz, g_z_y_0_xxzzz, g_z_y_0_xxzzzz, g_z_y_0_xyyyy, g_z_y_0_xyyyyy, g_z_y_0_xyyyyz, g_z_y_0_xyyyz, g_z_y_0_xyyyzz, g_z_y_0_xyyzz, g_z_y_0_xyyzzz, g_z_y_0_xyzzz, g_z_y_0_xyzzzz, g_z_y_0_xzzzz, g_z_y_0_xzzzzz, g_z_y_0_yyyyy, g_z_y_0_yyyyz, g_z_y_0_yyyzz, g_z_y_0_yyzzz, g_z_y_0_yzzzz, g_z_y_0_zzzzz, g_z_y_x_xxxxx, g_z_y_x_xxxxy, g_z_y_x_xxxxz, g_z_y_x_xxxyy, g_z_y_x_xxxyz, g_z_y_x_xxxzz, g_z_y_x_xxyyy, g_z_y_x_xxyyz, g_z_y_x_xxyzz, g_z_y_x_xxzzz, g_z_y_x_xyyyy, g_z_y_x_xyyyz, g_z_y_x_xyyzz, g_z_y_x_xyzzz, g_z_y_x_xzzzz, g_z_y_x_yyyyy, g_z_y_x_yyyyz, g_z_y_x_yyyzz, g_z_y_x_yyzzz, g_z_y_x_yzzzz, g_z_y_x_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_x_xxxxx[k] = -g_z_y_0_xxxxx[k] * ab_x + g_z_y_0_xxxxxx[k];

                g_z_y_x_xxxxy[k] = -g_z_y_0_xxxxy[k] * ab_x + g_z_y_0_xxxxxy[k];

                g_z_y_x_xxxxz[k] = -g_z_y_0_xxxxz[k] * ab_x + g_z_y_0_xxxxxz[k];

                g_z_y_x_xxxyy[k] = -g_z_y_0_xxxyy[k] * ab_x + g_z_y_0_xxxxyy[k];

                g_z_y_x_xxxyz[k] = -g_z_y_0_xxxyz[k] * ab_x + g_z_y_0_xxxxyz[k];

                g_z_y_x_xxxzz[k] = -g_z_y_0_xxxzz[k] * ab_x + g_z_y_0_xxxxzz[k];

                g_z_y_x_xxyyy[k] = -g_z_y_0_xxyyy[k] * ab_x + g_z_y_0_xxxyyy[k];

                g_z_y_x_xxyyz[k] = -g_z_y_0_xxyyz[k] * ab_x + g_z_y_0_xxxyyz[k];

                g_z_y_x_xxyzz[k] = -g_z_y_0_xxyzz[k] * ab_x + g_z_y_0_xxxyzz[k];

                g_z_y_x_xxzzz[k] = -g_z_y_0_xxzzz[k] * ab_x + g_z_y_0_xxxzzz[k];

                g_z_y_x_xyyyy[k] = -g_z_y_0_xyyyy[k] * ab_x + g_z_y_0_xxyyyy[k];

                g_z_y_x_xyyyz[k] = -g_z_y_0_xyyyz[k] * ab_x + g_z_y_0_xxyyyz[k];

                g_z_y_x_xyyzz[k] = -g_z_y_0_xyyzz[k] * ab_x + g_z_y_0_xxyyzz[k];

                g_z_y_x_xyzzz[k] = -g_z_y_0_xyzzz[k] * ab_x + g_z_y_0_xxyzzz[k];

                g_z_y_x_xzzzz[k] = -g_z_y_0_xzzzz[k] * ab_x + g_z_y_0_xxzzzz[k];

                g_z_y_x_yyyyy[k] = -g_z_y_0_yyyyy[k] * ab_x + g_z_y_0_xyyyyy[k];

                g_z_y_x_yyyyz[k] = -g_z_y_0_yyyyz[k] * ab_x + g_z_y_0_xyyyyz[k];

                g_z_y_x_yyyzz[k] = -g_z_y_0_yyyzz[k] * ab_x + g_z_y_0_xyyyzz[k];

                g_z_y_x_yyzzz[k] = -g_z_y_0_yyzzz[k] * ab_x + g_z_y_0_xyyzzz[k];

                g_z_y_x_yzzzz[k] = -g_z_y_0_yzzzz[k] * ab_x + g_z_y_0_xyzzzz[k];

                g_z_y_x_zzzzz[k] = -g_z_y_0_zzzzz[k] * ab_x + g_z_y_0_xzzzzz[k];
            }

            /// Set up 462-483 components of targeted buffer : cbuffer.data(

            auto g_z_y_y_xxxxx = cbuffer.data(ph_geom_11_off + 462 * ccomps * dcomps);

            auto g_z_y_y_xxxxy = cbuffer.data(ph_geom_11_off + 463 * ccomps * dcomps);

            auto g_z_y_y_xxxxz = cbuffer.data(ph_geom_11_off + 464 * ccomps * dcomps);

            auto g_z_y_y_xxxyy = cbuffer.data(ph_geom_11_off + 465 * ccomps * dcomps);

            auto g_z_y_y_xxxyz = cbuffer.data(ph_geom_11_off + 466 * ccomps * dcomps);

            auto g_z_y_y_xxxzz = cbuffer.data(ph_geom_11_off + 467 * ccomps * dcomps);

            auto g_z_y_y_xxyyy = cbuffer.data(ph_geom_11_off + 468 * ccomps * dcomps);

            auto g_z_y_y_xxyyz = cbuffer.data(ph_geom_11_off + 469 * ccomps * dcomps);

            auto g_z_y_y_xxyzz = cbuffer.data(ph_geom_11_off + 470 * ccomps * dcomps);

            auto g_z_y_y_xxzzz = cbuffer.data(ph_geom_11_off + 471 * ccomps * dcomps);

            auto g_z_y_y_xyyyy = cbuffer.data(ph_geom_11_off + 472 * ccomps * dcomps);

            auto g_z_y_y_xyyyz = cbuffer.data(ph_geom_11_off + 473 * ccomps * dcomps);

            auto g_z_y_y_xyyzz = cbuffer.data(ph_geom_11_off + 474 * ccomps * dcomps);

            auto g_z_y_y_xyzzz = cbuffer.data(ph_geom_11_off + 475 * ccomps * dcomps);

            auto g_z_y_y_xzzzz = cbuffer.data(ph_geom_11_off + 476 * ccomps * dcomps);

            auto g_z_y_y_yyyyy = cbuffer.data(ph_geom_11_off + 477 * ccomps * dcomps);

            auto g_z_y_y_yyyyz = cbuffer.data(ph_geom_11_off + 478 * ccomps * dcomps);

            auto g_z_y_y_yyyzz = cbuffer.data(ph_geom_11_off + 479 * ccomps * dcomps);

            auto g_z_y_y_yyzzz = cbuffer.data(ph_geom_11_off + 480 * ccomps * dcomps);

            auto g_z_y_y_yzzzz = cbuffer.data(ph_geom_11_off + 481 * ccomps * dcomps);

            auto g_z_y_y_zzzzz = cbuffer.data(ph_geom_11_off + 482 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_0_xxxxx, g_z_0_0_xxxxy, g_z_0_0_xxxxz, g_z_0_0_xxxyy, g_z_0_0_xxxyz, g_z_0_0_xxxzz, g_z_0_0_xxyyy, g_z_0_0_xxyyz, g_z_0_0_xxyzz, g_z_0_0_xxzzz, g_z_0_0_xyyyy, g_z_0_0_xyyyz, g_z_0_0_xyyzz, g_z_0_0_xyzzz, g_z_0_0_xzzzz, g_z_0_0_yyyyy, g_z_0_0_yyyyz, g_z_0_0_yyyzz, g_z_0_0_yyzzz, g_z_0_0_yzzzz, g_z_0_0_zzzzz, g_z_y_0_xxxxx, g_z_y_0_xxxxxy, g_z_y_0_xxxxy, g_z_y_0_xxxxyy, g_z_y_0_xxxxyz, g_z_y_0_xxxxz, g_z_y_0_xxxyy, g_z_y_0_xxxyyy, g_z_y_0_xxxyyz, g_z_y_0_xxxyz, g_z_y_0_xxxyzz, g_z_y_0_xxxzz, g_z_y_0_xxyyy, g_z_y_0_xxyyyy, g_z_y_0_xxyyyz, g_z_y_0_xxyyz, g_z_y_0_xxyyzz, g_z_y_0_xxyzz, g_z_y_0_xxyzzz, g_z_y_0_xxzzz, g_z_y_0_xyyyy, g_z_y_0_xyyyyy, g_z_y_0_xyyyyz, g_z_y_0_xyyyz, g_z_y_0_xyyyzz, g_z_y_0_xyyzz, g_z_y_0_xyyzzz, g_z_y_0_xyzzz, g_z_y_0_xyzzzz, g_z_y_0_xzzzz, g_z_y_0_yyyyy, g_z_y_0_yyyyyy, g_z_y_0_yyyyyz, g_z_y_0_yyyyz, g_z_y_0_yyyyzz, g_z_y_0_yyyzz, g_z_y_0_yyyzzz, g_z_y_0_yyzzz, g_z_y_0_yyzzzz, g_z_y_0_yzzzz, g_z_y_0_yzzzzz, g_z_y_0_zzzzz, g_z_y_y_xxxxx, g_z_y_y_xxxxy, g_z_y_y_xxxxz, g_z_y_y_xxxyy, g_z_y_y_xxxyz, g_z_y_y_xxxzz, g_z_y_y_xxyyy, g_z_y_y_xxyyz, g_z_y_y_xxyzz, g_z_y_y_xxzzz, g_z_y_y_xyyyy, g_z_y_y_xyyyz, g_z_y_y_xyyzz, g_z_y_y_xyzzz, g_z_y_y_xzzzz, g_z_y_y_yyyyy, g_z_y_y_yyyyz, g_z_y_y_yyyzz, g_z_y_y_yyzzz, g_z_y_y_yzzzz, g_z_y_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_y_xxxxx[k] = g_z_0_0_xxxxx[k] - g_z_y_0_xxxxx[k] * ab_y + g_z_y_0_xxxxxy[k];

                g_z_y_y_xxxxy[k] = g_z_0_0_xxxxy[k] - g_z_y_0_xxxxy[k] * ab_y + g_z_y_0_xxxxyy[k];

                g_z_y_y_xxxxz[k] = g_z_0_0_xxxxz[k] - g_z_y_0_xxxxz[k] * ab_y + g_z_y_0_xxxxyz[k];

                g_z_y_y_xxxyy[k] = g_z_0_0_xxxyy[k] - g_z_y_0_xxxyy[k] * ab_y + g_z_y_0_xxxyyy[k];

                g_z_y_y_xxxyz[k] = g_z_0_0_xxxyz[k] - g_z_y_0_xxxyz[k] * ab_y + g_z_y_0_xxxyyz[k];

                g_z_y_y_xxxzz[k] = g_z_0_0_xxxzz[k] - g_z_y_0_xxxzz[k] * ab_y + g_z_y_0_xxxyzz[k];

                g_z_y_y_xxyyy[k] = g_z_0_0_xxyyy[k] - g_z_y_0_xxyyy[k] * ab_y + g_z_y_0_xxyyyy[k];

                g_z_y_y_xxyyz[k] = g_z_0_0_xxyyz[k] - g_z_y_0_xxyyz[k] * ab_y + g_z_y_0_xxyyyz[k];

                g_z_y_y_xxyzz[k] = g_z_0_0_xxyzz[k] - g_z_y_0_xxyzz[k] * ab_y + g_z_y_0_xxyyzz[k];

                g_z_y_y_xxzzz[k] = g_z_0_0_xxzzz[k] - g_z_y_0_xxzzz[k] * ab_y + g_z_y_0_xxyzzz[k];

                g_z_y_y_xyyyy[k] = g_z_0_0_xyyyy[k] - g_z_y_0_xyyyy[k] * ab_y + g_z_y_0_xyyyyy[k];

                g_z_y_y_xyyyz[k] = g_z_0_0_xyyyz[k] - g_z_y_0_xyyyz[k] * ab_y + g_z_y_0_xyyyyz[k];

                g_z_y_y_xyyzz[k] = g_z_0_0_xyyzz[k] - g_z_y_0_xyyzz[k] * ab_y + g_z_y_0_xyyyzz[k];

                g_z_y_y_xyzzz[k] = g_z_0_0_xyzzz[k] - g_z_y_0_xyzzz[k] * ab_y + g_z_y_0_xyyzzz[k];

                g_z_y_y_xzzzz[k] = g_z_0_0_xzzzz[k] - g_z_y_0_xzzzz[k] * ab_y + g_z_y_0_xyzzzz[k];

                g_z_y_y_yyyyy[k] = g_z_0_0_yyyyy[k] - g_z_y_0_yyyyy[k] * ab_y + g_z_y_0_yyyyyy[k];

                g_z_y_y_yyyyz[k] = g_z_0_0_yyyyz[k] - g_z_y_0_yyyyz[k] * ab_y + g_z_y_0_yyyyyz[k];

                g_z_y_y_yyyzz[k] = g_z_0_0_yyyzz[k] - g_z_y_0_yyyzz[k] * ab_y + g_z_y_0_yyyyzz[k];

                g_z_y_y_yyzzz[k] = g_z_0_0_yyzzz[k] - g_z_y_0_yyzzz[k] * ab_y + g_z_y_0_yyyzzz[k];

                g_z_y_y_yzzzz[k] = g_z_0_0_yzzzz[k] - g_z_y_0_yzzzz[k] * ab_y + g_z_y_0_yyzzzz[k];

                g_z_y_y_zzzzz[k] = g_z_0_0_zzzzz[k] - g_z_y_0_zzzzz[k] * ab_y + g_z_y_0_yzzzzz[k];
            }

            /// Set up 483-504 components of targeted buffer : cbuffer.data(

            auto g_z_y_z_xxxxx = cbuffer.data(ph_geom_11_off + 483 * ccomps * dcomps);

            auto g_z_y_z_xxxxy = cbuffer.data(ph_geom_11_off + 484 * ccomps * dcomps);

            auto g_z_y_z_xxxxz = cbuffer.data(ph_geom_11_off + 485 * ccomps * dcomps);

            auto g_z_y_z_xxxyy = cbuffer.data(ph_geom_11_off + 486 * ccomps * dcomps);

            auto g_z_y_z_xxxyz = cbuffer.data(ph_geom_11_off + 487 * ccomps * dcomps);

            auto g_z_y_z_xxxzz = cbuffer.data(ph_geom_11_off + 488 * ccomps * dcomps);

            auto g_z_y_z_xxyyy = cbuffer.data(ph_geom_11_off + 489 * ccomps * dcomps);

            auto g_z_y_z_xxyyz = cbuffer.data(ph_geom_11_off + 490 * ccomps * dcomps);

            auto g_z_y_z_xxyzz = cbuffer.data(ph_geom_11_off + 491 * ccomps * dcomps);

            auto g_z_y_z_xxzzz = cbuffer.data(ph_geom_11_off + 492 * ccomps * dcomps);

            auto g_z_y_z_xyyyy = cbuffer.data(ph_geom_11_off + 493 * ccomps * dcomps);

            auto g_z_y_z_xyyyz = cbuffer.data(ph_geom_11_off + 494 * ccomps * dcomps);

            auto g_z_y_z_xyyzz = cbuffer.data(ph_geom_11_off + 495 * ccomps * dcomps);

            auto g_z_y_z_xyzzz = cbuffer.data(ph_geom_11_off + 496 * ccomps * dcomps);

            auto g_z_y_z_xzzzz = cbuffer.data(ph_geom_11_off + 497 * ccomps * dcomps);

            auto g_z_y_z_yyyyy = cbuffer.data(ph_geom_11_off + 498 * ccomps * dcomps);

            auto g_z_y_z_yyyyz = cbuffer.data(ph_geom_11_off + 499 * ccomps * dcomps);

            auto g_z_y_z_yyyzz = cbuffer.data(ph_geom_11_off + 500 * ccomps * dcomps);

            auto g_z_y_z_yyzzz = cbuffer.data(ph_geom_11_off + 501 * ccomps * dcomps);

            auto g_z_y_z_yzzzz = cbuffer.data(ph_geom_11_off + 502 * ccomps * dcomps);

            auto g_z_y_z_zzzzz = cbuffer.data(ph_geom_11_off + 503 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_0_xxxxx, g_0_y_0_xxxxy, g_0_y_0_xxxxz, g_0_y_0_xxxyy, g_0_y_0_xxxyz, g_0_y_0_xxxzz, g_0_y_0_xxyyy, g_0_y_0_xxyyz, g_0_y_0_xxyzz, g_0_y_0_xxzzz, g_0_y_0_xyyyy, g_0_y_0_xyyyz, g_0_y_0_xyyzz, g_0_y_0_xyzzz, g_0_y_0_xzzzz, g_0_y_0_yyyyy, g_0_y_0_yyyyz, g_0_y_0_yyyzz, g_0_y_0_yyzzz, g_0_y_0_yzzzz, g_0_y_0_zzzzz, g_z_y_0_xxxxx, g_z_y_0_xxxxxz, g_z_y_0_xxxxy, g_z_y_0_xxxxyz, g_z_y_0_xxxxz, g_z_y_0_xxxxzz, g_z_y_0_xxxyy, g_z_y_0_xxxyyz, g_z_y_0_xxxyz, g_z_y_0_xxxyzz, g_z_y_0_xxxzz, g_z_y_0_xxxzzz, g_z_y_0_xxyyy, g_z_y_0_xxyyyz, g_z_y_0_xxyyz, g_z_y_0_xxyyzz, g_z_y_0_xxyzz, g_z_y_0_xxyzzz, g_z_y_0_xxzzz, g_z_y_0_xxzzzz, g_z_y_0_xyyyy, g_z_y_0_xyyyyz, g_z_y_0_xyyyz, g_z_y_0_xyyyzz, g_z_y_0_xyyzz, g_z_y_0_xyyzzz, g_z_y_0_xyzzz, g_z_y_0_xyzzzz, g_z_y_0_xzzzz, g_z_y_0_xzzzzz, g_z_y_0_yyyyy, g_z_y_0_yyyyyz, g_z_y_0_yyyyz, g_z_y_0_yyyyzz, g_z_y_0_yyyzz, g_z_y_0_yyyzzz, g_z_y_0_yyzzz, g_z_y_0_yyzzzz, g_z_y_0_yzzzz, g_z_y_0_yzzzzz, g_z_y_0_zzzzz, g_z_y_0_zzzzzz, g_z_y_z_xxxxx, g_z_y_z_xxxxy, g_z_y_z_xxxxz, g_z_y_z_xxxyy, g_z_y_z_xxxyz, g_z_y_z_xxxzz, g_z_y_z_xxyyy, g_z_y_z_xxyyz, g_z_y_z_xxyzz, g_z_y_z_xxzzz, g_z_y_z_xyyyy, g_z_y_z_xyyyz, g_z_y_z_xyyzz, g_z_y_z_xyzzz, g_z_y_z_xzzzz, g_z_y_z_yyyyy, g_z_y_z_yyyyz, g_z_y_z_yyyzz, g_z_y_z_yyzzz, g_z_y_z_yzzzz, g_z_y_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_z_xxxxx[k] = -g_0_y_0_xxxxx[k] - g_z_y_0_xxxxx[k] * ab_z + g_z_y_0_xxxxxz[k];

                g_z_y_z_xxxxy[k] = -g_0_y_0_xxxxy[k] - g_z_y_0_xxxxy[k] * ab_z + g_z_y_0_xxxxyz[k];

                g_z_y_z_xxxxz[k] = -g_0_y_0_xxxxz[k] - g_z_y_0_xxxxz[k] * ab_z + g_z_y_0_xxxxzz[k];

                g_z_y_z_xxxyy[k] = -g_0_y_0_xxxyy[k] - g_z_y_0_xxxyy[k] * ab_z + g_z_y_0_xxxyyz[k];

                g_z_y_z_xxxyz[k] = -g_0_y_0_xxxyz[k] - g_z_y_0_xxxyz[k] * ab_z + g_z_y_0_xxxyzz[k];

                g_z_y_z_xxxzz[k] = -g_0_y_0_xxxzz[k] - g_z_y_0_xxxzz[k] * ab_z + g_z_y_0_xxxzzz[k];

                g_z_y_z_xxyyy[k] = -g_0_y_0_xxyyy[k] - g_z_y_0_xxyyy[k] * ab_z + g_z_y_0_xxyyyz[k];

                g_z_y_z_xxyyz[k] = -g_0_y_0_xxyyz[k] - g_z_y_0_xxyyz[k] * ab_z + g_z_y_0_xxyyzz[k];

                g_z_y_z_xxyzz[k] = -g_0_y_0_xxyzz[k] - g_z_y_0_xxyzz[k] * ab_z + g_z_y_0_xxyzzz[k];

                g_z_y_z_xxzzz[k] = -g_0_y_0_xxzzz[k] - g_z_y_0_xxzzz[k] * ab_z + g_z_y_0_xxzzzz[k];

                g_z_y_z_xyyyy[k] = -g_0_y_0_xyyyy[k] - g_z_y_0_xyyyy[k] * ab_z + g_z_y_0_xyyyyz[k];

                g_z_y_z_xyyyz[k] = -g_0_y_0_xyyyz[k] - g_z_y_0_xyyyz[k] * ab_z + g_z_y_0_xyyyzz[k];

                g_z_y_z_xyyzz[k] = -g_0_y_0_xyyzz[k] - g_z_y_0_xyyzz[k] * ab_z + g_z_y_0_xyyzzz[k];

                g_z_y_z_xyzzz[k] = -g_0_y_0_xyzzz[k] - g_z_y_0_xyzzz[k] * ab_z + g_z_y_0_xyzzzz[k];

                g_z_y_z_xzzzz[k] = -g_0_y_0_xzzzz[k] - g_z_y_0_xzzzz[k] * ab_z + g_z_y_0_xzzzzz[k];

                g_z_y_z_yyyyy[k] = -g_0_y_0_yyyyy[k] - g_z_y_0_yyyyy[k] * ab_z + g_z_y_0_yyyyyz[k];

                g_z_y_z_yyyyz[k] = -g_0_y_0_yyyyz[k] - g_z_y_0_yyyyz[k] * ab_z + g_z_y_0_yyyyzz[k];

                g_z_y_z_yyyzz[k] = -g_0_y_0_yyyzz[k] - g_z_y_0_yyyzz[k] * ab_z + g_z_y_0_yyyzzz[k];

                g_z_y_z_yyzzz[k] = -g_0_y_0_yyzzz[k] - g_z_y_0_yyzzz[k] * ab_z + g_z_y_0_yyzzzz[k];

                g_z_y_z_yzzzz[k] = -g_0_y_0_yzzzz[k] - g_z_y_0_yzzzz[k] * ab_z + g_z_y_0_yzzzzz[k];

                g_z_y_z_zzzzz[k] = -g_0_y_0_zzzzz[k] - g_z_y_0_zzzzz[k] * ab_z + g_z_y_0_zzzzzz[k];
            }

            /// Set up 504-525 components of targeted buffer : cbuffer.data(

            auto g_z_z_x_xxxxx = cbuffer.data(ph_geom_11_off + 504 * ccomps * dcomps);

            auto g_z_z_x_xxxxy = cbuffer.data(ph_geom_11_off + 505 * ccomps * dcomps);

            auto g_z_z_x_xxxxz = cbuffer.data(ph_geom_11_off + 506 * ccomps * dcomps);

            auto g_z_z_x_xxxyy = cbuffer.data(ph_geom_11_off + 507 * ccomps * dcomps);

            auto g_z_z_x_xxxyz = cbuffer.data(ph_geom_11_off + 508 * ccomps * dcomps);

            auto g_z_z_x_xxxzz = cbuffer.data(ph_geom_11_off + 509 * ccomps * dcomps);

            auto g_z_z_x_xxyyy = cbuffer.data(ph_geom_11_off + 510 * ccomps * dcomps);

            auto g_z_z_x_xxyyz = cbuffer.data(ph_geom_11_off + 511 * ccomps * dcomps);

            auto g_z_z_x_xxyzz = cbuffer.data(ph_geom_11_off + 512 * ccomps * dcomps);

            auto g_z_z_x_xxzzz = cbuffer.data(ph_geom_11_off + 513 * ccomps * dcomps);

            auto g_z_z_x_xyyyy = cbuffer.data(ph_geom_11_off + 514 * ccomps * dcomps);

            auto g_z_z_x_xyyyz = cbuffer.data(ph_geom_11_off + 515 * ccomps * dcomps);

            auto g_z_z_x_xyyzz = cbuffer.data(ph_geom_11_off + 516 * ccomps * dcomps);

            auto g_z_z_x_xyzzz = cbuffer.data(ph_geom_11_off + 517 * ccomps * dcomps);

            auto g_z_z_x_xzzzz = cbuffer.data(ph_geom_11_off + 518 * ccomps * dcomps);

            auto g_z_z_x_yyyyy = cbuffer.data(ph_geom_11_off + 519 * ccomps * dcomps);

            auto g_z_z_x_yyyyz = cbuffer.data(ph_geom_11_off + 520 * ccomps * dcomps);

            auto g_z_z_x_yyyzz = cbuffer.data(ph_geom_11_off + 521 * ccomps * dcomps);

            auto g_z_z_x_yyzzz = cbuffer.data(ph_geom_11_off + 522 * ccomps * dcomps);

            auto g_z_z_x_yzzzz = cbuffer.data(ph_geom_11_off + 523 * ccomps * dcomps);

            auto g_z_z_x_zzzzz = cbuffer.data(ph_geom_11_off + 524 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_0_xxxxx, g_z_z_0_xxxxxx, g_z_z_0_xxxxxy, g_z_z_0_xxxxxz, g_z_z_0_xxxxy, g_z_z_0_xxxxyy, g_z_z_0_xxxxyz, g_z_z_0_xxxxz, g_z_z_0_xxxxzz, g_z_z_0_xxxyy, g_z_z_0_xxxyyy, g_z_z_0_xxxyyz, g_z_z_0_xxxyz, g_z_z_0_xxxyzz, g_z_z_0_xxxzz, g_z_z_0_xxxzzz, g_z_z_0_xxyyy, g_z_z_0_xxyyyy, g_z_z_0_xxyyyz, g_z_z_0_xxyyz, g_z_z_0_xxyyzz, g_z_z_0_xxyzz, g_z_z_0_xxyzzz, g_z_z_0_xxzzz, g_z_z_0_xxzzzz, g_z_z_0_xyyyy, g_z_z_0_xyyyyy, g_z_z_0_xyyyyz, g_z_z_0_xyyyz, g_z_z_0_xyyyzz, g_z_z_0_xyyzz, g_z_z_0_xyyzzz, g_z_z_0_xyzzz, g_z_z_0_xyzzzz, g_z_z_0_xzzzz, g_z_z_0_xzzzzz, g_z_z_0_yyyyy, g_z_z_0_yyyyz, g_z_z_0_yyyzz, g_z_z_0_yyzzz, g_z_z_0_yzzzz, g_z_z_0_zzzzz, g_z_z_x_xxxxx, g_z_z_x_xxxxy, g_z_z_x_xxxxz, g_z_z_x_xxxyy, g_z_z_x_xxxyz, g_z_z_x_xxxzz, g_z_z_x_xxyyy, g_z_z_x_xxyyz, g_z_z_x_xxyzz, g_z_z_x_xxzzz, g_z_z_x_xyyyy, g_z_z_x_xyyyz, g_z_z_x_xyyzz, g_z_z_x_xyzzz, g_z_z_x_xzzzz, g_z_z_x_yyyyy, g_z_z_x_yyyyz, g_z_z_x_yyyzz, g_z_z_x_yyzzz, g_z_z_x_yzzzz, g_z_z_x_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_x_xxxxx[k] = -g_z_z_0_xxxxx[k] * ab_x + g_z_z_0_xxxxxx[k];

                g_z_z_x_xxxxy[k] = -g_z_z_0_xxxxy[k] * ab_x + g_z_z_0_xxxxxy[k];

                g_z_z_x_xxxxz[k] = -g_z_z_0_xxxxz[k] * ab_x + g_z_z_0_xxxxxz[k];

                g_z_z_x_xxxyy[k] = -g_z_z_0_xxxyy[k] * ab_x + g_z_z_0_xxxxyy[k];

                g_z_z_x_xxxyz[k] = -g_z_z_0_xxxyz[k] * ab_x + g_z_z_0_xxxxyz[k];

                g_z_z_x_xxxzz[k] = -g_z_z_0_xxxzz[k] * ab_x + g_z_z_0_xxxxzz[k];

                g_z_z_x_xxyyy[k] = -g_z_z_0_xxyyy[k] * ab_x + g_z_z_0_xxxyyy[k];

                g_z_z_x_xxyyz[k] = -g_z_z_0_xxyyz[k] * ab_x + g_z_z_0_xxxyyz[k];

                g_z_z_x_xxyzz[k] = -g_z_z_0_xxyzz[k] * ab_x + g_z_z_0_xxxyzz[k];

                g_z_z_x_xxzzz[k] = -g_z_z_0_xxzzz[k] * ab_x + g_z_z_0_xxxzzz[k];

                g_z_z_x_xyyyy[k] = -g_z_z_0_xyyyy[k] * ab_x + g_z_z_0_xxyyyy[k];

                g_z_z_x_xyyyz[k] = -g_z_z_0_xyyyz[k] * ab_x + g_z_z_0_xxyyyz[k];

                g_z_z_x_xyyzz[k] = -g_z_z_0_xyyzz[k] * ab_x + g_z_z_0_xxyyzz[k];

                g_z_z_x_xyzzz[k] = -g_z_z_0_xyzzz[k] * ab_x + g_z_z_0_xxyzzz[k];

                g_z_z_x_xzzzz[k] = -g_z_z_0_xzzzz[k] * ab_x + g_z_z_0_xxzzzz[k];

                g_z_z_x_yyyyy[k] = -g_z_z_0_yyyyy[k] * ab_x + g_z_z_0_xyyyyy[k];

                g_z_z_x_yyyyz[k] = -g_z_z_0_yyyyz[k] * ab_x + g_z_z_0_xyyyyz[k];

                g_z_z_x_yyyzz[k] = -g_z_z_0_yyyzz[k] * ab_x + g_z_z_0_xyyyzz[k];

                g_z_z_x_yyzzz[k] = -g_z_z_0_yyzzz[k] * ab_x + g_z_z_0_xyyzzz[k];

                g_z_z_x_yzzzz[k] = -g_z_z_0_yzzzz[k] * ab_x + g_z_z_0_xyzzzz[k];

                g_z_z_x_zzzzz[k] = -g_z_z_0_zzzzz[k] * ab_x + g_z_z_0_xzzzzz[k];
            }

            /// Set up 525-546 components of targeted buffer : cbuffer.data(

            auto g_z_z_y_xxxxx = cbuffer.data(ph_geom_11_off + 525 * ccomps * dcomps);

            auto g_z_z_y_xxxxy = cbuffer.data(ph_geom_11_off + 526 * ccomps * dcomps);

            auto g_z_z_y_xxxxz = cbuffer.data(ph_geom_11_off + 527 * ccomps * dcomps);

            auto g_z_z_y_xxxyy = cbuffer.data(ph_geom_11_off + 528 * ccomps * dcomps);

            auto g_z_z_y_xxxyz = cbuffer.data(ph_geom_11_off + 529 * ccomps * dcomps);

            auto g_z_z_y_xxxzz = cbuffer.data(ph_geom_11_off + 530 * ccomps * dcomps);

            auto g_z_z_y_xxyyy = cbuffer.data(ph_geom_11_off + 531 * ccomps * dcomps);

            auto g_z_z_y_xxyyz = cbuffer.data(ph_geom_11_off + 532 * ccomps * dcomps);

            auto g_z_z_y_xxyzz = cbuffer.data(ph_geom_11_off + 533 * ccomps * dcomps);

            auto g_z_z_y_xxzzz = cbuffer.data(ph_geom_11_off + 534 * ccomps * dcomps);

            auto g_z_z_y_xyyyy = cbuffer.data(ph_geom_11_off + 535 * ccomps * dcomps);

            auto g_z_z_y_xyyyz = cbuffer.data(ph_geom_11_off + 536 * ccomps * dcomps);

            auto g_z_z_y_xyyzz = cbuffer.data(ph_geom_11_off + 537 * ccomps * dcomps);

            auto g_z_z_y_xyzzz = cbuffer.data(ph_geom_11_off + 538 * ccomps * dcomps);

            auto g_z_z_y_xzzzz = cbuffer.data(ph_geom_11_off + 539 * ccomps * dcomps);

            auto g_z_z_y_yyyyy = cbuffer.data(ph_geom_11_off + 540 * ccomps * dcomps);

            auto g_z_z_y_yyyyz = cbuffer.data(ph_geom_11_off + 541 * ccomps * dcomps);

            auto g_z_z_y_yyyzz = cbuffer.data(ph_geom_11_off + 542 * ccomps * dcomps);

            auto g_z_z_y_yyzzz = cbuffer.data(ph_geom_11_off + 543 * ccomps * dcomps);

            auto g_z_z_y_yzzzz = cbuffer.data(ph_geom_11_off + 544 * ccomps * dcomps);

            auto g_z_z_y_zzzzz = cbuffer.data(ph_geom_11_off + 545 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_0_xxxxx, g_z_z_0_xxxxxy, g_z_z_0_xxxxy, g_z_z_0_xxxxyy, g_z_z_0_xxxxyz, g_z_z_0_xxxxz, g_z_z_0_xxxyy, g_z_z_0_xxxyyy, g_z_z_0_xxxyyz, g_z_z_0_xxxyz, g_z_z_0_xxxyzz, g_z_z_0_xxxzz, g_z_z_0_xxyyy, g_z_z_0_xxyyyy, g_z_z_0_xxyyyz, g_z_z_0_xxyyz, g_z_z_0_xxyyzz, g_z_z_0_xxyzz, g_z_z_0_xxyzzz, g_z_z_0_xxzzz, g_z_z_0_xyyyy, g_z_z_0_xyyyyy, g_z_z_0_xyyyyz, g_z_z_0_xyyyz, g_z_z_0_xyyyzz, g_z_z_0_xyyzz, g_z_z_0_xyyzzz, g_z_z_0_xyzzz, g_z_z_0_xyzzzz, g_z_z_0_xzzzz, g_z_z_0_yyyyy, g_z_z_0_yyyyyy, g_z_z_0_yyyyyz, g_z_z_0_yyyyz, g_z_z_0_yyyyzz, g_z_z_0_yyyzz, g_z_z_0_yyyzzz, g_z_z_0_yyzzz, g_z_z_0_yyzzzz, g_z_z_0_yzzzz, g_z_z_0_yzzzzz, g_z_z_0_zzzzz, g_z_z_y_xxxxx, g_z_z_y_xxxxy, g_z_z_y_xxxxz, g_z_z_y_xxxyy, g_z_z_y_xxxyz, g_z_z_y_xxxzz, g_z_z_y_xxyyy, g_z_z_y_xxyyz, g_z_z_y_xxyzz, g_z_z_y_xxzzz, g_z_z_y_xyyyy, g_z_z_y_xyyyz, g_z_z_y_xyyzz, g_z_z_y_xyzzz, g_z_z_y_xzzzz, g_z_z_y_yyyyy, g_z_z_y_yyyyz, g_z_z_y_yyyzz, g_z_z_y_yyzzz, g_z_z_y_yzzzz, g_z_z_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_y_xxxxx[k] = -g_z_z_0_xxxxx[k] * ab_y + g_z_z_0_xxxxxy[k];

                g_z_z_y_xxxxy[k] = -g_z_z_0_xxxxy[k] * ab_y + g_z_z_0_xxxxyy[k];

                g_z_z_y_xxxxz[k] = -g_z_z_0_xxxxz[k] * ab_y + g_z_z_0_xxxxyz[k];

                g_z_z_y_xxxyy[k] = -g_z_z_0_xxxyy[k] * ab_y + g_z_z_0_xxxyyy[k];

                g_z_z_y_xxxyz[k] = -g_z_z_0_xxxyz[k] * ab_y + g_z_z_0_xxxyyz[k];

                g_z_z_y_xxxzz[k] = -g_z_z_0_xxxzz[k] * ab_y + g_z_z_0_xxxyzz[k];

                g_z_z_y_xxyyy[k] = -g_z_z_0_xxyyy[k] * ab_y + g_z_z_0_xxyyyy[k];

                g_z_z_y_xxyyz[k] = -g_z_z_0_xxyyz[k] * ab_y + g_z_z_0_xxyyyz[k];

                g_z_z_y_xxyzz[k] = -g_z_z_0_xxyzz[k] * ab_y + g_z_z_0_xxyyzz[k];

                g_z_z_y_xxzzz[k] = -g_z_z_0_xxzzz[k] * ab_y + g_z_z_0_xxyzzz[k];

                g_z_z_y_xyyyy[k] = -g_z_z_0_xyyyy[k] * ab_y + g_z_z_0_xyyyyy[k];

                g_z_z_y_xyyyz[k] = -g_z_z_0_xyyyz[k] * ab_y + g_z_z_0_xyyyyz[k];

                g_z_z_y_xyyzz[k] = -g_z_z_0_xyyzz[k] * ab_y + g_z_z_0_xyyyzz[k];

                g_z_z_y_xyzzz[k] = -g_z_z_0_xyzzz[k] * ab_y + g_z_z_0_xyyzzz[k];

                g_z_z_y_xzzzz[k] = -g_z_z_0_xzzzz[k] * ab_y + g_z_z_0_xyzzzz[k];

                g_z_z_y_yyyyy[k] = -g_z_z_0_yyyyy[k] * ab_y + g_z_z_0_yyyyyy[k];

                g_z_z_y_yyyyz[k] = -g_z_z_0_yyyyz[k] * ab_y + g_z_z_0_yyyyyz[k];

                g_z_z_y_yyyzz[k] = -g_z_z_0_yyyzz[k] * ab_y + g_z_z_0_yyyyzz[k];

                g_z_z_y_yyzzz[k] = -g_z_z_0_yyzzz[k] * ab_y + g_z_z_0_yyyzzz[k];

                g_z_z_y_yzzzz[k] = -g_z_z_0_yzzzz[k] * ab_y + g_z_z_0_yyzzzz[k];

                g_z_z_y_zzzzz[k] = -g_z_z_0_zzzzz[k] * ab_y + g_z_z_0_yzzzzz[k];
            }

            /// Set up 546-567 components of targeted buffer : cbuffer.data(

            auto g_z_z_z_xxxxx = cbuffer.data(ph_geom_11_off + 546 * ccomps * dcomps);

            auto g_z_z_z_xxxxy = cbuffer.data(ph_geom_11_off + 547 * ccomps * dcomps);

            auto g_z_z_z_xxxxz = cbuffer.data(ph_geom_11_off + 548 * ccomps * dcomps);

            auto g_z_z_z_xxxyy = cbuffer.data(ph_geom_11_off + 549 * ccomps * dcomps);

            auto g_z_z_z_xxxyz = cbuffer.data(ph_geom_11_off + 550 * ccomps * dcomps);

            auto g_z_z_z_xxxzz = cbuffer.data(ph_geom_11_off + 551 * ccomps * dcomps);

            auto g_z_z_z_xxyyy = cbuffer.data(ph_geom_11_off + 552 * ccomps * dcomps);

            auto g_z_z_z_xxyyz = cbuffer.data(ph_geom_11_off + 553 * ccomps * dcomps);

            auto g_z_z_z_xxyzz = cbuffer.data(ph_geom_11_off + 554 * ccomps * dcomps);

            auto g_z_z_z_xxzzz = cbuffer.data(ph_geom_11_off + 555 * ccomps * dcomps);

            auto g_z_z_z_xyyyy = cbuffer.data(ph_geom_11_off + 556 * ccomps * dcomps);

            auto g_z_z_z_xyyyz = cbuffer.data(ph_geom_11_off + 557 * ccomps * dcomps);

            auto g_z_z_z_xyyzz = cbuffer.data(ph_geom_11_off + 558 * ccomps * dcomps);

            auto g_z_z_z_xyzzz = cbuffer.data(ph_geom_11_off + 559 * ccomps * dcomps);

            auto g_z_z_z_xzzzz = cbuffer.data(ph_geom_11_off + 560 * ccomps * dcomps);

            auto g_z_z_z_yyyyy = cbuffer.data(ph_geom_11_off + 561 * ccomps * dcomps);

            auto g_z_z_z_yyyyz = cbuffer.data(ph_geom_11_off + 562 * ccomps * dcomps);

            auto g_z_z_z_yyyzz = cbuffer.data(ph_geom_11_off + 563 * ccomps * dcomps);

            auto g_z_z_z_yyzzz = cbuffer.data(ph_geom_11_off + 564 * ccomps * dcomps);

            auto g_z_z_z_yzzzz = cbuffer.data(ph_geom_11_off + 565 * ccomps * dcomps);

            auto g_z_z_z_zzzzz = cbuffer.data(ph_geom_11_off + 566 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_0_xxxxx, g_0_z_0_xxxxy, g_0_z_0_xxxxz, g_0_z_0_xxxyy, g_0_z_0_xxxyz, g_0_z_0_xxxzz, g_0_z_0_xxyyy, g_0_z_0_xxyyz, g_0_z_0_xxyzz, g_0_z_0_xxzzz, g_0_z_0_xyyyy, g_0_z_0_xyyyz, g_0_z_0_xyyzz, g_0_z_0_xyzzz, g_0_z_0_xzzzz, g_0_z_0_yyyyy, g_0_z_0_yyyyz, g_0_z_0_yyyzz, g_0_z_0_yyzzz, g_0_z_0_yzzzz, g_0_z_0_zzzzz, g_z_0_0_xxxxx, g_z_0_0_xxxxy, g_z_0_0_xxxxz, g_z_0_0_xxxyy, g_z_0_0_xxxyz, g_z_0_0_xxxzz, g_z_0_0_xxyyy, g_z_0_0_xxyyz, g_z_0_0_xxyzz, g_z_0_0_xxzzz, g_z_0_0_xyyyy, g_z_0_0_xyyyz, g_z_0_0_xyyzz, g_z_0_0_xyzzz, g_z_0_0_xzzzz, g_z_0_0_yyyyy, g_z_0_0_yyyyz, g_z_0_0_yyyzz, g_z_0_0_yyzzz, g_z_0_0_yzzzz, g_z_0_0_zzzzz, g_z_z_0_xxxxx, g_z_z_0_xxxxxz, g_z_z_0_xxxxy, g_z_z_0_xxxxyz, g_z_z_0_xxxxz, g_z_z_0_xxxxzz, g_z_z_0_xxxyy, g_z_z_0_xxxyyz, g_z_z_0_xxxyz, g_z_z_0_xxxyzz, g_z_z_0_xxxzz, g_z_z_0_xxxzzz, g_z_z_0_xxyyy, g_z_z_0_xxyyyz, g_z_z_0_xxyyz, g_z_z_0_xxyyzz, g_z_z_0_xxyzz, g_z_z_0_xxyzzz, g_z_z_0_xxzzz, g_z_z_0_xxzzzz, g_z_z_0_xyyyy, g_z_z_0_xyyyyz, g_z_z_0_xyyyz, g_z_z_0_xyyyzz, g_z_z_0_xyyzz, g_z_z_0_xyyzzz, g_z_z_0_xyzzz, g_z_z_0_xyzzzz, g_z_z_0_xzzzz, g_z_z_0_xzzzzz, g_z_z_0_yyyyy, g_z_z_0_yyyyyz, g_z_z_0_yyyyz, g_z_z_0_yyyyzz, g_z_z_0_yyyzz, g_z_z_0_yyyzzz, g_z_z_0_yyzzz, g_z_z_0_yyzzzz, g_z_z_0_yzzzz, g_z_z_0_yzzzzz, g_z_z_0_zzzzz, g_z_z_0_zzzzzz, g_z_z_z_xxxxx, g_z_z_z_xxxxy, g_z_z_z_xxxxz, g_z_z_z_xxxyy, g_z_z_z_xxxyz, g_z_z_z_xxxzz, g_z_z_z_xxyyy, g_z_z_z_xxyyz, g_z_z_z_xxyzz, g_z_z_z_xxzzz, g_z_z_z_xyyyy, g_z_z_z_xyyyz, g_z_z_z_xyyzz, g_z_z_z_xyzzz, g_z_z_z_xzzzz, g_z_z_z_yyyyy, g_z_z_z_yyyyz, g_z_z_z_yyyzz, g_z_z_z_yyzzz, g_z_z_z_yzzzz, g_z_z_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_z_xxxxx[k] = -g_0_z_0_xxxxx[k] + g_z_0_0_xxxxx[k] - g_z_z_0_xxxxx[k] * ab_z + g_z_z_0_xxxxxz[k];

                g_z_z_z_xxxxy[k] = -g_0_z_0_xxxxy[k] + g_z_0_0_xxxxy[k] - g_z_z_0_xxxxy[k] * ab_z + g_z_z_0_xxxxyz[k];

                g_z_z_z_xxxxz[k] = -g_0_z_0_xxxxz[k] + g_z_0_0_xxxxz[k] - g_z_z_0_xxxxz[k] * ab_z + g_z_z_0_xxxxzz[k];

                g_z_z_z_xxxyy[k] = -g_0_z_0_xxxyy[k] + g_z_0_0_xxxyy[k] - g_z_z_0_xxxyy[k] * ab_z + g_z_z_0_xxxyyz[k];

                g_z_z_z_xxxyz[k] = -g_0_z_0_xxxyz[k] + g_z_0_0_xxxyz[k] - g_z_z_0_xxxyz[k] * ab_z + g_z_z_0_xxxyzz[k];

                g_z_z_z_xxxzz[k] = -g_0_z_0_xxxzz[k] + g_z_0_0_xxxzz[k] - g_z_z_0_xxxzz[k] * ab_z + g_z_z_0_xxxzzz[k];

                g_z_z_z_xxyyy[k] = -g_0_z_0_xxyyy[k] + g_z_0_0_xxyyy[k] - g_z_z_0_xxyyy[k] * ab_z + g_z_z_0_xxyyyz[k];

                g_z_z_z_xxyyz[k] = -g_0_z_0_xxyyz[k] + g_z_0_0_xxyyz[k] - g_z_z_0_xxyyz[k] * ab_z + g_z_z_0_xxyyzz[k];

                g_z_z_z_xxyzz[k] = -g_0_z_0_xxyzz[k] + g_z_0_0_xxyzz[k] - g_z_z_0_xxyzz[k] * ab_z + g_z_z_0_xxyzzz[k];

                g_z_z_z_xxzzz[k] = -g_0_z_0_xxzzz[k] + g_z_0_0_xxzzz[k] - g_z_z_0_xxzzz[k] * ab_z + g_z_z_0_xxzzzz[k];

                g_z_z_z_xyyyy[k] = -g_0_z_0_xyyyy[k] + g_z_0_0_xyyyy[k] - g_z_z_0_xyyyy[k] * ab_z + g_z_z_0_xyyyyz[k];

                g_z_z_z_xyyyz[k] = -g_0_z_0_xyyyz[k] + g_z_0_0_xyyyz[k] - g_z_z_0_xyyyz[k] * ab_z + g_z_z_0_xyyyzz[k];

                g_z_z_z_xyyzz[k] = -g_0_z_0_xyyzz[k] + g_z_0_0_xyyzz[k] - g_z_z_0_xyyzz[k] * ab_z + g_z_z_0_xyyzzz[k];

                g_z_z_z_xyzzz[k] = -g_0_z_0_xyzzz[k] + g_z_0_0_xyzzz[k] - g_z_z_0_xyzzz[k] * ab_z + g_z_z_0_xyzzzz[k];

                g_z_z_z_xzzzz[k] = -g_0_z_0_xzzzz[k] + g_z_0_0_xzzzz[k] - g_z_z_0_xzzzz[k] * ab_z + g_z_z_0_xzzzzz[k];

                g_z_z_z_yyyyy[k] = -g_0_z_0_yyyyy[k] + g_z_0_0_yyyyy[k] - g_z_z_0_yyyyy[k] * ab_z + g_z_z_0_yyyyyz[k];

                g_z_z_z_yyyyz[k] = -g_0_z_0_yyyyz[k] + g_z_0_0_yyyyz[k] - g_z_z_0_yyyyz[k] * ab_z + g_z_z_0_yyyyzz[k];

                g_z_z_z_yyyzz[k] = -g_0_z_0_yyyzz[k] + g_z_0_0_yyyzz[k] - g_z_z_0_yyyzz[k] * ab_z + g_z_z_0_yyyzzz[k];

                g_z_z_z_yyzzz[k] = -g_0_z_0_yyzzz[k] + g_z_0_0_yyzzz[k] - g_z_z_0_yyzzz[k] * ab_z + g_z_z_0_yyzzzz[k];

                g_z_z_z_yzzzz[k] = -g_0_z_0_yzzzz[k] + g_z_0_0_yzzzz[k] - g_z_z_0_yzzzz[k] * ab_z + g_z_z_0_yzzzzz[k];

                g_z_z_z_zzzzz[k] = -g_0_z_0_zzzzz[k] + g_z_0_0_zzzzz[k] - g_z_z_0_zzzzz[k] * ab_z + g_z_z_0_zzzzzz[k];
            }
        }
    }
}

} // erirec namespace

