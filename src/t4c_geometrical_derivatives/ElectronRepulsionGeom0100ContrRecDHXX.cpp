#include "ElectronRepulsionGeom0100ContrRecDHXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_dhxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_dhxx,
                                            const size_t idx_phxx,
                                            const size_t idx_geom_01_phxx,
                                            const size_t idx_geom_01_pixx,
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
            /// Set up components of auxilary buffer : PHSS

            const auto ph_off = idx_phxx + i * dcomps + j;

            auto g_x_xxxxx = cbuffer.data(ph_off + 0 * ccomps * dcomps);

            auto g_x_xxxxy = cbuffer.data(ph_off + 1 * ccomps * dcomps);

            auto g_x_xxxxz = cbuffer.data(ph_off + 2 * ccomps * dcomps);

            auto g_x_xxxyy = cbuffer.data(ph_off + 3 * ccomps * dcomps);

            auto g_x_xxxyz = cbuffer.data(ph_off + 4 * ccomps * dcomps);

            auto g_x_xxxzz = cbuffer.data(ph_off + 5 * ccomps * dcomps);

            auto g_x_xxyyy = cbuffer.data(ph_off + 6 * ccomps * dcomps);

            auto g_x_xxyyz = cbuffer.data(ph_off + 7 * ccomps * dcomps);

            auto g_x_xxyzz = cbuffer.data(ph_off + 8 * ccomps * dcomps);

            auto g_x_xxzzz = cbuffer.data(ph_off + 9 * ccomps * dcomps);

            auto g_x_xyyyy = cbuffer.data(ph_off + 10 * ccomps * dcomps);

            auto g_x_xyyyz = cbuffer.data(ph_off + 11 * ccomps * dcomps);

            auto g_x_xyyzz = cbuffer.data(ph_off + 12 * ccomps * dcomps);

            auto g_x_xyzzz = cbuffer.data(ph_off + 13 * ccomps * dcomps);

            auto g_x_xzzzz = cbuffer.data(ph_off + 14 * ccomps * dcomps);

            auto g_x_yyyyy = cbuffer.data(ph_off + 15 * ccomps * dcomps);

            auto g_x_yyyyz = cbuffer.data(ph_off + 16 * ccomps * dcomps);

            auto g_x_yyyzz = cbuffer.data(ph_off + 17 * ccomps * dcomps);

            auto g_x_yyzzz = cbuffer.data(ph_off + 18 * ccomps * dcomps);

            auto g_x_yzzzz = cbuffer.data(ph_off + 19 * ccomps * dcomps);

            auto g_x_zzzzz = cbuffer.data(ph_off + 20 * ccomps * dcomps);

            auto g_y_xxxxx = cbuffer.data(ph_off + 21 * ccomps * dcomps);

            auto g_y_xxxxy = cbuffer.data(ph_off + 22 * ccomps * dcomps);

            auto g_y_xxxxz = cbuffer.data(ph_off + 23 * ccomps * dcomps);

            auto g_y_xxxyy = cbuffer.data(ph_off + 24 * ccomps * dcomps);

            auto g_y_xxxyz = cbuffer.data(ph_off + 25 * ccomps * dcomps);

            auto g_y_xxxzz = cbuffer.data(ph_off + 26 * ccomps * dcomps);

            auto g_y_xxyyy = cbuffer.data(ph_off + 27 * ccomps * dcomps);

            auto g_y_xxyyz = cbuffer.data(ph_off + 28 * ccomps * dcomps);

            auto g_y_xxyzz = cbuffer.data(ph_off + 29 * ccomps * dcomps);

            auto g_y_xxzzz = cbuffer.data(ph_off + 30 * ccomps * dcomps);

            auto g_y_xyyyy = cbuffer.data(ph_off + 31 * ccomps * dcomps);

            auto g_y_xyyyz = cbuffer.data(ph_off + 32 * ccomps * dcomps);

            auto g_y_xyyzz = cbuffer.data(ph_off + 33 * ccomps * dcomps);

            auto g_y_xyzzz = cbuffer.data(ph_off + 34 * ccomps * dcomps);

            auto g_y_xzzzz = cbuffer.data(ph_off + 35 * ccomps * dcomps);

            auto g_y_yyyyy = cbuffer.data(ph_off + 36 * ccomps * dcomps);

            auto g_y_yyyyz = cbuffer.data(ph_off + 37 * ccomps * dcomps);

            auto g_y_yyyzz = cbuffer.data(ph_off + 38 * ccomps * dcomps);

            auto g_y_yyzzz = cbuffer.data(ph_off + 39 * ccomps * dcomps);

            auto g_y_yzzzz = cbuffer.data(ph_off + 40 * ccomps * dcomps);

            auto g_y_zzzzz = cbuffer.data(ph_off + 41 * ccomps * dcomps);

            auto g_z_xxxxx = cbuffer.data(ph_off + 42 * ccomps * dcomps);

            auto g_z_xxxxy = cbuffer.data(ph_off + 43 * ccomps * dcomps);

            auto g_z_xxxxz = cbuffer.data(ph_off + 44 * ccomps * dcomps);

            auto g_z_xxxyy = cbuffer.data(ph_off + 45 * ccomps * dcomps);

            auto g_z_xxxyz = cbuffer.data(ph_off + 46 * ccomps * dcomps);

            auto g_z_xxxzz = cbuffer.data(ph_off + 47 * ccomps * dcomps);

            auto g_z_xxyyy = cbuffer.data(ph_off + 48 * ccomps * dcomps);

            auto g_z_xxyyz = cbuffer.data(ph_off + 49 * ccomps * dcomps);

            auto g_z_xxyzz = cbuffer.data(ph_off + 50 * ccomps * dcomps);

            auto g_z_xxzzz = cbuffer.data(ph_off + 51 * ccomps * dcomps);

            auto g_z_xyyyy = cbuffer.data(ph_off + 52 * ccomps * dcomps);

            auto g_z_xyyyz = cbuffer.data(ph_off + 53 * ccomps * dcomps);

            auto g_z_xyyzz = cbuffer.data(ph_off + 54 * ccomps * dcomps);

            auto g_z_xyzzz = cbuffer.data(ph_off + 55 * ccomps * dcomps);

            auto g_z_xzzzz = cbuffer.data(ph_off + 56 * ccomps * dcomps);

            auto g_z_yyyyy = cbuffer.data(ph_off + 57 * ccomps * dcomps);

            auto g_z_yyyyz = cbuffer.data(ph_off + 58 * ccomps * dcomps);

            auto g_z_yyyzz = cbuffer.data(ph_off + 59 * ccomps * dcomps);

            auto g_z_yyzzz = cbuffer.data(ph_off + 60 * ccomps * dcomps);

            auto g_z_yzzzz = cbuffer.data(ph_off + 61 * ccomps * dcomps);

            auto g_z_zzzzz = cbuffer.data(ph_off + 62 * ccomps * dcomps);

            /// Set up components of auxilary buffer : PHSS

            const auto ph_geom_01_off = idx_geom_01_phxx + i * dcomps + j;

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

            /// Set up components of auxilary buffer : PISS

            const auto pi_geom_01_off = idx_geom_01_pixx + i * dcomps + j;

            auto g_0_x_x_xxxxxx = cbuffer.data(pi_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_x_xxxxxy = cbuffer.data(pi_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_x_xxxxxz = cbuffer.data(pi_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_x_xxxxyy = cbuffer.data(pi_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_x_xxxxyz = cbuffer.data(pi_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_x_xxxxzz = cbuffer.data(pi_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_x_xxxyyy = cbuffer.data(pi_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_x_xxxyyz = cbuffer.data(pi_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_x_xxxyzz = cbuffer.data(pi_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_x_xxxzzz = cbuffer.data(pi_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_x_xxyyyy = cbuffer.data(pi_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_x_xxyyyz = cbuffer.data(pi_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_x_xxyyzz = cbuffer.data(pi_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_x_xxyzzz = cbuffer.data(pi_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_x_xxzzzz = cbuffer.data(pi_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_x_xyyyyy = cbuffer.data(pi_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_x_xyyyyz = cbuffer.data(pi_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_x_xyyyzz = cbuffer.data(pi_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_x_xyyzzz = cbuffer.data(pi_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_x_xyzzzz = cbuffer.data(pi_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_x_xzzzzz = cbuffer.data(pi_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_x_yyyyyy = cbuffer.data(pi_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_x_yyyyyz = cbuffer.data(pi_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_x_yyyyzz = cbuffer.data(pi_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_x_yyyzzz = cbuffer.data(pi_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_x_yyzzzz = cbuffer.data(pi_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_x_yzzzzz = cbuffer.data(pi_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_x_zzzzzz = cbuffer.data(pi_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_y_xxxxxx = cbuffer.data(pi_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_y_xxxxxy = cbuffer.data(pi_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_y_xxxxxz = cbuffer.data(pi_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_y_xxxxyy = cbuffer.data(pi_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_y_xxxxyz = cbuffer.data(pi_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_y_xxxxzz = cbuffer.data(pi_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_y_xxxyyy = cbuffer.data(pi_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_y_xxxyyz = cbuffer.data(pi_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_y_xxxyzz = cbuffer.data(pi_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_y_xxxzzz = cbuffer.data(pi_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_y_xxyyyy = cbuffer.data(pi_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_y_xxyyyz = cbuffer.data(pi_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_y_xxyyzz = cbuffer.data(pi_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_y_xxyzzz = cbuffer.data(pi_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_y_xxzzzz = cbuffer.data(pi_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_y_xyyyyy = cbuffer.data(pi_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_y_xyyyyz = cbuffer.data(pi_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_y_xyyyzz = cbuffer.data(pi_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_y_xyyzzz = cbuffer.data(pi_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_y_xyzzzz = cbuffer.data(pi_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_y_xzzzzz = cbuffer.data(pi_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_y_yyyyyy = cbuffer.data(pi_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_y_yyyyyz = cbuffer.data(pi_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_y_yyyyzz = cbuffer.data(pi_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_y_yyyzzz = cbuffer.data(pi_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_y_yyzzzz = cbuffer.data(pi_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_y_yzzzzz = cbuffer.data(pi_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_y_zzzzzz = cbuffer.data(pi_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_z_xxxxxx = cbuffer.data(pi_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_z_xxxxxy = cbuffer.data(pi_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_z_xxxxxz = cbuffer.data(pi_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_z_xxxxyy = cbuffer.data(pi_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_z_xxxxyz = cbuffer.data(pi_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_z_xxxxzz = cbuffer.data(pi_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_z_xxxyyy = cbuffer.data(pi_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_z_xxxyyz = cbuffer.data(pi_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_z_xxxyzz = cbuffer.data(pi_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_z_xxxzzz = cbuffer.data(pi_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_z_xxyyyy = cbuffer.data(pi_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_z_xxyyyz = cbuffer.data(pi_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_z_xxyyzz = cbuffer.data(pi_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_z_xxyzzz = cbuffer.data(pi_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_z_xxzzzz = cbuffer.data(pi_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_z_xyyyyy = cbuffer.data(pi_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_z_xyyyyz = cbuffer.data(pi_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_z_xyyyzz = cbuffer.data(pi_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_z_xyyzzz = cbuffer.data(pi_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_z_xyzzzz = cbuffer.data(pi_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_z_xzzzzz = cbuffer.data(pi_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_z_yyyyyy = cbuffer.data(pi_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_z_yyyyyz = cbuffer.data(pi_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_z_yyyyzz = cbuffer.data(pi_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_z_yyyzzz = cbuffer.data(pi_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_z_yyzzzz = cbuffer.data(pi_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_z_yzzzzz = cbuffer.data(pi_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_z_zzzzzz = cbuffer.data(pi_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_y_x_xxxxxx = cbuffer.data(pi_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_y_x_xxxxxy = cbuffer.data(pi_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_y_x_xxxxxz = cbuffer.data(pi_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_y_x_xxxxyy = cbuffer.data(pi_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_y_x_xxxxyz = cbuffer.data(pi_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_y_x_xxxxzz = cbuffer.data(pi_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_y_x_xxxyyy = cbuffer.data(pi_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_y_x_xxxyyz = cbuffer.data(pi_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_y_x_xxxyzz = cbuffer.data(pi_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_y_x_xxxzzz = cbuffer.data(pi_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_y_x_xxyyyy = cbuffer.data(pi_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_y_x_xxyyyz = cbuffer.data(pi_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_y_x_xxyyzz = cbuffer.data(pi_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_y_x_xxyzzz = cbuffer.data(pi_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_y_x_xxzzzz = cbuffer.data(pi_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_y_x_xyyyyy = cbuffer.data(pi_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_y_x_xyyyyz = cbuffer.data(pi_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_y_x_xyyyzz = cbuffer.data(pi_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_y_x_xyyzzz = cbuffer.data(pi_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_y_x_xyzzzz = cbuffer.data(pi_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_y_x_xzzzzz = cbuffer.data(pi_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_y_x_yyyyyy = cbuffer.data(pi_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_y_x_yyyyyz = cbuffer.data(pi_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_y_x_yyyyzz = cbuffer.data(pi_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_y_x_yyyzzz = cbuffer.data(pi_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_y_x_yyzzzz = cbuffer.data(pi_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_y_x_yzzzzz = cbuffer.data(pi_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_y_x_zzzzzz = cbuffer.data(pi_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_y_y_xxxxxx = cbuffer.data(pi_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_y_y_xxxxxy = cbuffer.data(pi_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_y_y_xxxxxz = cbuffer.data(pi_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_y_y_xxxxyy = cbuffer.data(pi_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_y_y_xxxxyz = cbuffer.data(pi_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_y_y_xxxxzz = cbuffer.data(pi_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_y_y_xxxyyy = cbuffer.data(pi_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_y_y_xxxyyz = cbuffer.data(pi_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_y_y_xxxyzz = cbuffer.data(pi_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_y_y_xxxzzz = cbuffer.data(pi_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_y_y_xxyyyy = cbuffer.data(pi_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_y_y_xxyyyz = cbuffer.data(pi_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_y_y_xxyyzz = cbuffer.data(pi_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_y_y_xxyzzz = cbuffer.data(pi_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_y_y_xxzzzz = cbuffer.data(pi_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_y_y_xyyyyy = cbuffer.data(pi_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_y_y_xyyyyz = cbuffer.data(pi_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_y_y_xyyyzz = cbuffer.data(pi_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_y_y_xyyzzz = cbuffer.data(pi_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_y_y_xyzzzz = cbuffer.data(pi_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_y_y_xzzzzz = cbuffer.data(pi_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_y_y_yyyyyy = cbuffer.data(pi_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_y_y_yyyyyz = cbuffer.data(pi_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_y_y_yyyyzz = cbuffer.data(pi_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_y_y_yyyzzz = cbuffer.data(pi_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_y_y_yyzzzz = cbuffer.data(pi_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_y_y_yzzzzz = cbuffer.data(pi_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_y_y_zzzzzz = cbuffer.data(pi_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_y_z_xxxxxx = cbuffer.data(pi_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_y_z_xxxxxy = cbuffer.data(pi_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_y_z_xxxxxz = cbuffer.data(pi_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_y_z_xxxxyy = cbuffer.data(pi_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_y_z_xxxxyz = cbuffer.data(pi_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_y_z_xxxxzz = cbuffer.data(pi_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_y_z_xxxyyy = cbuffer.data(pi_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_y_z_xxxyyz = cbuffer.data(pi_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_y_z_xxxyzz = cbuffer.data(pi_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_y_z_xxxzzz = cbuffer.data(pi_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_y_z_xxyyyy = cbuffer.data(pi_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_y_z_xxyyyz = cbuffer.data(pi_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_y_z_xxyyzz = cbuffer.data(pi_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_y_z_xxyzzz = cbuffer.data(pi_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_y_z_xxzzzz = cbuffer.data(pi_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_y_z_xyyyyy = cbuffer.data(pi_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_y_z_xyyyyz = cbuffer.data(pi_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_y_z_xyyyzz = cbuffer.data(pi_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_y_z_xyyzzz = cbuffer.data(pi_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_y_z_xyzzzz = cbuffer.data(pi_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_y_z_xzzzzz = cbuffer.data(pi_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_y_z_yyyyyy = cbuffer.data(pi_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_y_z_yyyyyz = cbuffer.data(pi_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_y_z_yyyyzz = cbuffer.data(pi_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_y_z_yyyzzz = cbuffer.data(pi_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_y_z_yyzzzz = cbuffer.data(pi_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_y_z_yzzzzz = cbuffer.data(pi_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_y_z_zzzzzz = cbuffer.data(pi_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_z_x_xxxxxx = cbuffer.data(pi_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_z_x_xxxxxy = cbuffer.data(pi_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_z_x_xxxxxz = cbuffer.data(pi_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_z_x_xxxxyy = cbuffer.data(pi_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_z_x_xxxxyz = cbuffer.data(pi_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_z_x_xxxxzz = cbuffer.data(pi_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_z_x_xxxyyy = cbuffer.data(pi_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_z_x_xxxyyz = cbuffer.data(pi_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_z_x_xxxyzz = cbuffer.data(pi_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_z_x_xxxzzz = cbuffer.data(pi_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_z_x_xxyyyy = cbuffer.data(pi_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_z_x_xxyyyz = cbuffer.data(pi_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_z_x_xxyyzz = cbuffer.data(pi_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_z_x_xxyzzz = cbuffer.data(pi_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_z_x_xxzzzz = cbuffer.data(pi_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_z_x_xyyyyy = cbuffer.data(pi_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_z_x_xyyyyz = cbuffer.data(pi_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_z_x_xyyyzz = cbuffer.data(pi_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_z_x_xyyzzz = cbuffer.data(pi_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_z_x_xyzzzz = cbuffer.data(pi_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_z_x_xzzzzz = cbuffer.data(pi_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_z_x_yyyyyy = cbuffer.data(pi_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_z_x_yyyyyz = cbuffer.data(pi_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_z_x_yyyyzz = cbuffer.data(pi_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_z_x_yyyzzz = cbuffer.data(pi_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_z_x_yyzzzz = cbuffer.data(pi_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_z_x_yzzzzz = cbuffer.data(pi_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_z_x_zzzzzz = cbuffer.data(pi_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_z_y_xxxxxx = cbuffer.data(pi_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_z_y_xxxxxy = cbuffer.data(pi_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_z_y_xxxxxz = cbuffer.data(pi_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_z_y_xxxxyy = cbuffer.data(pi_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_z_y_xxxxyz = cbuffer.data(pi_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_z_y_xxxxzz = cbuffer.data(pi_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_z_y_xxxyyy = cbuffer.data(pi_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_z_y_xxxyyz = cbuffer.data(pi_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_z_y_xxxyzz = cbuffer.data(pi_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_z_y_xxxzzz = cbuffer.data(pi_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_z_y_xxyyyy = cbuffer.data(pi_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_z_y_xxyyyz = cbuffer.data(pi_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_z_y_xxyyzz = cbuffer.data(pi_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_z_y_xxyzzz = cbuffer.data(pi_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_z_y_xxzzzz = cbuffer.data(pi_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_z_y_xyyyyy = cbuffer.data(pi_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_z_y_xyyyyz = cbuffer.data(pi_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_z_y_xyyyzz = cbuffer.data(pi_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_z_y_xyyzzz = cbuffer.data(pi_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_z_y_xyzzzz = cbuffer.data(pi_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_z_y_xzzzzz = cbuffer.data(pi_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_z_y_yyyyyy = cbuffer.data(pi_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_z_y_yyyyyz = cbuffer.data(pi_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_z_y_yyyyzz = cbuffer.data(pi_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_z_y_yyyzzz = cbuffer.data(pi_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_z_y_yyzzzz = cbuffer.data(pi_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_z_y_yzzzzz = cbuffer.data(pi_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_z_y_zzzzzz = cbuffer.data(pi_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_z_z_xxxxxx = cbuffer.data(pi_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_z_z_xxxxxy = cbuffer.data(pi_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_z_z_xxxxxz = cbuffer.data(pi_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_z_z_xxxxyy = cbuffer.data(pi_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_z_z_xxxxyz = cbuffer.data(pi_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_z_z_xxxxzz = cbuffer.data(pi_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_z_z_xxxyyy = cbuffer.data(pi_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_z_z_xxxyyz = cbuffer.data(pi_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_z_z_xxxyzz = cbuffer.data(pi_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_z_z_xxxzzz = cbuffer.data(pi_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_z_z_xxyyyy = cbuffer.data(pi_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_z_z_xxyyyz = cbuffer.data(pi_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_z_z_xxyyzz = cbuffer.data(pi_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_z_z_xxyzzz = cbuffer.data(pi_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_z_z_xxzzzz = cbuffer.data(pi_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_z_z_xyyyyy = cbuffer.data(pi_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_z_z_xyyyyz = cbuffer.data(pi_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_z_z_xyyyzz = cbuffer.data(pi_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_z_z_xyyzzz = cbuffer.data(pi_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_z_z_xyzzzz = cbuffer.data(pi_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_z_z_xzzzzz = cbuffer.data(pi_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_z_z_yyyyyy = cbuffer.data(pi_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_z_z_yyyyyz = cbuffer.data(pi_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_z_z_yyyyzz = cbuffer.data(pi_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_z_z_yyyzzz = cbuffer.data(pi_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_z_z_yyzzzz = cbuffer.data(pi_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_z_z_yzzzzz = cbuffer.data(pi_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_z_z_zzzzzz = cbuffer.data(pi_geom_01_off + 251 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_dhxx

            const auto dh_geom_01_off = idx_geom_01_dhxx + i * dcomps + j;

            /// Set up 0-21 components of targeted buffer : cbuffer.data(

            auto g_0_x_xx_xxxxx = cbuffer.data(dh_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xx_xxxxy = cbuffer.data(dh_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xx_xxxxz = cbuffer.data(dh_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xx_xxxyy = cbuffer.data(dh_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xx_xxxyz = cbuffer.data(dh_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xx_xxxzz = cbuffer.data(dh_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xx_xxyyy = cbuffer.data(dh_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xx_xxyyz = cbuffer.data(dh_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xx_xxyzz = cbuffer.data(dh_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xx_xxzzz = cbuffer.data(dh_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xx_xyyyy = cbuffer.data(dh_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xx_xyyyz = cbuffer.data(dh_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xx_xyyzz = cbuffer.data(dh_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xx_xyzzz = cbuffer.data(dh_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xx_xzzzz = cbuffer.data(dh_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xx_yyyyy = cbuffer.data(dh_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xx_yyyyz = cbuffer.data(dh_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xx_yyyzz = cbuffer.data(dh_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xx_yyzzz = cbuffer.data(dh_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xx_yzzzz = cbuffer.data(dh_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xx_zzzzz = cbuffer.data(dh_geom_01_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_x_xxxxx, g_0_x_x_xxxxxx, g_0_x_x_xxxxxy, g_0_x_x_xxxxxz, g_0_x_x_xxxxy, g_0_x_x_xxxxyy, g_0_x_x_xxxxyz, g_0_x_x_xxxxz, g_0_x_x_xxxxzz, g_0_x_x_xxxyy, g_0_x_x_xxxyyy, g_0_x_x_xxxyyz, g_0_x_x_xxxyz, g_0_x_x_xxxyzz, g_0_x_x_xxxzz, g_0_x_x_xxxzzz, g_0_x_x_xxyyy, g_0_x_x_xxyyyy, g_0_x_x_xxyyyz, g_0_x_x_xxyyz, g_0_x_x_xxyyzz, g_0_x_x_xxyzz, g_0_x_x_xxyzzz, g_0_x_x_xxzzz, g_0_x_x_xxzzzz, g_0_x_x_xyyyy, g_0_x_x_xyyyyy, g_0_x_x_xyyyyz, g_0_x_x_xyyyz, g_0_x_x_xyyyzz, g_0_x_x_xyyzz, g_0_x_x_xyyzzz, g_0_x_x_xyzzz, g_0_x_x_xyzzzz, g_0_x_x_xzzzz, g_0_x_x_xzzzzz, g_0_x_x_yyyyy, g_0_x_x_yyyyz, g_0_x_x_yyyzz, g_0_x_x_yyzzz, g_0_x_x_yzzzz, g_0_x_x_zzzzz, g_0_x_xx_xxxxx, g_0_x_xx_xxxxy, g_0_x_xx_xxxxz, g_0_x_xx_xxxyy, g_0_x_xx_xxxyz, g_0_x_xx_xxxzz, g_0_x_xx_xxyyy, g_0_x_xx_xxyyz, g_0_x_xx_xxyzz, g_0_x_xx_xxzzz, g_0_x_xx_xyyyy, g_0_x_xx_xyyyz, g_0_x_xx_xyyzz, g_0_x_xx_xyzzz, g_0_x_xx_xzzzz, g_0_x_xx_yyyyy, g_0_x_xx_yyyyz, g_0_x_xx_yyyzz, g_0_x_xx_yyzzz, g_0_x_xx_yzzzz, g_0_x_xx_zzzzz, g_x_xxxxx, g_x_xxxxy, g_x_xxxxz, g_x_xxxyy, g_x_xxxyz, g_x_xxxzz, g_x_xxyyy, g_x_xxyyz, g_x_xxyzz, g_x_xxzzz, g_x_xyyyy, g_x_xyyyz, g_x_xyyzz, g_x_xyzzz, g_x_xzzzz, g_x_yyyyy, g_x_yyyyz, g_x_yyyzz, g_x_yyzzz, g_x_yzzzz, g_x_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xx_xxxxx[k] = g_x_xxxxx[k] - g_0_x_x_xxxxx[k] * ab_x + g_0_x_x_xxxxxx[k];

                g_0_x_xx_xxxxy[k] = g_x_xxxxy[k] - g_0_x_x_xxxxy[k] * ab_x + g_0_x_x_xxxxxy[k];

                g_0_x_xx_xxxxz[k] = g_x_xxxxz[k] - g_0_x_x_xxxxz[k] * ab_x + g_0_x_x_xxxxxz[k];

                g_0_x_xx_xxxyy[k] = g_x_xxxyy[k] - g_0_x_x_xxxyy[k] * ab_x + g_0_x_x_xxxxyy[k];

                g_0_x_xx_xxxyz[k] = g_x_xxxyz[k] - g_0_x_x_xxxyz[k] * ab_x + g_0_x_x_xxxxyz[k];

                g_0_x_xx_xxxzz[k] = g_x_xxxzz[k] - g_0_x_x_xxxzz[k] * ab_x + g_0_x_x_xxxxzz[k];

                g_0_x_xx_xxyyy[k] = g_x_xxyyy[k] - g_0_x_x_xxyyy[k] * ab_x + g_0_x_x_xxxyyy[k];

                g_0_x_xx_xxyyz[k] = g_x_xxyyz[k] - g_0_x_x_xxyyz[k] * ab_x + g_0_x_x_xxxyyz[k];

                g_0_x_xx_xxyzz[k] = g_x_xxyzz[k] - g_0_x_x_xxyzz[k] * ab_x + g_0_x_x_xxxyzz[k];

                g_0_x_xx_xxzzz[k] = g_x_xxzzz[k] - g_0_x_x_xxzzz[k] * ab_x + g_0_x_x_xxxzzz[k];

                g_0_x_xx_xyyyy[k] = g_x_xyyyy[k] - g_0_x_x_xyyyy[k] * ab_x + g_0_x_x_xxyyyy[k];

                g_0_x_xx_xyyyz[k] = g_x_xyyyz[k] - g_0_x_x_xyyyz[k] * ab_x + g_0_x_x_xxyyyz[k];

                g_0_x_xx_xyyzz[k] = g_x_xyyzz[k] - g_0_x_x_xyyzz[k] * ab_x + g_0_x_x_xxyyzz[k];

                g_0_x_xx_xyzzz[k] = g_x_xyzzz[k] - g_0_x_x_xyzzz[k] * ab_x + g_0_x_x_xxyzzz[k];

                g_0_x_xx_xzzzz[k] = g_x_xzzzz[k] - g_0_x_x_xzzzz[k] * ab_x + g_0_x_x_xxzzzz[k];

                g_0_x_xx_yyyyy[k] = g_x_yyyyy[k] - g_0_x_x_yyyyy[k] * ab_x + g_0_x_x_xyyyyy[k];

                g_0_x_xx_yyyyz[k] = g_x_yyyyz[k] - g_0_x_x_yyyyz[k] * ab_x + g_0_x_x_xyyyyz[k];

                g_0_x_xx_yyyzz[k] = g_x_yyyzz[k] - g_0_x_x_yyyzz[k] * ab_x + g_0_x_x_xyyyzz[k];

                g_0_x_xx_yyzzz[k] = g_x_yyzzz[k] - g_0_x_x_yyzzz[k] * ab_x + g_0_x_x_xyyzzz[k];

                g_0_x_xx_yzzzz[k] = g_x_yzzzz[k] - g_0_x_x_yzzzz[k] * ab_x + g_0_x_x_xyzzzz[k];

                g_0_x_xx_zzzzz[k] = g_x_zzzzz[k] - g_0_x_x_zzzzz[k] * ab_x + g_0_x_x_xzzzzz[k];
            }

            /// Set up 21-42 components of targeted buffer : cbuffer.data(

            auto g_0_x_xy_xxxxx = cbuffer.data(dh_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xy_xxxxy = cbuffer.data(dh_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xy_xxxxz = cbuffer.data(dh_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xy_xxxyy = cbuffer.data(dh_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xy_xxxyz = cbuffer.data(dh_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xy_xxxzz = cbuffer.data(dh_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xy_xxyyy = cbuffer.data(dh_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xy_xxyyz = cbuffer.data(dh_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xy_xxyzz = cbuffer.data(dh_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xy_xxzzz = cbuffer.data(dh_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xy_xyyyy = cbuffer.data(dh_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xy_xyyyz = cbuffer.data(dh_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xy_xyyzz = cbuffer.data(dh_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xy_xyzzz = cbuffer.data(dh_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xy_xzzzz = cbuffer.data(dh_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xy_yyyyy = cbuffer.data(dh_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xy_yyyyz = cbuffer.data(dh_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xy_yyyzz = cbuffer.data(dh_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xy_yyzzz = cbuffer.data(dh_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xy_yzzzz = cbuffer.data(dh_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xy_zzzzz = cbuffer.data(dh_geom_01_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_x_xxxxx, g_0_x_x_xxxxxy, g_0_x_x_xxxxy, g_0_x_x_xxxxyy, g_0_x_x_xxxxyz, g_0_x_x_xxxxz, g_0_x_x_xxxyy, g_0_x_x_xxxyyy, g_0_x_x_xxxyyz, g_0_x_x_xxxyz, g_0_x_x_xxxyzz, g_0_x_x_xxxzz, g_0_x_x_xxyyy, g_0_x_x_xxyyyy, g_0_x_x_xxyyyz, g_0_x_x_xxyyz, g_0_x_x_xxyyzz, g_0_x_x_xxyzz, g_0_x_x_xxyzzz, g_0_x_x_xxzzz, g_0_x_x_xyyyy, g_0_x_x_xyyyyy, g_0_x_x_xyyyyz, g_0_x_x_xyyyz, g_0_x_x_xyyyzz, g_0_x_x_xyyzz, g_0_x_x_xyyzzz, g_0_x_x_xyzzz, g_0_x_x_xyzzzz, g_0_x_x_xzzzz, g_0_x_x_yyyyy, g_0_x_x_yyyyyy, g_0_x_x_yyyyyz, g_0_x_x_yyyyz, g_0_x_x_yyyyzz, g_0_x_x_yyyzz, g_0_x_x_yyyzzz, g_0_x_x_yyzzz, g_0_x_x_yyzzzz, g_0_x_x_yzzzz, g_0_x_x_yzzzzz, g_0_x_x_zzzzz, g_0_x_xy_xxxxx, g_0_x_xy_xxxxy, g_0_x_xy_xxxxz, g_0_x_xy_xxxyy, g_0_x_xy_xxxyz, g_0_x_xy_xxxzz, g_0_x_xy_xxyyy, g_0_x_xy_xxyyz, g_0_x_xy_xxyzz, g_0_x_xy_xxzzz, g_0_x_xy_xyyyy, g_0_x_xy_xyyyz, g_0_x_xy_xyyzz, g_0_x_xy_xyzzz, g_0_x_xy_xzzzz, g_0_x_xy_yyyyy, g_0_x_xy_yyyyz, g_0_x_xy_yyyzz, g_0_x_xy_yyzzz, g_0_x_xy_yzzzz, g_0_x_xy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xy_xxxxx[k] = -g_0_x_x_xxxxx[k] * ab_y + g_0_x_x_xxxxxy[k];

                g_0_x_xy_xxxxy[k] = -g_0_x_x_xxxxy[k] * ab_y + g_0_x_x_xxxxyy[k];

                g_0_x_xy_xxxxz[k] = -g_0_x_x_xxxxz[k] * ab_y + g_0_x_x_xxxxyz[k];

                g_0_x_xy_xxxyy[k] = -g_0_x_x_xxxyy[k] * ab_y + g_0_x_x_xxxyyy[k];

                g_0_x_xy_xxxyz[k] = -g_0_x_x_xxxyz[k] * ab_y + g_0_x_x_xxxyyz[k];

                g_0_x_xy_xxxzz[k] = -g_0_x_x_xxxzz[k] * ab_y + g_0_x_x_xxxyzz[k];

                g_0_x_xy_xxyyy[k] = -g_0_x_x_xxyyy[k] * ab_y + g_0_x_x_xxyyyy[k];

                g_0_x_xy_xxyyz[k] = -g_0_x_x_xxyyz[k] * ab_y + g_0_x_x_xxyyyz[k];

                g_0_x_xy_xxyzz[k] = -g_0_x_x_xxyzz[k] * ab_y + g_0_x_x_xxyyzz[k];

                g_0_x_xy_xxzzz[k] = -g_0_x_x_xxzzz[k] * ab_y + g_0_x_x_xxyzzz[k];

                g_0_x_xy_xyyyy[k] = -g_0_x_x_xyyyy[k] * ab_y + g_0_x_x_xyyyyy[k];

                g_0_x_xy_xyyyz[k] = -g_0_x_x_xyyyz[k] * ab_y + g_0_x_x_xyyyyz[k];

                g_0_x_xy_xyyzz[k] = -g_0_x_x_xyyzz[k] * ab_y + g_0_x_x_xyyyzz[k];

                g_0_x_xy_xyzzz[k] = -g_0_x_x_xyzzz[k] * ab_y + g_0_x_x_xyyzzz[k];

                g_0_x_xy_xzzzz[k] = -g_0_x_x_xzzzz[k] * ab_y + g_0_x_x_xyzzzz[k];

                g_0_x_xy_yyyyy[k] = -g_0_x_x_yyyyy[k] * ab_y + g_0_x_x_yyyyyy[k];

                g_0_x_xy_yyyyz[k] = -g_0_x_x_yyyyz[k] * ab_y + g_0_x_x_yyyyyz[k];

                g_0_x_xy_yyyzz[k] = -g_0_x_x_yyyzz[k] * ab_y + g_0_x_x_yyyyzz[k];

                g_0_x_xy_yyzzz[k] = -g_0_x_x_yyzzz[k] * ab_y + g_0_x_x_yyyzzz[k];

                g_0_x_xy_yzzzz[k] = -g_0_x_x_yzzzz[k] * ab_y + g_0_x_x_yyzzzz[k];

                g_0_x_xy_zzzzz[k] = -g_0_x_x_zzzzz[k] * ab_y + g_0_x_x_yzzzzz[k];
            }

            /// Set up 42-63 components of targeted buffer : cbuffer.data(

            auto g_0_x_xz_xxxxx = cbuffer.data(dh_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xz_xxxxy = cbuffer.data(dh_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xz_xxxxz = cbuffer.data(dh_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xz_xxxyy = cbuffer.data(dh_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xz_xxxyz = cbuffer.data(dh_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xz_xxxzz = cbuffer.data(dh_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xz_xxyyy = cbuffer.data(dh_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xz_xxyyz = cbuffer.data(dh_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xz_xxyzz = cbuffer.data(dh_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xz_xxzzz = cbuffer.data(dh_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xz_xyyyy = cbuffer.data(dh_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xz_xyyyz = cbuffer.data(dh_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xz_xyyzz = cbuffer.data(dh_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xz_xyzzz = cbuffer.data(dh_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xz_xzzzz = cbuffer.data(dh_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xz_yyyyy = cbuffer.data(dh_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xz_yyyyz = cbuffer.data(dh_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xz_yyyzz = cbuffer.data(dh_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_xz_yyzzz = cbuffer.data(dh_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xz_yzzzz = cbuffer.data(dh_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xz_zzzzz = cbuffer.data(dh_geom_01_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_x_xxxxx, g_0_x_x_xxxxxz, g_0_x_x_xxxxy, g_0_x_x_xxxxyz, g_0_x_x_xxxxz, g_0_x_x_xxxxzz, g_0_x_x_xxxyy, g_0_x_x_xxxyyz, g_0_x_x_xxxyz, g_0_x_x_xxxyzz, g_0_x_x_xxxzz, g_0_x_x_xxxzzz, g_0_x_x_xxyyy, g_0_x_x_xxyyyz, g_0_x_x_xxyyz, g_0_x_x_xxyyzz, g_0_x_x_xxyzz, g_0_x_x_xxyzzz, g_0_x_x_xxzzz, g_0_x_x_xxzzzz, g_0_x_x_xyyyy, g_0_x_x_xyyyyz, g_0_x_x_xyyyz, g_0_x_x_xyyyzz, g_0_x_x_xyyzz, g_0_x_x_xyyzzz, g_0_x_x_xyzzz, g_0_x_x_xyzzzz, g_0_x_x_xzzzz, g_0_x_x_xzzzzz, g_0_x_x_yyyyy, g_0_x_x_yyyyyz, g_0_x_x_yyyyz, g_0_x_x_yyyyzz, g_0_x_x_yyyzz, g_0_x_x_yyyzzz, g_0_x_x_yyzzz, g_0_x_x_yyzzzz, g_0_x_x_yzzzz, g_0_x_x_yzzzzz, g_0_x_x_zzzzz, g_0_x_x_zzzzzz, g_0_x_xz_xxxxx, g_0_x_xz_xxxxy, g_0_x_xz_xxxxz, g_0_x_xz_xxxyy, g_0_x_xz_xxxyz, g_0_x_xz_xxxzz, g_0_x_xz_xxyyy, g_0_x_xz_xxyyz, g_0_x_xz_xxyzz, g_0_x_xz_xxzzz, g_0_x_xz_xyyyy, g_0_x_xz_xyyyz, g_0_x_xz_xyyzz, g_0_x_xz_xyzzz, g_0_x_xz_xzzzz, g_0_x_xz_yyyyy, g_0_x_xz_yyyyz, g_0_x_xz_yyyzz, g_0_x_xz_yyzzz, g_0_x_xz_yzzzz, g_0_x_xz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xz_xxxxx[k] = -g_0_x_x_xxxxx[k] * ab_z + g_0_x_x_xxxxxz[k];

                g_0_x_xz_xxxxy[k] = -g_0_x_x_xxxxy[k] * ab_z + g_0_x_x_xxxxyz[k];

                g_0_x_xz_xxxxz[k] = -g_0_x_x_xxxxz[k] * ab_z + g_0_x_x_xxxxzz[k];

                g_0_x_xz_xxxyy[k] = -g_0_x_x_xxxyy[k] * ab_z + g_0_x_x_xxxyyz[k];

                g_0_x_xz_xxxyz[k] = -g_0_x_x_xxxyz[k] * ab_z + g_0_x_x_xxxyzz[k];

                g_0_x_xz_xxxzz[k] = -g_0_x_x_xxxzz[k] * ab_z + g_0_x_x_xxxzzz[k];

                g_0_x_xz_xxyyy[k] = -g_0_x_x_xxyyy[k] * ab_z + g_0_x_x_xxyyyz[k];

                g_0_x_xz_xxyyz[k] = -g_0_x_x_xxyyz[k] * ab_z + g_0_x_x_xxyyzz[k];

                g_0_x_xz_xxyzz[k] = -g_0_x_x_xxyzz[k] * ab_z + g_0_x_x_xxyzzz[k];

                g_0_x_xz_xxzzz[k] = -g_0_x_x_xxzzz[k] * ab_z + g_0_x_x_xxzzzz[k];

                g_0_x_xz_xyyyy[k] = -g_0_x_x_xyyyy[k] * ab_z + g_0_x_x_xyyyyz[k];

                g_0_x_xz_xyyyz[k] = -g_0_x_x_xyyyz[k] * ab_z + g_0_x_x_xyyyzz[k];

                g_0_x_xz_xyyzz[k] = -g_0_x_x_xyyzz[k] * ab_z + g_0_x_x_xyyzzz[k];

                g_0_x_xz_xyzzz[k] = -g_0_x_x_xyzzz[k] * ab_z + g_0_x_x_xyzzzz[k];

                g_0_x_xz_xzzzz[k] = -g_0_x_x_xzzzz[k] * ab_z + g_0_x_x_xzzzzz[k];

                g_0_x_xz_yyyyy[k] = -g_0_x_x_yyyyy[k] * ab_z + g_0_x_x_yyyyyz[k];

                g_0_x_xz_yyyyz[k] = -g_0_x_x_yyyyz[k] * ab_z + g_0_x_x_yyyyzz[k];

                g_0_x_xz_yyyzz[k] = -g_0_x_x_yyyzz[k] * ab_z + g_0_x_x_yyyzzz[k];

                g_0_x_xz_yyzzz[k] = -g_0_x_x_yyzzz[k] * ab_z + g_0_x_x_yyzzzz[k];

                g_0_x_xz_yzzzz[k] = -g_0_x_x_yzzzz[k] * ab_z + g_0_x_x_yzzzzz[k];

                g_0_x_xz_zzzzz[k] = -g_0_x_x_zzzzz[k] * ab_z + g_0_x_x_zzzzzz[k];
            }

            /// Set up 63-84 components of targeted buffer : cbuffer.data(

            auto g_0_x_yy_xxxxx = cbuffer.data(dh_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_yy_xxxxy = cbuffer.data(dh_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_yy_xxxxz = cbuffer.data(dh_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_yy_xxxyy = cbuffer.data(dh_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_yy_xxxyz = cbuffer.data(dh_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_yy_xxxzz = cbuffer.data(dh_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_yy_xxyyy = cbuffer.data(dh_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_yy_xxyyz = cbuffer.data(dh_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_yy_xxyzz = cbuffer.data(dh_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_yy_xxzzz = cbuffer.data(dh_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_yy_xyyyy = cbuffer.data(dh_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_yy_xyyyz = cbuffer.data(dh_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_yy_xyyzz = cbuffer.data(dh_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_yy_xyzzz = cbuffer.data(dh_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_yy_xzzzz = cbuffer.data(dh_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_yy_yyyyy = cbuffer.data(dh_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_yy_yyyyz = cbuffer.data(dh_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_yy_yyyzz = cbuffer.data(dh_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_yy_yyzzz = cbuffer.data(dh_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_yy_yzzzz = cbuffer.data(dh_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_yy_zzzzz = cbuffer.data(dh_geom_01_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_y_xxxxx, g_0_x_y_xxxxxy, g_0_x_y_xxxxy, g_0_x_y_xxxxyy, g_0_x_y_xxxxyz, g_0_x_y_xxxxz, g_0_x_y_xxxyy, g_0_x_y_xxxyyy, g_0_x_y_xxxyyz, g_0_x_y_xxxyz, g_0_x_y_xxxyzz, g_0_x_y_xxxzz, g_0_x_y_xxyyy, g_0_x_y_xxyyyy, g_0_x_y_xxyyyz, g_0_x_y_xxyyz, g_0_x_y_xxyyzz, g_0_x_y_xxyzz, g_0_x_y_xxyzzz, g_0_x_y_xxzzz, g_0_x_y_xyyyy, g_0_x_y_xyyyyy, g_0_x_y_xyyyyz, g_0_x_y_xyyyz, g_0_x_y_xyyyzz, g_0_x_y_xyyzz, g_0_x_y_xyyzzz, g_0_x_y_xyzzz, g_0_x_y_xyzzzz, g_0_x_y_xzzzz, g_0_x_y_yyyyy, g_0_x_y_yyyyyy, g_0_x_y_yyyyyz, g_0_x_y_yyyyz, g_0_x_y_yyyyzz, g_0_x_y_yyyzz, g_0_x_y_yyyzzz, g_0_x_y_yyzzz, g_0_x_y_yyzzzz, g_0_x_y_yzzzz, g_0_x_y_yzzzzz, g_0_x_y_zzzzz, g_0_x_yy_xxxxx, g_0_x_yy_xxxxy, g_0_x_yy_xxxxz, g_0_x_yy_xxxyy, g_0_x_yy_xxxyz, g_0_x_yy_xxxzz, g_0_x_yy_xxyyy, g_0_x_yy_xxyyz, g_0_x_yy_xxyzz, g_0_x_yy_xxzzz, g_0_x_yy_xyyyy, g_0_x_yy_xyyyz, g_0_x_yy_xyyzz, g_0_x_yy_xyzzz, g_0_x_yy_xzzzz, g_0_x_yy_yyyyy, g_0_x_yy_yyyyz, g_0_x_yy_yyyzz, g_0_x_yy_yyzzz, g_0_x_yy_yzzzz, g_0_x_yy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yy_xxxxx[k] = -g_0_x_y_xxxxx[k] * ab_y + g_0_x_y_xxxxxy[k];

                g_0_x_yy_xxxxy[k] = -g_0_x_y_xxxxy[k] * ab_y + g_0_x_y_xxxxyy[k];

                g_0_x_yy_xxxxz[k] = -g_0_x_y_xxxxz[k] * ab_y + g_0_x_y_xxxxyz[k];

                g_0_x_yy_xxxyy[k] = -g_0_x_y_xxxyy[k] * ab_y + g_0_x_y_xxxyyy[k];

                g_0_x_yy_xxxyz[k] = -g_0_x_y_xxxyz[k] * ab_y + g_0_x_y_xxxyyz[k];

                g_0_x_yy_xxxzz[k] = -g_0_x_y_xxxzz[k] * ab_y + g_0_x_y_xxxyzz[k];

                g_0_x_yy_xxyyy[k] = -g_0_x_y_xxyyy[k] * ab_y + g_0_x_y_xxyyyy[k];

                g_0_x_yy_xxyyz[k] = -g_0_x_y_xxyyz[k] * ab_y + g_0_x_y_xxyyyz[k];

                g_0_x_yy_xxyzz[k] = -g_0_x_y_xxyzz[k] * ab_y + g_0_x_y_xxyyzz[k];

                g_0_x_yy_xxzzz[k] = -g_0_x_y_xxzzz[k] * ab_y + g_0_x_y_xxyzzz[k];

                g_0_x_yy_xyyyy[k] = -g_0_x_y_xyyyy[k] * ab_y + g_0_x_y_xyyyyy[k];

                g_0_x_yy_xyyyz[k] = -g_0_x_y_xyyyz[k] * ab_y + g_0_x_y_xyyyyz[k];

                g_0_x_yy_xyyzz[k] = -g_0_x_y_xyyzz[k] * ab_y + g_0_x_y_xyyyzz[k];

                g_0_x_yy_xyzzz[k] = -g_0_x_y_xyzzz[k] * ab_y + g_0_x_y_xyyzzz[k];

                g_0_x_yy_xzzzz[k] = -g_0_x_y_xzzzz[k] * ab_y + g_0_x_y_xyzzzz[k];

                g_0_x_yy_yyyyy[k] = -g_0_x_y_yyyyy[k] * ab_y + g_0_x_y_yyyyyy[k];

                g_0_x_yy_yyyyz[k] = -g_0_x_y_yyyyz[k] * ab_y + g_0_x_y_yyyyyz[k];

                g_0_x_yy_yyyzz[k] = -g_0_x_y_yyyzz[k] * ab_y + g_0_x_y_yyyyzz[k];

                g_0_x_yy_yyzzz[k] = -g_0_x_y_yyzzz[k] * ab_y + g_0_x_y_yyyzzz[k];

                g_0_x_yy_yzzzz[k] = -g_0_x_y_yzzzz[k] * ab_y + g_0_x_y_yyzzzz[k];

                g_0_x_yy_zzzzz[k] = -g_0_x_y_zzzzz[k] * ab_y + g_0_x_y_yzzzzz[k];
            }

            /// Set up 84-105 components of targeted buffer : cbuffer.data(

            auto g_0_x_yz_xxxxx = cbuffer.data(dh_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_yz_xxxxy = cbuffer.data(dh_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_yz_xxxxz = cbuffer.data(dh_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_yz_xxxyy = cbuffer.data(dh_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_yz_xxxyz = cbuffer.data(dh_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_yz_xxxzz = cbuffer.data(dh_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_yz_xxyyy = cbuffer.data(dh_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_yz_xxyyz = cbuffer.data(dh_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_yz_xxyzz = cbuffer.data(dh_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_yz_xxzzz = cbuffer.data(dh_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_yz_xyyyy = cbuffer.data(dh_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_yz_xyyyz = cbuffer.data(dh_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_yz_xyyzz = cbuffer.data(dh_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_yz_xyzzz = cbuffer.data(dh_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_yz_xzzzz = cbuffer.data(dh_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_yz_yyyyy = cbuffer.data(dh_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_yz_yyyyz = cbuffer.data(dh_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_yz_yyyzz = cbuffer.data(dh_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_yz_yyzzz = cbuffer.data(dh_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_yz_yzzzz = cbuffer.data(dh_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_yz_zzzzz = cbuffer.data(dh_geom_01_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yz_xxxxx, g_0_x_yz_xxxxy, g_0_x_yz_xxxxz, g_0_x_yz_xxxyy, g_0_x_yz_xxxyz, g_0_x_yz_xxxzz, g_0_x_yz_xxyyy, g_0_x_yz_xxyyz, g_0_x_yz_xxyzz, g_0_x_yz_xxzzz, g_0_x_yz_xyyyy, g_0_x_yz_xyyyz, g_0_x_yz_xyyzz, g_0_x_yz_xyzzz, g_0_x_yz_xzzzz, g_0_x_yz_yyyyy, g_0_x_yz_yyyyz, g_0_x_yz_yyyzz, g_0_x_yz_yyzzz, g_0_x_yz_yzzzz, g_0_x_yz_zzzzz, g_0_x_z_xxxxx, g_0_x_z_xxxxxy, g_0_x_z_xxxxy, g_0_x_z_xxxxyy, g_0_x_z_xxxxyz, g_0_x_z_xxxxz, g_0_x_z_xxxyy, g_0_x_z_xxxyyy, g_0_x_z_xxxyyz, g_0_x_z_xxxyz, g_0_x_z_xxxyzz, g_0_x_z_xxxzz, g_0_x_z_xxyyy, g_0_x_z_xxyyyy, g_0_x_z_xxyyyz, g_0_x_z_xxyyz, g_0_x_z_xxyyzz, g_0_x_z_xxyzz, g_0_x_z_xxyzzz, g_0_x_z_xxzzz, g_0_x_z_xyyyy, g_0_x_z_xyyyyy, g_0_x_z_xyyyyz, g_0_x_z_xyyyz, g_0_x_z_xyyyzz, g_0_x_z_xyyzz, g_0_x_z_xyyzzz, g_0_x_z_xyzzz, g_0_x_z_xyzzzz, g_0_x_z_xzzzz, g_0_x_z_yyyyy, g_0_x_z_yyyyyy, g_0_x_z_yyyyyz, g_0_x_z_yyyyz, g_0_x_z_yyyyzz, g_0_x_z_yyyzz, g_0_x_z_yyyzzz, g_0_x_z_yyzzz, g_0_x_z_yyzzzz, g_0_x_z_yzzzz, g_0_x_z_yzzzzz, g_0_x_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yz_xxxxx[k] = -g_0_x_z_xxxxx[k] * ab_y + g_0_x_z_xxxxxy[k];

                g_0_x_yz_xxxxy[k] = -g_0_x_z_xxxxy[k] * ab_y + g_0_x_z_xxxxyy[k];

                g_0_x_yz_xxxxz[k] = -g_0_x_z_xxxxz[k] * ab_y + g_0_x_z_xxxxyz[k];

                g_0_x_yz_xxxyy[k] = -g_0_x_z_xxxyy[k] * ab_y + g_0_x_z_xxxyyy[k];

                g_0_x_yz_xxxyz[k] = -g_0_x_z_xxxyz[k] * ab_y + g_0_x_z_xxxyyz[k];

                g_0_x_yz_xxxzz[k] = -g_0_x_z_xxxzz[k] * ab_y + g_0_x_z_xxxyzz[k];

                g_0_x_yz_xxyyy[k] = -g_0_x_z_xxyyy[k] * ab_y + g_0_x_z_xxyyyy[k];

                g_0_x_yz_xxyyz[k] = -g_0_x_z_xxyyz[k] * ab_y + g_0_x_z_xxyyyz[k];

                g_0_x_yz_xxyzz[k] = -g_0_x_z_xxyzz[k] * ab_y + g_0_x_z_xxyyzz[k];

                g_0_x_yz_xxzzz[k] = -g_0_x_z_xxzzz[k] * ab_y + g_0_x_z_xxyzzz[k];

                g_0_x_yz_xyyyy[k] = -g_0_x_z_xyyyy[k] * ab_y + g_0_x_z_xyyyyy[k];

                g_0_x_yz_xyyyz[k] = -g_0_x_z_xyyyz[k] * ab_y + g_0_x_z_xyyyyz[k];

                g_0_x_yz_xyyzz[k] = -g_0_x_z_xyyzz[k] * ab_y + g_0_x_z_xyyyzz[k];

                g_0_x_yz_xyzzz[k] = -g_0_x_z_xyzzz[k] * ab_y + g_0_x_z_xyyzzz[k];

                g_0_x_yz_xzzzz[k] = -g_0_x_z_xzzzz[k] * ab_y + g_0_x_z_xyzzzz[k];

                g_0_x_yz_yyyyy[k] = -g_0_x_z_yyyyy[k] * ab_y + g_0_x_z_yyyyyy[k];

                g_0_x_yz_yyyyz[k] = -g_0_x_z_yyyyz[k] * ab_y + g_0_x_z_yyyyyz[k];

                g_0_x_yz_yyyzz[k] = -g_0_x_z_yyyzz[k] * ab_y + g_0_x_z_yyyyzz[k];

                g_0_x_yz_yyzzz[k] = -g_0_x_z_yyzzz[k] * ab_y + g_0_x_z_yyyzzz[k];

                g_0_x_yz_yzzzz[k] = -g_0_x_z_yzzzz[k] * ab_y + g_0_x_z_yyzzzz[k];

                g_0_x_yz_zzzzz[k] = -g_0_x_z_zzzzz[k] * ab_y + g_0_x_z_yzzzzz[k];
            }

            /// Set up 105-126 components of targeted buffer : cbuffer.data(

            auto g_0_x_zz_xxxxx = cbuffer.data(dh_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_zz_xxxxy = cbuffer.data(dh_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_zz_xxxxz = cbuffer.data(dh_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_zz_xxxyy = cbuffer.data(dh_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_zz_xxxyz = cbuffer.data(dh_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_zz_xxxzz = cbuffer.data(dh_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_zz_xxyyy = cbuffer.data(dh_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_zz_xxyyz = cbuffer.data(dh_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_zz_xxyzz = cbuffer.data(dh_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_zz_xxzzz = cbuffer.data(dh_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_zz_xyyyy = cbuffer.data(dh_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_zz_xyyyz = cbuffer.data(dh_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_zz_xyyzz = cbuffer.data(dh_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_zz_xyzzz = cbuffer.data(dh_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_zz_xzzzz = cbuffer.data(dh_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_x_zz_yyyyy = cbuffer.data(dh_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_zz_yyyyz = cbuffer.data(dh_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_zz_yyyzz = cbuffer.data(dh_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_zz_yyzzz = cbuffer.data(dh_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_zz_yzzzz = cbuffer.data(dh_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_zz_zzzzz = cbuffer.data(dh_geom_01_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_z_xxxxx, g_0_x_z_xxxxxz, g_0_x_z_xxxxy, g_0_x_z_xxxxyz, g_0_x_z_xxxxz, g_0_x_z_xxxxzz, g_0_x_z_xxxyy, g_0_x_z_xxxyyz, g_0_x_z_xxxyz, g_0_x_z_xxxyzz, g_0_x_z_xxxzz, g_0_x_z_xxxzzz, g_0_x_z_xxyyy, g_0_x_z_xxyyyz, g_0_x_z_xxyyz, g_0_x_z_xxyyzz, g_0_x_z_xxyzz, g_0_x_z_xxyzzz, g_0_x_z_xxzzz, g_0_x_z_xxzzzz, g_0_x_z_xyyyy, g_0_x_z_xyyyyz, g_0_x_z_xyyyz, g_0_x_z_xyyyzz, g_0_x_z_xyyzz, g_0_x_z_xyyzzz, g_0_x_z_xyzzz, g_0_x_z_xyzzzz, g_0_x_z_xzzzz, g_0_x_z_xzzzzz, g_0_x_z_yyyyy, g_0_x_z_yyyyyz, g_0_x_z_yyyyz, g_0_x_z_yyyyzz, g_0_x_z_yyyzz, g_0_x_z_yyyzzz, g_0_x_z_yyzzz, g_0_x_z_yyzzzz, g_0_x_z_yzzzz, g_0_x_z_yzzzzz, g_0_x_z_zzzzz, g_0_x_z_zzzzzz, g_0_x_zz_xxxxx, g_0_x_zz_xxxxy, g_0_x_zz_xxxxz, g_0_x_zz_xxxyy, g_0_x_zz_xxxyz, g_0_x_zz_xxxzz, g_0_x_zz_xxyyy, g_0_x_zz_xxyyz, g_0_x_zz_xxyzz, g_0_x_zz_xxzzz, g_0_x_zz_xyyyy, g_0_x_zz_xyyyz, g_0_x_zz_xyyzz, g_0_x_zz_xyzzz, g_0_x_zz_xzzzz, g_0_x_zz_yyyyy, g_0_x_zz_yyyyz, g_0_x_zz_yyyzz, g_0_x_zz_yyzzz, g_0_x_zz_yzzzz, g_0_x_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_zz_xxxxx[k] = -g_0_x_z_xxxxx[k] * ab_z + g_0_x_z_xxxxxz[k];

                g_0_x_zz_xxxxy[k] = -g_0_x_z_xxxxy[k] * ab_z + g_0_x_z_xxxxyz[k];

                g_0_x_zz_xxxxz[k] = -g_0_x_z_xxxxz[k] * ab_z + g_0_x_z_xxxxzz[k];

                g_0_x_zz_xxxyy[k] = -g_0_x_z_xxxyy[k] * ab_z + g_0_x_z_xxxyyz[k];

                g_0_x_zz_xxxyz[k] = -g_0_x_z_xxxyz[k] * ab_z + g_0_x_z_xxxyzz[k];

                g_0_x_zz_xxxzz[k] = -g_0_x_z_xxxzz[k] * ab_z + g_0_x_z_xxxzzz[k];

                g_0_x_zz_xxyyy[k] = -g_0_x_z_xxyyy[k] * ab_z + g_0_x_z_xxyyyz[k];

                g_0_x_zz_xxyyz[k] = -g_0_x_z_xxyyz[k] * ab_z + g_0_x_z_xxyyzz[k];

                g_0_x_zz_xxyzz[k] = -g_0_x_z_xxyzz[k] * ab_z + g_0_x_z_xxyzzz[k];

                g_0_x_zz_xxzzz[k] = -g_0_x_z_xxzzz[k] * ab_z + g_0_x_z_xxzzzz[k];

                g_0_x_zz_xyyyy[k] = -g_0_x_z_xyyyy[k] * ab_z + g_0_x_z_xyyyyz[k];

                g_0_x_zz_xyyyz[k] = -g_0_x_z_xyyyz[k] * ab_z + g_0_x_z_xyyyzz[k];

                g_0_x_zz_xyyzz[k] = -g_0_x_z_xyyzz[k] * ab_z + g_0_x_z_xyyzzz[k];

                g_0_x_zz_xyzzz[k] = -g_0_x_z_xyzzz[k] * ab_z + g_0_x_z_xyzzzz[k];

                g_0_x_zz_xzzzz[k] = -g_0_x_z_xzzzz[k] * ab_z + g_0_x_z_xzzzzz[k];

                g_0_x_zz_yyyyy[k] = -g_0_x_z_yyyyy[k] * ab_z + g_0_x_z_yyyyyz[k];

                g_0_x_zz_yyyyz[k] = -g_0_x_z_yyyyz[k] * ab_z + g_0_x_z_yyyyzz[k];

                g_0_x_zz_yyyzz[k] = -g_0_x_z_yyyzz[k] * ab_z + g_0_x_z_yyyzzz[k];

                g_0_x_zz_yyzzz[k] = -g_0_x_z_yyzzz[k] * ab_z + g_0_x_z_yyzzzz[k];

                g_0_x_zz_yzzzz[k] = -g_0_x_z_yzzzz[k] * ab_z + g_0_x_z_yzzzzz[k];

                g_0_x_zz_zzzzz[k] = -g_0_x_z_zzzzz[k] * ab_z + g_0_x_z_zzzzzz[k];
            }

            /// Set up 126-147 components of targeted buffer : cbuffer.data(

            auto g_0_y_xx_xxxxx = cbuffer.data(dh_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_y_xx_xxxxy = cbuffer.data(dh_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_y_xx_xxxxz = cbuffer.data(dh_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_y_xx_xxxyy = cbuffer.data(dh_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_y_xx_xxxyz = cbuffer.data(dh_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_y_xx_xxxzz = cbuffer.data(dh_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_y_xx_xxyyy = cbuffer.data(dh_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_y_xx_xxyyz = cbuffer.data(dh_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_y_xx_xxyzz = cbuffer.data(dh_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_y_xx_xxzzz = cbuffer.data(dh_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_y_xx_xyyyy = cbuffer.data(dh_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_y_xx_xyyyz = cbuffer.data(dh_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_y_xx_xyyzz = cbuffer.data(dh_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_y_xx_xyzzz = cbuffer.data(dh_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_y_xx_xzzzz = cbuffer.data(dh_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_y_xx_yyyyy = cbuffer.data(dh_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_y_xx_yyyyz = cbuffer.data(dh_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_y_xx_yyyzz = cbuffer.data(dh_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_y_xx_yyzzz = cbuffer.data(dh_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_y_xx_yzzzz = cbuffer.data(dh_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_y_xx_zzzzz = cbuffer.data(dh_geom_01_off + 146 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_x_xxxxx, g_0_y_x_xxxxxx, g_0_y_x_xxxxxy, g_0_y_x_xxxxxz, g_0_y_x_xxxxy, g_0_y_x_xxxxyy, g_0_y_x_xxxxyz, g_0_y_x_xxxxz, g_0_y_x_xxxxzz, g_0_y_x_xxxyy, g_0_y_x_xxxyyy, g_0_y_x_xxxyyz, g_0_y_x_xxxyz, g_0_y_x_xxxyzz, g_0_y_x_xxxzz, g_0_y_x_xxxzzz, g_0_y_x_xxyyy, g_0_y_x_xxyyyy, g_0_y_x_xxyyyz, g_0_y_x_xxyyz, g_0_y_x_xxyyzz, g_0_y_x_xxyzz, g_0_y_x_xxyzzz, g_0_y_x_xxzzz, g_0_y_x_xxzzzz, g_0_y_x_xyyyy, g_0_y_x_xyyyyy, g_0_y_x_xyyyyz, g_0_y_x_xyyyz, g_0_y_x_xyyyzz, g_0_y_x_xyyzz, g_0_y_x_xyyzzz, g_0_y_x_xyzzz, g_0_y_x_xyzzzz, g_0_y_x_xzzzz, g_0_y_x_xzzzzz, g_0_y_x_yyyyy, g_0_y_x_yyyyz, g_0_y_x_yyyzz, g_0_y_x_yyzzz, g_0_y_x_yzzzz, g_0_y_x_zzzzz, g_0_y_xx_xxxxx, g_0_y_xx_xxxxy, g_0_y_xx_xxxxz, g_0_y_xx_xxxyy, g_0_y_xx_xxxyz, g_0_y_xx_xxxzz, g_0_y_xx_xxyyy, g_0_y_xx_xxyyz, g_0_y_xx_xxyzz, g_0_y_xx_xxzzz, g_0_y_xx_xyyyy, g_0_y_xx_xyyyz, g_0_y_xx_xyyzz, g_0_y_xx_xyzzz, g_0_y_xx_xzzzz, g_0_y_xx_yyyyy, g_0_y_xx_yyyyz, g_0_y_xx_yyyzz, g_0_y_xx_yyzzz, g_0_y_xx_yzzzz, g_0_y_xx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xx_xxxxx[k] = -g_0_y_x_xxxxx[k] * ab_x + g_0_y_x_xxxxxx[k];

                g_0_y_xx_xxxxy[k] = -g_0_y_x_xxxxy[k] * ab_x + g_0_y_x_xxxxxy[k];

                g_0_y_xx_xxxxz[k] = -g_0_y_x_xxxxz[k] * ab_x + g_0_y_x_xxxxxz[k];

                g_0_y_xx_xxxyy[k] = -g_0_y_x_xxxyy[k] * ab_x + g_0_y_x_xxxxyy[k];

                g_0_y_xx_xxxyz[k] = -g_0_y_x_xxxyz[k] * ab_x + g_0_y_x_xxxxyz[k];

                g_0_y_xx_xxxzz[k] = -g_0_y_x_xxxzz[k] * ab_x + g_0_y_x_xxxxzz[k];

                g_0_y_xx_xxyyy[k] = -g_0_y_x_xxyyy[k] * ab_x + g_0_y_x_xxxyyy[k];

                g_0_y_xx_xxyyz[k] = -g_0_y_x_xxyyz[k] * ab_x + g_0_y_x_xxxyyz[k];

                g_0_y_xx_xxyzz[k] = -g_0_y_x_xxyzz[k] * ab_x + g_0_y_x_xxxyzz[k];

                g_0_y_xx_xxzzz[k] = -g_0_y_x_xxzzz[k] * ab_x + g_0_y_x_xxxzzz[k];

                g_0_y_xx_xyyyy[k] = -g_0_y_x_xyyyy[k] * ab_x + g_0_y_x_xxyyyy[k];

                g_0_y_xx_xyyyz[k] = -g_0_y_x_xyyyz[k] * ab_x + g_0_y_x_xxyyyz[k];

                g_0_y_xx_xyyzz[k] = -g_0_y_x_xyyzz[k] * ab_x + g_0_y_x_xxyyzz[k];

                g_0_y_xx_xyzzz[k] = -g_0_y_x_xyzzz[k] * ab_x + g_0_y_x_xxyzzz[k];

                g_0_y_xx_xzzzz[k] = -g_0_y_x_xzzzz[k] * ab_x + g_0_y_x_xxzzzz[k];

                g_0_y_xx_yyyyy[k] = -g_0_y_x_yyyyy[k] * ab_x + g_0_y_x_xyyyyy[k];

                g_0_y_xx_yyyyz[k] = -g_0_y_x_yyyyz[k] * ab_x + g_0_y_x_xyyyyz[k];

                g_0_y_xx_yyyzz[k] = -g_0_y_x_yyyzz[k] * ab_x + g_0_y_x_xyyyzz[k];

                g_0_y_xx_yyzzz[k] = -g_0_y_x_yyzzz[k] * ab_x + g_0_y_x_xyyzzz[k];

                g_0_y_xx_yzzzz[k] = -g_0_y_x_yzzzz[k] * ab_x + g_0_y_x_xyzzzz[k];

                g_0_y_xx_zzzzz[k] = -g_0_y_x_zzzzz[k] * ab_x + g_0_y_x_xzzzzz[k];
            }

            /// Set up 147-168 components of targeted buffer : cbuffer.data(

            auto g_0_y_xy_xxxxx = cbuffer.data(dh_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_y_xy_xxxxy = cbuffer.data(dh_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_y_xy_xxxxz = cbuffer.data(dh_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_y_xy_xxxyy = cbuffer.data(dh_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_y_xy_xxxyz = cbuffer.data(dh_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_y_xy_xxxzz = cbuffer.data(dh_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_y_xy_xxyyy = cbuffer.data(dh_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_y_xy_xxyyz = cbuffer.data(dh_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_y_xy_xxyzz = cbuffer.data(dh_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_y_xy_xxzzz = cbuffer.data(dh_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_y_xy_xyyyy = cbuffer.data(dh_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_y_xy_xyyyz = cbuffer.data(dh_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_y_xy_xyyzz = cbuffer.data(dh_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_y_xy_xyzzz = cbuffer.data(dh_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_y_xy_xzzzz = cbuffer.data(dh_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_y_xy_yyyyy = cbuffer.data(dh_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_y_xy_yyyyz = cbuffer.data(dh_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_y_xy_yyyzz = cbuffer.data(dh_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_y_xy_yyzzz = cbuffer.data(dh_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_y_xy_yzzzz = cbuffer.data(dh_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_y_xy_zzzzz = cbuffer.data(dh_geom_01_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xy_xxxxx, g_0_y_xy_xxxxy, g_0_y_xy_xxxxz, g_0_y_xy_xxxyy, g_0_y_xy_xxxyz, g_0_y_xy_xxxzz, g_0_y_xy_xxyyy, g_0_y_xy_xxyyz, g_0_y_xy_xxyzz, g_0_y_xy_xxzzz, g_0_y_xy_xyyyy, g_0_y_xy_xyyyz, g_0_y_xy_xyyzz, g_0_y_xy_xyzzz, g_0_y_xy_xzzzz, g_0_y_xy_yyyyy, g_0_y_xy_yyyyz, g_0_y_xy_yyyzz, g_0_y_xy_yyzzz, g_0_y_xy_yzzzz, g_0_y_xy_zzzzz, g_0_y_y_xxxxx, g_0_y_y_xxxxxx, g_0_y_y_xxxxxy, g_0_y_y_xxxxxz, g_0_y_y_xxxxy, g_0_y_y_xxxxyy, g_0_y_y_xxxxyz, g_0_y_y_xxxxz, g_0_y_y_xxxxzz, g_0_y_y_xxxyy, g_0_y_y_xxxyyy, g_0_y_y_xxxyyz, g_0_y_y_xxxyz, g_0_y_y_xxxyzz, g_0_y_y_xxxzz, g_0_y_y_xxxzzz, g_0_y_y_xxyyy, g_0_y_y_xxyyyy, g_0_y_y_xxyyyz, g_0_y_y_xxyyz, g_0_y_y_xxyyzz, g_0_y_y_xxyzz, g_0_y_y_xxyzzz, g_0_y_y_xxzzz, g_0_y_y_xxzzzz, g_0_y_y_xyyyy, g_0_y_y_xyyyyy, g_0_y_y_xyyyyz, g_0_y_y_xyyyz, g_0_y_y_xyyyzz, g_0_y_y_xyyzz, g_0_y_y_xyyzzz, g_0_y_y_xyzzz, g_0_y_y_xyzzzz, g_0_y_y_xzzzz, g_0_y_y_xzzzzz, g_0_y_y_yyyyy, g_0_y_y_yyyyz, g_0_y_y_yyyzz, g_0_y_y_yyzzz, g_0_y_y_yzzzz, g_0_y_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xy_xxxxx[k] = -g_0_y_y_xxxxx[k] * ab_x + g_0_y_y_xxxxxx[k];

                g_0_y_xy_xxxxy[k] = -g_0_y_y_xxxxy[k] * ab_x + g_0_y_y_xxxxxy[k];

                g_0_y_xy_xxxxz[k] = -g_0_y_y_xxxxz[k] * ab_x + g_0_y_y_xxxxxz[k];

                g_0_y_xy_xxxyy[k] = -g_0_y_y_xxxyy[k] * ab_x + g_0_y_y_xxxxyy[k];

                g_0_y_xy_xxxyz[k] = -g_0_y_y_xxxyz[k] * ab_x + g_0_y_y_xxxxyz[k];

                g_0_y_xy_xxxzz[k] = -g_0_y_y_xxxzz[k] * ab_x + g_0_y_y_xxxxzz[k];

                g_0_y_xy_xxyyy[k] = -g_0_y_y_xxyyy[k] * ab_x + g_0_y_y_xxxyyy[k];

                g_0_y_xy_xxyyz[k] = -g_0_y_y_xxyyz[k] * ab_x + g_0_y_y_xxxyyz[k];

                g_0_y_xy_xxyzz[k] = -g_0_y_y_xxyzz[k] * ab_x + g_0_y_y_xxxyzz[k];

                g_0_y_xy_xxzzz[k] = -g_0_y_y_xxzzz[k] * ab_x + g_0_y_y_xxxzzz[k];

                g_0_y_xy_xyyyy[k] = -g_0_y_y_xyyyy[k] * ab_x + g_0_y_y_xxyyyy[k];

                g_0_y_xy_xyyyz[k] = -g_0_y_y_xyyyz[k] * ab_x + g_0_y_y_xxyyyz[k];

                g_0_y_xy_xyyzz[k] = -g_0_y_y_xyyzz[k] * ab_x + g_0_y_y_xxyyzz[k];

                g_0_y_xy_xyzzz[k] = -g_0_y_y_xyzzz[k] * ab_x + g_0_y_y_xxyzzz[k];

                g_0_y_xy_xzzzz[k] = -g_0_y_y_xzzzz[k] * ab_x + g_0_y_y_xxzzzz[k];

                g_0_y_xy_yyyyy[k] = -g_0_y_y_yyyyy[k] * ab_x + g_0_y_y_xyyyyy[k];

                g_0_y_xy_yyyyz[k] = -g_0_y_y_yyyyz[k] * ab_x + g_0_y_y_xyyyyz[k];

                g_0_y_xy_yyyzz[k] = -g_0_y_y_yyyzz[k] * ab_x + g_0_y_y_xyyyzz[k];

                g_0_y_xy_yyzzz[k] = -g_0_y_y_yyzzz[k] * ab_x + g_0_y_y_xyyzzz[k];

                g_0_y_xy_yzzzz[k] = -g_0_y_y_yzzzz[k] * ab_x + g_0_y_y_xyzzzz[k];

                g_0_y_xy_zzzzz[k] = -g_0_y_y_zzzzz[k] * ab_x + g_0_y_y_xzzzzz[k];
            }

            /// Set up 168-189 components of targeted buffer : cbuffer.data(

            auto g_0_y_xz_xxxxx = cbuffer.data(dh_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_y_xz_xxxxy = cbuffer.data(dh_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_y_xz_xxxxz = cbuffer.data(dh_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_y_xz_xxxyy = cbuffer.data(dh_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_y_xz_xxxyz = cbuffer.data(dh_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_y_xz_xxxzz = cbuffer.data(dh_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_y_xz_xxyyy = cbuffer.data(dh_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_y_xz_xxyyz = cbuffer.data(dh_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_y_xz_xxyzz = cbuffer.data(dh_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_y_xz_xxzzz = cbuffer.data(dh_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_y_xz_xyyyy = cbuffer.data(dh_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_y_xz_xyyyz = cbuffer.data(dh_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_y_xz_xyyzz = cbuffer.data(dh_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_y_xz_xyzzz = cbuffer.data(dh_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_y_xz_xzzzz = cbuffer.data(dh_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_y_xz_yyyyy = cbuffer.data(dh_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_y_xz_yyyyz = cbuffer.data(dh_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_y_xz_yyyzz = cbuffer.data(dh_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_y_xz_yyzzz = cbuffer.data(dh_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_y_xz_yzzzz = cbuffer.data(dh_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_y_xz_zzzzz = cbuffer.data(dh_geom_01_off + 188 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xz_xxxxx, g_0_y_xz_xxxxy, g_0_y_xz_xxxxz, g_0_y_xz_xxxyy, g_0_y_xz_xxxyz, g_0_y_xz_xxxzz, g_0_y_xz_xxyyy, g_0_y_xz_xxyyz, g_0_y_xz_xxyzz, g_0_y_xz_xxzzz, g_0_y_xz_xyyyy, g_0_y_xz_xyyyz, g_0_y_xz_xyyzz, g_0_y_xz_xyzzz, g_0_y_xz_xzzzz, g_0_y_xz_yyyyy, g_0_y_xz_yyyyz, g_0_y_xz_yyyzz, g_0_y_xz_yyzzz, g_0_y_xz_yzzzz, g_0_y_xz_zzzzz, g_0_y_z_xxxxx, g_0_y_z_xxxxxx, g_0_y_z_xxxxxy, g_0_y_z_xxxxxz, g_0_y_z_xxxxy, g_0_y_z_xxxxyy, g_0_y_z_xxxxyz, g_0_y_z_xxxxz, g_0_y_z_xxxxzz, g_0_y_z_xxxyy, g_0_y_z_xxxyyy, g_0_y_z_xxxyyz, g_0_y_z_xxxyz, g_0_y_z_xxxyzz, g_0_y_z_xxxzz, g_0_y_z_xxxzzz, g_0_y_z_xxyyy, g_0_y_z_xxyyyy, g_0_y_z_xxyyyz, g_0_y_z_xxyyz, g_0_y_z_xxyyzz, g_0_y_z_xxyzz, g_0_y_z_xxyzzz, g_0_y_z_xxzzz, g_0_y_z_xxzzzz, g_0_y_z_xyyyy, g_0_y_z_xyyyyy, g_0_y_z_xyyyyz, g_0_y_z_xyyyz, g_0_y_z_xyyyzz, g_0_y_z_xyyzz, g_0_y_z_xyyzzz, g_0_y_z_xyzzz, g_0_y_z_xyzzzz, g_0_y_z_xzzzz, g_0_y_z_xzzzzz, g_0_y_z_yyyyy, g_0_y_z_yyyyz, g_0_y_z_yyyzz, g_0_y_z_yyzzz, g_0_y_z_yzzzz, g_0_y_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xz_xxxxx[k] = -g_0_y_z_xxxxx[k] * ab_x + g_0_y_z_xxxxxx[k];

                g_0_y_xz_xxxxy[k] = -g_0_y_z_xxxxy[k] * ab_x + g_0_y_z_xxxxxy[k];

                g_0_y_xz_xxxxz[k] = -g_0_y_z_xxxxz[k] * ab_x + g_0_y_z_xxxxxz[k];

                g_0_y_xz_xxxyy[k] = -g_0_y_z_xxxyy[k] * ab_x + g_0_y_z_xxxxyy[k];

                g_0_y_xz_xxxyz[k] = -g_0_y_z_xxxyz[k] * ab_x + g_0_y_z_xxxxyz[k];

                g_0_y_xz_xxxzz[k] = -g_0_y_z_xxxzz[k] * ab_x + g_0_y_z_xxxxzz[k];

                g_0_y_xz_xxyyy[k] = -g_0_y_z_xxyyy[k] * ab_x + g_0_y_z_xxxyyy[k];

                g_0_y_xz_xxyyz[k] = -g_0_y_z_xxyyz[k] * ab_x + g_0_y_z_xxxyyz[k];

                g_0_y_xz_xxyzz[k] = -g_0_y_z_xxyzz[k] * ab_x + g_0_y_z_xxxyzz[k];

                g_0_y_xz_xxzzz[k] = -g_0_y_z_xxzzz[k] * ab_x + g_0_y_z_xxxzzz[k];

                g_0_y_xz_xyyyy[k] = -g_0_y_z_xyyyy[k] * ab_x + g_0_y_z_xxyyyy[k];

                g_0_y_xz_xyyyz[k] = -g_0_y_z_xyyyz[k] * ab_x + g_0_y_z_xxyyyz[k];

                g_0_y_xz_xyyzz[k] = -g_0_y_z_xyyzz[k] * ab_x + g_0_y_z_xxyyzz[k];

                g_0_y_xz_xyzzz[k] = -g_0_y_z_xyzzz[k] * ab_x + g_0_y_z_xxyzzz[k];

                g_0_y_xz_xzzzz[k] = -g_0_y_z_xzzzz[k] * ab_x + g_0_y_z_xxzzzz[k];

                g_0_y_xz_yyyyy[k] = -g_0_y_z_yyyyy[k] * ab_x + g_0_y_z_xyyyyy[k];

                g_0_y_xz_yyyyz[k] = -g_0_y_z_yyyyz[k] * ab_x + g_0_y_z_xyyyyz[k];

                g_0_y_xz_yyyzz[k] = -g_0_y_z_yyyzz[k] * ab_x + g_0_y_z_xyyyzz[k];

                g_0_y_xz_yyzzz[k] = -g_0_y_z_yyzzz[k] * ab_x + g_0_y_z_xyyzzz[k];

                g_0_y_xz_yzzzz[k] = -g_0_y_z_yzzzz[k] * ab_x + g_0_y_z_xyzzzz[k];

                g_0_y_xz_zzzzz[k] = -g_0_y_z_zzzzz[k] * ab_x + g_0_y_z_xzzzzz[k];
            }

            /// Set up 189-210 components of targeted buffer : cbuffer.data(

            auto g_0_y_yy_xxxxx = cbuffer.data(dh_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_y_yy_xxxxy = cbuffer.data(dh_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_y_yy_xxxxz = cbuffer.data(dh_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_y_yy_xxxyy = cbuffer.data(dh_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_y_yy_xxxyz = cbuffer.data(dh_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_y_yy_xxxzz = cbuffer.data(dh_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_y_yy_xxyyy = cbuffer.data(dh_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_y_yy_xxyyz = cbuffer.data(dh_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_y_yy_xxyzz = cbuffer.data(dh_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_y_yy_xxzzz = cbuffer.data(dh_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_y_yy_xyyyy = cbuffer.data(dh_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_y_yy_xyyyz = cbuffer.data(dh_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_y_yy_xyyzz = cbuffer.data(dh_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_y_yy_xyzzz = cbuffer.data(dh_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_y_yy_xzzzz = cbuffer.data(dh_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_y_yy_yyyyy = cbuffer.data(dh_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_y_yy_yyyyz = cbuffer.data(dh_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_y_yy_yyyzz = cbuffer.data(dh_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_y_yy_yyzzz = cbuffer.data(dh_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_y_yy_yzzzz = cbuffer.data(dh_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_y_yy_zzzzz = cbuffer.data(dh_geom_01_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_y_xxxxx, g_0_y_y_xxxxxy, g_0_y_y_xxxxy, g_0_y_y_xxxxyy, g_0_y_y_xxxxyz, g_0_y_y_xxxxz, g_0_y_y_xxxyy, g_0_y_y_xxxyyy, g_0_y_y_xxxyyz, g_0_y_y_xxxyz, g_0_y_y_xxxyzz, g_0_y_y_xxxzz, g_0_y_y_xxyyy, g_0_y_y_xxyyyy, g_0_y_y_xxyyyz, g_0_y_y_xxyyz, g_0_y_y_xxyyzz, g_0_y_y_xxyzz, g_0_y_y_xxyzzz, g_0_y_y_xxzzz, g_0_y_y_xyyyy, g_0_y_y_xyyyyy, g_0_y_y_xyyyyz, g_0_y_y_xyyyz, g_0_y_y_xyyyzz, g_0_y_y_xyyzz, g_0_y_y_xyyzzz, g_0_y_y_xyzzz, g_0_y_y_xyzzzz, g_0_y_y_xzzzz, g_0_y_y_yyyyy, g_0_y_y_yyyyyy, g_0_y_y_yyyyyz, g_0_y_y_yyyyz, g_0_y_y_yyyyzz, g_0_y_y_yyyzz, g_0_y_y_yyyzzz, g_0_y_y_yyzzz, g_0_y_y_yyzzzz, g_0_y_y_yzzzz, g_0_y_y_yzzzzz, g_0_y_y_zzzzz, g_0_y_yy_xxxxx, g_0_y_yy_xxxxy, g_0_y_yy_xxxxz, g_0_y_yy_xxxyy, g_0_y_yy_xxxyz, g_0_y_yy_xxxzz, g_0_y_yy_xxyyy, g_0_y_yy_xxyyz, g_0_y_yy_xxyzz, g_0_y_yy_xxzzz, g_0_y_yy_xyyyy, g_0_y_yy_xyyyz, g_0_y_yy_xyyzz, g_0_y_yy_xyzzz, g_0_y_yy_xzzzz, g_0_y_yy_yyyyy, g_0_y_yy_yyyyz, g_0_y_yy_yyyzz, g_0_y_yy_yyzzz, g_0_y_yy_yzzzz, g_0_y_yy_zzzzz, g_y_xxxxx, g_y_xxxxy, g_y_xxxxz, g_y_xxxyy, g_y_xxxyz, g_y_xxxzz, g_y_xxyyy, g_y_xxyyz, g_y_xxyzz, g_y_xxzzz, g_y_xyyyy, g_y_xyyyz, g_y_xyyzz, g_y_xyzzz, g_y_xzzzz, g_y_yyyyy, g_y_yyyyz, g_y_yyyzz, g_y_yyzzz, g_y_yzzzz, g_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yy_xxxxx[k] = g_y_xxxxx[k] - g_0_y_y_xxxxx[k] * ab_y + g_0_y_y_xxxxxy[k];

                g_0_y_yy_xxxxy[k] = g_y_xxxxy[k] - g_0_y_y_xxxxy[k] * ab_y + g_0_y_y_xxxxyy[k];

                g_0_y_yy_xxxxz[k] = g_y_xxxxz[k] - g_0_y_y_xxxxz[k] * ab_y + g_0_y_y_xxxxyz[k];

                g_0_y_yy_xxxyy[k] = g_y_xxxyy[k] - g_0_y_y_xxxyy[k] * ab_y + g_0_y_y_xxxyyy[k];

                g_0_y_yy_xxxyz[k] = g_y_xxxyz[k] - g_0_y_y_xxxyz[k] * ab_y + g_0_y_y_xxxyyz[k];

                g_0_y_yy_xxxzz[k] = g_y_xxxzz[k] - g_0_y_y_xxxzz[k] * ab_y + g_0_y_y_xxxyzz[k];

                g_0_y_yy_xxyyy[k] = g_y_xxyyy[k] - g_0_y_y_xxyyy[k] * ab_y + g_0_y_y_xxyyyy[k];

                g_0_y_yy_xxyyz[k] = g_y_xxyyz[k] - g_0_y_y_xxyyz[k] * ab_y + g_0_y_y_xxyyyz[k];

                g_0_y_yy_xxyzz[k] = g_y_xxyzz[k] - g_0_y_y_xxyzz[k] * ab_y + g_0_y_y_xxyyzz[k];

                g_0_y_yy_xxzzz[k] = g_y_xxzzz[k] - g_0_y_y_xxzzz[k] * ab_y + g_0_y_y_xxyzzz[k];

                g_0_y_yy_xyyyy[k] = g_y_xyyyy[k] - g_0_y_y_xyyyy[k] * ab_y + g_0_y_y_xyyyyy[k];

                g_0_y_yy_xyyyz[k] = g_y_xyyyz[k] - g_0_y_y_xyyyz[k] * ab_y + g_0_y_y_xyyyyz[k];

                g_0_y_yy_xyyzz[k] = g_y_xyyzz[k] - g_0_y_y_xyyzz[k] * ab_y + g_0_y_y_xyyyzz[k];

                g_0_y_yy_xyzzz[k] = g_y_xyzzz[k] - g_0_y_y_xyzzz[k] * ab_y + g_0_y_y_xyyzzz[k];

                g_0_y_yy_xzzzz[k] = g_y_xzzzz[k] - g_0_y_y_xzzzz[k] * ab_y + g_0_y_y_xyzzzz[k];

                g_0_y_yy_yyyyy[k] = g_y_yyyyy[k] - g_0_y_y_yyyyy[k] * ab_y + g_0_y_y_yyyyyy[k];

                g_0_y_yy_yyyyz[k] = g_y_yyyyz[k] - g_0_y_y_yyyyz[k] * ab_y + g_0_y_y_yyyyyz[k];

                g_0_y_yy_yyyzz[k] = g_y_yyyzz[k] - g_0_y_y_yyyzz[k] * ab_y + g_0_y_y_yyyyzz[k];

                g_0_y_yy_yyzzz[k] = g_y_yyzzz[k] - g_0_y_y_yyzzz[k] * ab_y + g_0_y_y_yyyzzz[k];

                g_0_y_yy_yzzzz[k] = g_y_yzzzz[k] - g_0_y_y_yzzzz[k] * ab_y + g_0_y_y_yyzzzz[k];

                g_0_y_yy_zzzzz[k] = g_y_zzzzz[k] - g_0_y_y_zzzzz[k] * ab_y + g_0_y_y_yzzzzz[k];
            }

            /// Set up 210-231 components of targeted buffer : cbuffer.data(

            auto g_0_y_yz_xxxxx = cbuffer.data(dh_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_y_yz_xxxxy = cbuffer.data(dh_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_y_yz_xxxxz = cbuffer.data(dh_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_y_yz_xxxyy = cbuffer.data(dh_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_y_yz_xxxyz = cbuffer.data(dh_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_y_yz_xxxzz = cbuffer.data(dh_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_y_yz_xxyyy = cbuffer.data(dh_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_y_yz_xxyyz = cbuffer.data(dh_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_y_yz_xxyzz = cbuffer.data(dh_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_y_yz_xxzzz = cbuffer.data(dh_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_y_yz_xyyyy = cbuffer.data(dh_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_y_yz_xyyyz = cbuffer.data(dh_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_y_yz_xyyzz = cbuffer.data(dh_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_y_yz_xyzzz = cbuffer.data(dh_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_y_yz_xzzzz = cbuffer.data(dh_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_y_yz_yyyyy = cbuffer.data(dh_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_y_yz_yyyyz = cbuffer.data(dh_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_y_yz_yyyzz = cbuffer.data(dh_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_y_yz_yyzzz = cbuffer.data(dh_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_y_yz_yzzzz = cbuffer.data(dh_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_y_yz_zzzzz = cbuffer.data(dh_geom_01_off + 230 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_y_xxxxx, g_0_y_y_xxxxxz, g_0_y_y_xxxxy, g_0_y_y_xxxxyz, g_0_y_y_xxxxz, g_0_y_y_xxxxzz, g_0_y_y_xxxyy, g_0_y_y_xxxyyz, g_0_y_y_xxxyz, g_0_y_y_xxxyzz, g_0_y_y_xxxzz, g_0_y_y_xxxzzz, g_0_y_y_xxyyy, g_0_y_y_xxyyyz, g_0_y_y_xxyyz, g_0_y_y_xxyyzz, g_0_y_y_xxyzz, g_0_y_y_xxyzzz, g_0_y_y_xxzzz, g_0_y_y_xxzzzz, g_0_y_y_xyyyy, g_0_y_y_xyyyyz, g_0_y_y_xyyyz, g_0_y_y_xyyyzz, g_0_y_y_xyyzz, g_0_y_y_xyyzzz, g_0_y_y_xyzzz, g_0_y_y_xyzzzz, g_0_y_y_xzzzz, g_0_y_y_xzzzzz, g_0_y_y_yyyyy, g_0_y_y_yyyyyz, g_0_y_y_yyyyz, g_0_y_y_yyyyzz, g_0_y_y_yyyzz, g_0_y_y_yyyzzz, g_0_y_y_yyzzz, g_0_y_y_yyzzzz, g_0_y_y_yzzzz, g_0_y_y_yzzzzz, g_0_y_y_zzzzz, g_0_y_y_zzzzzz, g_0_y_yz_xxxxx, g_0_y_yz_xxxxy, g_0_y_yz_xxxxz, g_0_y_yz_xxxyy, g_0_y_yz_xxxyz, g_0_y_yz_xxxzz, g_0_y_yz_xxyyy, g_0_y_yz_xxyyz, g_0_y_yz_xxyzz, g_0_y_yz_xxzzz, g_0_y_yz_xyyyy, g_0_y_yz_xyyyz, g_0_y_yz_xyyzz, g_0_y_yz_xyzzz, g_0_y_yz_xzzzz, g_0_y_yz_yyyyy, g_0_y_yz_yyyyz, g_0_y_yz_yyyzz, g_0_y_yz_yyzzz, g_0_y_yz_yzzzz, g_0_y_yz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yz_xxxxx[k] = -g_0_y_y_xxxxx[k] * ab_z + g_0_y_y_xxxxxz[k];

                g_0_y_yz_xxxxy[k] = -g_0_y_y_xxxxy[k] * ab_z + g_0_y_y_xxxxyz[k];

                g_0_y_yz_xxxxz[k] = -g_0_y_y_xxxxz[k] * ab_z + g_0_y_y_xxxxzz[k];

                g_0_y_yz_xxxyy[k] = -g_0_y_y_xxxyy[k] * ab_z + g_0_y_y_xxxyyz[k];

                g_0_y_yz_xxxyz[k] = -g_0_y_y_xxxyz[k] * ab_z + g_0_y_y_xxxyzz[k];

                g_0_y_yz_xxxzz[k] = -g_0_y_y_xxxzz[k] * ab_z + g_0_y_y_xxxzzz[k];

                g_0_y_yz_xxyyy[k] = -g_0_y_y_xxyyy[k] * ab_z + g_0_y_y_xxyyyz[k];

                g_0_y_yz_xxyyz[k] = -g_0_y_y_xxyyz[k] * ab_z + g_0_y_y_xxyyzz[k];

                g_0_y_yz_xxyzz[k] = -g_0_y_y_xxyzz[k] * ab_z + g_0_y_y_xxyzzz[k];

                g_0_y_yz_xxzzz[k] = -g_0_y_y_xxzzz[k] * ab_z + g_0_y_y_xxzzzz[k];

                g_0_y_yz_xyyyy[k] = -g_0_y_y_xyyyy[k] * ab_z + g_0_y_y_xyyyyz[k];

                g_0_y_yz_xyyyz[k] = -g_0_y_y_xyyyz[k] * ab_z + g_0_y_y_xyyyzz[k];

                g_0_y_yz_xyyzz[k] = -g_0_y_y_xyyzz[k] * ab_z + g_0_y_y_xyyzzz[k];

                g_0_y_yz_xyzzz[k] = -g_0_y_y_xyzzz[k] * ab_z + g_0_y_y_xyzzzz[k];

                g_0_y_yz_xzzzz[k] = -g_0_y_y_xzzzz[k] * ab_z + g_0_y_y_xzzzzz[k];

                g_0_y_yz_yyyyy[k] = -g_0_y_y_yyyyy[k] * ab_z + g_0_y_y_yyyyyz[k];

                g_0_y_yz_yyyyz[k] = -g_0_y_y_yyyyz[k] * ab_z + g_0_y_y_yyyyzz[k];

                g_0_y_yz_yyyzz[k] = -g_0_y_y_yyyzz[k] * ab_z + g_0_y_y_yyyzzz[k];

                g_0_y_yz_yyzzz[k] = -g_0_y_y_yyzzz[k] * ab_z + g_0_y_y_yyzzzz[k];

                g_0_y_yz_yzzzz[k] = -g_0_y_y_yzzzz[k] * ab_z + g_0_y_y_yzzzzz[k];

                g_0_y_yz_zzzzz[k] = -g_0_y_y_zzzzz[k] * ab_z + g_0_y_y_zzzzzz[k];
            }

            /// Set up 231-252 components of targeted buffer : cbuffer.data(

            auto g_0_y_zz_xxxxx = cbuffer.data(dh_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_y_zz_xxxxy = cbuffer.data(dh_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_y_zz_xxxxz = cbuffer.data(dh_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_y_zz_xxxyy = cbuffer.data(dh_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_y_zz_xxxyz = cbuffer.data(dh_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_y_zz_xxxzz = cbuffer.data(dh_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_y_zz_xxyyy = cbuffer.data(dh_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_y_zz_xxyyz = cbuffer.data(dh_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_y_zz_xxyzz = cbuffer.data(dh_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_y_zz_xxzzz = cbuffer.data(dh_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_y_zz_xyyyy = cbuffer.data(dh_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_y_zz_xyyyz = cbuffer.data(dh_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_y_zz_xyyzz = cbuffer.data(dh_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_y_zz_xyzzz = cbuffer.data(dh_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_y_zz_xzzzz = cbuffer.data(dh_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_y_zz_yyyyy = cbuffer.data(dh_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_y_zz_yyyyz = cbuffer.data(dh_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_y_zz_yyyzz = cbuffer.data(dh_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_y_zz_yyzzz = cbuffer.data(dh_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_y_zz_yzzzz = cbuffer.data(dh_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_y_zz_zzzzz = cbuffer.data(dh_geom_01_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_z_xxxxx, g_0_y_z_xxxxxz, g_0_y_z_xxxxy, g_0_y_z_xxxxyz, g_0_y_z_xxxxz, g_0_y_z_xxxxzz, g_0_y_z_xxxyy, g_0_y_z_xxxyyz, g_0_y_z_xxxyz, g_0_y_z_xxxyzz, g_0_y_z_xxxzz, g_0_y_z_xxxzzz, g_0_y_z_xxyyy, g_0_y_z_xxyyyz, g_0_y_z_xxyyz, g_0_y_z_xxyyzz, g_0_y_z_xxyzz, g_0_y_z_xxyzzz, g_0_y_z_xxzzz, g_0_y_z_xxzzzz, g_0_y_z_xyyyy, g_0_y_z_xyyyyz, g_0_y_z_xyyyz, g_0_y_z_xyyyzz, g_0_y_z_xyyzz, g_0_y_z_xyyzzz, g_0_y_z_xyzzz, g_0_y_z_xyzzzz, g_0_y_z_xzzzz, g_0_y_z_xzzzzz, g_0_y_z_yyyyy, g_0_y_z_yyyyyz, g_0_y_z_yyyyz, g_0_y_z_yyyyzz, g_0_y_z_yyyzz, g_0_y_z_yyyzzz, g_0_y_z_yyzzz, g_0_y_z_yyzzzz, g_0_y_z_yzzzz, g_0_y_z_yzzzzz, g_0_y_z_zzzzz, g_0_y_z_zzzzzz, g_0_y_zz_xxxxx, g_0_y_zz_xxxxy, g_0_y_zz_xxxxz, g_0_y_zz_xxxyy, g_0_y_zz_xxxyz, g_0_y_zz_xxxzz, g_0_y_zz_xxyyy, g_0_y_zz_xxyyz, g_0_y_zz_xxyzz, g_0_y_zz_xxzzz, g_0_y_zz_xyyyy, g_0_y_zz_xyyyz, g_0_y_zz_xyyzz, g_0_y_zz_xyzzz, g_0_y_zz_xzzzz, g_0_y_zz_yyyyy, g_0_y_zz_yyyyz, g_0_y_zz_yyyzz, g_0_y_zz_yyzzz, g_0_y_zz_yzzzz, g_0_y_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_zz_xxxxx[k] = -g_0_y_z_xxxxx[k] * ab_z + g_0_y_z_xxxxxz[k];

                g_0_y_zz_xxxxy[k] = -g_0_y_z_xxxxy[k] * ab_z + g_0_y_z_xxxxyz[k];

                g_0_y_zz_xxxxz[k] = -g_0_y_z_xxxxz[k] * ab_z + g_0_y_z_xxxxzz[k];

                g_0_y_zz_xxxyy[k] = -g_0_y_z_xxxyy[k] * ab_z + g_0_y_z_xxxyyz[k];

                g_0_y_zz_xxxyz[k] = -g_0_y_z_xxxyz[k] * ab_z + g_0_y_z_xxxyzz[k];

                g_0_y_zz_xxxzz[k] = -g_0_y_z_xxxzz[k] * ab_z + g_0_y_z_xxxzzz[k];

                g_0_y_zz_xxyyy[k] = -g_0_y_z_xxyyy[k] * ab_z + g_0_y_z_xxyyyz[k];

                g_0_y_zz_xxyyz[k] = -g_0_y_z_xxyyz[k] * ab_z + g_0_y_z_xxyyzz[k];

                g_0_y_zz_xxyzz[k] = -g_0_y_z_xxyzz[k] * ab_z + g_0_y_z_xxyzzz[k];

                g_0_y_zz_xxzzz[k] = -g_0_y_z_xxzzz[k] * ab_z + g_0_y_z_xxzzzz[k];

                g_0_y_zz_xyyyy[k] = -g_0_y_z_xyyyy[k] * ab_z + g_0_y_z_xyyyyz[k];

                g_0_y_zz_xyyyz[k] = -g_0_y_z_xyyyz[k] * ab_z + g_0_y_z_xyyyzz[k];

                g_0_y_zz_xyyzz[k] = -g_0_y_z_xyyzz[k] * ab_z + g_0_y_z_xyyzzz[k];

                g_0_y_zz_xyzzz[k] = -g_0_y_z_xyzzz[k] * ab_z + g_0_y_z_xyzzzz[k];

                g_0_y_zz_xzzzz[k] = -g_0_y_z_xzzzz[k] * ab_z + g_0_y_z_xzzzzz[k];

                g_0_y_zz_yyyyy[k] = -g_0_y_z_yyyyy[k] * ab_z + g_0_y_z_yyyyyz[k];

                g_0_y_zz_yyyyz[k] = -g_0_y_z_yyyyz[k] * ab_z + g_0_y_z_yyyyzz[k];

                g_0_y_zz_yyyzz[k] = -g_0_y_z_yyyzz[k] * ab_z + g_0_y_z_yyyzzz[k];

                g_0_y_zz_yyzzz[k] = -g_0_y_z_yyzzz[k] * ab_z + g_0_y_z_yyzzzz[k];

                g_0_y_zz_yzzzz[k] = -g_0_y_z_yzzzz[k] * ab_z + g_0_y_z_yzzzzz[k];

                g_0_y_zz_zzzzz[k] = -g_0_y_z_zzzzz[k] * ab_z + g_0_y_z_zzzzzz[k];
            }

            /// Set up 252-273 components of targeted buffer : cbuffer.data(

            auto g_0_z_xx_xxxxx = cbuffer.data(dh_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_z_xx_xxxxy = cbuffer.data(dh_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_z_xx_xxxxz = cbuffer.data(dh_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_z_xx_xxxyy = cbuffer.data(dh_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_z_xx_xxxyz = cbuffer.data(dh_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_z_xx_xxxzz = cbuffer.data(dh_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_z_xx_xxyyy = cbuffer.data(dh_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_z_xx_xxyyz = cbuffer.data(dh_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_z_xx_xxyzz = cbuffer.data(dh_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_z_xx_xxzzz = cbuffer.data(dh_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_z_xx_xyyyy = cbuffer.data(dh_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_z_xx_xyyyz = cbuffer.data(dh_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_z_xx_xyyzz = cbuffer.data(dh_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_z_xx_xyzzz = cbuffer.data(dh_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_z_xx_xzzzz = cbuffer.data(dh_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_z_xx_yyyyy = cbuffer.data(dh_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_z_xx_yyyyz = cbuffer.data(dh_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_z_xx_yyyzz = cbuffer.data(dh_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_z_xx_yyzzz = cbuffer.data(dh_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_z_xx_yzzzz = cbuffer.data(dh_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_z_xx_zzzzz = cbuffer.data(dh_geom_01_off + 272 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_x_xxxxx, g_0_z_x_xxxxxx, g_0_z_x_xxxxxy, g_0_z_x_xxxxxz, g_0_z_x_xxxxy, g_0_z_x_xxxxyy, g_0_z_x_xxxxyz, g_0_z_x_xxxxz, g_0_z_x_xxxxzz, g_0_z_x_xxxyy, g_0_z_x_xxxyyy, g_0_z_x_xxxyyz, g_0_z_x_xxxyz, g_0_z_x_xxxyzz, g_0_z_x_xxxzz, g_0_z_x_xxxzzz, g_0_z_x_xxyyy, g_0_z_x_xxyyyy, g_0_z_x_xxyyyz, g_0_z_x_xxyyz, g_0_z_x_xxyyzz, g_0_z_x_xxyzz, g_0_z_x_xxyzzz, g_0_z_x_xxzzz, g_0_z_x_xxzzzz, g_0_z_x_xyyyy, g_0_z_x_xyyyyy, g_0_z_x_xyyyyz, g_0_z_x_xyyyz, g_0_z_x_xyyyzz, g_0_z_x_xyyzz, g_0_z_x_xyyzzz, g_0_z_x_xyzzz, g_0_z_x_xyzzzz, g_0_z_x_xzzzz, g_0_z_x_xzzzzz, g_0_z_x_yyyyy, g_0_z_x_yyyyz, g_0_z_x_yyyzz, g_0_z_x_yyzzz, g_0_z_x_yzzzz, g_0_z_x_zzzzz, g_0_z_xx_xxxxx, g_0_z_xx_xxxxy, g_0_z_xx_xxxxz, g_0_z_xx_xxxyy, g_0_z_xx_xxxyz, g_0_z_xx_xxxzz, g_0_z_xx_xxyyy, g_0_z_xx_xxyyz, g_0_z_xx_xxyzz, g_0_z_xx_xxzzz, g_0_z_xx_xyyyy, g_0_z_xx_xyyyz, g_0_z_xx_xyyzz, g_0_z_xx_xyzzz, g_0_z_xx_xzzzz, g_0_z_xx_yyyyy, g_0_z_xx_yyyyz, g_0_z_xx_yyyzz, g_0_z_xx_yyzzz, g_0_z_xx_yzzzz, g_0_z_xx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xx_xxxxx[k] = -g_0_z_x_xxxxx[k] * ab_x + g_0_z_x_xxxxxx[k];

                g_0_z_xx_xxxxy[k] = -g_0_z_x_xxxxy[k] * ab_x + g_0_z_x_xxxxxy[k];

                g_0_z_xx_xxxxz[k] = -g_0_z_x_xxxxz[k] * ab_x + g_0_z_x_xxxxxz[k];

                g_0_z_xx_xxxyy[k] = -g_0_z_x_xxxyy[k] * ab_x + g_0_z_x_xxxxyy[k];

                g_0_z_xx_xxxyz[k] = -g_0_z_x_xxxyz[k] * ab_x + g_0_z_x_xxxxyz[k];

                g_0_z_xx_xxxzz[k] = -g_0_z_x_xxxzz[k] * ab_x + g_0_z_x_xxxxzz[k];

                g_0_z_xx_xxyyy[k] = -g_0_z_x_xxyyy[k] * ab_x + g_0_z_x_xxxyyy[k];

                g_0_z_xx_xxyyz[k] = -g_0_z_x_xxyyz[k] * ab_x + g_0_z_x_xxxyyz[k];

                g_0_z_xx_xxyzz[k] = -g_0_z_x_xxyzz[k] * ab_x + g_0_z_x_xxxyzz[k];

                g_0_z_xx_xxzzz[k] = -g_0_z_x_xxzzz[k] * ab_x + g_0_z_x_xxxzzz[k];

                g_0_z_xx_xyyyy[k] = -g_0_z_x_xyyyy[k] * ab_x + g_0_z_x_xxyyyy[k];

                g_0_z_xx_xyyyz[k] = -g_0_z_x_xyyyz[k] * ab_x + g_0_z_x_xxyyyz[k];

                g_0_z_xx_xyyzz[k] = -g_0_z_x_xyyzz[k] * ab_x + g_0_z_x_xxyyzz[k];

                g_0_z_xx_xyzzz[k] = -g_0_z_x_xyzzz[k] * ab_x + g_0_z_x_xxyzzz[k];

                g_0_z_xx_xzzzz[k] = -g_0_z_x_xzzzz[k] * ab_x + g_0_z_x_xxzzzz[k];

                g_0_z_xx_yyyyy[k] = -g_0_z_x_yyyyy[k] * ab_x + g_0_z_x_xyyyyy[k];

                g_0_z_xx_yyyyz[k] = -g_0_z_x_yyyyz[k] * ab_x + g_0_z_x_xyyyyz[k];

                g_0_z_xx_yyyzz[k] = -g_0_z_x_yyyzz[k] * ab_x + g_0_z_x_xyyyzz[k];

                g_0_z_xx_yyzzz[k] = -g_0_z_x_yyzzz[k] * ab_x + g_0_z_x_xyyzzz[k];

                g_0_z_xx_yzzzz[k] = -g_0_z_x_yzzzz[k] * ab_x + g_0_z_x_xyzzzz[k];

                g_0_z_xx_zzzzz[k] = -g_0_z_x_zzzzz[k] * ab_x + g_0_z_x_xzzzzz[k];
            }

            /// Set up 273-294 components of targeted buffer : cbuffer.data(

            auto g_0_z_xy_xxxxx = cbuffer.data(dh_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_z_xy_xxxxy = cbuffer.data(dh_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_z_xy_xxxxz = cbuffer.data(dh_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_z_xy_xxxyy = cbuffer.data(dh_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_z_xy_xxxyz = cbuffer.data(dh_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_z_xy_xxxzz = cbuffer.data(dh_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_z_xy_xxyyy = cbuffer.data(dh_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_z_xy_xxyyz = cbuffer.data(dh_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_z_xy_xxyzz = cbuffer.data(dh_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_z_xy_xxzzz = cbuffer.data(dh_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_z_xy_xyyyy = cbuffer.data(dh_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_z_xy_xyyyz = cbuffer.data(dh_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_z_xy_xyyzz = cbuffer.data(dh_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_z_xy_xyzzz = cbuffer.data(dh_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_z_xy_xzzzz = cbuffer.data(dh_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_z_xy_yyyyy = cbuffer.data(dh_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_z_xy_yyyyz = cbuffer.data(dh_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_z_xy_yyyzz = cbuffer.data(dh_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_z_xy_yyzzz = cbuffer.data(dh_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_z_xy_yzzzz = cbuffer.data(dh_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_z_xy_zzzzz = cbuffer.data(dh_geom_01_off + 293 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xy_xxxxx, g_0_z_xy_xxxxy, g_0_z_xy_xxxxz, g_0_z_xy_xxxyy, g_0_z_xy_xxxyz, g_0_z_xy_xxxzz, g_0_z_xy_xxyyy, g_0_z_xy_xxyyz, g_0_z_xy_xxyzz, g_0_z_xy_xxzzz, g_0_z_xy_xyyyy, g_0_z_xy_xyyyz, g_0_z_xy_xyyzz, g_0_z_xy_xyzzz, g_0_z_xy_xzzzz, g_0_z_xy_yyyyy, g_0_z_xy_yyyyz, g_0_z_xy_yyyzz, g_0_z_xy_yyzzz, g_0_z_xy_yzzzz, g_0_z_xy_zzzzz, g_0_z_y_xxxxx, g_0_z_y_xxxxxx, g_0_z_y_xxxxxy, g_0_z_y_xxxxxz, g_0_z_y_xxxxy, g_0_z_y_xxxxyy, g_0_z_y_xxxxyz, g_0_z_y_xxxxz, g_0_z_y_xxxxzz, g_0_z_y_xxxyy, g_0_z_y_xxxyyy, g_0_z_y_xxxyyz, g_0_z_y_xxxyz, g_0_z_y_xxxyzz, g_0_z_y_xxxzz, g_0_z_y_xxxzzz, g_0_z_y_xxyyy, g_0_z_y_xxyyyy, g_0_z_y_xxyyyz, g_0_z_y_xxyyz, g_0_z_y_xxyyzz, g_0_z_y_xxyzz, g_0_z_y_xxyzzz, g_0_z_y_xxzzz, g_0_z_y_xxzzzz, g_0_z_y_xyyyy, g_0_z_y_xyyyyy, g_0_z_y_xyyyyz, g_0_z_y_xyyyz, g_0_z_y_xyyyzz, g_0_z_y_xyyzz, g_0_z_y_xyyzzz, g_0_z_y_xyzzz, g_0_z_y_xyzzzz, g_0_z_y_xzzzz, g_0_z_y_xzzzzz, g_0_z_y_yyyyy, g_0_z_y_yyyyz, g_0_z_y_yyyzz, g_0_z_y_yyzzz, g_0_z_y_yzzzz, g_0_z_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xy_xxxxx[k] = -g_0_z_y_xxxxx[k] * ab_x + g_0_z_y_xxxxxx[k];

                g_0_z_xy_xxxxy[k] = -g_0_z_y_xxxxy[k] * ab_x + g_0_z_y_xxxxxy[k];

                g_0_z_xy_xxxxz[k] = -g_0_z_y_xxxxz[k] * ab_x + g_0_z_y_xxxxxz[k];

                g_0_z_xy_xxxyy[k] = -g_0_z_y_xxxyy[k] * ab_x + g_0_z_y_xxxxyy[k];

                g_0_z_xy_xxxyz[k] = -g_0_z_y_xxxyz[k] * ab_x + g_0_z_y_xxxxyz[k];

                g_0_z_xy_xxxzz[k] = -g_0_z_y_xxxzz[k] * ab_x + g_0_z_y_xxxxzz[k];

                g_0_z_xy_xxyyy[k] = -g_0_z_y_xxyyy[k] * ab_x + g_0_z_y_xxxyyy[k];

                g_0_z_xy_xxyyz[k] = -g_0_z_y_xxyyz[k] * ab_x + g_0_z_y_xxxyyz[k];

                g_0_z_xy_xxyzz[k] = -g_0_z_y_xxyzz[k] * ab_x + g_0_z_y_xxxyzz[k];

                g_0_z_xy_xxzzz[k] = -g_0_z_y_xxzzz[k] * ab_x + g_0_z_y_xxxzzz[k];

                g_0_z_xy_xyyyy[k] = -g_0_z_y_xyyyy[k] * ab_x + g_0_z_y_xxyyyy[k];

                g_0_z_xy_xyyyz[k] = -g_0_z_y_xyyyz[k] * ab_x + g_0_z_y_xxyyyz[k];

                g_0_z_xy_xyyzz[k] = -g_0_z_y_xyyzz[k] * ab_x + g_0_z_y_xxyyzz[k];

                g_0_z_xy_xyzzz[k] = -g_0_z_y_xyzzz[k] * ab_x + g_0_z_y_xxyzzz[k];

                g_0_z_xy_xzzzz[k] = -g_0_z_y_xzzzz[k] * ab_x + g_0_z_y_xxzzzz[k];

                g_0_z_xy_yyyyy[k] = -g_0_z_y_yyyyy[k] * ab_x + g_0_z_y_xyyyyy[k];

                g_0_z_xy_yyyyz[k] = -g_0_z_y_yyyyz[k] * ab_x + g_0_z_y_xyyyyz[k];

                g_0_z_xy_yyyzz[k] = -g_0_z_y_yyyzz[k] * ab_x + g_0_z_y_xyyyzz[k];

                g_0_z_xy_yyzzz[k] = -g_0_z_y_yyzzz[k] * ab_x + g_0_z_y_xyyzzz[k];

                g_0_z_xy_yzzzz[k] = -g_0_z_y_yzzzz[k] * ab_x + g_0_z_y_xyzzzz[k];

                g_0_z_xy_zzzzz[k] = -g_0_z_y_zzzzz[k] * ab_x + g_0_z_y_xzzzzz[k];
            }

            /// Set up 294-315 components of targeted buffer : cbuffer.data(

            auto g_0_z_xz_xxxxx = cbuffer.data(dh_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_z_xz_xxxxy = cbuffer.data(dh_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_z_xz_xxxxz = cbuffer.data(dh_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_z_xz_xxxyy = cbuffer.data(dh_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_z_xz_xxxyz = cbuffer.data(dh_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_z_xz_xxxzz = cbuffer.data(dh_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_z_xz_xxyyy = cbuffer.data(dh_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_z_xz_xxyyz = cbuffer.data(dh_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_z_xz_xxyzz = cbuffer.data(dh_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_z_xz_xxzzz = cbuffer.data(dh_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_z_xz_xyyyy = cbuffer.data(dh_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_z_xz_xyyyz = cbuffer.data(dh_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_z_xz_xyyzz = cbuffer.data(dh_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_z_xz_xyzzz = cbuffer.data(dh_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_z_xz_xzzzz = cbuffer.data(dh_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_z_xz_yyyyy = cbuffer.data(dh_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_z_xz_yyyyz = cbuffer.data(dh_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_z_xz_yyyzz = cbuffer.data(dh_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_z_xz_yyzzz = cbuffer.data(dh_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_z_xz_yzzzz = cbuffer.data(dh_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_z_xz_zzzzz = cbuffer.data(dh_geom_01_off + 314 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xz_xxxxx, g_0_z_xz_xxxxy, g_0_z_xz_xxxxz, g_0_z_xz_xxxyy, g_0_z_xz_xxxyz, g_0_z_xz_xxxzz, g_0_z_xz_xxyyy, g_0_z_xz_xxyyz, g_0_z_xz_xxyzz, g_0_z_xz_xxzzz, g_0_z_xz_xyyyy, g_0_z_xz_xyyyz, g_0_z_xz_xyyzz, g_0_z_xz_xyzzz, g_0_z_xz_xzzzz, g_0_z_xz_yyyyy, g_0_z_xz_yyyyz, g_0_z_xz_yyyzz, g_0_z_xz_yyzzz, g_0_z_xz_yzzzz, g_0_z_xz_zzzzz, g_0_z_z_xxxxx, g_0_z_z_xxxxxx, g_0_z_z_xxxxxy, g_0_z_z_xxxxxz, g_0_z_z_xxxxy, g_0_z_z_xxxxyy, g_0_z_z_xxxxyz, g_0_z_z_xxxxz, g_0_z_z_xxxxzz, g_0_z_z_xxxyy, g_0_z_z_xxxyyy, g_0_z_z_xxxyyz, g_0_z_z_xxxyz, g_0_z_z_xxxyzz, g_0_z_z_xxxzz, g_0_z_z_xxxzzz, g_0_z_z_xxyyy, g_0_z_z_xxyyyy, g_0_z_z_xxyyyz, g_0_z_z_xxyyz, g_0_z_z_xxyyzz, g_0_z_z_xxyzz, g_0_z_z_xxyzzz, g_0_z_z_xxzzz, g_0_z_z_xxzzzz, g_0_z_z_xyyyy, g_0_z_z_xyyyyy, g_0_z_z_xyyyyz, g_0_z_z_xyyyz, g_0_z_z_xyyyzz, g_0_z_z_xyyzz, g_0_z_z_xyyzzz, g_0_z_z_xyzzz, g_0_z_z_xyzzzz, g_0_z_z_xzzzz, g_0_z_z_xzzzzz, g_0_z_z_yyyyy, g_0_z_z_yyyyz, g_0_z_z_yyyzz, g_0_z_z_yyzzz, g_0_z_z_yzzzz, g_0_z_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xz_xxxxx[k] = -g_0_z_z_xxxxx[k] * ab_x + g_0_z_z_xxxxxx[k];

                g_0_z_xz_xxxxy[k] = -g_0_z_z_xxxxy[k] * ab_x + g_0_z_z_xxxxxy[k];

                g_0_z_xz_xxxxz[k] = -g_0_z_z_xxxxz[k] * ab_x + g_0_z_z_xxxxxz[k];

                g_0_z_xz_xxxyy[k] = -g_0_z_z_xxxyy[k] * ab_x + g_0_z_z_xxxxyy[k];

                g_0_z_xz_xxxyz[k] = -g_0_z_z_xxxyz[k] * ab_x + g_0_z_z_xxxxyz[k];

                g_0_z_xz_xxxzz[k] = -g_0_z_z_xxxzz[k] * ab_x + g_0_z_z_xxxxzz[k];

                g_0_z_xz_xxyyy[k] = -g_0_z_z_xxyyy[k] * ab_x + g_0_z_z_xxxyyy[k];

                g_0_z_xz_xxyyz[k] = -g_0_z_z_xxyyz[k] * ab_x + g_0_z_z_xxxyyz[k];

                g_0_z_xz_xxyzz[k] = -g_0_z_z_xxyzz[k] * ab_x + g_0_z_z_xxxyzz[k];

                g_0_z_xz_xxzzz[k] = -g_0_z_z_xxzzz[k] * ab_x + g_0_z_z_xxxzzz[k];

                g_0_z_xz_xyyyy[k] = -g_0_z_z_xyyyy[k] * ab_x + g_0_z_z_xxyyyy[k];

                g_0_z_xz_xyyyz[k] = -g_0_z_z_xyyyz[k] * ab_x + g_0_z_z_xxyyyz[k];

                g_0_z_xz_xyyzz[k] = -g_0_z_z_xyyzz[k] * ab_x + g_0_z_z_xxyyzz[k];

                g_0_z_xz_xyzzz[k] = -g_0_z_z_xyzzz[k] * ab_x + g_0_z_z_xxyzzz[k];

                g_0_z_xz_xzzzz[k] = -g_0_z_z_xzzzz[k] * ab_x + g_0_z_z_xxzzzz[k];

                g_0_z_xz_yyyyy[k] = -g_0_z_z_yyyyy[k] * ab_x + g_0_z_z_xyyyyy[k];

                g_0_z_xz_yyyyz[k] = -g_0_z_z_yyyyz[k] * ab_x + g_0_z_z_xyyyyz[k];

                g_0_z_xz_yyyzz[k] = -g_0_z_z_yyyzz[k] * ab_x + g_0_z_z_xyyyzz[k];

                g_0_z_xz_yyzzz[k] = -g_0_z_z_yyzzz[k] * ab_x + g_0_z_z_xyyzzz[k];

                g_0_z_xz_yzzzz[k] = -g_0_z_z_yzzzz[k] * ab_x + g_0_z_z_xyzzzz[k];

                g_0_z_xz_zzzzz[k] = -g_0_z_z_zzzzz[k] * ab_x + g_0_z_z_xzzzzz[k];
            }

            /// Set up 315-336 components of targeted buffer : cbuffer.data(

            auto g_0_z_yy_xxxxx = cbuffer.data(dh_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_z_yy_xxxxy = cbuffer.data(dh_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_z_yy_xxxxz = cbuffer.data(dh_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_z_yy_xxxyy = cbuffer.data(dh_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_z_yy_xxxyz = cbuffer.data(dh_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_z_yy_xxxzz = cbuffer.data(dh_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_z_yy_xxyyy = cbuffer.data(dh_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_z_yy_xxyyz = cbuffer.data(dh_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_z_yy_xxyzz = cbuffer.data(dh_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_z_yy_xxzzz = cbuffer.data(dh_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_z_yy_xyyyy = cbuffer.data(dh_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_z_yy_xyyyz = cbuffer.data(dh_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_z_yy_xyyzz = cbuffer.data(dh_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_z_yy_xyzzz = cbuffer.data(dh_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_z_yy_xzzzz = cbuffer.data(dh_geom_01_off + 329 * ccomps * dcomps);

            auto g_0_z_yy_yyyyy = cbuffer.data(dh_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_z_yy_yyyyz = cbuffer.data(dh_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_z_yy_yyyzz = cbuffer.data(dh_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_z_yy_yyzzz = cbuffer.data(dh_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_z_yy_yzzzz = cbuffer.data(dh_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_z_yy_zzzzz = cbuffer.data(dh_geom_01_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_y_xxxxx, g_0_z_y_xxxxxy, g_0_z_y_xxxxy, g_0_z_y_xxxxyy, g_0_z_y_xxxxyz, g_0_z_y_xxxxz, g_0_z_y_xxxyy, g_0_z_y_xxxyyy, g_0_z_y_xxxyyz, g_0_z_y_xxxyz, g_0_z_y_xxxyzz, g_0_z_y_xxxzz, g_0_z_y_xxyyy, g_0_z_y_xxyyyy, g_0_z_y_xxyyyz, g_0_z_y_xxyyz, g_0_z_y_xxyyzz, g_0_z_y_xxyzz, g_0_z_y_xxyzzz, g_0_z_y_xxzzz, g_0_z_y_xyyyy, g_0_z_y_xyyyyy, g_0_z_y_xyyyyz, g_0_z_y_xyyyz, g_0_z_y_xyyyzz, g_0_z_y_xyyzz, g_0_z_y_xyyzzz, g_0_z_y_xyzzz, g_0_z_y_xyzzzz, g_0_z_y_xzzzz, g_0_z_y_yyyyy, g_0_z_y_yyyyyy, g_0_z_y_yyyyyz, g_0_z_y_yyyyz, g_0_z_y_yyyyzz, g_0_z_y_yyyzz, g_0_z_y_yyyzzz, g_0_z_y_yyzzz, g_0_z_y_yyzzzz, g_0_z_y_yzzzz, g_0_z_y_yzzzzz, g_0_z_y_zzzzz, g_0_z_yy_xxxxx, g_0_z_yy_xxxxy, g_0_z_yy_xxxxz, g_0_z_yy_xxxyy, g_0_z_yy_xxxyz, g_0_z_yy_xxxzz, g_0_z_yy_xxyyy, g_0_z_yy_xxyyz, g_0_z_yy_xxyzz, g_0_z_yy_xxzzz, g_0_z_yy_xyyyy, g_0_z_yy_xyyyz, g_0_z_yy_xyyzz, g_0_z_yy_xyzzz, g_0_z_yy_xzzzz, g_0_z_yy_yyyyy, g_0_z_yy_yyyyz, g_0_z_yy_yyyzz, g_0_z_yy_yyzzz, g_0_z_yy_yzzzz, g_0_z_yy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yy_xxxxx[k] = -g_0_z_y_xxxxx[k] * ab_y + g_0_z_y_xxxxxy[k];

                g_0_z_yy_xxxxy[k] = -g_0_z_y_xxxxy[k] * ab_y + g_0_z_y_xxxxyy[k];

                g_0_z_yy_xxxxz[k] = -g_0_z_y_xxxxz[k] * ab_y + g_0_z_y_xxxxyz[k];

                g_0_z_yy_xxxyy[k] = -g_0_z_y_xxxyy[k] * ab_y + g_0_z_y_xxxyyy[k];

                g_0_z_yy_xxxyz[k] = -g_0_z_y_xxxyz[k] * ab_y + g_0_z_y_xxxyyz[k];

                g_0_z_yy_xxxzz[k] = -g_0_z_y_xxxzz[k] * ab_y + g_0_z_y_xxxyzz[k];

                g_0_z_yy_xxyyy[k] = -g_0_z_y_xxyyy[k] * ab_y + g_0_z_y_xxyyyy[k];

                g_0_z_yy_xxyyz[k] = -g_0_z_y_xxyyz[k] * ab_y + g_0_z_y_xxyyyz[k];

                g_0_z_yy_xxyzz[k] = -g_0_z_y_xxyzz[k] * ab_y + g_0_z_y_xxyyzz[k];

                g_0_z_yy_xxzzz[k] = -g_0_z_y_xxzzz[k] * ab_y + g_0_z_y_xxyzzz[k];

                g_0_z_yy_xyyyy[k] = -g_0_z_y_xyyyy[k] * ab_y + g_0_z_y_xyyyyy[k];

                g_0_z_yy_xyyyz[k] = -g_0_z_y_xyyyz[k] * ab_y + g_0_z_y_xyyyyz[k];

                g_0_z_yy_xyyzz[k] = -g_0_z_y_xyyzz[k] * ab_y + g_0_z_y_xyyyzz[k];

                g_0_z_yy_xyzzz[k] = -g_0_z_y_xyzzz[k] * ab_y + g_0_z_y_xyyzzz[k];

                g_0_z_yy_xzzzz[k] = -g_0_z_y_xzzzz[k] * ab_y + g_0_z_y_xyzzzz[k];

                g_0_z_yy_yyyyy[k] = -g_0_z_y_yyyyy[k] * ab_y + g_0_z_y_yyyyyy[k];

                g_0_z_yy_yyyyz[k] = -g_0_z_y_yyyyz[k] * ab_y + g_0_z_y_yyyyyz[k];

                g_0_z_yy_yyyzz[k] = -g_0_z_y_yyyzz[k] * ab_y + g_0_z_y_yyyyzz[k];

                g_0_z_yy_yyzzz[k] = -g_0_z_y_yyzzz[k] * ab_y + g_0_z_y_yyyzzz[k];

                g_0_z_yy_yzzzz[k] = -g_0_z_y_yzzzz[k] * ab_y + g_0_z_y_yyzzzz[k];

                g_0_z_yy_zzzzz[k] = -g_0_z_y_zzzzz[k] * ab_y + g_0_z_y_yzzzzz[k];
            }

            /// Set up 336-357 components of targeted buffer : cbuffer.data(

            auto g_0_z_yz_xxxxx = cbuffer.data(dh_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_z_yz_xxxxy = cbuffer.data(dh_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_z_yz_xxxxz = cbuffer.data(dh_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_z_yz_xxxyy = cbuffer.data(dh_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_z_yz_xxxyz = cbuffer.data(dh_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_z_yz_xxxzz = cbuffer.data(dh_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_z_yz_xxyyy = cbuffer.data(dh_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_z_yz_xxyyz = cbuffer.data(dh_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_z_yz_xxyzz = cbuffer.data(dh_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_z_yz_xxzzz = cbuffer.data(dh_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_z_yz_xyyyy = cbuffer.data(dh_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_z_yz_xyyyz = cbuffer.data(dh_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_z_yz_xyyzz = cbuffer.data(dh_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_z_yz_xyzzz = cbuffer.data(dh_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_z_yz_xzzzz = cbuffer.data(dh_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_z_yz_yyyyy = cbuffer.data(dh_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_z_yz_yyyyz = cbuffer.data(dh_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_z_yz_yyyzz = cbuffer.data(dh_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_z_yz_yyzzz = cbuffer.data(dh_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_z_yz_yzzzz = cbuffer.data(dh_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_z_yz_zzzzz = cbuffer.data(dh_geom_01_off + 356 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yz_xxxxx, g_0_z_yz_xxxxy, g_0_z_yz_xxxxz, g_0_z_yz_xxxyy, g_0_z_yz_xxxyz, g_0_z_yz_xxxzz, g_0_z_yz_xxyyy, g_0_z_yz_xxyyz, g_0_z_yz_xxyzz, g_0_z_yz_xxzzz, g_0_z_yz_xyyyy, g_0_z_yz_xyyyz, g_0_z_yz_xyyzz, g_0_z_yz_xyzzz, g_0_z_yz_xzzzz, g_0_z_yz_yyyyy, g_0_z_yz_yyyyz, g_0_z_yz_yyyzz, g_0_z_yz_yyzzz, g_0_z_yz_yzzzz, g_0_z_yz_zzzzz, g_0_z_z_xxxxx, g_0_z_z_xxxxxy, g_0_z_z_xxxxy, g_0_z_z_xxxxyy, g_0_z_z_xxxxyz, g_0_z_z_xxxxz, g_0_z_z_xxxyy, g_0_z_z_xxxyyy, g_0_z_z_xxxyyz, g_0_z_z_xxxyz, g_0_z_z_xxxyzz, g_0_z_z_xxxzz, g_0_z_z_xxyyy, g_0_z_z_xxyyyy, g_0_z_z_xxyyyz, g_0_z_z_xxyyz, g_0_z_z_xxyyzz, g_0_z_z_xxyzz, g_0_z_z_xxyzzz, g_0_z_z_xxzzz, g_0_z_z_xyyyy, g_0_z_z_xyyyyy, g_0_z_z_xyyyyz, g_0_z_z_xyyyz, g_0_z_z_xyyyzz, g_0_z_z_xyyzz, g_0_z_z_xyyzzz, g_0_z_z_xyzzz, g_0_z_z_xyzzzz, g_0_z_z_xzzzz, g_0_z_z_yyyyy, g_0_z_z_yyyyyy, g_0_z_z_yyyyyz, g_0_z_z_yyyyz, g_0_z_z_yyyyzz, g_0_z_z_yyyzz, g_0_z_z_yyyzzz, g_0_z_z_yyzzz, g_0_z_z_yyzzzz, g_0_z_z_yzzzz, g_0_z_z_yzzzzz, g_0_z_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yz_xxxxx[k] = -g_0_z_z_xxxxx[k] * ab_y + g_0_z_z_xxxxxy[k];

                g_0_z_yz_xxxxy[k] = -g_0_z_z_xxxxy[k] * ab_y + g_0_z_z_xxxxyy[k];

                g_0_z_yz_xxxxz[k] = -g_0_z_z_xxxxz[k] * ab_y + g_0_z_z_xxxxyz[k];

                g_0_z_yz_xxxyy[k] = -g_0_z_z_xxxyy[k] * ab_y + g_0_z_z_xxxyyy[k];

                g_0_z_yz_xxxyz[k] = -g_0_z_z_xxxyz[k] * ab_y + g_0_z_z_xxxyyz[k];

                g_0_z_yz_xxxzz[k] = -g_0_z_z_xxxzz[k] * ab_y + g_0_z_z_xxxyzz[k];

                g_0_z_yz_xxyyy[k] = -g_0_z_z_xxyyy[k] * ab_y + g_0_z_z_xxyyyy[k];

                g_0_z_yz_xxyyz[k] = -g_0_z_z_xxyyz[k] * ab_y + g_0_z_z_xxyyyz[k];

                g_0_z_yz_xxyzz[k] = -g_0_z_z_xxyzz[k] * ab_y + g_0_z_z_xxyyzz[k];

                g_0_z_yz_xxzzz[k] = -g_0_z_z_xxzzz[k] * ab_y + g_0_z_z_xxyzzz[k];

                g_0_z_yz_xyyyy[k] = -g_0_z_z_xyyyy[k] * ab_y + g_0_z_z_xyyyyy[k];

                g_0_z_yz_xyyyz[k] = -g_0_z_z_xyyyz[k] * ab_y + g_0_z_z_xyyyyz[k];

                g_0_z_yz_xyyzz[k] = -g_0_z_z_xyyzz[k] * ab_y + g_0_z_z_xyyyzz[k];

                g_0_z_yz_xyzzz[k] = -g_0_z_z_xyzzz[k] * ab_y + g_0_z_z_xyyzzz[k];

                g_0_z_yz_xzzzz[k] = -g_0_z_z_xzzzz[k] * ab_y + g_0_z_z_xyzzzz[k];

                g_0_z_yz_yyyyy[k] = -g_0_z_z_yyyyy[k] * ab_y + g_0_z_z_yyyyyy[k];

                g_0_z_yz_yyyyz[k] = -g_0_z_z_yyyyz[k] * ab_y + g_0_z_z_yyyyyz[k];

                g_0_z_yz_yyyzz[k] = -g_0_z_z_yyyzz[k] * ab_y + g_0_z_z_yyyyzz[k];

                g_0_z_yz_yyzzz[k] = -g_0_z_z_yyzzz[k] * ab_y + g_0_z_z_yyyzzz[k];

                g_0_z_yz_yzzzz[k] = -g_0_z_z_yzzzz[k] * ab_y + g_0_z_z_yyzzzz[k];

                g_0_z_yz_zzzzz[k] = -g_0_z_z_zzzzz[k] * ab_y + g_0_z_z_yzzzzz[k];
            }

            /// Set up 357-378 components of targeted buffer : cbuffer.data(

            auto g_0_z_zz_xxxxx = cbuffer.data(dh_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_z_zz_xxxxy = cbuffer.data(dh_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_z_zz_xxxxz = cbuffer.data(dh_geom_01_off + 359 * ccomps * dcomps);

            auto g_0_z_zz_xxxyy = cbuffer.data(dh_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_z_zz_xxxyz = cbuffer.data(dh_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_z_zz_xxxzz = cbuffer.data(dh_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_z_zz_xxyyy = cbuffer.data(dh_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_z_zz_xxyyz = cbuffer.data(dh_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_z_zz_xxyzz = cbuffer.data(dh_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_z_zz_xxzzz = cbuffer.data(dh_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_z_zz_xyyyy = cbuffer.data(dh_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_z_zz_xyyyz = cbuffer.data(dh_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_z_zz_xyyzz = cbuffer.data(dh_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_z_zz_xyzzz = cbuffer.data(dh_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_z_zz_xzzzz = cbuffer.data(dh_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_z_zz_yyyyy = cbuffer.data(dh_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_z_zz_yyyyz = cbuffer.data(dh_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_z_zz_yyyzz = cbuffer.data(dh_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_z_zz_yyzzz = cbuffer.data(dh_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_z_zz_yzzzz = cbuffer.data(dh_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_z_zz_zzzzz = cbuffer.data(dh_geom_01_off + 377 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_z_xxxxx, g_0_z_z_xxxxxz, g_0_z_z_xxxxy, g_0_z_z_xxxxyz, g_0_z_z_xxxxz, g_0_z_z_xxxxzz, g_0_z_z_xxxyy, g_0_z_z_xxxyyz, g_0_z_z_xxxyz, g_0_z_z_xxxyzz, g_0_z_z_xxxzz, g_0_z_z_xxxzzz, g_0_z_z_xxyyy, g_0_z_z_xxyyyz, g_0_z_z_xxyyz, g_0_z_z_xxyyzz, g_0_z_z_xxyzz, g_0_z_z_xxyzzz, g_0_z_z_xxzzz, g_0_z_z_xxzzzz, g_0_z_z_xyyyy, g_0_z_z_xyyyyz, g_0_z_z_xyyyz, g_0_z_z_xyyyzz, g_0_z_z_xyyzz, g_0_z_z_xyyzzz, g_0_z_z_xyzzz, g_0_z_z_xyzzzz, g_0_z_z_xzzzz, g_0_z_z_xzzzzz, g_0_z_z_yyyyy, g_0_z_z_yyyyyz, g_0_z_z_yyyyz, g_0_z_z_yyyyzz, g_0_z_z_yyyzz, g_0_z_z_yyyzzz, g_0_z_z_yyzzz, g_0_z_z_yyzzzz, g_0_z_z_yzzzz, g_0_z_z_yzzzzz, g_0_z_z_zzzzz, g_0_z_z_zzzzzz, g_0_z_zz_xxxxx, g_0_z_zz_xxxxy, g_0_z_zz_xxxxz, g_0_z_zz_xxxyy, g_0_z_zz_xxxyz, g_0_z_zz_xxxzz, g_0_z_zz_xxyyy, g_0_z_zz_xxyyz, g_0_z_zz_xxyzz, g_0_z_zz_xxzzz, g_0_z_zz_xyyyy, g_0_z_zz_xyyyz, g_0_z_zz_xyyzz, g_0_z_zz_xyzzz, g_0_z_zz_xzzzz, g_0_z_zz_yyyyy, g_0_z_zz_yyyyz, g_0_z_zz_yyyzz, g_0_z_zz_yyzzz, g_0_z_zz_yzzzz, g_0_z_zz_zzzzz, g_z_xxxxx, g_z_xxxxy, g_z_xxxxz, g_z_xxxyy, g_z_xxxyz, g_z_xxxzz, g_z_xxyyy, g_z_xxyyz, g_z_xxyzz, g_z_xxzzz, g_z_xyyyy, g_z_xyyyz, g_z_xyyzz, g_z_xyzzz, g_z_xzzzz, g_z_yyyyy, g_z_yyyyz, g_z_yyyzz, g_z_yyzzz, g_z_yzzzz, g_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_zz_xxxxx[k] = g_z_xxxxx[k] - g_0_z_z_xxxxx[k] * ab_z + g_0_z_z_xxxxxz[k];

                g_0_z_zz_xxxxy[k] = g_z_xxxxy[k] - g_0_z_z_xxxxy[k] * ab_z + g_0_z_z_xxxxyz[k];

                g_0_z_zz_xxxxz[k] = g_z_xxxxz[k] - g_0_z_z_xxxxz[k] * ab_z + g_0_z_z_xxxxzz[k];

                g_0_z_zz_xxxyy[k] = g_z_xxxyy[k] - g_0_z_z_xxxyy[k] * ab_z + g_0_z_z_xxxyyz[k];

                g_0_z_zz_xxxyz[k] = g_z_xxxyz[k] - g_0_z_z_xxxyz[k] * ab_z + g_0_z_z_xxxyzz[k];

                g_0_z_zz_xxxzz[k] = g_z_xxxzz[k] - g_0_z_z_xxxzz[k] * ab_z + g_0_z_z_xxxzzz[k];

                g_0_z_zz_xxyyy[k] = g_z_xxyyy[k] - g_0_z_z_xxyyy[k] * ab_z + g_0_z_z_xxyyyz[k];

                g_0_z_zz_xxyyz[k] = g_z_xxyyz[k] - g_0_z_z_xxyyz[k] * ab_z + g_0_z_z_xxyyzz[k];

                g_0_z_zz_xxyzz[k] = g_z_xxyzz[k] - g_0_z_z_xxyzz[k] * ab_z + g_0_z_z_xxyzzz[k];

                g_0_z_zz_xxzzz[k] = g_z_xxzzz[k] - g_0_z_z_xxzzz[k] * ab_z + g_0_z_z_xxzzzz[k];

                g_0_z_zz_xyyyy[k] = g_z_xyyyy[k] - g_0_z_z_xyyyy[k] * ab_z + g_0_z_z_xyyyyz[k];

                g_0_z_zz_xyyyz[k] = g_z_xyyyz[k] - g_0_z_z_xyyyz[k] * ab_z + g_0_z_z_xyyyzz[k];

                g_0_z_zz_xyyzz[k] = g_z_xyyzz[k] - g_0_z_z_xyyzz[k] * ab_z + g_0_z_z_xyyzzz[k];

                g_0_z_zz_xyzzz[k] = g_z_xyzzz[k] - g_0_z_z_xyzzz[k] * ab_z + g_0_z_z_xyzzzz[k];

                g_0_z_zz_xzzzz[k] = g_z_xzzzz[k] - g_0_z_z_xzzzz[k] * ab_z + g_0_z_z_xzzzzz[k];

                g_0_z_zz_yyyyy[k] = g_z_yyyyy[k] - g_0_z_z_yyyyy[k] * ab_z + g_0_z_z_yyyyyz[k];

                g_0_z_zz_yyyyz[k] = g_z_yyyyz[k] - g_0_z_z_yyyyz[k] * ab_z + g_0_z_z_yyyyzz[k];

                g_0_z_zz_yyyzz[k] = g_z_yyyzz[k] - g_0_z_z_yyyzz[k] * ab_z + g_0_z_z_yyyzzz[k];

                g_0_z_zz_yyzzz[k] = g_z_yyzzz[k] - g_0_z_z_yyzzz[k] * ab_z + g_0_z_z_yyzzzz[k];

                g_0_z_zz_yzzzz[k] = g_z_yzzzz[k] - g_0_z_z_yzzzz[k] * ab_z + g_0_z_z_yzzzzz[k];

                g_0_z_zz_zzzzz[k] = g_z_zzzzz[k] - g_0_z_z_zzzzz[k] * ab_z + g_0_z_z_zzzzzz[k];
            }
        }
    }
}

} // erirec namespace

