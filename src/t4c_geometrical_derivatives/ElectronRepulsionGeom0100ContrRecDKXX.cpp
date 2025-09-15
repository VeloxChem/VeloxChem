#include "ElectronRepulsionGeom0100ContrRecDKXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_dkxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_dkxx,
                                            const size_t idx_pkxx,
                                            const size_t idx_geom_01_pkxx,
                                            const size_t idx_geom_01_plxx,
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
            /// Set up components of auxilary buffer : PKSS

            const auto pk_off = idx_pkxx + i * dcomps + j;

            auto g_x_xxxxxxx = cbuffer.data(pk_off + 0 * ccomps * dcomps);

            auto g_x_xxxxxxy = cbuffer.data(pk_off + 1 * ccomps * dcomps);

            auto g_x_xxxxxxz = cbuffer.data(pk_off + 2 * ccomps * dcomps);

            auto g_x_xxxxxyy = cbuffer.data(pk_off + 3 * ccomps * dcomps);

            auto g_x_xxxxxyz = cbuffer.data(pk_off + 4 * ccomps * dcomps);

            auto g_x_xxxxxzz = cbuffer.data(pk_off + 5 * ccomps * dcomps);

            auto g_x_xxxxyyy = cbuffer.data(pk_off + 6 * ccomps * dcomps);

            auto g_x_xxxxyyz = cbuffer.data(pk_off + 7 * ccomps * dcomps);

            auto g_x_xxxxyzz = cbuffer.data(pk_off + 8 * ccomps * dcomps);

            auto g_x_xxxxzzz = cbuffer.data(pk_off + 9 * ccomps * dcomps);

            auto g_x_xxxyyyy = cbuffer.data(pk_off + 10 * ccomps * dcomps);

            auto g_x_xxxyyyz = cbuffer.data(pk_off + 11 * ccomps * dcomps);

            auto g_x_xxxyyzz = cbuffer.data(pk_off + 12 * ccomps * dcomps);

            auto g_x_xxxyzzz = cbuffer.data(pk_off + 13 * ccomps * dcomps);

            auto g_x_xxxzzzz = cbuffer.data(pk_off + 14 * ccomps * dcomps);

            auto g_x_xxyyyyy = cbuffer.data(pk_off + 15 * ccomps * dcomps);

            auto g_x_xxyyyyz = cbuffer.data(pk_off + 16 * ccomps * dcomps);

            auto g_x_xxyyyzz = cbuffer.data(pk_off + 17 * ccomps * dcomps);

            auto g_x_xxyyzzz = cbuffer.data(pk_off + 18 * ccomps * dcomps);

            auto g_x_xxyzzzz = cbuffer.data(pk_off + 19 * ccomps * dcomps);

            auto g_x_xxzzzzz = cbuffer.data(pk_off + 20 * ccomps * dcomps);

            auto g_x_xyyyyyy = cbuffer.data(pk_off + 21 * ccomps * dcomps);

            auto g_x_xyyyyyz = cbuffer.data(pk_off + 22 * ccomps * dcomps);

            auto g_x_xyyyyzz = cbuffer.data(pk_off + 23 * ccomps * dcomps);

            auto g_x_xyyyzzz = cbuffer.data(pk_off + 24 * ccomps * dcomps);

            auto g_x_xyyzzzz = cbuffer.data(pk_off + 25 * ccomps * dcomps);

            auto g_x_xyzzzzz = cbuffer.data(pk_off + 26 * ccomps * dcomps);

            auto g_x_xzzzzzz = cbuffer.data(pk_off + 27 * ccomps * dcomps);

            auto g_x_yyyyyyy = cbuffer.data(pk_off + 28 * ccomps * dcomps);

            auto g_x_yyyyyyz = cbuffer.data(pk_off + 29 * ccomps * dcomps);

            auto g_x_yyyyyzz = cbuffer.data(pk_off + 30 * ccomps * dcomps);

            auto g_x_yyyyzzz = cbuffer.data(pk_off + 31 * ccomps * dcomps);

            auto g_x_yyyzzzz = cbuffer.data(pk_off + 32 * ccomps * dcomps);

            auto g_x_yyzzzzz = cbuffer.data(pk_off + 33 * ccomps * dcomps);

            auto g_x_yzzzzzz = cbuffer.data(pk_off + 34 * ccomps * dcomps);

            auto g_x_zzzzzzz = cbuffer.data(pk_off + 35 * ccomps * dcomps);

            auto g_y_xxxxxxx = cbuffer.data(pk_off + 36 * ccomps * dcomps);

            auto g_y_xxxxxxy = cbuffer.data(pk_off + 37 * ccomps * dcomps);

            auto g_y_xxxxxxz = cbuffer.data(pk_off + 38 * ccomps * dcomps);

            auto g_y_xxxxxyy = cbuffer.data(pk_off + 39 * ccomps * dcomps);

            auto g_y_xxxxxyz = cbuffer.data(pk_off + 40 * ccomps * dcomps);

            auto g_y_xxxxxzz = cbuffer.data(pk_off + 41 * ccomps * dcomps);

            auto g_y_xxxxyyy = cbuffer.data(pk_off + 42 * ccomps * dcomps);

            auto g_y_xxxxyyz = cbuffer.data(pk_off + 43 * ccomps * dcomps);

            auto g_y_xxxxyzz = cbuffer.data(pk_off + 44 * ccomps * dcomps);

            auto g_y_xxxxzzz = cbuffer.data(pk_off + 45 * ccomps * dcomps);

            auto g_y_xxxyyyy = cbuffer.data(pk_off + 46 * ccomps * dcomps);

            auto g_y_xxxyyyz = cbuffer.data(pk_off + 47 * ccomps * dcomps);

            auto g_y_xxxyyzz = cbuffer.data(pk_off + 48 * ccomps * dcomps);

            auto g_y_xxxyzzz = cbuffer.data(pk_off + 49 * ccomps * dcomps);

            auto g_y_xxxzzzz = cbuffer.data(pk_off + 50 * ccomps * dcomps);

            auto g_y_xxyyyyy = cbuffer.data(pk_off + 51 * ccomps * dcomps);

            auto g_y_xxyyyyz = cbuffer.data(pk_off + 52 * ccomps * dcomps);

            auto g_y_xxyyyzz = cbuffer.data(pk_off + 53 * ccomps * dcomps);

            auto g_y_xxyyzzz = cbuffer.data(pk_off + 54 * ccomps * dcomps);

            auto g_y_xxyzzzz = cbuffer.data(pk_off + 55 * ccomps * dcomps);

            auto g_y_xxzzzzz = cbuffer.data(pk_off + 56 * ccomps * dcomps);

            auto g_y_xyyyyyy = cbuffer.data(pk_off + 57 * ccomps * dcomps);

            auto g_y_xyyyyyz = cbuffer.data(pk_off + 58 * ccomps * dcomps);

            auto g_y_xyyyyzz = cbuffer.data(pk_off + 59 * ccomps * dcomps);

            auto g_y_xyyyzzz = cbuffer.data(pk_off + 60 * ccomps * dcomps);

            auto g_y_xyyzzzz = cbuffer.data(pk_off + 61 * ccomps * dcomps);

            auto g_y_xyzzzzz = cbuffer.data(pk_off + 62 * ccomps * dcomps);

            auto g_y_xzzzzzz = cbuffer.data(pk_off + 63 * ccomps * dcomps);

            auto g_y_yyyyyyy = cbuffer.data(pk_off + 64 * ccomps * dcomps);

            auto g_y_yyyyyyz = cbuffer.data(pk_off + 65 * ccomps * dcomps);

            auto g_y_yyyyyzz = cbuffer.data(pk_off + 66 * ccomps * dcomps);

            auto g_y_yyyyzzz = cbuffer.data(pk_off + 67 * ccomps * dcomps);

            auto g_y_yyyzzzz = cbuffer.data(pk_off + 68 * ccomps * dcomps);

            auto g_y_yyzzzzz = cbuffer.data(pk_off + 69 * ccomps * dcomps);

            auto g_y_yzzzzzz = cbuffer.data(pk_off + 70 * ccomps * dcomps);

            auto g_y_zzzzzzz = cbuffer.data(pk_off + 71 * ccomps * dcomps);

            auto g_z_xxxxxxx = cbuffer.data(pk_off + 72 * ccomps * dcomps);

            auto g_z_xxxxxxy = cbuffer.data(pk_off + 73 * ccomps * dcomps);

            auto g_z_xxxxxxz = cbuffer.data(pk_off + 74 * ccomps * dcomps);

            auto g_z_xxxxxyy = cbuffer.data(pk_off + 75 * ccomps * dcomps);

            auto g_z_xxxxxyz = cbuffer.data(pk_off + 76 * ccomps * dcomps);

            auto g_z_xxxxxzz = cbuffer.data(pk_off + 77 * ccomps * dcomps);

            auto g_z_xxxxyyy = cbuffer.data(pk_off + 78 * ccomps * dcomps);

            auto g_z_xxxxyyz = cbuffer.data(pk_off + 79 * ccomps * dcomps);

            auto g_z_xxxxyzz = cbuffer.data(pk_off + 80 * ccomps * dcomps);

            auto g_z_xxxxzzz = cbuffer.data(pk_off + 81 * ccomps * dcomps);

            auto g_z_xxxyyyy = cbuffer.data(pk_off + 82 * ccomps * dcomps);

            auto g_z_xxxyyyz = cbuffer.data(pk_off + 83 * ccomps * dcomps);

            auto g_z_xxxyyzz = cbuffer.data(pk_off + 84 * ccomps * dcomps);

            auto g_z_xxxyzzz = cbuffer.data(pk_off + 85 * ccomps * dcomps);

            auto g_z_xxxzzzz = cbuffer.data(pk_off + 86 * ccomps * dcomps);

            auto g_z_xxyyyyy = cbuffer.data(pk_off + 87 * ccomps * dcomps);

            auto g_z_xxyyyyz = cbuffer.data(pk_off + 88 * ccomps * dcomps);

            auto g_z_xxyyyzz = cbuffer.data(pk_off + 89 * ccomps * dcomps);

            auto g_z_xxyyzzz = cbuffer.data(pk_off + 90 * ccomps * dcomps);

            auto g_z_xxyzzzz = cbuffer.data(pk_off + 91 * ccomps * dcomps);

            auto g_z_xxzzzzz = cbuffer.data(pk_off + 92 * ccomps * dcomps);

            auto g_z_xyyyyyy = cbuffer.data(pk_off + 93 * ccomps * dcomps);

            auto g_z_xyyyyyz = cbuffer.data(pk_off + 94 * ccomps * dcomps);

            auto g_z_xyyyyzz = cbuffer.data(pk_off + 95 * ccomps * dcomps);

            auto g_z_xyyyzzz = cbuffer.data(pk_off + 96 * ccomps * dcomps);

            auto g_z_xyyzzzz = cbuffer.data(pk_off + 97 * ccomps * dcomps);

            auto g_z_xyzzzzz = cbuffer.data(pk_off + 98 * ccomps * dcomps);

            auto g_z_xzzzzzz = cbuffer.data(pk_off + 99 * ccomps * dcomps);

            auto g_z_yyyyyyy = cbuffer.data(pk_off + 100 * ccomps * dcomps);

            auto g_z_yyyyyyz = cbuffer.data(pk_off + 101 * ccomps * dcomps);

            auto g_z_yyyyyzz = cbuffer.data(pk_off + 102 * ccomps * dcomps);

            auto g_z_yyyyzzz = cbuffer.data(pk_off + 103 * ccomps * dcomps);

            auto g_z_yyyzzzz = cbuffer.data(pk_off + 104 * ccomps * dcomps);

            auto g_z_yyzzzzz = cbuffer.data(pk_off + 105 * ccomps * dcomps);

            auto g_z_yzzzzzz = cbuffer.data(pk_off + 106 * ccomps * dcomps);

            auto g_z_zzzzzzz = cbuffer.data(pk_off + 107 * ccomps * dcomps);

            /// Set up components of auxilary buffer : PKSS

            const auto pk_geom_01_off = idx_geom_01_pkxx + i * dcomps + j;

            auto g_0_x_x_xxxxxxx = cbuffer.data(pk_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_x_xxxxxxy = cbuffer.data(pk_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_x_xxxxxxz = cbuffer.data(pk_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_x_xxxxxyy = cbuffer.data(pk_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_x_xxxxxyz = cbuffer.data(pk_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_x_xxxxxzz = cbuffer.data(pk_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_x_xxxxyyy = cbuffer.data(pk_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_x_xxxxyyz = cbuffer.data(pk_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_x_xxxxyzz = cbuffer.data(pk_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_x_xxxxzzz = cbuffer.data(pk_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_x_xxxyyyy = cbuffer.data(pk_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_x_xxxyyyz = cbuffer.data(pk_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_x_xxxyyzz = cbuffer.data(pk_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_x_xxxyzzz = cbuffer.data(pk_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_x_xxxzzzz = cbuffer.data(pk_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_x_xxyyyyy = cbuffer.data(pk_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_x_xxyyyyz = cbuffer.data(pk_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_x_xxyyyzz = cbuffer.data(pk_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_x_xxyyzzz = cbuffer.data(pk_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_x_xxyzzzz = cbuffer.data(pk_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_x_xxzzzzz = cbuffer.data(pk_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_x_xyyyyyy = cbuffer.data(pk_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_x_xyyyyyz = cbuffer.data(pk_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_x_xyyyyzz = cbuffer.data(pk_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_x_xyyyzzz = cbuffer.data(pk_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_x_xyyzzzz = cbuffer.data(pk_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_x_xyzzzzz = cbuffer.data(pk_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_x_xzzzzzz = cbuffer.data(pk_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_x_yyyyyyy = cbuffer.data(pk_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_x_yyyyyyz = cbuffer.data(pk_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_x_yyyyyzz = cbuffer.data(pk_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_x_yyyyzzz = cbuffer.data(pk_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_x_yyyzzzz = cbuffer.data(pk_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_x_yyzzzzz = cbuffer.data(pk_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_x_yzzzzzz = cbuffer.data(pk_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_x_zzzzzzz = cbuffer.data(pk_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_y_xxxxxxx = cbuffer.data(pk_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_y_xxxxxxy = cbuffer.data(pk_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_y_xxxxxxz = cbuffer.data(pk_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_y_xxxxxyy = cbuffer.data(pk_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_y_xxxxxyz = cbuffer.data(pk_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_y_xxxxxzz = cbuffer.data(pk_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_y_xxxxyyy = cbuffer.data(pk_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_y_xxxxyyz = cbuffer.data(pk_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_y_xxxxyzz = cbuffer.data(pk_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_y_xxxxzzz = cbuffer.data(pk_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_y_xxxyyyy = cbuffer.data(pk_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_y_xxxyyyz = cbuffer.data(pk_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_y_xxxyyzz = cbuffer.data(pk_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_y_xxxyzzz = cbuffer.data(pk_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_y_xxxzzzz = cbuffer.data(pk_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_y_xxyyyyy = cbuffer.data(pk_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_y_xxyyyyz = cbuffer.data(pk_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_y_xxyyyzz = cbuffer.data(pk_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_y_xxyyzzz = cbuffer.data(pk_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_y_xxyzzzz = cbuffer.data(pk_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_y_xxzzzzz = cbuffer.data(pk_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_y_xyyyyyy = cbuffer.data(pk_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_y_xyyyyyz = cbuffer.data(pk_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_y_xyyyyzz = cbuffer.data(pk_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_y_xyyyzzz = cbuffer.data(pk_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_y_xyyzzzz = cbuffer.data(pk_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_y_xyzzzzz = cbuffer.data(pk_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_y_xzzzzzz = cbuffer.data(pk_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_y_yyyyyyy = cbuffer.data(pk_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_y_yyyyyyz = cbuffer.data(pk_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_y_yyyyyzz = cbuffer.data(pk_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_y_yyyyzzz = cbuffer.data(pk_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_y_yyyzzzz = cbuffer.data(pk_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_y_yyzzzzz = cbuffer.data(pk_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_y_yzzzzzz = cbuffer.data(pk_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_y_zzzzzzz = cbuffer.data(pk_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_z_xxxxxxx = cbuffer.data(pk_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_z_xxxxxxy = cbuffer.data(pk_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_z_xxxxxxz = cbuffer.data(pk_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_z_xxxxxyy = cbuffer.data(pk_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_z_xxxxxyz = cbuffer.data(pk_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_z_xxxxxzz = cbuffer.data(pk_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_z_xxxxyyy = cbuffer.data(pk_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_z_xxxxyyz = cbuffer.data(pk_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_z_xxxxyzz = cbuffer.data(pk_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_z_xxxxzzz = cbuffer.data(pk_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_z_xxxyyyy = cbuffer.data(pk_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_z_xxxyyyz = cbuffer.data(pk_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_z_xxxyyzz = cbuffer.data(pk_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_z_xxxyzzz = cbuffer.data(pk_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_z_xxxzzzz = cbuffer.data(pk_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_z_xxyyyyy = cbuffer.data(pk_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_z_xxyyyyz = cbuffer.data(pk_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_z_xxyyyzz = cbuffer.data(pk_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_z_xxyyzzz = cbuffer.data(pk_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_z_xxyzzzz = cbuffer.data(pk_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_z_xxzzzzz = cbuffer.data(pk_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_z_xyyyyyy = cbuffer.data(pk_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_z_xyyyyyz = cbuffer.data(pk_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_z_xyyyyzz = cbuffer.data(pk_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_z_xyyyzzz = cbuffer.data(pk_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_z_xyyzzzz = cbuffer.data(pk_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_z_xyzzzzz = cbuffer.data(pk_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_z_xzzzzzz = cbuffer.data(pk_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_z_yyyyyyy = cbuffer.data(pk_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_z_yyyyyyz = cbuffer.data(pk_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_z_yyyyyzz = cbuffer.data(pk_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_z_yyyyzzz = cbuffer.data(pk_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_z_yyyzzzz = cbuffer.data(pk_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_z_yyzzzzz = cbuffer.data(pk_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_z_yzzzzzz = cbuffer.data(pk_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_z_zzzzzzz = cbuffer.data(pk_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_y_x_xxxxxxx = cbuffer.data(pk_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_y_x_xxxxxxy = cbuffer.data(pk_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_y_x_xxxxxxz = cbuffer.data(pk_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_y_x_xxxxxyy = cbuffer.data(pk_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_y_x_xxxxxyz = cbuffer.data(pk_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_y_x_xxxxxzz = cbuffer.data(pk_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_y_x_xxxxyyy = cbuffer.data(pk_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_y_x_xxxxyyz = cbuffer.data(pk_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_y_x_xxxxyzz = cbuffer.data(pk_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_y_x_xxxxzzz = cbuffer.data(pk_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_y_x_xxxyyyy = cbuffer.data(pk_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_y_x_xxxyyyz = cbuffer.data(pk_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_y_x_xxxyyzz = cbuffer.data(pk_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_y_x_xxxyzzz = cbuffer.data(pk_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_y_x_xxxzzzz = cbuffer.data(pk_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_y_x_xxyyyyy = cbuffer.data(pk_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_y_x_xxyyyyz = cbuffer.data(pk_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_y_x_xxyyyzz = cbuffer.data(pk_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_y_x_xxyyzzz = cbuffer.data(pk_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_y_x_xxyzzzz = cbuffer.data(pk_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_y_x_xxzzzzz = cbuffer.data(pk_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_y_x_xyyyyyy = cbuffer.data(pk_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_y_x_xyyyyyz = cbuffer.data(pk_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_y_x_xyyyyzz = cbuffer.data(pk_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_y_x_xyyyzzz = cbuffer.data(pk_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_y_x_xyyzzzz = cbuffer.data(pk_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_y_x_xyzzzzz = cbuffer.data(pk_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_y_x_xzzzzzz = cbuffer.data(pk_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_y_x_yyyyyyy = cbuffer.data(pk_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_y_x_yyyyyyz = cbuffer.data(pk_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_y_x_yyyyyzz = cbuffer.data(pk_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_y_x_yyyyzzz = cbuffer.data(pk_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_y_x_yyyzzzz = cbuffer.data(pk_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_y_x_yyzzzzz = cbuffer.data(pk_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_y_x_yzzzzzz = cbuffer.data(pk_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_y_x_zzzzzzz = cbuffer.data(pk_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_y_y_xxxxxxx = cbuffer.data(pk_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_y_y_xxxxxxy = cbuffer.data(pk_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_y_y_xxxxxxz = cbuffer.data(pk_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_y_y_xxxxxyy = cbuffer.data(pk_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_y_y_xxxxxyz = cbuffer.data(pk_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_y_y_xxxxxzz = cbuffer.data(pk_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_y_y_xxxxyyy = cbuffer.data(pk_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_y_y_xxxxyyz = cbuffer.data(pk_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_y_y_xxxxyzz = cbuffer.data(pk_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_y_y_xxxxzzz = cbuffer.data(pk_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_y_y_xxxyyyy = cbuffer.data(pk_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_y_y_xxxyyyz = cbuffer.data(pk_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_y_y_xxxyyzz = cbuffer.data(pk_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_y_y_xxxyzzz = cbuffer.data(pk_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_y_y_xxxzzzz = cbuffer.data(pk_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_y_y_xxyyyyy = cbuffer.data(pk_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_y_y_xxyyyyz = cbuffer.data(pk_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_y_y_xxyyyzz = cbuffer.data(pk_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_y_y_xxyyzzz = cbuffer.data(pk_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_y_y_xxyzzzz = cbuffer.data(pk_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_y_y_xxzzzzz = cbuffer.data(pk_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_y_y_xyyyyyy = cbuffer.data(pk_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_y_y_xyyyyyz = cbuffer.data(pk_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_y_y_xyyyyzz = cbuffer.data(pk_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_y_y_xyyyzzz = cbuffer.data(pk_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_y_y_xyyzzzz = cbuffer.data(pk_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_y_y_xyzzzzz = cbuffer.data(pk_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_y_y_xzzzzzz = cbuffer.data(pk_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_y_y_yyyyyyy = cbuffer.data(pk_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_y_y_yyyyyyz = cbuffer.data(pk_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_y_y_yyyyyzz = cbuffer.data(pk_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_y_y_yyyyzzz = cbuffer.data(pk_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_y_y_yyyzzzz = cbuffer.data(pk_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_y_y_yyzzzzz = cbuffer.data(pk_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_y_y_yzzzzzz = cbuffer.data(pk_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_y_y_zzzzzzz = cbuffer.data(pk_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_y_z_xxxxxxx = cbuffer.data(pk_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_y_z_xxxxxxy = cbuffer.data(pk_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_y_z_xxxxxxz = cbuffer.data(pk_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_y_z_xxxxxyy = cbuffer.data(pk_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_y_z_xxxxxyz = cbuffer.data(pk_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_y_z_xxxxxzz = cbuffer.data(pk_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_y_z_xxxxyyy = cbuffer.data(pk_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_y_z_xxxxyyz = cbuffer.data(pk_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_y_z_xxxxyzz = cbuffer.data(pk_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_y_z_xxxxzzz = cbuffer.data(pk_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_y_z_xxxyyyy = cbuffer.data(pk_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_y_z_xxxyyyz = cbuffer.data(pk_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_y_z_xxxyyzz = cbuffer.data(pk_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_y_z_xxxyzzz = cbuffer.data(pk_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_y_z_xxxzzzz = cbuffer.data(pk_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_y_z_xxyyyyy = cbuffer.data(pk_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_y_z_xxyyyyz = cbuffer.data(pk_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_y_z_xxyyyzz = cbuffer.data(pk_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_y_z_xxyyzzz = cbuffer.data(pk_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_y_z_xxyzzzz = cbuffer.data(pk_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_y_z_xxzzzzz = cbuffer.data(pk_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_y_z_xyyyyyy = cbuffer.data(pk_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_y_z_xyyyyyz = cbuffer.data(pk_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_y_z_xyyyyzz = cbuffer.data(pk_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_y_z_xyyyzzz = cbuffer.data(pk_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_y_z_xyyzzzz = cbuffer.data(pk_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_y_z_xyzzzzz = cbuffer.data(pk_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_y_z_xzzzzzz = cbuffer.data(pk_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_y_z_yyyyyyy = cbuffer.data(pk_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_y_z_yyyyyyz = cbuffer.data(pk_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_y_z_yyyyyzz = cbuffer.data(pk_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_y_z_yyyyzzz = cbuffer.data(pk_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_y_z_yyyzzzz = cbuffer.data(pk_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_y_z_yyzzzzz = cbuffer.data(pk_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_y_z_yzzzzzz = cbuffer.data(pk_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_y_z_zzzzzzz = cbuffer.data(pk_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_z_x_xxxxxxx = cbuffer.data(pk_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_z_x_xxxxxxy = cbuffer.data(pk_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_z_x_xxxxxxz = cbuffer.data(pk_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_z_x_xxxxxyy = cbuffer.data(pk_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_z_x_xxxxxyz = cbuffer.data(pk_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_z_x_xxxxxzz = cbuffer.data(pk_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_z_x_xxxxyyy = cbuffer.data(pk_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_z_x_xxxxyyz = cbuffer.data(pk_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_z_x_xxxxyzz = cbuffer.data(pk_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_z_x_xxxxzzz = cbuffer.data(pk_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_z_x_xxxyyyy = cbuffer.data(pk_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_z_x_xxxyyyz = cbuffer.data(pk_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_z_x_xxxyyzz = cbuffer.data(pk_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_z_x_xxxyzzz = cbuffer.data(pk_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_z_x_xxxzzzz = cbuffer.data(pk_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_z_x_xxyyyyy = cbuffer.data(pk_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_z_x_xxyyyyz = cbuffer.data(pk_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_z_x_xxyyyzz = cbuffer.data(pk_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_z_x_xxyyzzz = cbuffer.data(pk_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_z_x_xxyzzzz = cbuffer.data(pk_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_z_x_xxzzzzz = cbuffer.data(pk_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_z_x_xyyyyyy = cbuffer.data(pk_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_z_x_xyyyyyz = cbuffer.data(pk_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_z_x_xyyyyzz = cbuffer.data(pk_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_z_x_xyyyzzz = cbuffer.data(pk_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_z_x_xyyzzzz = cbuffer.data(pk_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_z_x_xyzzzzz = cbuffer.data(pk_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_z_x_xzzzzzz = cbuffer.data(pk_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_z_x_yyyyyyy = cbuffer.data(pk_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_z_x_yyyyyyz = cbuffer.data(pk_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_z_x_yyyyyzz = cbuffer.data(pk_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_z_x_yyyyzzz = cbuffer.data(pk_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_z_x_yyyzzzz = cbuffer.data(pk_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_z_x_yyzzzzz = cbuffer.data(pk_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_z_x_yzzzzzz = cbuffer.data(pk_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_z_x_zzzzzzz = cbuffer.data(pk_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_z_y_xxxxxxx = cbuffer.data(pk_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_z_y_xxxxxxy = cbuffer.data(pk_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_z_y_xxxxxxz = cbuffer.data(pk_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_z_y_xxxxxyy = cbuffer.data(pk_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_z_y_xxxxxyz = cbuffer.data(pk_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_z_y_xxxxxzz = cbuffer.data(pk_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_z_y_xxxxyyy = cbuffer.data(pk_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_z_y_xxxxyyz = cbuffer.data(pk_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_z_y_xxxxyzz = cbuffer.data(pk_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_z_y_xxxxzzz = cbuffer.data(pk_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_z_y_xxxyyyy = cbuffer.data(pk_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_z_y_xxxyyyz = cbuffer.data(pk_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_z_y_xxxyyzz = cbuffer.data(pk_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_z_y_xxxyzzz = cbuffer.data(pk_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_z_y_xxxzzzz = cbuffer.data(pk_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_z_y_xxyyyyy = cbuffer.data(pk_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_z_y_xxyyyyz = cbuffer.data(pk_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_z_y_xxyyyzz = cbuffer.data(pk_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_z_y_xxyyzzz = cbuffer.data(pk_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_z_y_xxyzzzz = cbuffer.data(pk_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_z_y_xxzzzzz = cbuffer.data(pk_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_z_y_xyyyyyy = cbuffer.data(pk_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_z_y_xyyyyyz = cbuffer.data(pk_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_z_y_xyyyyzz = cbuffer.data(pk_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_z_y_xyyyzzz = cbuffer.data(pk_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_z_y_xyyzzzz = cbuffer.data(pk_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_z_y_xyzzzzz = cbuffer.data(pk_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_z_y_xzzzzzz = cbuffer.data(pk_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_z_y_yyyyyyy = cbuffer.data(pk_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_z_y_yyyyyyz = cbuffer.data(pk_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_z_y_yyyyyzz = cbuffer.data(pk_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_z_y_yyyyzzz = cbuffer.data(pk_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_z_y_yyyzzzz = cbuffer.data(pk_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_z_y_yyzzzzz = cbuffer.data(pk_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_z_y_yzzzzzz = cbuffer.data(pk_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_z_y_zzzzzzz = cbuffer.data(pk_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_z_z_xxxxxxx = cbuffer.data(pk_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_z_z_xxxxxxy = cbuffer.data(pk_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_z_z_xxxxxxz = cbuffer.data(pk_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_z_z_xxxxxyy = cbuffer.data(pk_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_z_z_xxxxxyz = cbuffer.data(pk_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_z_z_xxxxxzz = cbuffer.data(pk_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_z_z_xxxxyyy = cbuffer.data(pk_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_z_z_xxxxyyz = cbuffer.data(pk_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_z_z_xxxxyzz = cbuffer.data(pk_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_z_z_xxxxzzz = cbuffer.data(pk_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_z_z_xxxyyyy = cbuffer.data(pk_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_z_z_xxxyyyz = cbuffer.data(pk_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_z_z_xxxyyzz = cbuffer.data(pk_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_z_z_xxxyzzz = cbuffer.data(pk_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_z_z_xxxzzzz = cbuffer.data(pk_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_z_z_xxyyyyy = cbuffer.data(pk_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_z_z_xxyyyyz = cbuffer.data(pk_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_z_z_xxyyyzz = cbuffer.data(pk_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_z_z_xxyyzzz = cbuffer.data(pk_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_z_z_xxyzzzz = cbuffer.data(pk_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_z_z_xxzzzzz = cbuffer.data(pk_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_z_z_xyyyyyy = cbuffer.data(pk_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_z_z_xyyyyyz = cbuffer.data(pk_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_z_z_xyyyyzz = cbuffer.data(pk_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_z_z_xyyyzzz = cbuffer.data(pk_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_z_z_xyyzzzz = cbuffer.data(pk_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_z_z_xyzzzzz = cbuffer.data(pk_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_z_z_xzzzzzz = cbuffer.data(pk_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_z_z_yyyyyyy = cbuffer.data(pk_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_z_z_yyyyyyz = cbuffer.data(pk_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_z_z_yyyyyzz = cbuffer.data(pk_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_z_z_yyyyzzz = cbuffer.data(pk_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_z_z_yyyzzzz = cbuffer.data(pk_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_z_z_yyzzzzz = cbuffer.data(pk_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_z_z_yzzzzzz = cbuffer.data(pk_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_z_z_zzzzzzz = cbuffer.data(pk_geom_01_off + 323 * ccomps * dcomps);

            /// Set up components of auxilary buffer : PLSS

            const auto pl_geom_01_off = idx_geom_01_plxx + i * dcomps + j;

            auto g_0_x_x_xxxxxxxx = cbuffer.data(pl_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_x_xxxxxxxy = cbuffer.data(pl_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_x_xxxxxxxz = cbuffer.data(pl_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_x_xxxxxxyy = cbuffer.data(pl_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_x_xxxxxxyz = cbuffer.data(pl_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_x_xxxxxxzz = cbuffer.data(pl_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_x_xxxxxyyy = cbuffer.data(pl_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_x_xxxxxyyz = cbuffer.data(pl_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_x_xxxxxyzz = cbuffer.data(pl_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_x_xxxxxzzz = cbuffer.data(pl_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_x_xxxxyyyy = cbuffer.data(pl_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_x_xxxxyyyz = cbuffer.data(pl_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_x_xxxxyyzz = cbuffer.data(pl_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_x_xxxxyzzz = cbuffer.data(pl_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_x_xxxxzzzz = cbuffer.data(pl_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_x_xxxyyyyy = cbuffer.data(pl_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_x_xxxyyyyz = cbuffer.data(pl_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_x_xxxyyyzz = cbuffer.data(pl_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_x_xxxyyzzz = cbuffer.data(pl_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_x_xxxyzzzz = cbuffer.data(pl_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_x_xxxzzzzz = cbuffer.data(pl_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_x_xxyyyyyy = cbuffer.data(pl_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_x_xxyyyyyz = cbuffer.data(pl_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_x_xxyyyyzz = cbuffer.data(pl_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_x_xxyyyzzz = cbuffer.data(pl_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_x_xxyyzzzz = cbuffer.data(pl_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_x_xxyzzzzz = cbuffer.data(pl_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_x_xxzzzzzz = cbuffer.data(pl_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_x_xyyyyyyy = cbuffer.data(pl_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_x_xyyyyyyz = cbuffer.data(pl_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_x_xyyyyyzz = cbuffer.data(pl_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_x_xyyyyzzz = cbuffer.data(pl_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_x_xyyyzzzz = cbuffer.data(pl_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_x_xyyzzzzz = cbuffer.data(pl_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_x_xyzzzzzz = cbuffer.data(pl_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_x_xzzzzzzz = cbuffer.data(pl_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_x_yyyyyyyy = cbuffer.data(pl_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_x_yyyyyyyz = cbuffer.data(pl_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_x_yyyyyyzz = cbuffer.data(pl_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_x_yyyyyzzz = cbuffer.data(pl_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_x_yyyyzzzz = cbuffer.data(pl_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_x_yyyzzzzz = cbuffer.data(pl_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_x_yyzzzzzz = cbuffer.data(pl_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_x_yzzzzzzz = cbuffer.data(pl_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_x_zzzzzzzz = cbuffer.data(pl_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_y_xxxxxxxx = cbuffer.data(pl_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_y_xxxxxxxy = cbuffer.data(pl_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_y_xxxxxxxz = cbuffer.data(pl_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_y_xxxxxxyy = cbuffer.data(pl_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_y_xxxxxxyz = cbuffer.data(pl_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_y_xxxxxxzz = cbuffer.data(pl_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_y_xxxxxyyy = cbuffer.data(pl_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_y_xxxxxyyz = cbuffer.data(pl_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_y_xxxxxyzz = cbuffer.data(pl_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_y_xxxxxzzz = cbuffer.data(pl_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_y_xxxxyyyy = cbuffer.data(pl_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_y_xxxxyyyz = cbuffer.data(pl_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_y_xxxxyyzz = cbuffer.data(pl_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_y_xxxxyzzz = cbuffer.data(pl_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_y_xxxxzzzz = cbuffer.data(pl_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_y_xxxyyyyy = cbuffer.data(pl_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_y_xxxyyyyz = cbuffer.data(pl_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_y_xxxyyyzz = cbuffer.data(pl_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_y_xxxyyzzz = cbuffer.data(pl_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_y_xxxyzzzz = cbuffer.data(pl_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_y_xxxzzzzz = cbuffer.data(pl_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_y_xxyyyyyy = cbuffer.data(pl_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_y_xxyyyyyz = cbuffer.data(pl_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_y_xxyyyyzz = cbuffer.data(pl_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_y_xxyyyzzz = cbuffer.data(pl_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_y_xxyyzzzz = cbuffer.data(pl_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_y_xxyzzzzz = cbuffer.data(pl_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_y_xxzzzzzz = cbuffer.data(pl_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_y_xyyyyyyy = cbuffer.data(pl_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_y_xyyyyyyz = cbuffer.data(pl_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_y_xyyyyyzz = cbuffer.data(pl_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_y_xyyyyzzz = cbuffer.data(pl_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_y_xyyyzzzz = cbuffer.data(pl_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_y_xyyzzzzz = cbuffer.data(pl_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_y_xyzzzzzz = cbuffer.data(pl_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_y_xzzzzzzz = cbuffer.data(pl_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_y_yyyyyyyy = cbuffer.data(pl_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_y_yyyyyyyz = cbuffer.data(pl_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_y_yyyyyyzz = cbuffer.data(pl_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_y_yyyyyzzz = cbuffer.data(pl_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_y_yyyyzzzz = cbuffer.data(pl_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_y_yyyzzzzz = cbuffer.data(pl_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_y_yyzzzzzz = cbuffer.data(pl_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_y_yzzzzzzz = cbuffer.data(pl_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_y_zzzzzzzz = cbuffer.data(pl_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_z_xxxxxxxx = cbuffer.data(pl_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_z_xxxxxxxy = cbuffer.data(pl_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_z_xxxxxxxz = cbuffer.data(pl_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_z_xxxxxxyy = cbuffer.data(pl_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_z_xxxxxxyz = cbuffer.data(pl_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_z_xxxxxxzz = cbuffer.data(pl_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_z_xxxxxyyy = cbuffer.data(pl_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_z_xxxxxyyz = cbuffer.data(pl_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_z_xxxxxyzz = cbuffer.data(pl_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_z_xxxxxzzz = cbuffer.data(pl_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_z_xxxxyyyy = cbuffer.data(pl_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_z_xxxxyyyz = cbuffer.data(pl_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_z_xxxxyyzz = cbuffer.data(pl_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_z_xxxxyzzz = cbuffer.data(pl_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_z_xxxxzzzz = cbuffer.data(pl_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_z_xxxyyyyy = cbuffer.data(pl_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_z_xxxyyyyz = cbuffer.data(pl_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_z_xxxyyyzz = cbuffer.data(pl_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_z_xxxyyzzz = cbuffer.data(pl_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_z_xxxyzzzz = cbuffer.data(pl_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_z_xxxzzzzz = cbuffer.data(pl_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_z_xxyyyyyy = cbuffer.data(pl_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_z_xxyyyyyz = cbuffer.data(pl_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_z_xxyyyyzz = cbuffer.data(pl_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_z_xxyyyzzz = cbuffer.data(pl_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_z_xxyyzzzz = cbuffer.data(pl_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_z_xxyzzzzz = cbuffer.data(pl_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_z_xxzzzzzz = cbuffer.data(pl_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_z_xyyyyyyy = cbuffer.data(pl_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_z_xyyyyyyz = cbuffer.data(pl_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_x_z_xyyyyyzz = cbuffer.data(pl_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_z_xyyyyzzz = cbuffer.data(pl_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_z_xyyyzzzz = cbuffer.data(pl_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_z_xyyzzzzz = cbuffer.data(pl_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_z_xyzzzzzz = cbuffer.data(pl_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_z_xzzzzzzz = cbuffer.data(pl_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_x_z_yyyyyyyy = cbuffer.data(pl_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_z_yyyyyyyz = cbuffer.data(pl_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_z_yyyyyyzz = cbuffer.data(pl_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_z_yyyyyzzz = cbuffer.data(pl_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_z_yyyyzzzz = cbuffer.data(pl_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_z_yyyzzzzz = cbuffer.data(pl_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_z_yyzzzzzz = cbuffer.data(pl_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_z_yzzzzzzz = cbuffer.data(pl_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_z_zzzzzzzz = cbuffer.data(pl_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_y_x_xxxxxxxx = cbuffer.data(pl_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_y_x_xxxxxxxy = cbuffer.data(pl_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_y_x_xxxxxxxz = cbuffer.data(pl_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_y_x_xxxxxxyy = cbuffer.data(pl_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_y_x_xxxxxxyz = cbuffer.data(pl_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_y_x_xxxxxxzz = cbuffer.data(pl_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_y_x_xxxxxyyy = cbuffer.data(pl_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_y_x_xxxxxyyz = cbuffer.data(pl_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_y_x_xxxxxyzz = cbuffer.data(pl_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_y_x_xxxxxzzz = cbuffer.data(pl_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_y_x_xxxxyyyy = cbuffer.data(pl_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_y_x_xxxxyyyz = cbuffer.data(pl_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_y_x_xxxxyyzz = cbuffer.data(pl_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_y_x_xxxxyzzz = cbuffer.data(pl_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_y_x_xxxxzzzz = cbuffer.data(pl_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_y_x_xxxyyyyy = cbuffer.data(pl_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_y_x_xxxyyyyz = cbuffer.data(pl_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_y_x_xxxyyyzz = cbuffer.data(pl_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_y_x_xxxyyzzz = cbuffer.data(pl_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_y_x_xxxyzzzz = cbuffer.data(pl_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_y_x_xxxzzzzz = cbuffer.data(pl_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_y_x_xxyyyyyy = cbuffer.data(pl_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_y_x_xxyyyyyz = cbuffer.data(pl_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_y_x_xxyyyyzz = cbuffer.data(pl_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_y_x_xxyyyzzz = cbuffer.data(pl_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_y_x_xxyyzzzz = cbuffer.data(pl_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_y_x_xxyzzzzz = cbuffer.data(pl_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_y_x_xxzzzzzz = cbuffer.data(pl_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_y_x_xyyyyyyy = cbuffer.data(pl_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_y_x_xyyyyyyz = cbuffer.data(pl_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_y_x_xyyyyyzz = cbuffer.data(pl_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_y_x_xyyyyzzz = cbuffer.data(pl_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_y_x_xyyyzzzz = cbuffer.data(pl_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_y_x_xyyzzzzz = cbuffer.data(pl_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_y_x_xyzzzzzz = cbuffer.data(pl_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_y_x_xzzzzzzz = cbuffer.data(pl_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_y_x_yyyyyyyy = cbuffer.data(pl_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_y_x_yyyyyyyz = cbuffer.data(pl_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_y_x_yyyyyyzz = cbuffer.data(pl_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_y_x_yyyyyzzz = cbuffer.data(pl_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_y_x_yyyyzzzz = cbuffer.data(pl_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_y_x_yyyzzzzz = cbuffer.data(pl_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_y_x_yyzzzzzz = cbuffer.data(pl_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_y_x_yzzzzzzz = cbuffer.data(pl_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_y_x_zzzzzzzz = cbuffer.data(pl_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_y_y_xxxxxxxx = cbuffer.data(pl_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_y_y_xxxxxxxy = cbuffer.data(pl_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_y_y_xxxxxxxz = cbuffer.data(pl_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_y_y_xxxxxxyy = cbuffer.data(pl_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_y_y_xxxxxxyz = cbuffer.data(pl_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_y_y_xxxxxxzz = cbuffer.data(pl_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_y_y_xxxxxyyy = cbuffer.data(pl_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_y_y_xxxxxyyz = cbuffer.data(pl_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_y_y_xxxxxyzz = cbuffer.data(pl_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_y_y_xxxxxzzz = cbuffer.data(pl_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_y_y_xxxxyyyy = cbuffer.data(pl_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_y_y_xxxxyyyz = cbuffer.data(pl_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_y_y_xxxxyyzz = cbuffer.data(pl_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_y_y_xxxxyzzz = cbuffer.data(pl_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_y_y_xxxxzzzz = cbuffer.data(pl_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_y_y_xxxyyyyy = cbuffer.data(pl_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_y_y_xxxyyyyz = cbuffer.data(pl_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_y_y_xxxyyyzz = cbuffer.data(pl_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_y_y_xxxyyzzz = cbuffer.data(pl_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_y_y_xxxyzzzz = cbuffer.data(pl_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_y_y_xxxzzzzz = cbuffer.data(pl_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_y_y_xxyyyyyy = cbuffer.data(pl_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_y_y_xxyyyyyz = cbuffer.data(pl_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_y_y_xxyyyyzz = cbuffer.data(pl_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_y_y_xxyyyzzz = cbuffer.data(pl_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_y_y_xxyyzzzz = cbuffer.data(pl_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_y_y_xxyzzzzz = cbuffer.data(pl_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_y_y_xxzzzzzz = cbuffer.data(pl_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_y_y_xyyyyyyy = cbuffer.data(pl_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_y_y_xyyyyyyz = cbuffer.data(pl_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_y_y_xyyyyyzz = cbuffer.data(pl_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_y_y_xyyyyzzz = cbuffer.data(pl_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_y_y_xyyyzzzz = cbuffer.data(pl_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_y_y_xyyzzzzz = cbuffer.data(pl_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_y_y_xyzzzzzz = cbuffer.data(pl_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_y_y_xzzzzzzz = cbuffer.data(pl_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_y_y_yyyyyyyy = cbuffer.data(pl_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_y_y_yyyyyyyz = cbuffer.data(pl_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_y_y_yyyyyyzz = cbuffer.data(pl_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_y_y_yyyyyzzz = cbuffer.data(pl_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_y_y_yyyyzzzz = cbuffer.data(pl_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_y_y_yyyzzzzz = cbuffer.data(pl_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_y_y_yyzzzzzz = cbuffer.data(pl_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_y_y_yzzzzzzz = cbuffer.data(pl_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_y_y_zzzzzzzz = cbuffer.data(pl_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_y_z_xxxxxxxx = cbuffer.data(pl_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_y_z_xxxxxxxy = cbuffer.data(pl_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_y_z_xxxxxxxz = cbuffer.data(pl_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_y_z_xxxxxxyy = cbuffer.data(pl_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_y_z_xxxxxxyz = cbuffer.data(pl_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_y_z_xxxxxxzz = cbuffer.data(pl_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_y_z_xxxxxyyy = cbuffer.data(pl_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_y_z_xxxxxyyz = cbuffer.data(pl_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_y_z_xxxxxyzz = cbuffer.data(pl_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_y_z_xxxxxzzz = cbuffer.data(pl_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_y_z_xxxxyyyy = cbuffer.data(pl_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_y_z_xxxxyyyz = cbuffer.data(pl_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_y_z_xxxxyyzz = cbuffer.data(pl_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_y_z_xxxxyzzz = cbuffer.data(pl_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_y_z_xxxxzzzz = cbuffer.data(pl_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_y_z_xxxyyyyy = cbuffer.data(pl_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_y_z_xxxyyyyz = cbuffer.data(pl_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_y_z_xxxyyyzz = cbuffer.data(pl_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_y_z_xxxyyzzz = cbuffer.data(pl_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_y_z_xxxyzzzz = cbuffer.data(pl_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_y_z_xxxzzzzz = cbuffer.data(pl_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_y_z_xxyyyyyy = cbuffer.data(pl_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_y_z_xxyyyyyz = cbuffer.data(pl_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_y_z_xxyyyyzz = cbuffer.data(pl_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_y_z_xxyyyzzz = cbuffer.data(pl_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_y_z_xxyyzzzz = cbuffer.data(pl_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_y_z_xxyzzzzz = cbuffer.data(pl_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_y_z_xxzzzzzz = cbuffer.data(pl_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_y_z_xyyyyyyy = cbuffer.data(pl_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_y_z_xyyyyyyz = cbuffer.data(pl_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_y_z_xyyyyyzz = cbuffer.data(pl_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_y_z_xyyyyzzz = cbuffer.data(pl_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_y_z_xyyyzzzz = cbuffer.data(pl_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_y_z_xyyzzzzz = cbuffer.data(pl_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_y_z_xyzzzzzz = cbuffer.data(pl_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_y_z_xzzzzzzz = cbuffer.data(pl_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_y_z_yyyyyyyy = cbuffer.data(pl_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_y_z_yyyyyyyz = cbuffer.data(pl_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_y_z_yyyyyyzz = cbuffer.data(pl_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_y_z_yyyyyzzz = cbuffer.data(pl_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_y_z_yyyyzzzz = cbuffer.data(pl_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_y_z_yyyzzzzz = cbuffer.data(pl_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_y_z_yyzzzzzz = cbuffer.data(pl_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_y_z_yzzzzzzz = cbuffer.data(pl_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_y_z_zzzzzzzz = cbuffer.data(pl_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_z_x_xxxxxxxx = cbuffer.data(pl_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_z_x_xxxxxxxy = cbuffer.data(pl_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_z_x_xxxxxxxz = cbuffer.data(pl_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_z_x_xxxxxxyy = cbuffer.data(pl_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_z_x_xxxxxxyz = cbuffer.data(pl_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_z_x_xxxxxxzz = cbuffer.data(pl_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_z_x_xxxxxyyy = cbuffer.data(pl_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_z_x_xxxxxyyz = cbuffer.data(pl_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_z_x_xxxxxyzz = cbuffer.data(pl_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_z_x_xxxxxzzz = cbuffer.data(pl_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_z_x_xxxxyyyy = cbuffer.data(pl_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_z_x_xxxxyyyz = cbuffer.data(pl_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_z_x_xxxxyyzz = cbuffer.data(pl_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_z_x_xxxxyzzz = cbuffer.data(pl_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_z_x_xxxxzzzz = cbuffer.data(pl_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_z_x_xxxyyyyy = cbuffer.data(pl_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_z_x_xxxyyyyz = cbuffer.data(pl_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_z_x_xxxyyyzz = cbuffer.data(pl_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_z_x_xxxyyzzz = cbuffer.data(pl_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_z_x_xxxyzzzz = cbuffer.data(pl_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_z_x_xxxzzzzz = cbuffer.data(pl_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_z_x_xxyyyyyy = cbuffer.data(pl_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_z_x_xxyyyyyz = cbuffer.data(pl_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_z_x_xxyyyyzz = cbuffer.data(pl_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_z_x_xxyyyzzz = cbuffer.data(pl_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_z_x_xxyyzzzz = cbuffer.data(pl_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_z_x_xxyzzzzz = cbuffer.data(pl_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_z_x_xxzzzzzz = cbuffer.data(pl_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_z_x_xyyyyyyy = cbuffer.data(pl_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_z_x_xyyyyyyz = cbuffer.data(pl_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_z_x_xyyyyyzz = cbuffer.data(pl_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_z_x_xyyyyzzz = cbuffer.data(pl_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_z_x_xyyyzzzz = cbuffer.data(pl_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_z_x_xyyzzzzz = cbuffer.data(pl_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_z_x_xyzzzzzz = cbuffer.data(pl_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_z_x_xzzzzzzz = cbuffer.data(pl_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_z_x_yyyyyyyy = cbuffer.data(pl_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_z_x_yyyyyyyz = cbuffer.data(pl_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_z_x_yyyyyyzz = cbuffer.data(pl_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_z_x_yyyyyzzz = cbuffer.data(pl_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_z_x_yyyyzzzz = cbuffer.data(pl_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_z_x_yyyzzzzz = cbuffer.data(pl_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_z_x_yyzzzzzz = cbuffer.data(pl_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_z_x_yzzzzzzz = cbuffer.data(pl_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_z_x_zzzzzzzz = cbuffer.data(pl_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_z_y_xxxxxxxx = cbuffer.data(pl_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_z_y_xxxxxxxy = cbuffer.data(pl_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_z_y_xxxxxxxz = cbuffer.data(pl_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_z_y_xxxxxxyy = cbuffer.data(pl_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_z_y_xxxxxxyz = cbuffer.data(pl_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_z_y_xxxxxxzz = cbuffer.data(pl_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_z_y_xxxxxyyy = cbuffer.data(pl_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_z_y_xxxxxyyz = cbuffer.data(pl_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_z_y_xxxxxyzz = cbuffer.data(pl_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_z_y_xxxxxzzz = cbuffer.data(pl_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_z_y_xxxxyyyy = cbuffer.data(pl_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_z_y_xxxxyyyz = cbuffer.data(pl_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_z_y_xxxxyyzz = cbuffer.data(pl_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_z_y_xxxxyzzz = cbuffer.data(pl_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_z_y_xxxxzzzz = cbuffer.data(pl_geom_01_off + 329 * ccomps * dcomps);

            auto g_0_z_y_xxxyyyyy = cbuffer.data(pl_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_z_y_xxxyyyyz = cbuffer.data(pl_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_z_y_xxxyyyzz = cbuffer.data(pl_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_z_y_xxxyyzzz = cbuffer.data(pl_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_z_y_xxxyzzzz = cbuffer.data(pl_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_z_y_xxxzzzzz = cbuffer.data(pl_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_z_y_xxyyyyyy = cbuffer.data(pl_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_z_y_xxyyyyyz = cbuffer.data(pl_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_z_y_xxyyyyzz = cbuffer.data(pl_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_z_y_xxyyyzzz = cbuffer.data(pl_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_z_y_xxyyzzzz = cbuffer.data(pl_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_z_y_xxyzzzzz = cbuffer.data(pl_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_z_y_xxzzzzzz = cbuffer.data(pl_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_z_y_xyyyyyyy = cbuffer.data(pl_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_z_y_xyyyyyyz = cbuffer.data(pl_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_z_y_xyyyyyzz = cbuffer.data(pl_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_z_y_xyyyyzzz = cbuffer.data(pl_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_z_y_xyyyzzzz = cbuffer.data(pl_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_z_y_xyyzzzzz = cbuffer.data(pl_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_z_y_xyzzzzzz = cbuffer.data(pl_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_z_y_xzzzzzzz = cbuffer.data(pl_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_z_y_yyyyyyyy = cbuffer.data(pl_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_z_y_yyyyyyyz = cbuffer.data(pl_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_z_y_yyyyyyzz = cbuffer.data(pl_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_z_y_yyyyyzzz = cbuffer.data(pl_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_z_y_yyyyzzzz = cbuffer.data(pl_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_z_y_yyyzzzzz = cbuffer.data(pl_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_z_y_yyzzzzzz = cbuffer.data(pl_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_z_y_yzzzzzzz = cbuffer.data(pl_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_z_y_zzzzzzzz = cbuffer.data(pl_geom_01_off + 359 * ccomps * dcomps);

            auto g_0_z_z_xxxxxxxx = cbuffer.data(pl_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_z_z_xxxxxxxy = cbuffer.data(pl_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_z_z_xxxxxxxz = cbuffer.data(pl_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_z_z_xxxxxxyy = cbuffer.data(pl_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_z_z_xxxxxxyz = cbuffer.data(pl_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_z_z_xxxxxxzz = cbuffer.data(pl_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_z_z_xxxxxyyy = cbuffer.data(pl_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_z_z_xxxxxyyz = cbuffer.data(pl_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_z_z_xxxxxyzz = cbuffer.data(pl_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_z_z_xxxxxzzz = cbuffer.data(pl_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_z_z_xxxxyyyy = cbuffer.data(pl_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_z_z_xxxxyyyz = cbuffer.data(pl_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_z_z_xxxxyyzz = cbuffer.data(pl_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_z_z_xxxxyzzz = cbuffer.data(pl_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_z_z_xxxxzzzz = cbuffer.data(pl_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_z_z_xxxyyyyy = cbuffer.data(pl_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_z_z_xxxyyyyz = cbuffer.data(pl_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_z_z_xxxyyyzz = cbuffer.data(pl_geom_01_off + 377 * ccomps * dcomps);

            auto g_0_z_z_xxxyyzzz = cbuffer.data(pl_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_z_z_xxxyzzzz = cbuffer.data(pl_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_z_z_xxxzzzzz = cbuffer.data(pl_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_z_z_xxyyyyyy = cbuffer.data(pl_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_z_z_xxyyyyyz = cbuffer.data(pl_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_z_z_xxyyyyzz = cbuffer.data(pl_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_z_z_xxyyyzzz = cbuffer.data(pl_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_z_z_xxyyzzzz = cbuffer.data(pl_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_z_z_xxyzzzzz = cbuffer.data(pl_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_z_z_xxzzzzzz = cbuffer.data(pl_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_z_z_xyyyyyyy = cbuffer.data(pl_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_z_z_xyyyyyyz = cbuffer.data(pl_geom_01_off + 389 * ccomps * dcomps);

            auto g_0_z_z_xyyyyyzz = cbuffer.data(pl_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_z_z_xyyyyzzz = cbuffer.data(pl_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_z_z_xyyyzzzz = cbuffer.data(pl_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_z_z_xyyzzzzz = cbuffer.data(pl_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_z_z_xyzzzzzz = cbuffer.data(pl_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_z_z_xzzzzzzz = cbuffer.data(pl_geom_01_off + 395 * ccomps * dcomps);

            auto g_0_z_z_yyyyyyyy = cbuffer.data(pl_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_z_z_yyyyyyyz = cbuffer.data(pl_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_z_z_yyyyyyzz = cbuffer.data(pl_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_z_z_yyyyyzzz = cbuffer.data(pl_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_z_z_yyyyzzzz = cbuffer.data(pl_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_z_z_yyyzzzzz = cbuffer.data(pl_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_z_z_yyzzzzzz = cbuffer.data(pl_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_z_z_yzzzzzzz = cbuffer.data(pl_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_z_z_zzzzzzzz = cbuffer.data(pl_geom_01_off + 404 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_dkxx

            const auto dk_geom_01_off = idx_geom_01_dkxx + i * dcomps + j;

            /// Set up 0-36 components of targeted buffer : cbuffer.data(

            auto g_0_x_xx_xxxxxxx = cbuffer.data(dk_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xx_xxxxxxy = cbuffer.data(dk_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xx_xxxxxxz = cbuffer.data(dk_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xx_xxxxxyy = cbuffer.data(dk_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xx_xxxxxyz = cbuffer.data(dk_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xx_xxxxxzz = cbuffer.data(dk_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xx_xxxxyyy = cbuffer.data(dk_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xx_xxxxyyz = cbuffer.data(dk_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xx_xxxxyzz = cbuffer.data(dk_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xx_xxxxzzz = cbuffer.data(dk_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xx_xxxyyyy = cbuffer.data(dk_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xx_xxxyyyz = cbuffer.data(dk_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xx_xxxyyzz = cbuffer.data(dk_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xx_xxxyzzz = cbuffer.data(dk_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xx_xxxzzzz = cbuffer.data(dk_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xx_xxyyyyy = cbuffer.data(dk_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xx_xxyyyyz = cbuffer.data(dk_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xx_xxyyyzz = cbuffer.data(dk_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xx_xxyyzzz = cbuffer.data(dk_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xx_xxyzzzz = cbuffer.data(dk_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xx_xxzzzzz = cbuffer.data(dk_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xx_xyyyyyy = cbuffer.data(dk_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xx_xyyyyyz = cbuffer.data(dk_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xx_xyyyyzz = cbuffer.data(dk_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xx_xyyyzzz = cbuffer.data(dk_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xx_xyyzzzz = cbuffer.data(dk_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xx_xyzzzzz = cbuffer.data(dk_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xx_xzzzzzz = cbuffer.data(dk_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xx_yyyyyyy = cbuffer.data(dk_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xx_yyyyyyz = cbuffer.data(dk_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xx_yyyyyzz = cbuffer.data(dk_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xx_yyyyzzz = cbuffer.data(dk_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xx_yyyzzzz = cbuffer.data(dk_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xx_yyzzzzz = cbuffer.data(dk_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xx_yzzzzzz = cbuffer.data(dk_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xx_zzzzzzz = cbuffer.data(dk_geom_01_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_x_xxxxxxx, g_0_x_x_xxxxxxxx, g_0_x_x_xxxxxxxy, g_0_x_x_xxxxxxxz, g_0_x_x_xxxxxxy, g_0_x_x_xxxxxxyy, g_0_x_x_xxxxxxyz, g_0_x_x_xxxxxxz, g_0_x_x_xxxxxxzz, g_0_x_x_xxxxxyy, g_0_x_x_xxxxxyyy, g_0_x_x_xxxxxyyz, g_0_x_x_xxxxxyz, g_0_x_x_xxxxxyzz, g_0_x_x_xxxxxzz, g_0_x_x_xxxxxzzz, g_0_x_x_xxxxyyy, g_0_x_x_xxxxyyyy, g_0_x_x_xxxxyyyz, g_0_x_x_xxxxyyz, g_0_x_x_xxxxyyzz, g_0_x_x_xxxxyzz, g_0_x_x_xxxxyzzz, g_0_x_x_xxxxzzz, g_0_x_x_xxxxzzzz, g_0_x_x_xxxyyyy, g_0_x_x_xxxyyyyy, g_0_x_x_xxxyyyyz, g_0_x_x_xxxyyyz, g_0_x_x_xxxyyyzz, g_0_x_x_xxxyyzz, g_0_x_x_xxxyyzzz, g_0_x_x_xxxyzzz, g_0_x_x_xxxyzzzz, g_0_x_x_xxxzzzz, g_0_x_x_xxxzzzzz, g_0_x_x_xxyyyyy, g_0_x_x_xxyyyyyy, g_0_x_x_xxyyyyyz, g_0_x_x_xxyyyyz, g_0_x_x_xxyyyyzz, g_0_x_x_xxyyyzz, g_0_x_x_xxyyyzzz, g_0_x_x_xxyyzzz, g_0_x_x_xxyyzzzz, g_0_x_x_xxyzzzz, g_0_x_x_xxyzzzzz, g_0_x_x_xxzzzzz, g_0_x_x_xxzzzzzz, g_0_x_x_xyyyyyy, g_0_x_x_xyyyyyyy, g_0_x_x_xyyyyyyz, g_0_x_x_xyyyyyz, g_0_x_x_xyyyyyzz, g_0_x_x_xyyyyzz, g_0_x_x_xyyyyzzz, g_0_x_x_xyyyzzz, g_0_x_x_xyyyzzzz, g_0_x_x_xyyzzzz, g_0_x_x_xyyzzzzz, g_0_x_x_xyzzzzz, g_0_x_x_xyzzzzzz, g_0_x_x_xzzzzzz, g_0_x_x_xzzzzzzz, g_0_x_x_yyyyyyy, g_0_x_x_yyyyyyz, g_0_x_x_yyyyyzz, g_0_x_x_yyyyzzz, g_0_x_x_yyyzzzz, g_0_x_x_yyzzzzz, g_0_x_x_yzzzzzz, g_0_x_x_zzzzzzz, g_0_x_xx_xxxxxxx, g_0_x_xx_xxxxxxy, g_0_x_xx_xxxxxxz, g_0_x_xx_xxxxxyy, g_0_x_xx_xxxxxyz, g_0_x_xx_xxxxxzz, g_0_x_xx_xxxxyyy, g_0_x_xx_xxxxyyz, g_0_x_xx_xxxxyzz, g_0_x_xx_xxxxzzz, g_0_x_xx_xxxyyyy, g_0_x_xx_xxxyyyz, g_0_x_xx_xxxyyzz, g_0_x_xx_xxxyzzz, g_0_x_xx_xxxzzzz, g_0_x_xx_xxyyyyy, g_0_x_xx_xxyyyyz, g_0_x_xx_xxyyyzz, g_0_x_xx_xxyyzzz, g_0_x_xx_xxyzzzz, g_0_x_xx_xxzzzzz, g_0_x_xx_xyyyyyy, g_0_x_xx_xyyyyyz, g_0_x_xx_xyyyyzz, g_0_x_xx_xyyyzzz, g_0_x_xx_xyyzzzz, g_0_x_xx_xyzzzzz, g_0_x_xx_xzzzzzz, g_0_x_xx_yyyyyyy, g_0_x_xx_yyyyyyz, g_0_x_xx_yyyyyzz, g_0_x_xx_yyyyzzz, g_0_x_xx_yyyzzzz, g_0_x_xx_yyzzzzz, g_0_x_xx_yzzzzzz, g_0_x_xx_zzzzzzz, g_x_xxxxxxx, g_x_xxxxxxy, g_x_xxxxxxz, g_x_xxxxxyy, g_x_xxxxxyz, g_x_xxxxxzz, g_x_xxxxyyy, g_x_xxxxyyz, g_x_xxxxyzz, g_x_xxxxzzz, g_x_xxxyyyy, g_x_xxxyyyz, g_x_xxxyyzz, g_x_xxxyzzz, g_x_xxxzzzz, g_x_xxyyyyy, g_x_xxyyyyz, g_x_xxyyyzz, g_x_xxyyzzz, g_x_xxyzzzz, g_x_xxzzzzz, g_x_xyyyyyy, g_x_xyyyyyz, g_x_xyyyyzz, g_x_xyyyzzz, g_x_xyyzzzz, g_x_xyzzzzz, g_x_xzzzzzz, g_x_yyyyyyy, g_x_yyyyyyz, g_x_yyyyyzz, g_x_yyyyzzz, g_x_yyyzzzz, g_x_yyzzzzz, g_x_yzzzzzz, g_x_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xx_xxxxxxx[k] = g_x_xxxxxxx[k] - g_0_x_x_xxxxxxx[k] * ab_x + g_0_x_x_xxxxxxxx[k];

                g_0_x_xx_xxxxxxy[k] = g_x_xxxxxxy[k] - g_0_x_x_xxxxxxy[k] * ab_x + g_0_x_x_xxxxxxxy[k];

                g_0_x_xx_xxxxxxz[k] = g_x_xxxxxxz[k] - g_0_x_x_xxxxxxz[k] * ab_x + g_0_x_x_xxxxxxxz[k];

                g_0_x_xx_xxxxxyy[k] = g_x_xxxxxyy[k] - g_0_x_x_xxxxxyy[k] * ab_x + g_0_x_x_xxxxxxyy[k];

                g_0_x_xx_xxxxxyz[k] = g_x_xxxxxyz[k] - g_0_x_x_xxxxxyz[k] * ab_x + g_0_x_x_xxxxxxyz[k];

                g_0_x_xx_xxxxxzz[k] = g_x_xxxxxzz[k] - g_0_x_x_xxxxxzz[k] * ab_x + g_0_x_x_xxxxxxzz[k];

                g_0_x_xx_xxxxyyy[k] = g_x_xxxxyyy[k] - g_0_x_x_xxxxyyy[k] * ab_x + g_0_x_x_xxxxxyyy[k];

                g_0_x_xx_xxxxyyz[k] = g_x_xxxxyyz[k] - g_0_x_x_xxxxyyz[k] * ab_x + g_0_x_x_xxxxxyyz[k];

                g_0_x_xx_xxxxyzz[k] = g_x_xxxxyzz[k] - g_0_x_x_xxxxyzz[k] * ab_x + g_0_x_x_xxxxxyzz[k];

                g_0_x_xx_xxxxzzz[k] = g_x_xxxxzzz[k] - g_0_x_x_xxxxzzz[k] * ab_x + g_0_x_x_xxxxxzzz[k];

                g_0_x_xx_xxxyyyy[k] = g_x_xxxyyyy[k] - g_0_x_x_xxxyyyy[k] * ab_x + g_0_x_x_xxxxyyyy[k];

                g_0_x_xx_xxxyyyz[k] = g_x_xxxyyyz[k] - g_0_x_x_xxxyyyz[k] * ab_x + g_0_x_x_xxxxyyyz[k];

                g_0_x_xx_xxxyyzz[k] = g_x_xxxyyzz[k] - g_0_x_x_xxxyyzz[k] * ab_x + g_0_x_x_xxxxyyzz[k];

                g_0_x_xx_xxxyzzz[k] = g_x_xxxyzzz[k] - g_0_x_x_xxxyzzz[k] * ab_x + g_0_x_x_xxxxyzzz[k];

                g_0_x_xx_xxxzzzz[k] = g_x_xxxzzzz[k] - g_0_x_x_xxxzzzz[k] * ab_x + g_0_x_x_xxxxzzzz[k];

                g_0_x_xx_xxyyyyy[k] = g_x_xxyyyyy[k] - g_0_x_x_xxyyyyy[k] * ab_x + g_0_x_x_xxxyyyyy[k];

                g_0_x_xx_xxyyyyz[k] = g_x_xxyyyyz[k] - g_0_x_x_xxyyyyz[k] * ab_x + g_0_x_x_xxxyyyyz[k];

                g_0_x_xx_xxyyyzz[k] = g_x_xxyyyzz[k] - g_0_x_x_xxyyyzz[k] * ab_x + g_0_x_x_xxxyyyzz[k];

                g_0_x_xx_xxyyzzz[k] = g_x_xxyyzzz[k] - g_0_x_x_xxyyzzz[k] * ab_x + g_0_x_x_xxxyyzzz[k];

                g_0_x_xx_xxyzzzz[k] = g_x_xxyzzzz[k] - g_0_x_x_xxyzzzz[k] * ab_x + g_0_x_x_xxxyzzzz[k];

                g_0_x_xx_xxzzzzz[k] = g_x_xxzzzzz[k] - g_0_x_x_xxzzzzz[k] * ab_x + g_0_x_x_xxxzzzzz[k];

                g_0_x_xx_xyyyyyy[k] = g_x_xyyyyyy[k] - g_0_x_x_xyyyyyy[k] * ab_x + g_0_x_x_xxyyyyyy[k];

                g_0_x_xx_xyyyyyz[k] = g_x_xyyyyyz[k] - g_0_x_x_xyyyyyz[k] * ab_x + g_0_x_x_xxyyyyyz[k];

                g_0_x_xx_xyyyyzz[k] = g_x_xyyyyzz[k] - g_0_x_x_xyyyyzz[k] * ab_x + g_0_x_x_xxyyyyzz[k];

                g_0_x_xx_xyyyzzz[k] = g_x_xyyyzzz[k] - g_0_x_x_xyyyzzz[k] * ab_x + g_0_x_x_xxyyyzzz[k];

                g_0_x_xx_xyyzzzz[k] = g_x_xyyzzzz[k] - g_0_x_x_xyyzzzz[k] * ab_x + g_0_x_x_xxyyzzzz[k];

                g_0_x_xx_xyzzzzz[k] = g_x_xyzzzzz[k] - g_0_x_x_xyzzzzz[k] * ab_x + g_0_x_x_xxyzzzzz[k];

                g_0_x_xx_xzzzzzz[k] = g_x_xzzzzzz[k] - g_0_x_x_xzzzzzz[k] * ab_x + g_0_x_x_xxzzzzzz[k];

                g_0_x_xx_yyyyyyy[k] = g_x_yyyyyyy[k] - g_0_x_x_yyyyyyy[k] * ab_x + g_0_x_x_xyyyyyyy[k];

                g_0_x_xx_yyyyyyz[k] = g_x_yyyyyyz[k] - g_0_x_x_yyyyyyz[k] * ab_x + g_0_x_x_xyyyyyyz[k];

                g_0_x_xx_yyyyyzz[k] = g_x_yyyyyzz[k] - g_0_x_x_yyyyyzz[k] * ab_x + g_0_x_x_xyyyyyzz[k];

                g_0_x_xx_yyyyzzz[k] = g_x_yyyyzzz[k] - g_0_x_x_yyyyzzz[k] * ab_x + g_0_x_x_xyyyyzzz[k];

                g_0_x_xx_yyyzzzz[k] = g_x_yyyzzzz[k] - g_0_x_x_yyyzzzz[k] * ab_x + g_0_x_x_xyyyzzzz[k];

                g_0_x_xx_yyzzzzz[k] = g_x_yyzzzzz[k] - g_0_x_x_yyzzzzz[k] * ab_x + g_0_x_x_xyyzzzzz[k];

                g_0_x_xx_yzzzzzz[k] = g_x_yzzzzzz[k] - g_0_x_x_yzzzzzz[k] * ab_x + g_0_x_x_xyzzzzzz[k];

                g_0_x_xx_zzzzzzz[k] = g_x_zzzzzzz[k] - g_0_x_x_zzzzzzz[k] * ab_x + g_0_x_x_xzzzzzzz[k];
            }

            /// Set up 36-72 components of targeted buffer : cbuffer.data(

            auto g_0_x_xy_xxxxxxx = cbuffer.data(dk_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xy_xxxxxxy = cbuffer.data(dk_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xy_xxxxxxz = cbuffer.data(dk_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xy_xxxxxyy = cbuffer.data(dk_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xy_xxxxxyz = cbuffer.data(dk_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xy_xxxxxzz = cbuffer.data(dk_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xy_xxxxyyy = cbuffer.data(dk_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xy_xxxxyyz = cbuffer.data(dk_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xy_xxxxyzz = cbuffer.data(dk_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xy_xxxxzzz = cbuffer.data(dk_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xy_xxxyyyy = cbuffer.data(dk_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xy_xxxyyyz = cbuffer.data(dk_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xy_xxxyyzz = cbuffer.data(dk_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xy_xxxyzzz = cbuffer.data(dk_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xy_xxxzzzz = cbuffer.data(dk_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xy_xxyyyyy = cbuffer.data(dk_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xy_xxyyyyz = cbuffer.data(dk_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xy_xxyyyzz = cbuffer.data(dk_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xy_xxyyzzz = cbuffer.data(dk_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xy_xxyzzzz = cbuffer.data(dk_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xy_xxzzzzz = cbuffer.data(dk_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xy_xyyyyyy = cbuffer.data(dk_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xy_xyyyyyz = cbuffer.data(dk_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xy_xyyyyzz = cbuffer.data(dk_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_xy_xyyyzzz = cbuffer.data(dk_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xy_xyyzzzz = cbuffer.data(dk_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xy_xyzzzzz = cbuffer.data(dk_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xy_xzzzzzz = cbuffer.data(dk_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xy_yyyyyyy = cbuffer.data(dk_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xy_yyyyyyz = cbuffer.data(dk_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xy_yyyyyzz = cbuffer.data(dk_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xy_yyyyzzz = cbuffer.data(dk_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xy_yyyzzzz = cbuffer.data(dk_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xy_yyzzzzz = cbuffer.data(dk_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xy_yzzzzzz = cbuffer.data(dk_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xy_zzzzzzz = cbuffer.data(dk_geom_01_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_x_xxxxxxx, g_0_x_x_xxxxxxxy, g_0_x_x_xxxxxxy, g_0_x_x_xxxxxxyy, g_0_x_x_xxxxxxyz, g_0_x_x_xxxxxxz, g_0_x_x_xxxxxyy, g_0_x_x_xxxxxyyy, g_0_x_x_xxxxxyyz, g_0_x_x_xxxxxyz, g_0_x_x_xxxxxyzz, g_0_x_x_xxxxxzz, g_0_x_x_xxxxyyy, g_0_x_x_xxxxyyyy, g_0_x_x_xxxxyyyz, g_0_x_x_xxxxyyz, g_0_x_x_xxxxyyzz, g_0_x_x_xxxxyzz, g_0_x_x_xxxxyzzz, g_0_x_x_xxxxzzz, g_0_x_x_xxxyyyy, g_0_x_x_xxxyyyyy, g_0_x_x_xxxyyyyz, g_0_x_x_xxxyyyz, g_0_x_x_xxxyyyzz, g_0_x_x_xxxyyzz, g_0_x_x_xxxyyzzz, g_0_x_x_xxxyzzz, g_0_x_x_xxxyzzzz, g_0_x_x_xxxzzzz, g_0_x_x_xxyyyyy, g_0_x_x_xxyyyyyy, g_0_x_x_xxyyyyyz, g_0_x_x_xxyyyyz, g_0_x_x_xxyyyyzz, g_0_x_x_xxyyyzz, g_0_x_x_xxyyyzzz, g_0_x_x_xxyyzzz, g_0_x_x_xxyyzzzz, g_0_x_x_xxyzzzz, g_0_x_x_xxyzzzzz, g_0_x_x_xxzzzzz, g_0_x_x_xyyyyyy, g_0_x_x_xyyyyyyy, g_0_x_x_xyyyyyyz, g_0_x_x_xyyyyyz, g_0_x_x_xyyyyyzz, g_0_x_x_xyyyyzz, g_0_x_x_xyyyyzzz, g_0_x_x_xyyyzzz, g_0_x_x_xyyyzzzz, g_0_x_x_xyyzzzz, g_0_x_x_xyyzzzzz, g_0_x_x_xyzzzzz, g_0_x_x_xyzzzzzz, g_0_x_x_xzzzzzz, g_0_x_x_yyyyyyy, g_0_x_x_yyyyyyyy, g_0_x_x_yyyyyyyz, g_0_x_x_yyyyyyz, g_0_x_x_yyyyyyzz, g_0_x_x_yyyyyzz, g_0_x_x_yyyyyzzz, g_0_x_x_yyyyzzz, g_0_x_x_yyyyzzzz, g_0_x_x_yyyzzzz, g_0_x_x_yyyzzzzz, g_0_x_x_yyzzzzz, g_0_x_x_yyzzzzzz, g_0_x_x_yzzzzzz, g_0_x_x_yzzzzzzz, g_0_x_x_zzzzzzz, g_0_x_xy_xxxxxxx, g_0_x_xy_xxxxxxy, g_0_x_xy_xxxxxxz, g_0_x_xy_xxxxxyy, g_0_x_xy_xxxxxyz, g_0_x_xy_xxxxxzz, g_0_x_xy_xxxxyyy, g_0_x_xy_xxxxyyz, g_0_x_xy_xxxxyzz, g_0_x_xy_xxxxzzz, g_0_x_xy_xxxyyyy, g_0_x_xy_xxxyyyz, g_0_x_xy_xxxyyzz, g_0_x_xy_xxxyzzz, g_0_x_xy_xxxzzzz, g_0_x_xy_xxyyyyy, g_0_x_xy_xxyyyyz, g_0_x_xy_xxyyyzz, g_0_x_xy_xxyyzzz, g_0_x_xy_xxyzzzz, g_0_x_xy_xxzzzzz, g_0_x_xy_xyyyyyy, g_0_x_xy_xyyyyyz, g_0_x_xy_xyyyyzz, g_0_x_xy_xyyyzzz, g_0_x_xy_xyyzzzz, g_0_x_xy_xyzzzzz, g_0_x_xy_xzzzzzz, g_0_x_xy_yyyyyyy, g_0_x_xy_yyyyyyz, g_0_x_xy_yyyyyzz, g_0_x_xy_yyyyzzz, g_0_x_xy_yyyzzzz, g_0_x_xy_yyzzzzz, g_0_x_xy_yzzzzzz, g_0_x_xy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xy_xxxxxxx[k] = -g_0_x_x_xxxxxxx[k] * ab_y + g_0_x_x_xxxxxxxy[k];

                g_0_x_xy_xxxxxxy[k] = -g_0_x_x_xxxxxxy[k] * ab_y + g_0_x_x_xxxxxxyy[k];

                g_0_x_xy_xxxxxxz[k] = -g_0_x_x_xxxxxxz[k] * ab_y + g_0_x_x_xxxxxxyz[k];

                g_0_x_xy_xxxxxyy[k] = -g_0_x_x_xxxxxyy[k] * ab_y + g_0_x_x_xxxxxyyy[k];

                g_0_x_xy_xxxxxyz[k] = -g_0_x_x_xxxxxyz[k] * ab_y + g_0_x_x_xxxxxyyz[k];

                g_0_x_xy_xxxxxzz[k] = -g_0_x_x_xxxxxzz[k] * ab_y + g_0_x_x_xxxxxyzz[k];

                g_0_x_xy_xxxxyyy[k] = -g_0_x_x_xxxxyyy[k] * ab_y + g_0_x_x_xxxxyyyy[k];

                g_0_x_xy_xxxxyyz[k] = -g_0_x_x_xxxxyyz[k] * ab_y + g_0_x_x_xxxxyyyz[k];

                g_0_x_xy_xxxxyzz[k] = -g_0_x_x_xxxxyzz[k] * ab_y + g_0_x_x_xxxxyyzz[k];

                g_0_x_xy_xxxxzzz[k] = -g_0_x_x_xxxxzzz[k] * ab_y + g_0_x_x_xxxxyzzz[k];

                g_0_x_xy_xxxyyyy[k] = -g_0_x_x_xxxyyyy[k] * ab_y + g_0_x_x_xxxyyyyy[k];

                g_0_x_xy_xxxyyyz[k] = -g_0_x_x_xxxyyyz[k] * ab_y + g_0_x_x_xxxyyyyz[k];

                g_0_x_xy_xxxyyzz[k] = -g_0_x_x_xxxyyzz[k] * ab_y + g_0_x_x_xxxyyyzz[k];

                g_0_x_xy_xxxyzzz[k] = -g_0_x_x_xxxyzzz[k] * ab_y + g_0_x_x_xxxyyzzz[k];

                g_0_x_xy_xxxzzzz[k] = -g_0_x_x_xxxzzzz[k] * ab_y + g_0_x_x_xxxyzzzz[k];

                g_0_x_xy_xxyyyyy[k] = -g_0_x_x_xxyyyyy[k] * ab_y + g_0_x_x_xxyyyyyy[k];

                g_0_x_xy_xxyyyyz[k] = -g_0_x_x_xxyyyyz[k] * ab_y + g_0_x_x_xxyyyyyz[k];

                g_0_x_xy_xxyyyzz[k] = -g_0_x_x_xxyyyzz[k] * ab_y + g_0_x_x_xxyyyyzz[k];

                g_0_x_xy_xxyyzzz[k] = -g_0_x_x_xxyyzzz[k] * ab_y + g_0_x_x_xxyyyzzz[k];

                g_0_x_xy_xxyzzzz[k] = -g_0_x_x_xxyzzzz[k] * ab_y + g_0_x_x_xxyyzzzz[k];

                g_0_x_xy_xxzzzzz[k] = -g_0_x_x_xxzzzzz[k] * ab_y + g_0_x_x_xxyzzzzz[k];

                g_0_x_xy_xyyyyyy[k] = -g_0_x_x_xyyyyyy[k] * ab_y + g_0_x_x_xyyyyyyy[k];

                g_0_x_xy_xyyyyyz[k] = -g_0_x_x_xyyyyyz[k] * ab_y + g_0_x_x_xyyyyyyz[k];

                g_0_x_xy_xyyyyzz[k] = -g_0_x_x_xyyyyzz[k] * ab_y + g_0_x_x_xyyyyyzz[k];

                g_0_x_xy_xyyyzzz[k] = -g_0_x_x_xyyyzzz[k] * ab_y + g_0_x_x_xyyyyzzz[k];

                g_0_x_xy_xyyzzzz[k] = -g_0_x_x_xyyzzzz[k] * ab_y + g_0_x_x_xyyyzzzz[k];

                g_0_x_xy_xyzzzzz[k] = -g_0_x_x_xyzzzzz[k] * ab_y + g_0_x_x_xyyzzzzz[k];

                g_0_x_xy_xzzzzzz[k] = -g_0_x_x_xzzzzzz[k] * ab_y + g_0_x_x_xyzzzzzz[k];

                g_0_x_xy_yyyyyyy[k] = -g_0_x_x_yyyyyyy[k] * ab_y + g_0_x_x_yyyyyyyy[k];

                g_0_x_xy_yyyyyyz[k] = -g_0_x_x_yyyyyyz[k] * ab_y + g_0_x_x_yyyyyyyz[k];

                g_0_x_xy_yyyyyzz[k] = -g_0_x_x_yyyyyzz[k] * ab_y + g_0_x_x_yyyyyyzz[k];

                g_0_x_xy_yyyyzzz[k] = -g_0_x_x_yyyyzzz[k] * ab_y + g_0_x_x_yyyyyzzz[k];

                g_0_x_xy_yyyzzzz[k] = -g_0_x_x_yyyzzzz[k] * ab_y + g_0_x_x_yyyyzzzz[k];

                g_0_x_xy_yyzzzzz[k] = -g_0_x_x_yyzzzzz[k] * ab_y + g_0_x_x_yyyzzzzz[k];

                g_0_x_xy_yzzzzzz[k] = -g_0_x_x_yzzzzzz[k] * ab_y + g_0_x_x_yyzzzzzz[k];

                g_0_x_xy_zzzzzzz[k] = -g_0_x_x_zzzzzzz[k] * ab_y + g_0_x_x_yzzzzzzz[k];
            }

            /// Set up 72-108 components of targeted buffer : cbuffer.data(

            auto g_0_x_xz_xxxxxxx = cbuffer.data(dk_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xz_xxxxxxy = cbuffer.data(dk_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xz_xxxxxxz = cbuffer.data(dk_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xz_xxxxxyy = cbuffer.data(dk_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xz_xxxxxyz = cbuffer.data(dk_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xz_xxxxxzz = cbuffer.data(dk_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xz_xxxxyyy = cbuffer.data(dk_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xz_xxxxyyz = cbuffer.data(dk_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xz_xxxxyzz = cbuffer.data(dk_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xz_xxxxzzz = cbuffer.data(dk_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xz_xxxyyyy = cbuffer.data(dk_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xz_xxxyyyz = cbuffer.data(dk_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xz_xxxyyzz = cbuffer.data(dk_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xz_xxxyzzz = cbuffer.data(dk_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xz_xxxzzzz = cbuffer.data(dk_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xz_xxyyyyy = cbuffer.data(dk_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xz_xxyyyyz = cbuffer.data(dk_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xz_xxyyyzz = cbuffer.data(dk_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_xz_xxyyzzz = cbuffer.data(dk_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_xz_xxyzzzz = cbuffer.data(dk_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_xz_xxzzzzz = cbuffer.data(dk_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_xz_xyyyyyy = cbuffer.data(dk_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_xz_xyyyyyz = cbuffer.data(dk_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_xz_xyyyyzz = cbuffer.data(dk_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_xz_xyyyzzz = cbuffer.data(dk_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_xz_xyyzzzz = cbuffer.data(dk_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_xz_xyzzzzz = cbuffer.data(dk_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_xz_xzzzzzz = cbuffer.data(dk_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_xz_yyyyyyy = cbuffer.data(dk_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_xz_yyyyyyz = cbuffer.data(dk_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_xz_yyyyyzz = cbuffer.data(dk_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_xz_yyyyzzz = cbuffer.data(dk_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_xz_yyyzzzz = cbuffer.data(dk_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_xz_yyzzzzz = cbuffer.data(dk_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_xz_yzzzzzz = cbuffer.data(dk_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_xz_zzzzzzz = cbuffer.data(dk_geom_01_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_x_xxxxxxx, g_0_x_x_xxxxxxxz, g_0_x_x_xxxxxxy, g_0_x_x_xxxxxxyz, g_0_x_x_xxxxxxz, g_0_x_x_xxxxxxzz, g_0_x_x_xxxxxyy, g_0_x_x_xxxxxyyz, g_0_x_x_xxxxxyz, g_0_x_x_xxxxxyzz, g_0_x_x_xxxxxzz, g_0_x_x_xxxxxzzz, g_0_x_x_xxxxyyy, g_0_x_x_xxxxyyyz, g_0_x_x_xxxxyyz, g_0_x_x_xxxxyyzz, g_0_x_x_xxxxyzz, g_0_x_x_xxxxyzzz, g_0_x_x_xxxxzzz, g_0_x_x_xxxxzzzz, g_0_x_x_xxxyyyy, g_0_x_x_xxxyyyyz, g_0_x_x_xxxyyyz, g_0_x_x_xxxyyyzz, g_0_x_x_xxxyyzz, g_0_x_x_xxxyyzzz, g_0_x_x_xxxyzzz, g_0_x_x_xxxyzzzz, g_0_x_x_xxxzzzz, g_0_x_x_xxxzzzzz, g_0_x_x_xxyyyyy, g_0_x_x_xxyyyyyz, g_0_x_x_xxyyyyz, g_0_x_x_xxyyyyzz, g_0_x_x_xxyyyzz, g_0_x_x_xxyyyzzz, g_0_x_x_xxyyzzz, g_0_x_x_xxyyzzzz, g_0_x_x_xxyzzzz, g_0_x_x_xxyzzzzz, g_0_x_x_xxzzzzz, g_0_x_x_xxzzzzzz, g_0_x_x_xyyyyyy, g_0_x_x_xyyyyyyz, g_0_x_x_xyyyyyz, g_0_x_x_xyyyyyzz, g_0_x_x_xyyyyzz, g_0_x_x_xyyyyzzz, g_0_x_x_xyyyzzz, g_0_x_x_xyyyzzzz, g_0_x_x_xyyzzzz, g_0_x_x_xyyzzzzz, g_0_x_x_xyzzzzz, g_0_x_x_xyzzzzzz, g_0_x_x_xzzzzzz, g_0_x_x_xzzzzzzz, g_0_x_x_yyyyyyy, g_0_x_x_yyyyyyyz, g_0_x_x_yyyyyyz, g_0_x_x_yyyyyyzz, g_0_x_x_yyyyyzz, g_0_x_x_yyyyyzzz, g_0_x_x_yyyyzzz, g_0_x_x_yyyyzzzz, g_0_x_x_yyyzzzz, g_0_x_x_yyyzzzzz, g_0_x_x_yyzzzzz, g_0_x_x_yyzzzzzz, g_0_x_x_yzzzzzz, g_0_x_x_yzzzzzzz, g_0_x_x_zzzzzzz, g_0_x_x_zzzzzzzz, g_0_x_xz_xxxxxxx, g_0_x_xz_xxxxxxy, g_0_x_xz_xxxxxxz, g_0_x_xz_xxxxxyy, g_0_x_xz_xxxxxyz, g_0_x_xz_xxxxxzz, g_0_x_xz_xxxxyyy, g_0_x_xz_xxxxyyz, g_0_x_xz_xxxxyzz, g_0_x_xz_xxxxzzz, g_0_x_xz_xxxyyyy, g_0_x_xz_xxxyyyz, g_0_x_xz_xxxyyzz, g_0_x_xz_xxxyzzz, g_0_x_xz_xxxzzzz, g_0_x_xz_xxyyyyy, g_0_x_xz_xxyyyyz, g_0_x_xz_xxyyyzz, g_0_x_xz_xxyyzzz, g_0_x_xz_xxyzzzz, g_0_x_xz_xxzzzzz, g_0_x_xz_xyyyyyy, g_0_x_xz_xyyyyyz, g_0_x_xz_xyyyyzz, g_0_x_xz_xyyyzzz, g_0_x_xz_xyyzzzz, g_0_x_xz_xyzzzzz, g_0_x_xz_xzzzzzz, g_0_x_xz_yyyyyyy, g_0_x_xz_yyyyyyz, g_0_x_xz_yyyyyzz, g_0_x_xz_yyyyzzz, g_0_x_xz_yyyzzzz, g_0_x_xz_yyzzzzz, g_0_x_xz_yzzzzzz, g_0_x_xz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xz_xxxxxxx[k] = -g_0_x_x_xxxxxxx[k] * ab_z + g_0_x_x_xxxxxxxz[k];

                g_0_x_xz_xxxxxxy[k] = -g_0_x_x_xxxxxxy[k] * ab_z + g_0_x_x_xxxxxxyz[k];

                g_0_x_xz_xxxxxxz[k] = -g_0_x_x_xxxxxxz[k] * ab_z + g_0_x_x_xxxxxxzz[k];

                g_0_x_xz_xxxxxyy[k] = -g_0_x_x_xxxxxyy[k] * ab_z + g_0_x_x_xxxxxyyz[k];

                g_0_x_xz_xxxxxyz[k] = -g_0_x_x_xxxxxyz[k] * ab_z + g_0_x_x_xxxxxyzz[k];

                g_0_x_xz_xxxxxzz[k] = -g_0_x_x_xxxxxzz[k] * ab_z + g_0_x_x_xxxxxzzz[k];

                g_0_x_xz_xxxxyyy[k] = -g_0_x_x_xxxxyyy[k] * ab_z + g_0_x_x_xxxxyyyz[k];

                g_0_x_xz_xxxxyyz[k] = -g_0_x_x_xxxxyyz[k] * ab_z + g_0_x_x_xxxxyyzz[k];

                g_0_x_xz_xxxxyzz[k] = -g_0_x_x_xxxxyzz[k] * ab_z + g_0_x_x_xxxxyzzz[k];

                g_0_x_xz_xxxxzzz[k] = -g_0_x_x_xxxxzzz[k] * ab_z + g_0_x_x_xxxxzzzz[k];

                g_0_x_xz_xxxyyyy[k] = -g_0_x_x_xxxyyyy[k] * ab_z + g_0_x_x_xxxyyyyz[k];

                g_0_x_xz_xxxyyyz[k] = -g_0_x_x_xxxyyyz[k] * ab_z + g_0_x_x_xxxyyyzz[k];

                g_0_x_xz_xxxyyzz[k] = -g_0_x_x_xxxyyzz[k] * ab_z + g_0_x_x_xxxyyzzz[k];

                g_0_x_xz_xxxyzzz[k] = -g_0_x_x_xxxyzzz[k] * ab_z + g_0_x_x_xxxyzzzz[k];

                g_0_x_xz_xxxzzzz[k] = -g_0_x_x_xxxzzzz[k] * ab_z + g_0_x_x_xxxzzzzz[k];

                g_0_x_xz_xxyyyyy[k] = -g_0_x_x_xxyyyyy[k] * ab_z + g_0_x_x_xxyyyyyz[k];

                g_0_x_xz_xxyyyyz[k] = -g_0_x_x_xxyyyyz[k] * ab_z + g_0_x_x_xxyyyyzz[k];

                g_0_x_xz_xxyyyzz[k] = -g_0_x_x_xxyyyzz[k] * ab_z + g_0_x_x_xxyyyzzz[k];

                g_0_x_xz_xxyyzzz[k] = -g_0_x_x_xxyyzzz[k] * ab_z + g_0_x_x_xxyyzzzz[k];

                g_0_x_xz_xxyzzzz[k] = -g_0_x_x_xxyzzzz[k] * ab_z + g_0_x_x_xxyzzzzz[k];

                g_0_x_xz_xxzzzzz[k] = -g_0_x_x_xxzzzzz[k] * ab_z + g_0_x_x_xxzzzzzz[k];

                g_0_x_xz_xyyyyyy[k] = -g_0_x_x_xyyyyyy[k] * ab_z + g_0_x_x_xyyyyyyz[k];

                g_0_x_xz_xyyyyyz[k] = -g_0_x_x_xyyyyyz[k] * ab_z + g_0_x_x_xyyyyyzz[k];

                g_0_x_xz_xyyyyzz[k] = -g_0_x_x_xyyyyzz[k] * ab_z + g_0_x_x_xyyyyzzz[k];

                g_0_x_xz_xyyyzzz[k] = -g_0_x_x_xyyyzzz[k] * ab_z + g_0_x_x_xyyyzzzz[k];

                g_0_x_xz_xyyzzzz[k] = -g_0_x_x_xyyzzzz[k] * ab_z + g_0_x_x_xyyzzzzz[k];

                g_0_x_xz_xyzzzzz[k] = -g_0_x_x_xyzzzzz[k] * ab_z + g_0_x_x_xyzzzzzz[k];

                g_0_x_xz_xzzzzzz[k] = -g_0_x_x_xzzzzzz[k] * ab_z + g_0_x_x_xzzzzzzz[k];

                g_0_x_xz_yyyyyyy[k] = -g_0_x_x_yyyyyyy[k] * ab_z + g_0_x_x_yyyyyyyz[k];

                g_0_x_xz_yyyyyyz[k] = -g_0_x_x_yyyyyyz[k] * ab_z + g_0_x_x_yyyyyyzz[k];

                g_0_x_xz_yyyyyzz[k] = -g_0_x_x_yyyyyzz[k] * ab_z + g_0_x_x_yyyyyzzz[k];

                g_0_x_xz_yyyyzzz[k] = -g_0_x_x_yyyyzzz[k] * ab_z + g_0_x_x_yyyyzzzz[k];

                g_0_x_xz_yyyzzzz[k] = -g_0_x_x_yyyzzzz[k] * ab_z + g_0_x_x_yyyzzzzz[k];

                g_0_x_xz_yyzzzzz[k] = -g_0_x_x_yyzzzzz[k] * ab_z + g_0_x_x_yyzzzzzz[k];

                g_0_x_xz_yzzzzzz[k] = -g_0_x_x_yzzzzzz[k] * ab_z + g_0_x_x_yzzzzzzz[k];

                g_0_x_xz_zzzzzzz[k] = -g_0_x_x_zzzzzzz[k] * ab_z + g_0_x_x_zzzzzzzz[k];
            }

            /// Set up 108-144 components of targeted buffer : cbuffer.data(

            auto g_0_x_yy_xxxxxxx = cbuffer.data(dk_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_yy_xxxxxxy = cbuffer.data(dk_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_yy_xxxxxxz = cbuffer.data(dk_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_yy_xxxxxyy = cbuffer.data(dk_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_yy_xxxxxyz = cbuffer.data(dk_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_yy_xxxxxzz = cbuffer.data(dk_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_yy_xxxxyyy = cbuffer.data(dk_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_yy_xxxxyyz = cbuffer.data(dk_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_yy_xxxxyzz = cbuffer.data(dk_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_yy_xxxxzzz = cbuffer.data(dk_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_yy_xxxyyyy = cbuffer.data(dk_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_yy_xxxyyyz = cbuffer.data(dk_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_x_yy_xxxyyzz = cbuffer.data(dk_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_yy_xxxyzzz = cbuffer.data(dk_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_yy_xxxzzzz = cbuffer.data(dk_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_yy_xxyyyyy = cbuffer.data(dk_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_yy_xxyyyyz = cbuffer.data(dk_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_yy_xxyyyzz = cbuffer.data(dk_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_x_yy_xxyyzzz = cbuffer.data(dk_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_yy_xxyzzzz = cbuffer.data(dk_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_yy_xxzzzzz = cbuffer.data(dk_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_yy_xyyyyyy = cbuffer.data(dk_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_yy_xyyyyyz = cbuffer.data(dk_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_yy_xyyyyzz = cbuffer.data(dk_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_yy_xyyyzzz = cbuffer.data(dk_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_yy_xyyzzzz = cbuffer.data(dk_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_yy_xyzzzzz = cbuffer.data(dk_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_x_yy_xzzzzzz = cbuffer.data(dk_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_yy_yyyyyyy = cbuffer.data(dk_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_yy_yyyyyyz = cbuffer.data(dk_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_x_yy_yyyyyzz = cbuffer.data(dk_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_yy_yyyyzzz = cbuffer.data(dk_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_x_yy_yyyzzzz = cbuffer.data(dk_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_yy_yyzzzzz = cbuffer.data(dk_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_yy_yzzzzzz = cbuffer.data(dk_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_yy_zzzzzzz = cbuffer.data(dk_geom_01_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_y_xxxxxxx, g_0_x_y_xxxxxxxy, g_0_x_y_xxxxxxy, g_0_x_y_xxxxxxyy, g_0_x_y_xxxxxxyz, g_0_x_y_xxxxxxz, g_0_x_y_xxxxxyy, g_0_x_y_xxxxxyyy, g_0_x_y_xxxxxyyz, g_0_x_y_xxxxxyz, g_0_x_y_xxxxxyzz, g_0_x_y_xxxxxzz, g_0_x_y_xxxxyyy, g_0_x_y_xxxxyyyy, g_0_x_y_xxxxyyyz, g_0_x_y_xxxxyyz, g_0_x_y_xxxxyyzz, g_0_x_y_xxxxyzz, g_0_x_y_xxxxyzzz, g_0_x_y_xxxxzzz, g_0_x_y_xxxyyyy, g_0_x_y_xxxyyyyy, g_0_x_y_xxxyyyyz, g_0_x_y_xxxyyyz, g_0_x_y_xxxyyyzz, g_0_x_y_xxxyyzz, g_0_x_y_xxxyyzzz, g_0_x_y_xxxyzzz, g_0_x_y_xxxyzzzz, g_0_x_y_xxxzzzz, g_0_x_y_xxyyyyy, g_0_x_y_xxyyyyyy, g_0_x_y_xxyyyyyz, g_0_x_y_xxyyyyz, g_0_x_y_xxyyyyzz, g_0_x_y_xxyyyzz, g_0_x_y_xxyyyzzz, g_0_x_y_xxyyzzz, g_0_x_y_xxyyzzzz, g_0_x_y_xxyzzzz, g_0_x_y_xxyzzzzz, g_0_x_y_xxzzzzz, g_0_x_y_xyyyyyy, g_0_x_y_xyyyyyyy, g_0_x_y_xyyyyyyz, g_0_x_y_xyyyyyz, g_0_x_y_xyyyyyzz, g_0_x_y_xyyyyzz, g_0_x_y_xyyyyzzz, g_0_x_y_xyyyzzz, g_0_x_y_xyyyzzzz, g_0_x_y_xyyzzzz, g_0_x_y_xyyzzzzz, g_0_x_y_xyzzzzz, g_0_x_y_xyzzzzzz, g_0_x_y_xzzzzzz, g_0_x_y_yyyyyyy, g_0_x_y_yyyyyyyy, g_0_x_y_yyyyyyyz, g_0_x_y_yyyyyyz, g_0_x_y_yyyyyyzz, g_0_x_y_yyyyyzz, g_0_x_y_yyyyyzzz, g_0_x_y_yyyyzzz, g_0_x_y_yyyyzzzz, g_0_x_y_yyyzzzz, g_0_x_y_yyyzzzzz, g_0_x_y_yyzzzzz, g_0_x_y_yyzzzzzz, g_0_x_y_yzzzzzz, g_0_x_y_yzzzzzzz, g_0_x_y_zzzzzzz, g_0_x_yy_xxxxxxx, g_0_x_yy_xxxxxxy, g_0_x_yy_xxxxxxz, g_0_x_yy_xxxxxyy, g_0_x_yy_xxxxxyz, g_0_x_yy_xxxxxzz, g_0_x_yy_xxxxyyy, g_0_x_yy_xxxxyyz, g_0_x_yy_xxxxyzz, g_0_x_yy_xxxxzzz, g_0_x_yy_xxxyyyy, g_0_x_yy_xxxyyyz, g_0_x_yy_xxxyyzz, g_0_x_yy_xxxyzzz, g_0_x_yy_xxxzzzz, g_0_x_yy_xxyyyyy, g_0_x_yy_xxyyyyz, g_0_x_yy_xxyyyzz, g_0_x_yy_xxyyzzz, g_0_x_yy_xxyzzzz, g_0_x_yy_xxzzzzz, g_0_x_yy_xyyyyyy, g_0_x_yy_xyyyyyz, g_0_x_yy_xyyyyzz, g_0_x_yy_xyyyzzz, g_0_x_yy_xyyzzzz, g_0_x_yy_xyzzzzz, g_0_x_yy_xzzzzzz, g_0_x_yy_yyyyyyy, g_0_x_yy_yyyyyyz, g_0_x_yy_yyyyyzz, g_0_x_yy_yyyyzzz, g_0_x_yy_yyyzzzz, g_0_x_yy_yyzzzzz, g_0_x_yy_yzzzzzz, g_0_x_yy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yy_xxxxxxx[k] = -g_0_x_y_xxxxxxx[k] * ab_y + g_0_x_y_xxxxxxxy[k];

                g_0_x_yy_xxxxxxy[k] = -g_0_x_y_xxxxxxy[k] * ab_y + g_0_x_y_xxxxxxyy[k];

                g_0_x_yy_xxxxxxz[k] = -g_0_x_y_xxxxxxz[k] * ab_y + g_0_x_y_xxxxxxyz[k];

                g_0_x_yy_xxxxxyy[k] = -g_0_x_y_xxxxxyy[k] * ab_y + g_0_x_y_xxxxxyyy[k];

                g_0_x_yy_xxxxxyz[k] = -g_0_x_y_xxxxxyz[k] * ab_y + g_0_x_y_xxxxxyyz[k];

                g_0_x_yy_xxxxxzz[k] = -g_0_x_y_xxxxxzz[k] * ab_y + g_0_x_y_xxxxxyzz[k];

                g_0_x_yy_xxxxyyy[k] = -g_0_x_y_xxxxyyy[k] * ab_y + g_0_x_y_xxxxyyyy[k];

                g_0_x_yy_xxxxyyz[k] = -g_0_x_y_xxxxyyz[k] * ab_y + g_0_x_y_xxxxyyyz[k];

                g_0_x_yy_xxxxyzz[k] = -g_0_x_y_xxxxyzz[k] * ab_y + g_0_x_y_xxxxyyzz[k];

                g_0_x_yy_xxxxzzz[k] = -g_0_x_y_xxxxzzz[k] * ab_y + g_0_x_y_xxxxyzzz[k];

                g_0_x_yy_xxxyyyy[k] = -g_0_x_y_xxxyyyy[k] * ab_y + g_0_x_y_xxxyyyyy[k];

                g_0_x_yy_xxxyyyz[k] = -g_0_x_y_xxxyyyz[k] * ab_y + g_0_x_y_xxxyyyyz[k];

                g_0_x_yy_xxxyyzz[k] = -g_0_x_y_xxxyyzz[k] * ab_y + g_0_x_y_xxxyyyzz[k];

                g_0_x_yy_xxxyzzz[k] = -g_0_x_y_xxxyzzz[k] * ab_y + g_0_x_y_xxxyyzzz[k];

                g_0_x_yy_xxxzzzz[k] = -g_0_x_y_xxxzzzz[k] * ab_y + g_0_x_y_xxxyzzzz[k];

                g_0_x_yy_xxyyyyy[k] = -g_0_x_y_xxyyyyy[k] * ab_y + g_0_x_y_xxyyyyyy[k];

                g_0_x_yy_xxyyyyz[k] = -g_0_x_y_xxyyyyz[k] * ab_y + g_0_x_y_xxyyyyyz[k];

                g_0_x_yy_xxyyyzz[k] = -g_0_x_y_xxyyyzz[k] * ab_y + g_0_x_y_xxyyyyzz[k];

                g_0_x_yy_xxyyzzz[k] = -g_0_x_y_xxyyzzz[k] * ab_y + g_0_x_y_xxyyyzzz[k];

                g_0_x_yy_xxyzzzz[k] = -g_0_x_y_xxyzzzz[k] * ab_y + g_0_x_y_xxyyzzzz[k];

                g_0_x_yy_xxzzzzz[k] = -g_0_x_y_xxzzzzz[k] * ab_y + g_0_x_y_xxyzzzzz[k];

                g_0_x_yy_xyyyyyy[k] = -g_0_x_y_xyyyyyy[k] * ab_y + g_0_x_y_xyyyyyyy[k];

                g_0_x_yy_xyyyyyz[k] = -g_0_x_y_xyyyyyz[k] * ab_y + g_0_x_y_xyyyyyyz[k];

                g_0_x_yy_xyyyyzz[k] = -g_0_x_y_xyyyyzz[k] * ab_y + g_0_x_y_xyyyyyzz[k];

                g_0_x_yy_xyyyzzz[k] = -g_0_x_y_xyyyzzz[k] * ab_y + g_0_x_y_xyyyyzzz[k];

                g_0_x_yy_xyyzzzz[k] = -g_0_x_y_xyyzzzz[k] * ab_y + g_0_x_y_xyyyzzzz[k];

                g_0_x_yy_xyzzzzz[k] = -g_0_x_y_xyzzzzz[k] * ab_y + g_0_x_y_xyyzzzzz[k];

                g_0_x_yy_xzzzzzz[k] = -g_0_x_y_xzzzzzz[k] * ab_y + g_0_x_y_xyzzzzzz[k];

                g_0_x_yy_yyyyyyy[k] = -g_0_x_y_yyyyyyy[k] * ab_y + g_0_x_y_yyyyyyyy[k];

                g_0_x_yy_yyyyyyz[k] = -g_0_x_y_yyyyyyz[k] * ab_y + g_0_x_y_yyyyyyyz[k];

                g_0_x_yy_yyyyyzz[k] = -g_0_x_y_yyyyyzz[k] * ab_y + g_0_x_y_yyyyyyzz[k];

                g_0_x_yy_yyyyzzz[k] = -g_0_x_y_yyyyzzz[k] * ab_y + g_0_x_y_yyyyyzzz[k];

                g_0_x_yy_yyyzzzz[k] = -g_0_x_y_yyyzzzz[k] * ab_y + g_0_x_y_yyyyzzzz[k];

                g_0_x_yy_yyzzzzz[k] = -g_0_x_y_yyzzzzz[k] * ab_y + g_0_x_y_yyyzzzzz[k];

                g_0_x_yy_yzzzzzz[k] = -g_0_x_y_yzzzzzz[k] * ab_y + g_0_x_y_yyzzzzzz[k];

                g_0_x_yy_zzzzzzz[k] = -g_0_x_y_zzzzzzz[k] * ab_y + g_0_x_y_yzzzzzzz[k];
            }

            /// Set up 144-180 components of targeted buffer : cbuffer.data(

            auto g_0_x_yz_xxxxxxx = cbuffer.data(dk_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_yz_xxxxxxy = cbuffer.data(dk_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_yz_xxxxxxz = cbuffer.data(dk_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_yz_xxxxxyy = cbuffer.data(dk_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_yz_xxxxxyz = cbuffer.data(dk_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_yz_xxxxxzz = cbuffer.data(dk_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_x_yz_xxxxyyy = cbuffer.data(dk_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_x_yz_xxxxyyz = cbuffer.data(dk_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_x_yz_xxxxyzz = cbuffer.data(dk_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_x_yz_xxxxzzz = cbuffer.data(dk_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_x_yz_xxxyyyy = cbuffer.data(dk_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_x_yz_xxxyyyz = cbuffer.data(dk_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_x_yz_xxxyyzz = cbuffer.data(dk_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_x_yz_xxxyzzz = cbuffer.data(dk_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_x_yz_xxxzzzz = cbuffer.data(dk_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_x_yz_xxyyyyy = cbuffer.data(dk_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_x_yz_xxyyyyz = cbuffer.data(dk_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_x_yz_xxyyyzz = cbuffer.data(dk_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_x_yz_xxyyzzz = cbuffer.data(dk_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_x_yz_xxyzzzz = cbuffer.data(dk_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_x_yz_xxzzzzz = cbuffer.data(dk_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_x_yz_xyyyyyy = cbuffer.data(dk_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_x_yz_xyyyyyz = cbuffer.data(dk_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_x_yz_xyyyyzz = cbuffer.data(dk_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_x_yz_xyyyzzz = cbuffer.data(dk_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_x_yz_xyyzzzz = cbuffer.data(dk_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_x_yz_xyzzzzz = cbuffer.data(dk_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_x_yz_xzzzzzz = cbuffer.data(dk_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_x_yz_yyyyyyy = cbuffer.data(dk_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_x_yz_yyyyyyz = cbuffer.data(dk_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_x_yz_yyyyyzz = cbuffer.data(dk_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_x_yz_yyyyzzz = cbuffer.data(dk_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_x_yz_yyyzzzz = cbuffer.data(dk_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_x_yz_yyzzzzz = cbuffer.data(dk_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_x_yz_yzzzzzz = cbuffer.data(dk_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_x_yz_zzzzzzz = cbuffer.data(dk_geom_01_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yz_xxxxxxx, g_0_x_yz_xxxxxxy, g_0_x_yz_xxxxxxz, g_0_x_yz_xxxxxyy, g_0_x_yz_xxxxxyz, g_0_x_yz_xxxxxzz, g_0_x_yz_xxxxyyy, g_0_x_yz_xxxxyyz, g_0_x_yz_xxxxyzz, g_0_x_yz_xxxxzzz, g_0_x_yz_xxxyyyy, g_0_x_yz_xxxyyyz, g_0_x_yz_xxxyyzz, g_0_x_yz_xxxyzzz, g_0_x_yz_xxxzzzz, g_0_x_yz_xxyyyyy, g_0_x_yz_xxyyyyz, g_0_x_yz_xxyyyzz, g_0_x_yz_xxyyzzz, g_0_x_yz_xxyzzzz, g_0_x_yz_xxzzzzz, g_0_x_yz_xyyyyyy, g_0_x_yz_xyyyyyz, g_0_x_yz_xyyyyzz, g_0_x_yz_xyyyzzz, g_0_x_yz_xyyzzzz, g_0_x_yz_xyzzzzz, g_0_x_yz_xzzzzzz, g_0_x_yz_yyyyyyy, g_0_x_yz_yyyyyyz, g_0_x_yz_yyyyyzz, g_0_x_yz_yyyyzzz, g_0_x_yz_yyyzzzz, g_0_x_yz_yyzzzzz, g_0_x_yz_yzzzzzz, g_0_x_yz_zzzzzzz, g_0_x_z_xxxxxxx, g_0_x_z_xxxxxxxy, g_0_x_z_xxxxxxy, g_0_x_z_xxxxxxyy, g_0_x_z_xxxxxxyz, g_0_x_z_xxxxxxz, g_0_x_z_xxxxxyy, g_0_x_z_xxxxxyyy, g_0_x_z_xxxxxyyz, g_0_x_z_xxxxxyz, g_0_x_z_xxxxxyzz, g_0_x_z_xxxxxzz, g_0_x_z_xxxxyyy, g_0_x_z_xxxxyyyy, g_0_x_z_xxxxyyyz, g_0_x_z_xxxxyyz, g_0_x_z_xxxxyyzz, g_0_x_z_xxxxyzz, g_0_x_z_xxxxyzzz, g_0_x_z_xxxxzzz, g_0_x_z_xxxyyyy, g_0_x_z_xxxyyyyy, g_0_x_z_xxxyyyyz, g_0_x_z_xxxyyyz, g_0_x_z_xxxyyyzz, g_0_x_z_xxxyyzz, g_0_x_z_xxxyyzzz, g_0_x_z_xxxyzzz, g_0_x_z_xxxyzzzz, g_0_x_z_xxxzzzz, g_0_x_z_xxyyyyy, g_0_x_z_xxyyyyyy, g_0_x_z_xxyyyyyz, g_0_x_z_xxyyyyz, g_0_x_z_xxyyyyzz, g_0_x_z_xxyyyzz, g_0_x_z_xxyyyzzz, g_0_x_z_xxyyzzz, g_0_x_z_xxyyzzzz, g_0_x_z_xxyzzzz, g_0_x_z_xxyzzzzz, g_0_x_z_xxzzzzz, g_0_x_z_xyyyyyy, g_0_x_z_xyyyyyyy, g_0_x_z_xyyyyyyz, g_0_x_z_xyyyyyz, g_0_x_z_xyyyyyzz, g_0_x_z_xyyyyzz, g_0_x_z_xyyyyzzz, g_0_x_z_xyyyzzz, g_0_x_z_xyyyzzzz, g_0_x_z_xyyzzzz, g_0_x_z_xyyzzzzz, g_0_x_z_xyzzzzz, g_0_x_z_xyzzzzzz, g_0_x_z_xzzzzzz, g_0_x_z_yyyyyyy, g_0_x_z_yyyyyyyy, g_0_x_z_yyyyyyyz, g_0_x_z_yyyyyyz, g_0_x_z_yyyyyyzz, g_0_x_z_yyyyyzz, g_0_x_z_yyyyyzzz, g_0_x_z_yyyyzzz, g_0_x_z_yyyyzzzz, g_0_x_z_yyyzzzz, g_0_x_z_yyyzzzzz, g_0_x_z_yyzzzzz, g_0_x_z_yyzzzzzz, g_0_x_z_yzzzzzz, g_0_x_z_yzzzzzzz, g_0_x_z_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yz_xxxxxxx[k] = -g_0_x_z_xxxxxxx[k] * ab_y + g_0_x_z_xxxxxxxy[k];

                g_0_x_yz_xxxxxxy[k] = -g_0_x_z_xxxxxxy[k] * ab_y + g_0_x_z_xxxxxxyy[k];

                g_0_x_yz_xxxxxxz[k] = -g_0_x_z_xxxxxxz[k] * ab_y + g_0_x_z_xxxxxxyz[k];

                g_0_x_yz_xxxxxyy[k] = -g_0_x_z_xxxxxyy[k] * ab_y + g_0_x_z_xxxxxyyy[k];

                g_0_x_yz_xxxxxyz[k] = -g_0_x_z_xxxxxyz[k] * ab_y + g_0_x_z_xxxxxyyz[k];

                g_0_x_yz_xxxxxzz[k] = -g_0_x_z_xxxxxzz[k] * ab_y + g_0_x_z_xxxxxyzz[k];

                g_0_x_yz_xxxxyyy[k] = -g_0_x_z_xxxxyyy[k] * ab_y + g_0_x_z_xxxxyyyy[k];

                g_0_x_yz_xxxxyyz[k] = -g_0_x_z_xxxxyyz[k] * ab_y + g_0_x_z_xxxxyyyz[k];

                g_0_x_yz_xxxxyzz[k] = -g_0_x_z_xxxxyzz[k] * ab_y + g_0_x_z_xxxxyyzz[k];

                g_0_x_yz_xxxxzzz[k] = -g_0_x_z_xxxxzzz[k] * ab_y + g_0_x_z_xxxxyzzz[k];

                g_0_x_yz_xxxyyyy[k] = -g_0_x_z_xxxyyyy[k] * ab_y + g_0_x_z_xxxyyyyy[k];

                g_0_x_yz_xxxyyyz[k] = -g_0_x_z_xxxyyyz[k] * ab_y + g_0_x_z_xxxyyyyz[k];

                g_0_x_yz_xxxyyzz[k] = -g_0_x_z_xxxyyzz[k] * ab_y + g_0_x_z_xxxyyyzz[k];

                g_0_x_yz_xxxyzzz[k] = -g_0_x_z_xxxyzzz[k] * ab_y + g_0_x_z_xxxyyzzz[k];

                g_0_x_yz_xxxzzzz[k] = -g_0_x_z_xxxzzzz[k] * ab_y + g_0_x_z_xxxyzzzz[k];

                g_0_x_yz_xxyyyyy[k] = -g_0_x_z_xxyyyyy[k] * ab_y + g_0_x_z_xxyyyyyy[k];

                g_0_x_yz_xxyyyyz[k] = -g_0_x_z_xxyyyyz[k] * ab_y + g_0_x_z_xxyyyyyz[k];

                g_0_x_yz_xxyyyzz[k] = -g_0_x_z_xxyyyzz[k] * ab_y + g_0_x_z_xxyyyyzz[k];

                g_0_x_yz_xxyyzzz[k] = -g_0_x_z_xxyyzzz[k] * ab_y + g_0_x_z_xxyyyzzz[k];

                g_0_x_yz_xxyzzzz[k] = -g_0_x_z_xxyzzzz[k] * ab_y + g_0_x_z_xxyyzzzz[k];

                g_0_x_yz_xxzzzzz[k] = -g_0_x_z_xxzzzzz[k] * ab_y + g_0_x_z_xxyzzzzz[k];

                g_0_x_yz_xyyyyyy[k] = -g_0_x_z_xyyyyyy[k] * ab_y + g_0_x_z_xyyyyyyy[k];

                g_0_x_yz_xyyyyyz[k] = -g_0_x_z_xyyyyyz[k] * ab_y + g_0_x_z_xyyyyyyz[k];

                g_0_x_yz_xyyyyzz[k] = -g_0_x_z_xyyyyzz[k] * ab_y + g_0_x_z_xyyyyyzz[k];

                g_0_x_yz_xyyyzzz[k] = -g_0_x_z_xyyyzzz[k] * ab_y + g_0_x_z_xyyyyzzz[k];

                g_0_x_yz_xyyzzzz[k] = -g_0_x_z_xyyzzzz[k] * ab_y + g_0_x_z_xyyyzzzz[k];

                g_0_x_yz_xyzzzzz[k] = -g_0_x_z_xyzzzzz[k] * ab_y + g_0_x_z_xyyzzzzz[k];

                g_0_x_yz_xzzzzzz[k] = -g_0_x_z_xzzzzzz[k] * ab_y + g_0_x_z_xyzzzzzz[k];

                g_0_x_yz_yyyyyyy[k] = -g_0_x_z_yyyyyyy[k] * ab_y + g_0_x_z_yyyyyyyy[k];

                g_0_x_yz_yyyyyyz[k] = -g_0_x_z_yyyyyyz[k] * ab_y + g_0_x_z_yyyyyyyz[k];

                g_0_x_yz_yyyyyzz[k] = -g_0_x_z_yyyyyzz[k] * ab_y + g_0_x_z_yyyyyyzz[k];

                g_0_x_yz_yyyyzzz[k] = -g_0_x_z_yyyyzzz[k] * ab_y + g_0_x_z_yyyyyzzz[k];

                g_0_x_yz_yyyzzzz[k] = -g_0_x_z_yyyzzzz[k] * ab_y + g_0_x_z_yyyyzzzz[k];

                g_0_x_yz_yyzzzzz[k] = -g_0_x_z_yyzzzzz[k] * ab_y + g_0_x_z_yyyzzzzz[k];

                g_0_x_yz_yzzzzzz[k] = -g_0_x_z_yzzzzzz[k] * ab_y + g_0_x_z_yyzzzzzz[k];

                g_0_x_yz_zzzzzzz[k] = -g_0_x_z_zzzzzzz[k] * ab_y + g_0_x_z_yzzzzzzz[k];
            }

            /// Set up 180-216 components of targeted buffer : cbuffer.data(

            auto g_0_x_zz_xxxxxxx = cbuffer.data(dk_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_x_zz_xxxxxxy = cbuffer.data(dk_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_x_zz_xxxxxxz = cbuffer.data(dk_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_x_zz_xxxxxyy = cbuffer.data(dk_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_x_zz_xxxxxyz = cbuffer.data(dk_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_x_zz_xxxxxzz = cbuffer.data(dk_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_x_zz_xxxxyyy = cbuffer.data(dk_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_x_zz_xxxxyyz = cbuffer.data(dk_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_x_zz_xxxxyzz = cbuffer.data(dk_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_x_zz_xxxxzzz = cbuffer.data(dk_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_x_zz_xxxyyyy = cbuffer.data(dk_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_x_zz_xxxyyyz = cbuffer.data(dk_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_x_zz_xxxyyzz = cbuffer.data(dk_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_x_zz_xxxyzzz = cbuffer.data(dk_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_x_zz_xxxzzzz = cbuffer.data(dk_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_x_zz_xxyyyyy = cbuffer.data(dk_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_x_zz_xxyyyyz = cbuffer.data(dk_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_x_zz_xxyyyzz = cbuffer.data(dk_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_x_zz_xxyyzzz = cbuffer.data(dk_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_x_zz_xxyzzzz = cbuffer.data(dk_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_x_zz_xxzzzzz = cbuffer.data(dk_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_x_zz_xyyyyyy = cbuffer.data(dk_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_x_zz_xyyyyyz = cbuffer.data(dk_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_x_zz_xyyyyzz = cbuffer.data(dk_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_x_zz_xyyyzzz = cbuffer.data(dk_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_x_zz_xyyzzzz = cbuffer.data(dk_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_x_zz_xyzzzzz = cbuffer.data(dk_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_x_zz_xzzzzzz = cbuffer.data(dk_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_x_zz_yyyyyyy = cbuffer.data(dk_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_x_zz_yyyyyyz = cbuffer.data(dk_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_x_zz_yyyyyzz = cbuffer.data(dk_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_x_zz_yyyyzzz = cbuffer.data(dk_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_x_zz_yyyzzzz = cbuffer.data(dk_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_x_zz_yyzzzzz = cbuffer.data(dk_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_x_zz_yzzzzzz = cbuffer.data(dk_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_x_zz_zzzzzzz = cbuffer.data(dk_geom_01_off + 215 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_z_xxxxxxx, g_0_x_z_xxxxxxxz, g_0_x_z_xxxxxxy, g_0_x_z_xxxxxxyz, g_0_x_z_xxxxxxz, g_0_x_z_xxxxxxzz, g_0_x_z_xxxxxyy, g_0_x_z_xxxxxyyz, g_0_x_z_xxxxxyz, g_0_x_z_xxxxxyzz, g_0_x_z_xxxxxzz, g_0_x_z_xxxxxzzz, g_0_x_z_xxxxyyy, g_0_x_z_xxxxyyyz, g_0_x_z_xxxxyyz, g_0_x_z_xxxxyyzz, g_0_x_z_xxxxyzz, g_0_x_z_xxxxyzzz, g_0_x_z_xxxxzzz, g_0_x_z_xxxxzzzz, g_0_x_z_xxxyyyy, g_0_x_z_xxxyyyyz, g_0_x_z_xxxyyyz, g_0_x_z_xxxyyyzz, g_0_x_z_xxxyyzz, g_0_x_z_xxxyyzzz, g_0_x_z_xxxyzzz, g_0_x_z_xxxyzzzz, g_0_x_z_xxxzzzz, g_0_x_z_xxxzzzzz, g_0_x_z_xxyyyyy, g_0_x_z_xxyyyyyz, g_0_x_z_xxyyyyz, g_0_x_z_xxyyyyzz, g_0_x_z_xxyyyzz, g_0_x_z_xxyyyzzz, g_0_x_z_xxyyzzz, g_0_x_z_xxyyzzzz, g_0_x_z_xxyzzzz, g_0_x_z_xxyzzzzz, g_0_x_z_xxzzzzz, g_0_x_z_xxzzzzzz, g_0_x_z_xyyyyyy, g_0_x_z_xyyyyyyz, g_0_x_z_xyyyyyz, g_0_x_z_xyyyyyzz, g_0_x_z_xyyyyzz, g_0_x_z_xyyyyzzz, g_0_x_z_xyyyzzz, g_0_x_z_xyyyzzzz, g_0_x_z_xyyzzzz, g_0_x_z_xyyzzzzz, g_0_x_z_xyzzzzz, g_0_x_z_xyzzzzzz, g_0_x_z_xzzzzzz, g_0_x_z_xzzzzzzz, g_0_x_z_yyyyyyy, g_0_x_z_yyyyyyyz, g_0_x_z_yyyyyyz, g_0_x_z_yyyyyyzz, g_0_x_z_yyyyyzz, g_0_x_z_yyyyyzzz, g_0_x_z_yyyyzzz, g_0_x_z_yyyyzzzz, g_0_x_z_yyyzzzz, g_0_x_z_yyyzzzzz, g_0_x_z_yyzzzzz, g_0_x_z_yyzzzzzz, g_0_x_z_yzzzzzz, g_0_x_z_yzzzzzzz, g_0_x_z_zzzzzzz, g_0_x_z_zzzzzzzz, g_0_x_zz_xxxxxxx, g_0_x_zz_xxxxxxy, g_0_x_zz_xxxxxxz, g_0_x_zz_xxxxxyy, g_0_x_zz_xxxxxyz, g_0_x_zz_xxxxxzz, g_0_x_zz_xxxxyyy, g_0_x_zz_xxxxyyz, g_0_x_zz_xxxxyzz, g_0_x_zz_xxxxzzz, g_0_x_zz_xxxyyyy, g_0_x_zz_xxxyyyz, g_0_x_zz_xxxyyzz, g_0_x_zz_xxxyzzz, g_0_x_zz_xxxzzzz, g_0_x_zz_xxyyyyy, g_0_x_zz_xxyyyyz, g_0_x_zz_xxyyyzz, g_0_x_zz_xxyyzzz, g_0_x_zz_xxyzzzz, g_0_x_zz_xxzzzzz, g_0_x_zz_xyyyyyy, g_0_x_zz_xyyyyyz, g_0_x_zz_xyyyyzz, g_0_x_zz_xyyyzzz, g_0_x_zz_xyyzzzz, g_0_x_zz_xyzzzzz, g_0_x_zz_xzzzzzz, g_0_x_zz_yyyyyyy, g_0_x_zz_yyyyyyz, g_0_x_zz_yyyyyzz, g_0_x_zz_yyyyzzz, g_0_x_zz_yyyzzzz, g_0_x_zz_yyzzzzz, g_0_x_zz_yzzzzzz, g_0_x_zz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_zz_xxxxxxx[k] = -g_0_x_z_xxxxxxx[k] * ab_z + g_0_x_z_xxxxxxxz[k];

                g_0_x_zz_xxxxxxy[k] = -g_0_x_z_xxxxxxy[k] * ab_z + g_0_x_z_xxxxxxyz[k];

                g_0_x_zz_xxxxxxz[k] = -g_0_x_z_xxxxxxz[k] * ab_z + g_0_x_z_xxxxxxzz[k];

                g_0_x_zz_xxxxxyy[k] = -g_0_x_z_xxxxxyy[k] * ab_z + g_0_x_z_xxxxxyyz[k];

                g_0_x_zz_xxxxxyz[k] = -g_0_x_z_xxxxxyz[k] * ab_z + g_0_x_z_xxxxxyzz[k];

                g_0_x_zz_xxxxxzz[k] = -g_0_x_z_xxxxxzz[k] * ab_z + g_0_x_z_xxxxxzzz[k];

                g_0_x_zz_xxxxyyy[k] = -g_0_x_z_xxxxyyy[k] * ab_z + g_0_x_z_xxxxyyyz[k];

                g_0_x_zz_xxxxyyz[k] = -g_0_x_z_xxxxyyz[k] * ab_z + g_0_x_z_xxxxyyzz[k];

                g_0_x_zz_xxxxyzz[k] = -g_0_x_z_xxxxyzz[k] * ab_z + g_0_x_z_xxxxyzzz[k];

                g_0_x_zz_xxxxzzz[k] = -g_0_x_z_xxxxzzz[k] * ab_z + g_0_x_z_xxxxzzzz[k];

                g_0_x_zz_xxxyyyy[k] = -g_0_x_z_xxxyyyy[k] * ab_z + g_0_x_z_xxxyyyyz[k];

                g_0_x_zz_xxxyyyz[k] = -g_0_x_z_xxxyyyz[k] * ab_z + g_0_x_z_xxxyyyzz[k];

                g_0_x_zz_xxxyyzz[k] = -g_0_x_z_xxxyyzz[k] * ab_z + g_0_x_z_xxxyyzzz[k];

                g_0_x_zz_xxxyzzz[k] = -g_0_x_z_xxxyzzz[k] * ab_z + g_0_x_z_xxxyzzzz[k];

                g_0_x_zz_xxxzzzz[k] = -g_0_x_z_xxxzzzz[k] * ab_z + g_0_x_z_xxxzzzzz[k];

                g_0_x_zz_xxyyyyy[k] = -g_0_x_z_xxyyyyy[k] * ab_z + g_0_x_z_xxyyyyyz[k];

                g_0_x_zz_xxyyyyz[k] = -g_0_x_z_xxyyyyz[k] * ab_z + g_0_x_z_xxyyyyzz[k];

                g_0_x_zz_xxyyyzz[k] = -g_0_x_z_xxyyyzz[k] * ab_z + g_0_x_z_xxyyyzzz[k];

                g_0_x_zz_xxyyzzz[k] = -g_0_x_z_xxyyzzz[k] * ab_z + g_0_x_z_xxyyzzzz[k];

                g_0_x_zz_xxyzzzz[k] = -g_0_x_z_xxyzzzz[k] * ab_z + g_0_x_z_xxyzzzzz[k];

                g_0_x_zz_xxzzzzz[k] = -g_0_x_z_xxzzzzz[k] * ab_z + g_0_x_z_xxzzzzzz[k];

                g_0_x_zz_xyyyyyy[k] = -g_0_x_z_xyyyyyy[k] * ab_z + g_0_x_z_xyyyyyyz[k];

                g_0_x_zz_xyyyyyz[k] = -g_0_x_z_xyyyyyz[k] * ab_z + g_0_x_z_xyyyyyzz[k];

                g_0_x_zz_xyyyyzz[k] = -g_0_x_z_xyyyyzz[k] * ab_z + g_0_x_z_xyyyyzzz[k];

                g_0_x_zz_xyyyzzz[k] = -g_0_x_z_xyyyzzz[k] * ab_z + g_0_x_z_xyyyzzzz[k];

                g_0_x_zz_xyyzzzz[k] = -g_0_x_z_xyyzzzz[k] * ab_z + g_0_x_z_xyyzzzzz[k];

                g_0_x_zz_xyzzzzz[k] = -g_0_x_z_xyzzzzz[k] * ab_z + g_0_x_z_xyzzzzzz[k];

                g_0_x_zz_xzzzzzz[k] = -g_0_x_z_xzzzzzz[k] * ab_z + g_0_x_z_xzzzzzzz[k];

                g_0_x_zz_yyyyyyy[k] = -g_0_x_z_yyyyyyy[k] * ab_z + g_0_x_z_yyyyyyyz[k];

                g_0_x_zz_yyyyyyz[k] = -g_0_x_z_yyyyyyz[k] * ab_z + g_0_x_z_yyyyyyzz[k];

                g_0_x_zz_yyyyyzz[k] = -g_0_x_z_yyyyyzz[k] * ab_z + g_0_x_z_yyyyyzzz[k];

                g_0_x_zz_yyyyzzz[k] = -g_0_x_z_yyyyzzz[k] * ab_z + g_0_x_z_yyyyzzzz[k];

                g_0_x_zz_yyyzzzz[k] = -g_0_x_z_yyyzzzz[k] * ab_z + g_0_x_z_yyyzzzzz[k];

                g_0_x_zz_yyzzzzz[k] = -g_0_x_z_yyzzzzz[k] * ab_z + g_0_x_z_yyzzzzzz[k];

                g_0_x_zz_yzzzzzz[k] = -g_0_x_z_yzzzzzz[k] * ab_z + g_0_x_z_yzzzzzzz[k];

                g_0_x_zz_zzzzzzz[k] = -g_0_x_z_zzzzzzz[k] * ab_z + g_0_x_z_zzzzzzzz[k];
            }

            /// Set up 216-252 components of targeted buffer : cbuffer.data(

            auto g_0_y_xx_xxxxxxx = cbuffer.data(dk_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_y_xx_xxxxxxy = cbuffer.data(dk_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_y_xx_xxxxxxz = cbuffer.data(dk_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_y_xx_xxxxxyy = cbuffer.data(dk_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_y_xx_xxxxxyz = cbuffer.data(dk_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_y_xx_xxxxxzz = cbuffer.data(dk_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_y_xx_xxxxyyy = cbuffer.data(dk_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_y_xx_xxxxyyz = cbuffer.data(dk_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_y_xx_xxxxyzz = cbuffer.data(dk_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_y_xx_xxxxzzz = cbuffer.data(dk_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_y_xx_xxxyyyy = cbuffer.data(dk_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_y_xx_xxxyyyz = cbuffer.data(dk_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_y_xx_xxxyyzz = cbuffer.data(dk_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_y_xx_xxxyzzz = cbuffer.data(dk_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_y_xx_xxxzzzz = cbuffer.data(dk_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_y_xx_xxyyyyy = cbuffer.data(dk_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_y_xx_xxyyyyz = cbuffer.data(dk_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_y_xx_xxyyyzz = cbuffer.data(dk_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_y_xx_xxyyzzz = cbuffer.data(dk_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_y_xx_xxyzzzz = cbuffer.data(dk_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_y_xx_xxzzzzz = cbuffer.data(dk_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_y_xx_xyyyyyy = cbuffer.data(dk_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_y_xx_xyyyyyz = cbuffer.data(dk_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_y_xx_xyyyyzz = cbuffer.data(dk_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_y_xx_xyyyzzz = cbuffer.data(dk_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_y_xx_xyyzzzz = cbuffer.data(dk_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_y_xx_xyzzzzz = cbuffer.data(dk_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_y_xx_xzzzzzz = cbuffer.data(dk_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_y_xx_yyyyyyy = cbuffer.data(dk_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_y_xx_yyyyyyz = cbuffer.data(dk_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_y_xx_yyyyyzz = cbuffer.data(dk_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_y_xx_yyyyzzz = cbuffer.data(dk_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_y_xx_yyyzzzz = cbuffer.data(dk_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_y_xx_yyzzzzz = cbuffer.data(dk_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_y_xx_yzzzzzz = cbuffer.data(dk_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_y_xx_zzzzzzz = cbuffer.data(dk_geom_01_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_x_xxxxxxx, g_0_y_x_xxxxxxxx, g_0_y_x_xxxxxxxy, g_0_y_x_xxxxxxxz, g_0_y_x_xxxxxxy, g_0_y_x_xxxxxxyy, g_0_y_x_xxxxxxyz, g_0_y_x_xxxxxxz, g_0_y_x_xxxxxxzz, g_0_y_x_xxxxxyy, g_0_y_x_xxxxxyyy, g_0_y_x_xxxxxyyz, g_0_y_x_xxxxxyz, g_0_y_x_xxxxxyzz, g_0_y_x_xxxxxzz, g_0_y_x_xxxxxzzz, g_0_y_x_xxxxyyy, g_0_y_x_xxxxyyyy, g_0_y_x_xxxxyyyz, g_0_y_x_xxxxyyz, g_0_y_x_xxxxyyzz, g_0_y_x_xxxxyzz, g_0_y_x_xxxxyzzz, g_0_y_x_xxxxzzz, g_0_y_x_xxxxzzzz, g_0_y_x_xxxyyyy, g_0_y_x_xxxyyyyy, g_0_y_x_xxxyyyyz, g_0_y_x_xxxyyyz, g_0_y_x_xxxyyyzz, g_0_y_x_xxxyyzz, g_0_y_x_xxxyyzzz, g_0_y_x_xxxyzzz, g_0_y_x_xxxyzzzz, g_0_y_x_xxxzzzz, g_0_y_x_xxxzzzzz, g_0_y_x_xxyyyyy, g_0_y_x_xxyyyyyy, g_0_y_x_xxyyyyyz, g_0_y_x_xxyyyyz, g_0_y_x_xxyyyyzz, g_0_y_x_xxyyyzz, g_0_y_x_xxyyyzzz, g_0_y_x_xxyyzzz, g_0_y_x_xxyyzzzz, g_0_y_x_xxyzzzz, g_0_y_x_xxyzzzzz, g_0_y_x_xxzzzzz, g_0_y_x_xxzzzzzz, g_0_y_x_xyyyyyy, g_0_y_x_xyyyyyyy, g_0_y_x_xyyyyyyz, g_0_y_x_xyyyyyz, g_0_y_x_xyyyyyzz, g_0_y_x_xyyyyzz, g_0_y_x_xyyyyzzz, g_0_y_x_xyyyzzz, g_0_y_x_xyyyzzzz, g_0_y_x_xyyzzzz, g_0_y_x_xyyzzzzz, g_0_y_x_xyzzzzz, g_0_y_x_xyzzzzzz, g_0_y_x_xzzzzzz, g_0_y_x_xzzzzzzz, g_0_y_x_yyyyyyy, g_0_y_x_yyyyyyz, g_0_y_x_yyyyyzz, g_0_y_x_yyyyzzz, g_0_y_x_yyyzzzz, g_0_y_x_yyzzzzz, g_0_y_x_yzzzzzz, g_0_y_x_zzzzzzz, g_0_y_xx_xxxxxxx, g_0_y_xx_xxxxxxy, g_0_y_xx_xxxxxxz, g_0_y_xx_xxxxxyy, g_0_y_xx_xxxxxyz, g_0_y_xx_xxxxxzz, g_0_y_xx_xxxxyyy, g_0_y_xx_xxxxyyz, g_0_y_xx_xxxxyzz, g_0_y_xx_xxxxzzz, g_0_y_xx_xxxyyyy, g_0_y_xx_xxxyyyz, g_0_y_xx_xxxyyzz, g_0_y_xx_xxxyzzz, g_0_y_xx_xxxzzzz, g_0_y_xx_xxyyyyy, g_0_y_xx_xxyyyyz, g_0_y_xx_xxyyyzz, g_0_y_xx_xxyyzzz, g_0_y_xx_xxyzzzz, g_0_y_xx_xxzzzzz, g_0_y_xx_xyyyyyy, g_0_y_xx_xyyyyyz, g_0_y_xx_xyyyyzz, g_0_y_xx_xyyyzzz, g_0_y_xx_xyyzzzz, g_0_y_xx_xyzzzzz, g_0_y_xx_xzzzzzz, g_0_y_xx_yyyyyyy, g_0_y_xx_yyyyyyz, g_0_y_xx_yyyyyzz, g_0_y_xx_yyyyzzz, g_0_y_xx_yyyzzzz, g_0_y_xx_yyzzzzz, g_0_y_xx_yzzzzzz, g_0_y_xx_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xx_xxxxxxx[k] = -g_0_y_x_xxxxxxx[k] * ab_x + g_0_y_x_xxxxxxxx[k];

                g_0_y_xx_xxxxxxy[k] = -g_0_y_x_xxxxxxy[k] * ab_x + g_0_y_x_xxxxxxxy[k];

                g_0_y_xx_xxxxxxz[k] = -g_0_y_x_xxxxxxz[k] * ab_x + g_0_y_x_xxxxxxxz[k];

                g_0_y_xx_xxxxxyy[k] = -g_0_y_x_xxxxxyy[k] * ab_x + g_0_y_x_xxxxxxyy[k];

                g_0_y_xx_xxxxxyz[k] = -g_0_y_x_xxxxxyz[k] * ab_x + g_0_y_x_xxxxxxyz[k];

                g_0_y_xx_xxxxxzz[k] = -g_0_y_x_xxxxxzz[k] * ab_x + g_0_y_x_xxxxxxzz[k];

                g_0_y_xx_xxxxyyy[k] = -g_0_y_x_xxxxyyy[k] * ab_x + g_0_y_x_xxxxxyyy[k];

                g_0_y_xx_xxxxyyz[k] = -g_0_y_x_xxxxyyz[k] * ab_x + g_0_y_x_xxxxxyyz[k];

                g_0_y_xx_xxxxyzz[k] = -g_0_y_x_xxxxyzz[k] * ab_x + g_0_y_x_xxxxxyzz[k];

                g_0_y_xx_xxxxzzz[k] = -g_0_y_x_xxxxzzz[k] * ab_x + g_0_y_x_xxxxxzzz[k];

                g_0_y_xx_xxxyyyy[k] = -g_0_y_x_xxxyyyy[k] * ab_x + g_0_y_x_xxxxyyyy[k];

                g_0_y_xx_xxxyyyz[k] = -g_0_y_x_xxxyyyz[k] * ab_x + g_0_y_x_xxxxyyyz[k];

                g_0_y_xx_xxxyyzz[k] = -g_0_y_x_xxxyyzz[k] * ab_x + g_0_y_x_xxxxyyzz[k];

                g_0_y_xx_xxxyzzz[k] = -g_0_y_x_xxxyzzz[k] * ab_x + g_0_y_x_xxxxyzzz[k];

                g_0_y_xx_xxxzzzz[k] = -g_0_y_x_xxxzzzz[k] * ab_x + g_0_y_x_xxxxzzzz[k];

                g_0_y_xx_xxyyyyy[k] = -g_0_y_x_xxyyyyy[k] * ab_x + g_0_y_x_xxxyyyyy[k];

                g_0_y_xx_xxyyyyz[k] = -g_0_y_x_xxyyyyz[k] * ab_x + g_0_y_x_xxxyyyyz[k];

                g_0_y_xx_xxyyyzz[k] = -g_0_y_x_xxyyyzz[k] * ab_x + g_0_y_x_xxxyyyzz[k];

                g_0_y_xx_xxyyzzz[k] = -g_0_y_x_xxyyzzz[k] * ab_x + g_0_y_x_xxxyyzzz[k];

                g_0_y_xx_xxyzzzz[k] = -g_0_y_x_xxyzzzz[k] * ab_x + g_0_y_x_xxxyzzzz[k];

                g_0_y_xx_xxzzzzz[k] = -g_0_y_x_xxzzzzz[k] * ab_x + g_0_y_x_xxxzzzzz[k];

                g_0_y_xx_xyyyyyy[k] = -g_0_y_x_xyyyyyy[k] * ab_x + g_0_y_x_xxyyyyyy[k];

                g_0_y_xx_xyyyyyz[k] = -g_0_y_x_xyyyyyz[k] * ab_x + g_0_y_x_xxyyyyyz[k];

                g_0_y_xx_xyyyyzz[k] = -g_0_y_x_xyyyyzz[k] * ab_x + g_0_y_x_xxyyyyzz[k];

                g_0_y_xx_xyyyzzz[k] = -g_0_y_x_xyyyzzz[k] * ab_x + g_0_y_x_xxyyyzzz[k];

                g_0_y_xx_xyyzzzz[k] = -g_0_y_x_xyyzzzz[k] * ab_x + g_0_y_x_xxyyzzzz[k];

                g_0_y_xx_xyzzzzz[k] = -g_0_y_x_xyzzzzz[k] * ab_x + g_0_y_x_xxyzzzzz[k];

                g_0_y_xx_xzzzzzz[k] = -g_0_y_x_xzzzzzz[k] * ab_x + g_0_y_x_xxzzzzzz[k];

                g_0_y_xx_yyyyyyy[k] = -g_0_y_x_yyyyyyy[k] * ab_x + g_0_y_x_xyyyyyyy[k];

                g_0_y_xx_yyyyyyz[k] = -g_0_y_x_yyyyyyz[k] * ab_x + g_0_y_x_xyyyyyyz[k];

                g_0_y_xx_yyyyyzz[k] = -g_0_y_x_yyyyyzz[k] * ab_x + g_0_y_x_xyyyyyzz[k];

                g_0_y_xx_yyyyzzz[k] = -g_0_y_x_yyyyzzz[k] * ab_x + g_0_y_x_xyyyyzzz[k];

                g_0_y_xx_yyyzzzz[k] = -g_0_y_x_yyyzzzz[k] * ab_x + g_0_y_x_xyyyzzzz[k];

                g_0_y_xx_yyzzzzz[k] = -g_0_y_x_yyzzzzz[k] * ab_x + g_0_y_x_xyyzzzzz[k];

                g_0_y_xx_yzzzzzz[k] = -g_0_y_x_yzzzzzz[k] * ab_x + g_0_y_x_xyzzzzzz[k];

                g_0_y_xx_zzzzzzz[k] = -g_0_y_x_zzzzzzz[k] * ab_x + g_0_y_x_xzzzzzzz[k];
            }

            /// Set up 252-288 components of targeted buffer : cbuffer.data(

            auto g_0_y_xy_xxxxxxx = cbuffer.data(dk_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_y_xy_xxxxxxy = cbuffer.data(dk_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_y_xy_xxxxxxz = cbuffer.data(dk_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_y_xy_xxxxxyy = cbuffer.data(dk_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_y_xy_xxxxxyz = cbuffer.data(dk_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_y_xy_xxxxxzz = cbuffer.data(dk_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_y_xy_xxxxyyy = cbuffer.data(dk_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_y_xy_xxxxyyz = cbuffer.data(dk_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_y_xy_xxxxyzz = cbuffer.data(dk_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_y_xy_xxxxzzz = cbuffer.data(dk_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_y_xy_xxxyyyy = cbuffer.data(dk_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_y_xy_xxxyyyz = cbuffer.data(dk_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_y_xy_xxxyyzz = cbuffer.data(dk_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_y_xy_xxxyzzz = cbuffer.data(dk_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_y_xy_xxxzzzz = cbuffer.data(dk_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_y_xy_xxyyyyy = cbuffer.data(dk_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_y_xy_xxyyyyz = cbuffer.data(dk_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_y_xy_xxyyyzz = cbuffer.data(dk_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_y_xy_xxyyzzz = cbuffer.data(dk_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_y_xy_xxyzzzz = cbuffer.data(dk_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_y_xy_xxzzzzz = cbuffer.data(dk_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_y_xy_xyyyyyy = cbuffer.data(dk_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_y_xy_xyyyyyz = cbuffer.data(dk_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_y_xy_xyyyyzz = cbuffer.data(dk_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_y_xy_xyyyzzz = cbuffer.data(dk_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_y_xy_xyyzzzz = cbuffer.data(dk_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_y_xy_xyzzzzz = cbuffer.data(dk_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_y_xy_xzzzzzz = cbuffer.data(dk_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_y_xy_yyyyyyy = cbuffer.data(dk_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_y_xy_yyyyyyz = cbuffer.data(dk_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_y_xy_yyyyyzz = cbuffer.data(dk_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_y_xy_yyyyzzz = cbuffer.data(dk_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_y_xy_yyyzzzz = cbuffer.data(dk_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_y_xy_yyzzzzz = cbuffer.data(dk_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_y_xy_yzzzzzz = cbuffer.data(dk_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_y_xy_zzzzzzz = cbuffer.data(dk_geom_01_off + 287 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xy_xxxxxxx, g_0_y_xy_xxxxxxy, g_0_y_xy_xxxxxxz, g_0_y_xy_xxxxxyy, g_0_y_xy_xxxxxyz, g_0_y_xy_xxxxxzz, g_0_y_xy_xxxxyyy, g_0_y_xy_xxxxyyz, g_0_y_xy_xxxxyzz, g_0_y_xy_xxxxzzz, g_0_y_xy_xxxyyyy, g_0_y_xy_xxxyyyz, g_0_y_xy_xxxyyzz, g_0_y_xy_xxxyzzz, g_0_y_xy_xxxzzzz, g_0_y_xy_xxyyyyy, g_0_y_xy_xxyyyyz, g_0_y_xy_xxyyyzz, g_0_y_xy_xxyyzzz, g_0_y_xy_xxyzzzz, g_0_y_xy_xxzzzzz, g_0_y_xy_xyyyyyy, g_0_y_xy_xyyyyyz, g_0_y_xy_xyyyyzz, g_0_y_xy_xyyyzzz, g_0_y_xy_xyyzzzz, g_0_y_xy_xyzzzzz, g_0_y_xy_xzzzzzz, g_0_y_xy_yyyyyyy, g_0_y_xy_yyyyyyz, g_0_y_xy_yyyyyzz, g_0_y_xy_yyyyzzz, g_0_y_xy_yyyzzzz, g_0_y_xy_yyzzzzz, g_0_y_xy_yzzzzzz, g_0_y_xy_zzzzzzz, g_0_y_y_xxxxxxx, g_0_y_y_xxxxxxxx, g_0_y_y_xxxxxxxy, g_0_y_y_xxxxxxxz, g_0_y_y_xxxxxxy, g_0_y_y_xxxxxxyy, g_0_y_y_xxxxxxyz, g_0_y_y_xxxxxxz, g_0_y_y_xxxxxxzz, g_0_y_y_xxxxxyy, g_0_y_y_xxxxxyyy, g_0_y_y_xxxxxyyz, g_0_y_y_xxxxxyz, g_0_y_y_xxxxxyzz, g_0_y_y_xxxxxzz, g_0_y_y_xxxxxzzz, g_0_y_y_xxxxyyy, g_0_y_y_xxxxyyyy, g_0_y_y_xxxxyyyz, g_0_y_y_xxxxyyz, g_0_y_y_xxxxyyzz, g_0_y_y_xxxxyzz, g_0_y_y_xxxxyzzz, g_0_y_y_xxxxzzz, g_0_y_y_xxxxzzzz, g_0_y_y_xxxyyyy, g_0_y_y_xxxyyyyy, g_0_y_y_xxxyyyyz, g_0_y_y_xxxyyyz, g_0_y_y_xxxyyyzz, g_0_y_y_xxxyyzz, g_0_y_y_xxxyyzzz, g_0_y_y_xxxyzzz, g_0_y_y_xxxyzzzz, g_0_y_y_xxxzzzz, g_0_y_y_xxxzzzzz, g_0_y_y_xxyyyyy, g_0_y_y_xxyyyyyy, g_0_y_y_xxyyyyyz, g_0_y_y_xxyyyyz, g_0_y_y_xxyyyyzz, g_0_y_y_xxyyyzz, g_0_y_y_xxyyyzzz, g_0_y_y_xxyyzzz, g_0_y_y_xxyyzzzz, g_0_y_y_xxyzzzz, g_0_y_y_xxyzzzzz, g_0_y_y_xxzzzzz, g_0_y_y_xxzzzzzz, g_0_y_y_xyyyyyy, g_0_y_y_xyyyyyyy, g_0_y_y_xyyyyyyz, g_0_y_y_xyyyyyz, g_0_y_y_xyyyyyzz, g_0_y_y_xyyyyzz, g_0_y_y_xyyyyzzz, g_0_y_y_xyyyzzz, g_0_y_y_xyyyzzzz, g_0_y_y_xyyzzzz, g_0_y_y_xyyzzzzz, g_0_y_y_xyzzzzz, g_0_y_y_xyzzzzzz, g_0_y_y_xzzzzzz, g_0_y_y_xzzzzzzz, g_0_y_y_yyyyyyy, g_0_y_y_yyyyyyz, g_0_y_y_yyyyyzz, g_0_y_y_yyyyzzz, g_0_y_y_yyyzzzz, g_0_y_y_yyzzzzz, g_0_y_y_yzzzzzz, g_0_y_y_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xy_xxxxxxx[k] = -g_0_y_y_xxxxxxx[k] * ab_x + g_0_y_y_xxxxxxxx[k];

                g_0_y_xy_xxxxxxy[k] = -g_0_y_y_xxxxxxy[k] * ab_x + g_0_y_y_xxxxxxxy[k];

                g_0_y_xy_xxxxxxz[k] = -g_0_y_y_xxxxxxz[k] * ab_x + g_0_y_y_xxxxxxxz[k];

                g_0_y_xy_xxxxxyy[k] = -g_0_y_y_xxxxxyy[k] * ab_x + g_0_y_y_xxxxxxyy[k];

                g_0_y_xy_xxxxxyz[k] = -g_0_y_y_xxxxxyz[k] * ab_x + g_0_y_y_xxxxxxyz[k];

                g_0_y_xy_xxxxxzz[k] = -g_0_y_y_xxxxxzz[k] * ab_x + g_0_y_y_xxxxxxzz[k];

                g_0_y_xy_xxxxyyy[k] = -g_0_y_y_xxxxyyy[k] * ab_x + g_0_y_y_xxxxxyyy[k];

                g_0_y_xy_xxxxyyz[k] = -g_0_y_y_xxxxyyz[k] * ab_x + g_0_y_y_xxxxxyyz[k];

                g_0_y_xy_xxxxyzz[k] = -g_0_y_y_xxxxyzz[k] * ab_x + g_0_y_y_xxxxxyzz[k];

                g_0_y_xy_xxxxzzz[k] = -g_0_y_y_xxxxzzz[k] * ab_x + g_0_y_y_xxxxxzzz[k];

                g_0_y_xy_xxxyyyy[k] = -g_0_y_y_xxxyyyy[k] * ab_x + g_0_y_y_xxxxyyyy[k];

                g_0_y_xy_xxxyyyz[k] = -g_0_y_y_xxxyyyz[k] * ab_x + g_0_y_y_xxxxyyyz[k];

                g_0_y_xy_xxxyyzz[k] = -g_0_y_y_xxxyyzz[k] * ab_x + g_0_y_y_xxxxyyzz[k];

                g_0_y_xy_xxxyzzz[k] = -g_0_y_y_xxxyzzz[k] * ab_x + g_0_y_y_xxxxyzzz[k];

                g_0_y_xy_xxxzzzz[k] = -g_0_y_y_xxxzzzz[k] * ab_x + g_0_y_y_xxxxzzzz[k];

                g_0_y_xy_xxyyyyy[k] = -g_0_y_y_xxyyyyy[k] * ab_x + g_0_y_y_xxxyyyyy[k];

                g_0_y_xy_xxyyyyz[k] = -g_0_y_y_xxyyyyz[k] * ab_x + g_0_y_y_xxxyyyyz[k];

                g_0_y_xy_xxyyyzz[k] = -g_0_y_y_xxyyyzz[k] * ab_x + g_0_y_y_xxxyyyzz[k];

                g_0_y_xy_xxyyzzz[k] = -g_0_y_y_xxyyzzz[k] * ab_x + g_0_y_y_xxxyyzzz[k];

                g_0_y_xy_xxyzzzz[k] = -g_0_y_y_xxyzzzz[k] * ab_x + g_0_y_y_xxxyzzzz[k];

                g_0_y_xy_xxzzzzz[k] = -g_0_y_y_xxzzzzz[k] * ab_x + g_0_y_y_xxxzzzzz[k];

                g_0_y_xy_xyyyyyy[k] = -g_0_y_y_xyyyyyy[k] * ab_x + g_0_y_y_xxyyyyyy[k];

                g_0_y_xy_xyyyyyz[k] = -g_0_y_y_xyyyyyz[k] * ab_x + g_0_y_y_xxyyyyyz[k];

                g_0_y_xy_xyyyyzz[k] = -g_0_y_y_xyyyyzz[k] * ab_x + g_0_y_y_xxyyyyzz[k];

                g_0_y_xy_xyyyzzz[k] = -g_0_y_y_xyyyzzz[k] * ab_x + g_0_y_y_xxyyyzzz[k];

                g_0_y_xy_xyyzzzz[k] = -g_0_y_y_xyyzzzz[k] * ab_x + g_0_y_y_xxyyzzzz[k];

                g_0_y_xy_xyzzzzz[k] = -g_0_y_y_xyzzzzz[k] * ab_x + g_0_y_y_xxyzzzzz[k];

                g_0_y_xy_xzzzzzz[k] = -g_0_y_y_xzzzzzz[k] * ab_x + g_0_y_y_xxzzzzzz[k];

                g_0_y_xy_yyyyyyy[k] = -g_0_y_y_yyyyyyy[k] * ab_x + g_0_y_y_xyyyyyyy[k];

                g_0_y_xy_yyyyyyz[k] = -g_0_y_y_yyyyyyz[k] * ab_x + g_0_y_y_xyyyyyyz[k];

                g_0_y_xy_yyyyyzz[k] = -g_0_y_y_yyyyyzz[k] * ab_x + g_0_y_y_xyyyyyzz[k];

                g_0_y_xy_yyyyzzz[k] = -g_0_y_y_yyyyzzz[k] * ab_x + g_0_y_y_xyyyyzzz[k];

                g_0_y_xy_yyyzzzz[k] = -g_0_y_y_yyyzzzz[k] * ab_x + g_0_y_y_xyyyzzzz[k];

                g_0_y_xy_yyzzzzz[k] = -g_0_y_y_yyzzzzz[k] * ab_x + g_0_y_y_xyyzzzzz[k];

                g_0_y_xy_yzzzzzz[k] = -g_0_y_y_yzzzzzz[k] * ab_x + g_0_y_y_xyzzzzzz[k];

                g_0_y_xy_zzzzzzz[k] = -g_0_y_y_zzzzzzz[k] * ab_x + g_0_y_y_xzzzzzzz[k];
            }

            /// Set up 288-324 components of targeted buffer : cbuffer.data(

            auto g_0_y_xz_xxxxxxx = cbuffer.data(dk_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_y_xz_xxxxxxy = cbuffer.data(dk_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_y_xz_xxxxxxz = cbuffer.data(dk_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_y_xz_xxxxxyy = cbuffer.data(dk_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_y_xz_xxxxxyz = cbuffer.data(dk_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_y_xz_xxxxxzz = cbuffer.data(dk_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_y_xz_xxxxyyy = cbuffer.data(dk_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_y_xz_xxxxyyz = cbuffer.data(dk_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_y_xz_xxxxyzz = cbuffer.data(dk_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_y_xz_xxxxzzz = cbuffer.data(dk_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_y_xz_xxxyyyy = cbuffer.data(dk_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_y_xz_xxxyyyz = cbuffer.data(dk_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_y_xz_xxxyyzz = cbuffer.data(dk_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_y_xz_xxxyzzz = cbuffer.data(dk_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_y_xz_xxxzzzz = cbuffer.data(dk_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_y_xz_xxyyyyy = cbuffer.data(dk_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_y_xz_xxyyyyz = cbuffer.data(dk_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_y_xz_xxyyyzz = cbuffer.data(dk_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_y_xz_xxyyzzz = cbuffer.data(dk_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_y_xz_xxyzzzz = cbuffer.data(dk_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_y_xz_xxzzzzz = cbuffer.data(dk_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_y_xz_xyyyyyy = cbuffer.data(dk_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_y_xz_xyyyyyz = cbuffer.data(dk_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_y_xz_xyyyyzz = cbuffer.data(dk_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_y_xz_xyyyzzz = cbuffer.data(dk_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_y_xz_xyyzzzz = cbuffer.data(dk_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_y_xz_xyzzzzz = cbuffer.data(dk_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_y_xz_xzzzzzz = cbuffer.data(dk_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_y_xz_yyyyyyy = cbuffer.data(dk_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_y_xz_yyyyyyz = cbuffer.data(dk_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_y_xz_yyyyyzz = cbuffer.data(dk_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_y_xz_yyyyzzz = cbuffer.data(dk_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_y_xz_yyyzzzz = cbuffer.data(dk_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_y_xz_yyzzzzz = cbuffer.data(dk_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_y_xz_yzzzzzz = cbuffer.data(dk_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_y_xz_zzzzzzz = cbuffer.data(dk_geom_01_off + 323 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xz_xxxxxxx, g_0_y_xz_xxxxxxy, g_0_y_xz_xxxxxxz, g_0_y_xz_xxxxxyy, g_0_y_xz_xxxxxyz, g_0_y_xz_xxxxxzz, g_0_y_xz_xxxxyyy, g_0_y_xz_xxxxyyz, g_0_y_xz_xxxxyzz, g_0_y_xz_xxxxzzz, g_0_y_xz_xxxyyyy, g_0_y_xz_xxxyyyz, g_0_y_xz_xxxyyzz, g_0_y_xz_xxxyzzz, g_0_y_xz_xxxzzzz, g_0_y_xz_xxyyyyy, g_0_y_xz_xxyyyyz, g_0_y_xz_xxyyyzz, g_0_y_xz_xxyyzzz, g_0_y_xz_xxyzzzz, g_0_y_xz_xxzzzzz, g_0_y_xz_xyyyyyy, g_0_y_xz_xyyyyyz, g_0_y_xz_xyyyyzz, g_0_y_xz_xyyyzzz, g_0_y_xz_xyyzzzz, g_0_y_xz_xyzzzzz, g_0_y_xz_xzzzzzz, g_0_y_xz_yyyyyyy, g_0_y_xz_yyyyyyz, g_0_y_xz_yyyyyzz, g_0_y_xz_yyyyzzz, g_0_y_xz_yyyzzzz, g_0_y_xz_yyzzzzz, g_0_y_xz_yzzzzzz, g_0_y_xz_zzzzzzz, g_0_y_z_xxxxxxx, g_0_y_z_xxxxxxxx, g_0_y_z_xxxxxxxy, g_0_y_z_xxxxxxxz, g_0_y_z_xxxxxxy, g_0_y_z_xxxxxxyy, g_0_y_z_xxxxxxyz, g_0_y_z_xxxxxxz, g_0_y_z_xxxxxxzz, g_0_y_z_xxxxxyy, g_0_y_z_xxxxxyyy, g_0_y_z_xxxxxyyz, g_0_y_z_xxxxxyz, g_0_y_z_xxxxxyzz, g_0_y_z_xxxxxzz, g_0_y_z_xxxxxzzz, g_0_y_z_xxxxyyy, g_0_y_z_xxxxyyyy, g_0_y_z_xxxxyyyz, g_0_y_z_xxxxyyz, g_0_y_z_xxxxyyzz, g_0_y_z_xxxxyzz, g_0_y_z_xxxxyzzz, g_0_y_z_xxxxzzz, g_0_y_z_xxxxzzzz, g_0_y_z_xxxyyyy, g_0_y_z_xxxyyyyy, g_0_y_z_xxxyyyyz, g_0_y_z_xxxyyyz, g_0_y_z_xxxyyyzz, g_0_y_z_xxxyyzz, g_0_y_z_xxxyyzzz, g_0_y_z_xxxyzzz, g_0_y_z_xxxyzzzz, g_0_y_z_xxxzzzz, g_0_y_z_xxxzzzzz, g_0_y_z_xxyyyyy, g_0_y_z_xxyyyyyy, g_0_y_z_xxyyyyyz, g_0_y_z_xxyyyyz, g_0_y_z_xxyyyyzz, g_0_y_z_xxyyyzz, g_0_y_z_xxyyyzzz, g_0_y_z_xxyyzzz, g_0_y_z_xxyyzzzz, g_0_y_z_xxyzzzz, g_0_y_z_xxyzzzzz, g_0_y_z_xxzzzzz, g_0_y_z_xxzzzzzz, g_0_y_z_xyyyyyy, g_0_y_z_xyyyyyyy, g_0_y_z_xyyyyyyz, g_0_y_z_xyyyyyz, g_0_y_z_xyyyyyzz, g_0_y_z_xyyyyzz, g_0_y_z_xyyyyzzz, g_0_y_z_xyyyzzz, g_0_y_z_xyyyzzzz, g_0_y_z_xyyzzzz, g_0_y_z_xyyzzzzz, g_0_y_z_xyzzzzz, g_0_y_z_xyzzzzzz, g_0_y_z_xzzzzzz, g_0_y_z_xzzzzzzz, g_0_y_z_yyyyyyy, g_0_y_z_yyyyyyz, g_0_y_z_yyyyyzz, g_0_y_z_yyyyzzz, g_0_y_z_yyyzzzz, g_0_y_z_yyzzzzz, g_0_y_z_yzzzzzz, g_0_y_z_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xz_xxxxxxx[k] = -g_0_y_z_xxxxxxx[k] * ab_x + g_0_y_z_xxxxxxxx[k];

                g_0_y_xz_xxxxxxy[k] = -g_0_y_z_xxxxxxy[k] * ab_x + g_0_y_z_xxxxxxxy[k];

                g_0_y_xz_xxxxxxz[k] = -g_0_y_z_xxxxxxz[k] * ab_x + g_0_y_z_xxxxxxxz[k];

                g_0_y_xz_xxxxxyy[k] = -g_0_y_z_xxxxxyy[k] * ab_x + g_0_y_z_xxxxxxyy[k];

                g_0_y_xz_xxxxxyz[k] = -g_0_y_z_xxxxxyz[k] * ab_x + g_0_y_z_xxxxxxyz[k];

                g_0_y_xz_xxxxxzz[k] = -g_0_y_z_xxxxxzz[k] * ab_x + g_0_y_z_xxxxxxzz[k];

                g_0_y_xz_xxxxyyy[k] = -g_0_y_z_xxxxyyy[k] * ab_x + g_0_y_z_xxxxxyyy[k];

                g_0_y_xz_xxxxyyz[k] = -g_0_y_z_xxxxyyz[k] * ab_x + g_0_y_z_xxxxxyyz[k];

                g_0_y_xz_xxxxyzz[k] = -g_0_y_z_xxxxyzz[k] * ab_x + g_0_y_z_xxxxxyzz[k];

                g_0_y_xz_xxxxzzz[k] = -g_0_y_z_xxxxzzz[k] * ab_x + g_0_y_z_xxxxxzzz[k];

                g_0_y_xz_xxxyyyy[k] = -g_0_y_z_xxxyyyy[k] * ab_x + g_0_y_z_xxxxyyyy[k];

                g_0_y_xz_xxxyyyz[k] = -g_0_y_z_xxxyyyz[k] * ab_x + g_0_y_z_xxxxyyyz[k];

                g_0_y_xz_xxxyyzz[k] = -g_0_y_z_xxxyyzz[k] * ab_x + g_0_y_z_xxxxyyzz[k];

                g_0_y_xz_xxxyzzz[k] = -g_0_y_z_xxxyzzz[k] * ab_x + g_0_y_z_xxxxyzzz[k];

                g_0_y_xz_xxxzzzz[k] = -g_0_y_z_xxxzzzz[k] * ab_x + g_0_y_z_xxxxzzzz[k];

                g_0_y_xz_xxyyyyy[k] = -g_0_y_z_xxyyyyy[k] * ab_x + g_0_y_z_xxxyyyyy[k];

                g_0_y_xz_xxyyyyz[k] = -g_0_y_z_xxyyyyz[k] * ab_x + g_0_y_z_xxxyyyyz[k];

                g_0_y_xz_xxyyyzz[k] = -g_0_y_z_xxyyyzz[k] * ab_x + g_0_y_z_xxxyyyzz[k];

                g_0_y_xz_xxyyzzz[k] = -g_0_y_z_xxyyzzz[k] * ab_x + g_0_y_z_xxxyyzzz[k];

                g_0_y_xz_xxyzzzz[k] = -g_0_y_z_xxyzzzz[k] * ab_x + g_0_y_z_xxxyzzzz[k];

                g_0_y_xz_xxzzzzz[k] = -g_0_y_z_xxzzzzz[k] * ab_x + g_0_y_z_xxxzzzzz[k];

                g_0_y_xz_xyyyyyy[k] = -g_0_y_z_xyyyyyy[k] * ab_x + g_0_y_z_xxyyyyyy[k];

                g_0_y_xz_xyyyyyz[k] = -g_0_y_z_xyyyyyz[k] * ab_x + g_0_y_z_xxyyyyyz[k];

                g_0_y_xz_xyyyyzz[k] = -g_0_y_z_xyyyyzz[k] * ab_x + g_0_y_z_xxyyyyzz[k];

                g_0_y_xz_xyyyzzz[k] = -g_0_y_z_xyyyzzz[k] * ab_x + g_0_y_z_xxyyyzzz[k];

                g_0_y_xz_xyyzzzz[k] = -g_0_y_z_xyyzzzz[k] * ab_x + g_0_y_z_xxyyzzzz[k];

                g_0_y_xz_xyzzzzz[k] = -g_0_y_z_xyzzzzz[k] * ab_x + g_0_y_z_xxyzzzzz[k];

                g_0_y_xz_xzzzzzz[k] = -g_0_y_z_xzzzzzz[k] * ab_x + g_0_y_z_xxzzzzzz[k];

                g_0_y_xz_yyyyyyy[k] = -g_0_y_z_yyyyyyy[k] * ab_x + g_0_y_z_xyyyyyyy[k];

                g_0_y_xz_yyyyyyz[k] = -g_0_y_z_yyyyyyz[k] * ab_x + g_0_y_z_xyyyyyyz[k];

                g_0_y_xz_yyyyyzz[k] = -g_0_y_z_yyyyyzz[k] * ab_x + g_0_y_z_xyyyyyzz[k];

                g_0_y_xz_yyyyzzz[k] = -g_0_y_z_yyyyzzz[k] * ab_x + g_0_y_z_xyyyyzzz[k];

                g_0_y_xz_yyyzzzz[k] = -g_0_y_z_yyyzzzz[k] * ab_x + g_0_y_z_xyyyzzzz[k];

                g_0_y_xz_yyzzzzz[k] = -g_0_y_z_yyzzzzz[k] * ab_x + g_0_y_z_xyyzzzzz[k];

                g_0_y_xz_yzzzzzz[k] = -g_0_y_z_yzzzzzz[k] * ab_x + g_0_y_z_xyzzzzzz[k];

                g_0_y_xz_zzzzzzz[k] = -g_0_y_z_zzzzzzz[k] * ab_x + g_0_y_z_xzzzzzzz[k];
            }

            /// Set up 324-360 components of targeted buffer : cbuffer.data(

            auto g_0_y_yy_xxxxxxx = cbuffer.data(dk_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_y_yy_xxxxxxy = cbuffer.data(dk_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_y_yy_xxxxxxz = cbuffer.data(dk_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_y_yy_xxxxxyy = cbuffer.data(dk_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_y_yy_xxxxxyz = cbuffer.data(dk_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_y_yy_xxxxxzz = cbuffer.data(dk_geom_01_off + 329 * ccomps * dcomps);

            auto g_0_y_yy_xxxxyyy = cbuffer.data(dk_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_y_yy_xxxxyyz = cbuffer.data(dk_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_y_yy_xxxxyzz = cbuffer.data(dk_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_y_yy_xxxxzzz = cbuffer.data(dk_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_y_yy_xxxyyyy = cbuffer.data(dk_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_y_yy_xxxyyyz = cbuffer.data(dk_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_y_yy_xxxyyzz = cbuffer.data(dk_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_y_yy_xxxyzzz = cbuffer.data(dk_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_y_yy_xxxzzzz = cbuffer.data(dk_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_y_yy_xxyyyyy = cbuffer.data(dk_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_y_yy_xxyyyyz = cbuffer.data(dk_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_y_yy_xxyyyzz = cbuffer.data(dk_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_y_yy_xxyyzzz = cbuffer.data(dk_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_y_yy_xxyzzzz = cbuffer.data(dk_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_y_yy_xxzzzzz = cbuffer.data(dk_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_y_yy_xyyyyyy = cbuffer.data(dk_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_y_yy_xyyyyyz = cbuffer.data(dk_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_y_yy_xyyyyzz = cbuffer.data(dk_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_y_yy_xyyyzzz = cbuffer.data(dk_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_y_yy_xyyzzzz = cbuffer.data(dk_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_y_yy_xyzzzzz = cbuffer.data(dk_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_y_yy_xzzzzzz = cbuffer.data(dk_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_y_yy_yyyyyyy = cbuffer.data(dk_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_y_yy_yyyyyyz = cbuffer.data(dk_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_y_yy_yyyyyzz = cbuffer.data(dk_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_y_yy_yyyyzzz = cbuffer.data(dk_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_y_yy_yyyzzzz = cbuffer.data(dk_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_y_yy_yyzzzzz = cbuffer.data(dk_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_y_yy_yzzzzzz = cbuffer.data(dk_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_y_yy_zzzzzzz = cbuffer.data(dk_geom_01_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_y_xxxxxxx, g_0_y_y_xxxxxxxy, g_0_y_y_xxxxxxy, g_0_y_y_xxxxxxyy, g_0_y_y_xxxxxxyz, g_0_y_y_xxxxxxz, g_0_y_y_xxxxxyy, g_0_y_y_xxxxxyyy, g_0_y_y_xxxxxyyz, g_0_y_y_xxxxxyz, g_0_y_y_xxxxxyzz, g_0_y_y_xxxxxzz, g_0_y_y_xxxxyyy, g_0_y_y_xxxxyyyy, g_0_y_y_xxxxyyyz, g_0_y_y_xxxxyyz, g_0_y_y_xxxxyyzz, g_0_y_y_xxxxyzz, g_0_y_y_xxxxyzzz, g_0_y_y_xxxxzzz, g_0_y_y_xxxyyyy, g_0_y_y_xxxyyyyy, g_0_y_y_xxxyyyyz, g_0_y_y_xxxyyyz, g_0_y_y_xxxyyyzz, g_0_y_y_xxxyyzz, g_0_y_y_xxxyyzzz, g_0_y_y_xxxyzzz, g_0_y_y_xxxyzzzz, g_0_y_y_xxxzzzz, g_0_y_y_xxyyyyy, g_0_y_y_xxyyyyyy, g_0_y_y_xxyyyyyz, g_0_y_y_xxyyyyz, g_0_y_y_xxyyyyzz, g_0_y_y_xxyyyzz, g_0_y_y_xxyyyzzz, g_0_y_y_xxyyzzz, g_0_y_y_xxyyzzzz, g_0_y_y_xxyzzzz, g_0_y_y_xxyzzzzz, g_0_y_y_xxzzzzz, g_0_y_y_xyyyyyy, g_0_y_y_xyyyyyyy, g_0_y_y_xyyyyyyz, g_0_y_y_xyyyyyz, g_0_y_y_xyyyyyzz, g_0_y_y_xyyyyzz, g_0_y_y_xyyyyzzz, g_0_y_y_xyyyzzz, g_0_y_y_xyyyzzzz, g_0_y_y_xyyzzzz, g_0_y_y_xyyzzzzz, g_0_y_y_xyzzzzz, g_0_y_y_xyzzzzzz, g_0_y_y_xzzzzzz, g_0_y_y_yyyyyyy, g_0_y_y_yyyyyyyy, g_0_y_y_yyyyyyyz, g_0_y_y_yyyyyyz, g_0_y_y_yyyyyyzz, g_0_y_y_yyyyyzz, g_0_y_y_yyyyyzzz, g_0_y_y_yyyyzzz, g_0_y_y_yyyyzzzz, g_0_y_y_yyyzzzz, g_0_y_y_yyyzzzzz, g_0_y_y_yyzzzzz, g_0_y_y_yyzzzzzz, g_0_y_y_yzzzzzz, g_0_y_y_yzzzzzzz, g_0_y_y_zzzzzzz, g_0_y_yy_xxxxxxx, g_0_y_yy_xxxxxxy, g_0_y_yy_xxxxxxz, g_0_y_yy_xxxxxyy, g_0_y_yy_xxxxxyz, g_0_y_yy_xxxxxzz, g_0_y_yy_xxxxyyy, g_0_y_yy_xxxxyyz, g_0_y_yy_xxxxyzz, g_0_y_yy_xxxxzzz, g_0_y_yy_xxxyyyy, g_0_y_yy_xxxyyyz, g_0_y_yy_xxxyyzz, g_0_y_yy_xxxyzzz, g_0_y_yy_xxxzzzz, g_0_y_yy_xxyyyyy, g_0_y_yy_xxyyyyz, g_0_y_yy_xxyyyzz, g_0_y_yy_xxyyzzz, g_0_y_yy_xxyzzzz, g_0_y_yy_xxzzzzz, g_0_y_yy_xyyyyyy, g_0_y_yy_xyyyyyz, g_0_y_yy_xyyyyzz, g_0_y_yy_xyyyzzz, g_0_y_yy_xyyzzzz, g_0_y_yy_xyzzzzz, g_0_y_yy_xzzzzzz, g_0_y_yy_yyyyyyy, g_0_y_yy_yyyyyyz, g_0_y_yy_yyyyyzz, g_0_y_yy_yyyyzzz, g_0_y_yy_yyyzzzz, g_0_y_yy_yyzzzzz, g_0_y_yy_yzzzzzz, g_0_y_yy_zzzzzzz, g_y_xxxxxxx, g_y_xxxxxxy, g_y_xxxxxxz, g_y_xxxxxyy, g_y_xxxxxyz, g_y_xxxxxzz, g_y_xxxxyyy, g_y_xxxxyyz, g_y_xxxxyzz, g_y_xxxxzzz, g_y_xxxyyyy, g_y_xxxyyyz, g_y_xxxyyzz, g_y_xxxyzzz, g_y_xxxzzzz, g_y_xxyyyyy, g_y_xxyyyyz, g_y_xxyyyzz, g_y_xxyyzzz, g_y_xxyzzzz, g_y_xxzzzzz, g_y_xyyyyyy, g_y_xyyyyyz, g_y_xyyyyzz, g_y_xyyyzzz, g_y_xyyzzzz, g_y_xyzzzzz, g_y_xzzzzzz, g_y_yyyyyyy, g_y_yyyyyyz, g_y_yyyyyzz, g_y_yyyyzzz, g_y_yyyzzzz, g_y_yyzzzzz, g_y_yzzzzzz, g_y_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yy_xxxxxxx[k] = g_y_xxxxxxx[k] - g_0_y_y_xxxxxxx[k] * ab_y + g_0_y_y_xxxxxxxy[k];

                g_0_y_yy_xxxxxxy[k] = g_y_xxxxxxy[k] - g_0_y_y_xxxxxxy[k] * ab_y + g_0_y_y_xxxxxxyy[k];

                g_0_y_yy_xxxxxxz[k] = g_y_xxxxxxz[k] - g_0_y_y_xxxxxxz[k] * ab_y + g_0_y_y_xxxxxxyz[k];

                g_0_y_yy_xxxxxyy[k] = g_y_xxxxxyy[k] - g_0_y_y_xxxxxyy[k] * ab_y + g_0_y_y_xxxxxyyy[k];

                g_0_y_yy_xxxxxyz[k] = g_y_xxxxxyz[k] - g_0_y_y_xxxxxyz[k] * ab_y + g_0_y_y_xxxxxyyz[k];

                g_0_y_yy_xxxxxzz[k] = g_y_xxxxxzz[k] - g_0_y_y_xxxxxzz[k] * ab_y + g_0_y_y_xxxxxyzz[k];

                g_0_y_yy_xxxxyyy[k] = g_y_xxxxyyy[k] - g_0_y_y_xxxxyyy[k] * ab_y + g_0_y_y_xxxxyyyy[k];

                g_0_y_yy_xxxxyyz[k] = g_y_xxxxyyz[k] - g_0_y_y_xxxxyyz[k] * ab_y + g_0_y_y_xxxxyyyz[k];

                g_0_y_yy_xxxxyzz[k] = g_y_xxxxyzz[k] - g_0_y_y_xxxxyzz[k] * ab_y + g_0_y_y_xxxxyyzz[k];

                g_0_y_yy_xxxxzzz[k] = g_y_xxxxzzz[k] - g_0_y_y_xxxxzzz[k] * ab_y + g_0_y_y_xxxxyzzz[k];

                g_0_y_yy_xxxyyyy[k] = g_y_xxxyyyy[k] - g_0_y_y_xxxyyyy[k] * ab_y + g_0_y_y_xxxyyyyy[k];

                g_0_y_yy_xxxyyyz[k] = g_y_xxxyyyz[k] - g_0_y_y_xxxyyyz[k] * ab_y + g_0_y_y_xxxyyyyz[k];

                g_0_y_yy_xxxyyzz[k] = g_y_xxxyyzz[k] - g_0_y_y_xxxyyzz[k] * ab_y + g_0_y_y_xxxyyyzz[k];

                g_0_y_yy_xxxyzzz[k] = g_y_xxxyzzz[k] - g_0_y_y_xxxyzzz[k] * ab_y + g_0_y_y_xxxyyzzz[k];

                g_0_y_yy_xxxzzzz[k] = g_y_xxxzzzz[k] - g_0_y_y_xxxzzzz[k] * ab_y + g_0_y_y_xxxyzzzz[k];

                g_0_y_yy_xxyyyyy[k] = g_y_xxyyyyy[k] - g_0_y_y_xxyyyyy[k] * ab_y + g_0_y_y_xxyyyyyy[k];

                g_0_y_yy_xxyyyyz[k] = g_y_xxyyyyz[k] - g_0_y_y_xxyyyyz[k] * ab_y + g_0_y_y_xxyyyyyz[k];

                g_0_y_yy_xxyyyzz[k] = g_y_xxyyyzz[k] - g_0_y_y_xxyyyzz[k] * ab_y + g_0_y_y_xxyyyyzz[k];

                g_0_y_yy_xxyyzzz[k] = g_y_xxyyzzz[k] - g_0_y_y_xxyyzzz[k] * ab_y + g_0_y_y_xxyyyzzz[k];

                g_0_y_yy_xxyzzzz[k] = g_y_xxyzzzz[k] - g_0_y_y_xxyzzzz[k] * ab_y + g_0_y_y_xxyyzzzz[k];

                g_0_y_yy_xxzzzzz[k] = g_y_xxzzzzz[k] - g_0_y_y_xxzzzzz[k] * ab_y + g_0_y_y_xxyzzzzz[k];

                g_0_y_yy_xyyyyyy[k] = g_y_xyyyyyy[k] - g_0_y_y_xyyyyyy[k] * ab_y + g_0_y_y_xyyyyyyy[k];

                g_0_y_yy_xyyyyyz[k] = g_y_xyyyyyz[k] - g_0_y_y_xyyyyyz[k] * ab_y + g_0_y_y_xyyyyyyz[k];

                g_0_y_yy_xyyyyzz[k] = g_y_xyyyyzz[k] - g_0_y_y_xyyyyzz[k] * ab_y + g_0_y_y_xyyyyyzz[k];

                g_0_y_yy_xyyyzzz[k] = g_y_xyyyzzz[k] - g_0_y_y_xyyyzzz[k] * ab_y + g_0_y_y_xyyyyzzz[k];

                g_0_y_yy_xyyzzzz[k] = g_y_xyyzzzz[k] - g_0_y_y_xyyzzzz[k] * ab_y + g_0_y_y_xyyyzzzz[k];

                g_0_y_yy_xyzzzzz[k] = g_y_xyzzzzz[k] - g_0_y_y_xyzzzzz[k] * ab_y + g_0_y_y_xyyzzzzz[k];

                g_0_y_yy_xzzzzzz[k] = g_y_xzzzzzz[k] - g_0_y_y_xzzzzzz[k] * ab_y + g_0_y_y_xyzzzzzz[k];

                g_0_y_yy_yyyyyyy[k] = g_y_yyyyyyy[k] - g_0_y_y_yyyyyyy[k] * ab_y + g_0_y_y_yyyyyyyy[k];

                g_0_y_yy_yyyyyyz[k] = g_y_yyyyyyz[k] - g_0_y_y_yyyyyyz[k] * ab_y + g_0_y_y_yyyyyyyz[k];

                g_0_y_yy_yyyyyzz[k] = g_y_yyyyyzz[k] - g_0_y_y_yyyyyzz[k] * ab_y + g_0_y_y_yyyyyyzz[k];

                g_0_y_yy_yyyyzzz[k] = g_y_yyyyzzz[k] - g_0_y_y_yyyyzzz[k] * ab_y + g_0_y_y_yyyyyzzz[k];

                g_0_y_yy_yyyzzzz[k] = g_y_yyyzzzz[k] - g_0_y_y_yyyzzzz[k] * ab_y + g_0_y_y_yyyyzzzz[k];

                g_0_y_yy_yyzzzzz[k] = g_y_yyzzzzz[k] - g_0_y_y_yyzzzzz[k] * ab_y + g_0_y_y_yyyzzzzz[k];

                g_0_y_yy_yzzzzzz[k] = g_y_yzzzzzz[k] - g_0_y_y_yzzzzzz[k] * ab_y + g_0_y_y_yyzzzzzz[k];

                g_0_y_yy_zzzzzzz[k] = g_y_zzzzzzz[k] - g_0_y_y_zzzzzzz[k] * ab_y + g_0_y_y_yzzzzzzz[k];
            }

            /// Set up 360-396 components of targeted buffer : cbuffer.data(

            auto g_0_y_yz_xxxxxxx = cbuffer.data(dk_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_y_yz_xxxxxxy = cbuffer.data(dk_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_y_yz_xxxxxxz = cbuffer.data(dk_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_y_yz_xxxxxyy = cbuffer.data(dk_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_y_yz_xxxxxyz = cbuffer.data(dk_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_y_yz_xxxxxzz = cbuffer.data(dk_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_y_yz_xxxxyyy = cbuffer.data(dk_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_y_yz_xxxxyyz = cbuffer.data(dk_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_y_yz_xxxxyzz = cbuffer.data(dk_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_y_yz_xxxxzzz = cbuffer.data(dk_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_y_yz_xxxyyyy = cbuffer.data(dk_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_y_yz_xxxyyyz = cbuffer.data(dk_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_y_yz_xxxyyzz = cbuffer.data(dk_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_y_yz_xxxyzzz = cbuffer.data(dk_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_y_yz_xxxzzzz = cbuffer.data(dk_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_y_yz_xxyyyyy = cbuffer.data(dk_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_y_yz_xxyyyyz = cbuffer.data(dk_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_y_yz_xxyyyzz = cbuffer.data(dk_geom_01_off + 377 * ccomps * dcomps);

            auto g_0_y_yz_xxyyzzz = cbuffer.data(dk_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_y_yz_xxyzzzz = cbuffer.data(dk_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_y_yz_xxzzzzz = cbuffer.data(dk_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_y_yz_xyyyyyy = cbuffer.data(dk_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_y_yz_xyyyyyz = cbuffer.data(dk_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_y_yz_xyyyyzz = cbuffer.data(dk_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_y_yz_xyyyzzz = cbuffer.data(dk_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_y_yz_xyyzzzz = cbuffer.data(dk_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_y_yz_xyzzzzz = cbuffer.data(dk_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_y_yz_xzzzzzz = cbuffer.data(dk_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_y_yz_yyyyyyy = cbuffer.data(dk_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_y_yz_yyyyyyz = cbuffer.data(dk_geom_01_off + 389 * ccomps * dcomps);

            auto g_0_y_yz_yyyyyzz = cbuffer.data(dk_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_y_yz_yyyyzzz = cbuffer.data(dk_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_y_yz_yyyzzzz = cbuffer.data(dk_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_y_yz_yyzzzzz = cbuffer.data(dk_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_y_yz_yzzzzzz = cbuffer.data(dk_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_y_yz_zzzzzzz = cbuffer.data(dk_geom_01_off + 395 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_y_xxxxxxx, g_0_y_y_xxxxxxxz, g_0_y_y_xxxxxxy, g_0_y_y_xxxxxxyz, g_0_y_y_xxxxxxz, g_0_y_y_xxxxxxzz, g_0_y_y_xxxxxyy, g_0_y_y_xxxxxyyz, g_0_y_y_xxxxxyz, g_0_y_y_xxxxxyzz, g_0_y_y_xxxxxzz, g_0_y_y_xxxxxzzz, g_0_y_y_xxxxyyy, g_0_y_y_xxxxyyyz, g_0_y_y_xxxxyyz, g_0_y_y_xxxxyyzz, g_0_y_y_xxxxyzz, g_0_y_y_xxxxyzzz, g_0_y_y_xxxxzzz, g_0_y_y_xxxxzzzz, g_0_y_y_xxxyyyy, g_0_y_y_xxxyyyyz, g_0_y_y_xxxyyyz, g_0_y_y_xxxyyyzz, g_0_y_y_xxxyyzz, g_0_y_y_xxxyyzzz, g_0_y_y_xxxyzzz, g_0_y_y_xxxyzzzz, g_0_y_y_xxxzzzz, g_0_y_y_xxxzzzzz, g_0_y_y_xxyyyyy, g_0_y_y_xxyyyyyz, g_0_y_y_xxyyyyz, g_0_y_y_xxyyyyzz, g_0_y_y_xxyyyzz, g_0_y_y_xxyyyzzz, g_0_y_y_xxyyzzz, g_0_y_y_xxyyzzzz, g_0_y_y_xxyzzzz, g_0_y_y_xxyzzzzz, g_0_y_y_xxzzzzz, g_0_y_y_xxzzzzzz, g_0_y_y_xyyyyyy, g_0_y_y_xyyyyyyz, g_0_y_y_xyyyyyz, g_0_y_y_xyyyyyzz, g_0_y_y_xyyyyzz, g_0_y_y_xyyyyzzz, g_0_y_y_xyyyzzz, g_0_y_y_xyyyzzzz, g_0_y_y_xyyzzzz, g_0_y_y_xyyzzzzz, g_0_y_y_xyzzzzz, g_0_y_y_xyzzzzzz, g_0_y_y_xzzzzzz, g_0_y_y_xzzzzzzz, g_0_y_y_yyyyyyy, g_0_y_y_yyyyyyyz, g_0_y_y_yyyyyyz, g_0_y_y_yyyyyyzz, g_0_y_y_yyyyyzz, g_0_y_y_yyyyyzzz, g_0_y_y_yyyyzzz, g_0_y_y_yyyyzzzz, g_0_y_y_yyyzzzz, g_0_y_y_yyyzzzzz, g_0_y_y_yyzzzzz, g_0_y_y_yyzzzzzz, g_0_y_y_yzzzzzz, g_0_y_y_yzzzzzzz, g_0_y_y_zzzzzzz, g_0_y_y_zzzzzzzz, g_0_y_yz_xxxxxxx, g_0_y_yz_xxxxxxy, g_0_y_yz_xxxxxxz, g_0_y_yz_xxxxxyy, g_0_y_yz_xxxxxyz, g_0_y_yz_xxxxxzz, g_0_y_yz_xxxxyyy, g_0_y_yz_xxxxyyz, g_0_y_yz_xxxxyzz, g_0_y_yz_xxxxzzz, g_0_y_yz_xxxyyyy, g_0_y_yz_xxxyyyz, g_0_y_yz_xxxyyzz, g_0_y_yz_xxxyzzz, g_0_y_yz_xxxzzzz, g_0_y_yz_xxyyyyy, g_0_y_yz_xxyyyyz, g_0_y_yz_xxyyyzz, g_0_y_yz_xxyyzzz, g_0_y_yz_xxyzzzz, g_0_y_yz_xxzzzzz, g_0_y_yz_xyyyyyy, g_0_y_yz_xyyyyyz, g_0_y_yz_xyyyyzz, g_0_y_yz_xyyyzzz, g_0_y_yz_xyyzzzz, g_0_y_yz_xyzzzzz, g_0_y_yz_xzzzzzz, g_0_y_yz_yyyyyyy, g_0_y_yz_yyyyyyz, g_0_y_yz_yyyyyzz, g_0_y_yz_yyyyzzz, g_0_y_yz_yyyzzzz, g_0_y_yz_yyzzzzz, g_0_y_yz_yzzzzzz, g_0_y_yz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yz_xxxxxxx[k] = -g_0_y_y_xxxxxxx[k] * ab_z + g_0_y_y_xxxxxxxz[k];

                g_0_y_yz_xxxxxxy[k] = -g_0_y_y_xxxxxxy[k] * ab_z + g_0_y_y_xxxxxxyz[k];

                g_0_y_yz_xxxxxxz[k] = -g_0_y_y_xxxxxxz[k] * ab_z + g_0_y_y_xxxxxxzz[k];

                g_0_y_yz_xxxxxyy[k] = -g_0_y_y_xxxxxyy[k] * ab_z + g_0_y_y_xxxxxyyz[k];

                g_0_y_yz_xxxxxyz[k] = -g_0_y_y_xxxxxyz[k] * ab_z + g_0_y_y_xxxxxyzz[k];

                g_0_y_yz_xxxxxzz[k] = -g_0_y_y_xxxxxzz[k] * ab_z + g_0_y_y_xxxxxzzz[k];

                g_0_y_yz_xxxxyyy[k] = -g_0_y_y_xxxxyyy[k] * ab_z + g_0_y_y_xxxxyyyz[k];

                g_0_y_yz_xxxxyyz[k] = -g_0_y_y_xxxxyyz[k] * ab_z + g_0_y_y_xxxxyyzz[k];

                g_0_y_yz_xxxxyzz[k] = -g_0_y_y_xxxxyzz[k] * ab_z + g_0_y_y_xxxxyzzz[k];

                g_0_y_yz_xxxxzzz[k] = -g_0_y_y_xxxxzzz[k] * ab_z + g_0_y_y_xxxxzzzz[k];

                g_0_y_yz_xxxyyyy[k] = -g_0_y_y_xxxyyyy[k] * ab_z + g_0_y_y_xxxyyyyz[k];

                g_0_y_yz_xxxyyyz[k] = -g_0_y_y_xxxyyyz[k] * ab_z + g_0_y_y_xxxyyyzz[k];

                g_0_y_yz_xxxyyzz[k] = -g_0_y_y_xxxyyzz[k] * ab_z + g_0_y_y_xxxyyzzz[k];

                g_0_y_yz_xxxyzzz[k] = -g_0_y_y_xxxyzzz[k] * ab_z + g_0_y_y_xxxyzzzz[k];

                g_0_y_yz_xxxzzzz[k] = -g_0_y_y_xxxzzzz[k] * ab_z + g_0_y_y_xxxzzzzz[k];

                g_0_y_yz_xxyyyyy[k] = -g_0_y_y_xxyyyyy[k] * ab_z + g_0_y_y_xxyyyyyz[k];

                g_0_y_yz_xxyyyyz[k] = -g_0_y_y_xxyyyyz[k] * ab_z + g_0_y_y_xxyyyyzz[k];

                g_0_y_yz_xxyyyzz[k] = -g_0_y_y_xxyyyzz[k] * ab_z + g_0_y_y_xxyyyzzz[k];

                g_0_y_yz_xxyyzzz[k] = -g_0_y_y_xxyyzzz[k] * ab_z + g_0_y_y_xxyyzzzz[k];

                g_0_y_yz_xxyzzzz[k] = -g_0_y_y_xxyzzzz[k] * ab_z + g_0_y_y_xxyzzzzz[k];

                g_0_y_yz_xxzzzzz[k] = -g_0_y_y_xxzzzzz[k] * ab_z + g_0_y_y_xxzzzzzz[k];

                g_0_y_yz_xyyyyyy[k] = -g_0_y_y_xyyyyyy[k] * ab_z + g_0_y_y_xyyyyyyz[k];

                g_0_y_yz_xyyyyyz[k] = -g_0_y_y_xyyyyyz[k] * ab_z + g_0_y_y_xyyyyyzz[k];

                g_0_y_yz_xyyyyzz[k] = -g_0_y_y_xyyyyzz[k] * ab_z + g_0_y_y_xyyyyzzz[k];

                g_0_y_yz_xyyyzzz[k] = -g_0_y_y_xyyyzzz[k] * ab_z + g_0_y_y_xyyyzzzz[k];

                g_0_y_yz_xyyzzzz[k] = -g_0_y_y_xyyzzzz[k] * ab_z + g_0_y_y_xyyzzzzz[k];

                g_0_y_yz_xyzzzzz[k] = -g_0_y_y_xyzzzzz[k] * ab_z + g_0_y_y_xyzzzzzz[k];

                g_0_y_yz_xzzzzzz[k] = -g_0_y_y_xzzzzzz[k] * ab_z + g_0_y_y_xzzzzzzz[k];

                g_0_y_yz_yyyyyyy[k] = -g_0_y_y_yyyyyyy[k] * ab_z + g_0_y_y_yyyyyyyz[k];

                g_0_y_yz_yyyyyyz[k] = -g_0_y_y_yyyyyyz[k] * ab_z + g_0_y_y_yyyyyyzz[k];

                g_0_y_yz_yyyyyzz[k] = -g_0_y_y_yyyyyzz[k] * ab_z + g_0_y_y_yyyyyzzz[k];

                g_0_y_yz_yyyyzzz[k] = -g_0_y_y_yyyyzzz[k] * ab_z + g_0_y_y_yyyyzzzz[k];

                g_0_y_yz_yyyzzzz[k] = -g_0_y_y_yyyzzzz[k] * ab_z + g_0_y_y_yyyzzzzz[k];

                g_0_y_yz_yyzzzzz[k] = -g_0_y_y_yyzzzzz[k] * ab_z + g_0_y_y_yyzzzzzz[k];

                g_0_y_yz_yzzzzzz[k] = -g_0_y_y_yzzzzzz[k] * ab_z + g_0_y_y_yzzzzzzz[k];

                g_0_y_yz_zzzzzzz[k] = -g_0_y_y_zzzzzzz[k] * ab_z + g_0_y_y_zzzzzzzz[k];
            }

            /// Set up 396-432 components of targeted buffer : cbuffer.data(

            auto g_0_y_zz_xxxxxxx = cbuffer.data(dk_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_y_zz_xxxxxxy = cbuffer.data(dk_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_y_zz_xxxxxxz = cbuffer.data(dk_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_y_zz_xxxxxyy = cbuffer.data(dk_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_y_zz_xxxxxyz = cbuffer.data(dk_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_y_zz_xxxxxzz = cbuffer.data(dk_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_y_zz_xxxxyyy = cbuffer.data(dk_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_y_zz_xxxxyyz = cbuffer.data(dk_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_y_zz_xxxxyzz = cbuffer.data(dk_geom_01_off + 404 * ccomps * dcomps);

            auto g_0_y_zz_xxxxzzz = cbuffer.data(dk_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_y_zz_xxxyyyy = cbuffer.data(dk_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_y_zz_xxxyyyz = cbuffer.data(dk_geom_01_off + 407 * ccomps * dcomps);

            auto g_0_y_zz_xxxyyzz = cbuffer.data(dk_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_y_zz_xxxyzzz = cbuffer.data(dk_geom_01_off + 409 * ccomps * dcomps);

            auto g_0_y_zz_xxxzzzz = cbuffer.data(dk_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_y_zz_xxyyyyy = cbuffer.data(dk_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_y_zz_xxyyyyz = cbuffer.data(dk_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_y_zz_xxyyyzz = cbuffer.data(dk_geom_01_off + 413 * ccomps * dcomps);

            auto g_0_y_zz_xxyyzzz = cbuffer.data(dk_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_y_zz_xxyzzzz = cbuffer.data(dk_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_y_zz_xxzzzzz = cbuffer.data(dk_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_y_zz_xyyyyyy = cbuffer.data(dk_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_y_zz_xyyyyyz = cbuffer.data(dk_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_y_zz_xyyyyzz = cbuffer.data(dk_geom_01_off + 419 * ccomps * dcomps);

            auto g_0_y_zz_xyyyzzz = cbuffer.data(dk_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_y_zz_xyyzzzz = cbuffer.data(dk_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_y_zz_xyzzzzz = cbuffer.data(dk_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_y_zz_xzzzzzz = cbuffer.data(dk_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_y_zz_yyyyyyy = cbuffer.data(dk_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_y_zz_yyyyyyz = cbuffer.data(dk_geom_01_off + 425 * ccomps * dcomps);

            auto g_0_y_zz_yyyyyzz = cbuffer.data(dk_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_y_zz_yyyyzzz = cbuffer.data(dk_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_y_zz_yyyzzzz = cbuffer.data(dk_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_y_zz_yyzzzzz = cbuffer.data(dk_geom_01_off + 429 * ccomps * dcomps);

            auto g_0_y_zz_yzzzzzz = cbuffer.data(dk_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_y_zz_zzzzzzz = cbuffer.data(dk_geom_01_off + 431 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_z_xxxxxxx, g_0_y_z_xxxxxxxz, g_0_y_z_xxxxxxy, g_0_y_z_xxxxxxyz, g_0_y_z_xxxxxxz, g_0_y_z_xxxxxxzz, g_0_y_z_xxxxxyy, g_0_y_z_xxxxxyyz, g_0_y_z_xxxxxyz, g_0_y_z_xxxxxyzz, g_0_y_z_xxxxxzz, g_0_y_z_xxxxxzzz, g_0_y_z_xxxxyyy, g_0_y_z_xxxxyyyz, g_0_y_z_xxxxyyz, g_0_y_z_xxxxyyzz, g_0_y_z_xxxxyzz, g_0_y_z_xxxxyzzz, g_0_y_z_xxxxzzz, g_0_y_z_xxxxzzzz, g_0_y_z_xxxyyyy, g_0_y_z_xxxyyyyz, g_0_y_z_xxxyyyz, g_0_y_z_xxxyyyzz, g_0_y_z_xxxyyzz, g_0_y_z_xxxyyzzz, g_0_y_z_xxxyzzz, g_0_y_z_xxxyzzzz, g_0_y_z_xxxzzzz, g_0_y_z_xxxzzzzz, g_0_y_z_xxyyyyy, g_0_y_z_xxyyyyyz, g_0_y_z_xxyyyyz, g_0_y_z_xxyyyyzz, g_0_y_z_xxyyyzz, g_0_y_z_xxyyyzzz, g_0_y_z_xxyyzzz, g_0_y_z_xxyyzzzz, g_0_y_z_xxyzzzz, g_0_y_z_xxyzzzzz, g_0_y_z_xxzzzzz, g_0_y_z_xxzzzzzz, g_0_y_z_xyyyyyy, g_0_y_z_xyyyyyyz, g_0_y_z_xyyyyyz, g_0_y_z_xyyyyyzz, g_0_y_z_xyyyyzz, g_0_y_z_xyyyyzzz, g_0_y_z_xyyyzzz, g_0_y_z_xyyyzzzz, g_0_y_z_xyyzzzz, g_0_y_z_xyyzzzzz, g_0_y_z_xyzzzzz, g_0_y_z_xyzzzzzz, g_0_y_z_xzzzzzz, g_0_y_z_xzzzzzzz, g_0_y_z_yyyyyyy, g_0_y_z_yyyyyyyz, g_0_y_z_yyyyyyz, g_0_y_z_yyyyyyzz, g_0_y_z_yyyyyzz, g_0_y_z_yyyyyzzz, g_0_y_z_yyyyzzz, g_0_y_z_yyyyzzzz, g_0_y_z_yyyzzzz, g_0_y_z_yyyzzzzz, g_0_y_z_yyzzzzz, g_0_y_z_yyzzzzzz, g_0_y_z_yzzzzzz, g_0_y_z_yzzzzzzz, g_0_y_z_zzzzzzz, g_0_y_z_zzzzzzzz, g_0_y_zz_xxxxxxx, g_0_y_zz_xxxxxxy, g_0_y_zz_xxxxxxz, g_0_y_zz_xxxxxyy, g_0_y_zz_xxxxxyz, g_0_y_zz_xxxxxzz, g_0_y_zz_xxxxyyy, g_0_y_zz_xxxxyyz, g_0_y_zz_xxxxyzz, g_0_y_zz_xxxxzzz, g_0_y_zz_xxxyyyy, g_0_y_zz_xxxyyyz, g_0_y_zz_xxxyyzz, g_0_y_zz_xxxyzzz, g_0_y_zz_xxxzzzz, g_0_y_zz_xxyyyyy, g_0_y_zz_xxyyyyz, g_0_y_zz_xxyyyzz, g_0_y_zz_xxyyzzz, g_0_y_zz_xxyzzzz, g_0_y_zz_xxzzzzz, g_0_y_zz_xyyyyyy, g_0_y_zz_xyyyyyz, g_0_y_zz_xyyyyzz, g_0_y_zz_xyyyzzz, g_0_y_zz_xyyzzzz, g_0_y_zz_xyzzzzz, g_0_y_zz_xzzzzzz, g_0_y_zz_yyyyyyy, g_0_y_zz_yyyyyyz, g_0_y_zz_yyyyyzz, g_0_y_zz_yyyyzzz, g_0_y_zz_yyyzzzz, g_0_y_zz_yyzzzzz, g_0_y_zz_yzzzzzz, g_0_y_zz_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_zz_xxxxxxx[k] = -g_0_y_z_xxxxxxx[k] * ab_z + g_0_y_z_xxxxxxxz[k];

                g_0_y_zz_xxxxxxy[k] = -g_0_y_z_xxxxxxy[k] * ab_z + g_0_y_z_xxxxxxyz[k];

                g_0_y_zz_xxxxxxz[k] = -g_0_y_z_xxxxxxz[k] * ab_z + g_0_y_z_xxxxxxzz[k];

                g_0_y_zz_xxxxxyy[k] = -g_0_y_z_xxxxxyy[k] * ab_z + g_0_y_z_xxxxxyyz[k];

                g_0_y_zz_xxxxxyz[k] = -g_0_y_z_xxxxxyz[k] * ab_z + g_0_y_z_xxxxxyzz[k];

                g_0_y_zz_xxxxxzz[k] = -g_0_y_z_xxxxxzz[k] * ab_z + g_0_y_z_xxxxxzzz[k];

                g_0_y_zz_xxxxyyy[k] = -g_0_y_z_xxxxyyy[k] * ab_z + g_0_y_z_xxxxyyyz[k];

                g_0_y_zz_xxxxyyz[k] = -g_0_y_z_xxxxyyz[k] * ab_z + g_0_y_z_xxxxyyzz[k];

                g_0_y_zz_xxxxyzz[k] = -g_0_y_z_xxxxyzz[k] * ab_z + g_0_y_z_xxxxyzzz[k];

                g_0_y_zz_xxxxzzz[k] = -g_0_y_z_xxxxzzz[k] * ab_z + g_0_y_z_xxxxzzzz[k];

                g_0_y_zz_xxxyyyy[k] = -g_0_y_z_xxxyyyy[k] * ab_z + g_0_y_z_xxxyyyyz[k];

                g_0_y_zz_xxxyyyz[k] = -g_0_y_z_xxxyyyz[k] * ab_z + g_0_y_z_xxxyyyzz[k];

                g_0_y_zz_xxxyyzz[k] = -g_0_y_z_xxxyyzz[k] * ab_z + g_0_y_z_xxxyyzzz[k];

                g_0_y_zz_xxxyzzz[k] = -g_0_y_z_xxxyzzz[k] * ab_z + g_0_y_z_xxxyzzzz[k];

                g_0_y_zz_xxxzzzz[k] = -g_0_y_z_xxxzzzz[k] * ab_z + g_0_y_z_xxxzzzzz[k];

                g_0_y_zz_xxyyyyy[k] = -g_0_y_z_xxyyyyy[k] * ab_z + g_0_y_z_xxyyyyyz[k];

                g_0_y_zz_xxyyyyz[k] = -g_0_y_z_xxyyyyz[k] * ab_z + g_0_y_z_xxyyyyzz[k];

                g_0_y_zz_xxyyyzz[k] = -g_0_y_z_xxyyyzz[k] * ab_z + g_0_y_z_xxyyyzzz[k];

                g_0_y_zz_xxyyzzz[k] = -g_0_y_z_xxyyzzz[k] * ab_z + g_0_y_z_xxyyzzzz[k];

                g_0_y_zz_xxyzzzz[k] = -g_0_y_z_xxyzzzz[k] * ab_z + g_0_y_z_xxyzzzzz[k];

                g_0_y_zz_xxzzzzz[k] = -g_0_y_z_xxzzzzz[k] * ab_z + g_0_y_z_xxzzzzzz[k];

                g_0_y_zz_xyyyyyy[k] = -g_0_y_z_xyyyyyy[k] * ab_z + g_0_y_z_xyyyyyyz[k];

                g_0_y_zz_xyyyyyz[k] = -g_0_y_z_xyyyyyz[k] * ab_z + g_0_y_z_xyyyyyzz[k];

                g_0_y_zz_xyyyyzz[k] = -g_0_y_z_xyyyyzz[k] * ab_z + g_0_y_z_xyyyyzzz[k];

                g_0_y_zz_xyyyzzz[k] = -g_0_y_z_xyyyzzz[k] * ab_z + g_0_y_z_xyyyzzzz[k];

                g_0_y_zz_xyyzzzz[k] = -g_0_y_z_xyyzzzz[k] * ab_z + g_0_y_z_xyyzzzzz[k];

                g_0_y_zz_xyzzzzz[k] = -g_0_y_z_xyzzzzz[k] * ab_z + g_0_y_z_xyzzzzzz[k];

                g_0_y_zz_xzzzzzz[k] = -g_0_y_z_xzzzzzz[k] * ab_z + g_0_y_z_xzzzzzzz[k];

                g_0_y_zz_yyyyyyy[k] = -g_0_y_z_yyyyyyy[k] * ab_z + g_0_y_z_yyyyyyyz[k];

                g_0_y_zz_yyyyyyz[k] = -g_0_y_z_yyyyyyz[k] * ab_z + g_0_y_z_yyyyyyzz[k];

                g_0_y_zz_yyyyyzz[k] = -g_0_y_z_yyyyyzz[k] * ab_z + g_0_y_z_yyyyyzzz[k];

                g_0_y_zz_yyyyzzz[k] = -g_0_y_z_yyyyzzz[k] * ab_z + g_0_y_z_yyyyzzzz[k];

                g_0_y_zz_yyyzzzz[k] = -g_0_y_z_yyyzzzz[k] * ab_z + g_0_y_z_yyyzzzzz[k];

                g_0_y_zz_yyzzzzz[k] = -g_0_y_z_yyzzzzz[k] * ab_z + g_0_y_z_yyzzzzzz[k];

                g_0_y_zz_yzzzzzz[k] = -g_0_y_z_yzzzzzz[k] * ab_z + g_0_y_z_yzzzzzzz[k];

                g_0_y_zz_zzzzzzz[k] = -g_0_y_z_zzzzzzz[k] * ab_z + g_0_y_z_zzzzzzzz[k];
            }

            /// Set up 432-468 components of targeted buffer : cbuffer.data(

            auto g_0_z_xx_xxxxxxx = cbuffer.data(dk_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_z_xx_xxxxxxy = cbuffer.data(dk_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_z_xx_xxxxxxz = cbuffer.data(dk_geom_01_off + 434 * ccomps * dcomps);

            auto g_0_z_xx_xxxxxyy = cbuffer.data(dk_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_z_xx_xxxxxyz = cbuffer.data(dk_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_z_xx_xxxxxzz = cbuffer.data(dk_geom_01_off + 437 * ccomps * dcomps);

            auto g_0_z_xx_xxxxyyy = cbuffer.data(dk_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_z_xx_xxxxyyz = cbuffer.data(dk_geom_01_off + 439 * ccomps * dcomps);

            auto g_0_z_xx_xxxxyzz = cbuffer.data(dk_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_z_xx_xxxxzzz = cbuffer.data(dk_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_z_xx_xxxyyyy = cbuffer.data(dk_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_z_xx_xxxyyyz = cbuffer.data(dk_geom_01_off + 443 * ccomps * dcomps);

            auto g_0_z_xx_xxxyyzz = cbuffer.data(dk_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_z_xx_xxxyzzz = cbuffer.data(dk_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_z_xx_xxxzzzz = cbuffer.data(dk_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_z_xx_xxyyyyy = cbuffer.data(dk_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_z_xx_xxyyyyz = cbuffer.data(dk_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_z_xx_xxyyyzz = cbuffer.data(dk_geom_01_off + 449 * ccomps * dcomps);

            auto g_0_z_xx_xxyyzzz = cbuffer.data(dk_geom_01_off + 450 * ccomps * dcomps);

            auto g_0_z_xx_xxyzzzz = cbuffer.data(dk_geom_01_off + 451 * ccomps * dcomps);

            auto g_0_z_xx_xxzzzzz = cbuffer.data(dk_geom_01_off + 452 * ccomps * dcomps);

            auto g_0_z_xx_xyyyyyy = cbuffer.data(dk_geom_01_off + 453 * ccomps * dcomps);

            auto g_0_z_xx_xyyyyyz = cbuffer.data(dk_geom_01_off + 454 * ccomps * dcomps);

            auto g_0_z_xx_xyyyyzz = cbuffer.data(dk_geom_01_off + 455 * ccomps * dcomps);

            auto g_0_z_xx_xyyyzzz = cbuffer.data(dk_geom_01_off + 456 * ccomps * dcomps);

            auto g_0_z_xx_xyyzzzz = cbuffer.data(dk_geom_01_off + 457 * ccomps * dcomps);

            auto g_0_z_xx_xyzzzzz = cbuffer.data(dk_geom_01_off + 458 * ccomps * dcomps);

            auto g_0_z_xx_xzzzzzz = cbuffer.data(dk_geom_01_off + 459 * ccomps * dcomps);

            auto g_0_z_xx_yyyyyyy = cbuffer.data(dk_geom_01_off + 460 * ccomps * dcomps);

            auto g_0_z_xx_yyyyyyz = cbuffer.data(dk_geom_01_off + 461 * ccomps * dcomps);

            auto g_0_z_xx_yyyyyzz = cbuffer.data(dk_geom_01_off + 462 * ccomps * dcomps);

            auto g_0_z_xx_yyyyzzz = cbuffer.data(dk_geom_01_off + 463 * ccomps * dcomps);

            auto g_0_z_xx_yyyzzzz = cbuffer.data(dk_geom_01_off + 464 * ccomps * dcomps);

            auto g_0_z_xx_yyzzzzz = cbuffer.data(dk_geom_01_off + 465 * ccomps * dcomps);

            auto g_0_z_xx_yzzzzzz = cbuffer.data(dk_geom_01_off + 466 * ccomps * dcomps);

            auto g_0_z_xx_zzzzzzz = cbuffer.data(dk_geom_01_off + 467 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_x_xxxxxxx, g_0_z_x_xxxxxxxx, g_0_z_x_xxxxxxxy, g_0_z_x_xxxxxxxz, g_0_z_x_xxxxxxy, g_0_z_x_xxxxxxyy, g_0_z_x_xxxxxxyz, g_0_z_x_xxxxxxz, g_0_z_x_xxxxxxzz, g_0_z_x_xxxxxyy, g_0_z_x_xxxxxyyy, g_0_z_x_xxxxxyyz, g_0_z_x_xxxxxyz, g_0_z_x_xxxxxyzz, g_0_z_x_xxxxxzz, g_0_z_x_xxxxxzzz, g_0_z_x_xxxxyyy, g_0_z_x_xxxxyyyy, g_0_z_x_xxxxyyyz, g_0_z_x_xxxxyyz, g_0_z_x_xxxxyyzz, g_0_z_x_xxxxyzz, g_0_z_x_xxxxyzzz, g_0_z_x_xxxxzzz, g_0_z_x_xxxxzzzz, g_0_z_x_xxxyyyy, g_0_z_x_xxxyyyyy, g_0_z_x_xxxyyyyz, g_0_z_x_xxxyyyz, g_0_z_x_xxxyyyzz, g_0_z_x_xxxyyzz, g_0_z_x_xxxyyzzz, g_0_z_x_xxxyzzz, g_0_z_x_xxxyzzzz, g_0_z_x_xxxzzzz, g_0_z_x_xxxzzzzz, g_0_z_x_xxyyyyy, g_0_z_x_xxyyyyyy, g_0_z_x_xxyyyyyz, g_0_z_x_xxyyyyz, g_0_z_x_xxyyyyzz, g_0_z_x_xxyyyzz, g_0_z_x_xxyyyzzz, g_0_z_x_xxyyzzz, g_0_z_x_xxyyzzzz, g_0_z_x_xxyzzzz, g_0_z_x_xxyzzzzz, g_0_z_x_xxzzzzz, g_0_z_x_xxzzzzzz, g_0_z_x_xyyyyyy, g_0_z_x_xyyyyyyy, g_0_z_x_xyyyyyyz, g_0_z_x_xyyyyyz, g_0_z_x_xyyyyyzz, g_0_z_x_xyyyyzz, g_0_z_x_xyyyyzzz, g_0_z_x_xyyyzzz, g_0_z_x_xyyyzzzz, g_0_z_x_xyyzzzz, g_0_z_x_xyyzzzzz, g_0_z_x_xyzzzzz, g_0_z_x_xyzzzzzz, g_0_z_x_xzzzzzz, g_0_z_x_xzzzzzzz, g_0_z_x_yyyyyyy, g_0_z_x_yyyyyyz, g_0_z_x_yyyyyzz, g_0_z_x_yyyyzzz, g_0_z_x_yyyzzzz, g_0_z_x_yyzzzzz, g_0_z_x_yzzzzzz, g_0_z_x_zzzzzzz, g_0_z_xx_xxxxxxx, g_0_z_xx_xxxxxxy, g_0_z_xx_xxxxxxz, g_0_z_xx_xxxxxyy, g_0_z_xx_xxxxxyz, g_0_z_xx_xxxxxzz, g_0_z_xx_xxxxyyy, g_0_z_xx_xxxxyyz, g_0_z_xx_xxxxyzz, g_0_z_xx_xxxxzzz, g_0_z_xx_xxxyyyy, g_0_z_xx_xxxyyyz, g_0_z_xx_xxxyyzz, g_0_z_xx_xxxyzzz, g_0_z_xx_xxxzzzz, g_0_z_xx_xxyyyyy, g_0_z_xx_xxyyyyz, g_0_z_xx_xxyyyzz, g_0_z_xx_xxyyzzz, g_0_z_xx_xxyzzzz, g_0_z_xx_xxzzzzz, g_0_z_xx_xyyyyyy, g_0_z_xx_xyyyyyz, g_0_z_xx_xyyyyzz, g_0_z_xx_xyyyzzz, g_0_z_xx_xyyzzzz, g_0_z_xx_xyzzzzz, g_0_z_xx_xzzzzzz, g_0_z_xx_yyyyyyy, g_0_z_xx_yyyyyyz, g_0_z_xx_yyyyyzz, g_0_z_xx_yyyyzzz, g_0_z_xx_yyyzzzz, g_0_z_xx_yyzzzzz, g_0_z_xx_yzzzzzz, g_0_z_xx_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xx_xxxxxxx[k] = -g_0_z_x_xxxxxxx[k] * ab_x + g_0_z_x_xxxxxxxx[k];

                g_0_z_xx_xxxxxxy[k] = -g_0_z_x_xxxxxxy[k] * ab_x + g_0_z_x_xxxxxxxy[k];

                g_0_z_xx_xxxxxxz[k] = -g_0_z_x_xxxxxxz[k] * ab_x + g_0_z_x_xxxxxxxz[k];

                g_0_z_xx_xxxxxyy[k] = -g_0_z_x_xxxxxyy[k] * ab_x + g_0_z_x_xxxxxxyy[k];

                g_0_z_xx_xxxxxyz[k] = -g_0_z_x_xxxxxyz[k] * ab_x + g_0_z_x_xxxxxxyz[k];

                g_0_z_xx_xxxxxzz[k] = -g_0_z_x_xxxxxzz[k] * ab_x + g_0_z_x_xxxxxxzz[k];

                g_0_z_xx_xxxxyyy[k] = -g_0_z_x_xxxxyyy[k] * ab_x + g_0_z_x_xxxxxyyy[k];

                g_0_z_xx_xxxxyyz[k] = -g_0_z_x_xxxxyyz[k] * ab_x + g_0_z_x_xxxxxyyz[k];

                g_0_z_xx_xxxxyzz[k] = -g_0_z_x_xxxxyzz[k] * ab_x + g_0_z_x_xxxxxyzz[k];

                g_0_z_xx_xxxxzzz[k] = -g_0_z_x_xxxxzzz[k] * ab_x + g_0_z_x_xxxxxzzz[k];

                g_0_z_xx_xxxyyyy[k] = -g_0_z_x_xxxyyyy[k] * ab_x + g_0_z_x_xxxxyyyy[k];

                g_0_z_xx_xxxyyyz[k] = -g_0_z_x_xxxyyyz[k] * ab_x + g_0_z_x_xxxxyyyz[k];

                g_0_z_xx_xxxyyzz[k] = -g_0_z_x_xxxyyzz[k] * ab_x + g_0_z_x_xxxxyyzz[k];

                g_0_z_xx_xxxyzzz[k] = -g_0_z_x_xxxyzzz[k] * ab_x + g_0_z_x_xxxxyzzz[k];

                g_0_z_xx_xxxzzzz[k] = -g_0_z_x_xxxzzzz[k] * ab_x + g_0_z_x_xxxxzzzz[k];

                g_0_z_xx_xxyyyyy[k] = -g_0_z_x_xxyyyyy[k] * ab_x + g_0_z_x_xxxyyyyy[k];

                g_0_z_xx_xxyyyyz[k] = -g_0_z_x_xxyyyyz[k] * ab_x + g_0_z_x_xxxyyyyz[k];

                g_0_z_xx_xxyyyzz[k] = -g_0_z_x_xxyyyzz[k] * ab_x + g_0_z_x_xxxyyyzz[k];

                g_0_z_xx_xxyyzzz[k] = -g_0_z_x_xxyyzzz[k] * ab_x + g_0_z_x_xxxyyzzz[k];

                g_0_z_xx_xxyzzzz[k] = -g_0_z_x_xxyzzzz[k] * ab_x + g_0_z_x_xxxyzzzz[k];

                g_0_z_xx_xxzzzzz[k] = -g_0_z_x_xxzzzzz[k] * ab_x + g_0_z_x_xxxzzzzz[k];

                g_0_z_xx_xyyyyyy[k] = -g_0_z_x_xyyyyyy[k] * ab_x + g_0_z_x_xxyyyyyy[k];

                g_0_z_xx_xyyyyyz[k] = -g_0_z_x_xyyyyyz[k] * ab_x + g_0_z_x_xxyyyyyz[k];

                g_0_z_xx_xyyyyzz[k] = -g_0_z_x_xyyyyzz[k] * ab_x + g_0_z_x_xxyyyyzz[k];

                g_0_z_xx_xyyyzzz[k] = -g_0_z_x_xyyyzzz[k] * ab_x + g_0_z_x_xxyyyzzz[k];

                g_0_z_xx_xyyzzzz[k] = -g_0_z_x_xyyzzzz[k] * ab_x + g_0_z_x_xxyyzzzz[k];

                g_0_z_xx_xyzzzzz[k] = -g_0_z_x_xyzzzzz[k] * ab_x + g_0_z_x_xxyzzzzz[k];

                g_0_z_xx_xzzzzzz[k] = -g_0_z_x_xzzzzzz[k] * ab_x + g_0_z_x_xxzzzzzz[k];

                g_0_z_xx_yyyyyyy[k] = -g_0_z_x_yyyyyyy[k] * ab_x + g_0_z_x_xyyyyyyy[k];

                g_0_z_xx_yyyyyyz[k] = -g_0_z_x_yyyyyyz[k] * ab_x + g_0_z_x_xyyyyyyz[k];

                g_0_z_xx_yyyyyzz[k] = -g_0_z_x_yyyyyzz[k] * ab_x + g_0_z_x_xyyyyyzz[k];

                g_0_z_xx_yyyyzzz[k] = -g_0_z_x_yyyyzzz[k] * ab_x + g_0_z_x_xyyyyzzz[k];

                g_0_z_xx_yyyzzzz[k] = -g_0_z_x_yyyzzzz[k] * ab_x + g_0_z_x_xyyyzzzz[k];

                g_0_z_xx_yyzzzzz[k] = -g_0_z_x_yyzzzzz[k] * ab_x + g_0_z_x_xyyzzzzz[k];

                g_0_z_xx_yzzzzzz[k] = -g_0_z_x_yzzzzzz[k] * ab_x + g_0_z_x_xyzzzzzz[k];

                g_0_z_xx_zzzzzzz[k] = -g_0_z_x_zzzzzzz[k] * ab_x + g_0_z_x_xzzzzzzz[k];
            }

            /// Set up 468-504 components of targeted buffer : cbuffer.data(

            auto g_0_z_xy_xxxxxxx = cbuffer.data(dk_geom_01_off + 468 * ccomps * dcomps);

            auto g_0_z_xy_xxxxxxy = cbuffer.data(dk_geom_01_off + 469 * ccomps * dcomps);

            auto g_0_z_xy_xxxxxxz = cbuffer.data(dk_geom_01_off + 470 * ccomps * dcomps);

            auto g_0_z_xy_xxxxxyy = cbuffer.data(dk_geom_01_off + 471 * ccomps * dcomps);

            auto g_0_z_xy_xxxxxyz = cbuffer.data(dk_geom_01_off + 472 * ccomps * dcomps);

            auto g_0_z_xy_xxxxxzz = cbuffer.data(dk_geom_01_off + 473 * ccomps * dcomps);

            auto g_0_z_xy_xxxxyyy = cbuffer.data(dk_geom_01_off + 474 * ccomps * dcomps);

            auto g_0_z_xy_xxxxyyz = cbuffer.data(dk_geom_01_off + 475 * ccomps * dcomps);

            auto g_0_z_xy_xxxxyzz = cbuffer.data(dk_geom_01_off + 476 * ccomps * dcomps);

            auto g_0_z_xy_xxxxzzz = cbuffer.data(dk_geom_01_off + 477 * ccomps * dcomps);

            auto g_0_z_xy_xxxyyyy = cbuffer.data(dk_geom_01_off + 478 * ccomps * dcomps);

            auto g_0_z_xy_xxxyyyz = cbuffer.data(dk_geom_01_off + 479 * ccomps * dcomps);

            auto g_0_z_xy_xxxyyzz = cbuffer.data(dk_geom_01_off + 480 * ccomps * dcomps);

            auto g_0_z_xy_xxxyzzz = cbuffer.data(dk_geom_01_off + 481 * ccomps * dcomps);

            auto g_0_z_xy_xxxzzzz = cbuffer.data(dk_geom_01_off + 482 * ccomps * dcomps);

            auto g_0_z_xy_xxyyyyy = cbuffer.data(dk_geom_01_off + 483 * ccomps * dcomps);

            auto g_0_z_xy_xxyyyyz = cbuffer.data(dk_geom_01_off + 484 * ccomps * dcomps);

            auto g_0_z_xy_xxyyyzz = cbuffer.data(dk_geom_01_off + 485 * ccomps * dcomps);

            auto g_0_z_xy_xxyyzzz = cbuffer.data(dk_geom_01_off + 486 * ccomps * dcomps);

            auto g_0_z_xy_xxyzzzz = cbuffer.data(dk_geom_01_off + 487 * ccomps * dcomps);

            auto g_0_z_xy_xxzzzzz = cbuffer.data(dk_geom_01_off + 488 * ccomps * dcomps);

            auto g_0_z_xy_xyyyyyy = cbuffer.data(dk_geom_01_off + 489 * ccomps * dcomps);

            auto g_0_z_xy_xyyyyyz = cbuffer.data(dk_geom_01_off + 490 * ccomps * dcomps);

            auto g_0_z_xy_xyyyyzz = cbuffer.data(dk_geom_01_off + 491 * ccomps * dcomps);

            auto g_0_z_xy_xyyyzzz = cbuffer.data(dk_geom_01_off + 492 * ccomps * dcomps);

            auto g_0_z_xy_xyyzzzz = cbuffer.data(dk_geom_01_off + 493 * ccomps * dcomps);

            auto g_0_z_xy_xyzzzzz = cbuffer.data(dk_geom_01_off + 494 * ccomps * dcomps);

            auto g_0_z_xy_xzzzzzz = cbuffer.data(dk_geom_01_off + 495 * ccomps * dcomps);

            auto g_0_z_xy_yyyyyyy = cbuffer.data(dk_geom_01_off + 496 * ccomps * dcomps);

            auto g_0_z_xy_yyyyyyz = cbuffer.data(dk_geom_01_off + 497 * ccomps * dcomps);

            auto g_0_z_xy_yyyyyzz = cbuffer.data(dk_geom_01_off + 498 * ccomps * dcomps);

            auto g_0_z_xy_yyyyzzz = cbuffer.data(dk_geom_01_off + 499 * ccomps * dcomps);

            auto g_0_z_xy_yyyzzzz = cbuffer.data(dk_geom_01_off + 500 * ccomps * dcomps);

            auto g_0_z_xy_yyzzzzz = cbuffer.data(dk_geom_01_off + 501 * ccomps * dcomps);

            auto g_0_z_xy_yzzzzzz = cbuffer.data(dk_geom_01_off + 502 * ccomps * dcomps);

            auto g_0_z_xy_zzzzzzz = cbuffer.data(dk_geom_01_off + 503 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xy_xxxxxxx, g_0_z_xy_xxxxxxy, g_0_z_xy_xxxxxxz, g_0_z_xy_xxxxxyy, g_0_z_xy_xxxxxyz, g_0_z_xy_xxxxxzz, g_0_z_xy_xxxxyyy, g_0_z_xy_xxxxyyz, g_0_z_xy_xxxxyzz, g_0_z_xy_xxxxzzz, g_0_z_xy_xxxyyyy, g_0_z_xy_xxxyyyz, g_0_z_xy_xxxyyzz, g_0_z_xy_xxxyzzz, g_0_z_xy_xxxzzzz, g_0_z_xy_xxyyyyy, g_0_z_xy_xxyyyyz, g_0_z_xy_xxyyyzz, g_0_z_xy_xxyyzzz, g_0_z_xy_xxyzzzz, g_0_z_xy_xxzzzzz, g_0_z_xy_xyyyyyy, g_0_z_xy_xyyyyyz, g_0_z_xy_xyyyyzz, g_0_z_xy_xyyyzzz, g_0_z_xy_xyyzzzz, g_0_z_xy_xyzzzzz, g_0_z_xy_xzzzzzz, g_0_z_xy_yyyyyyy, g_0_z_xy_yyyyyyz, g_0_z_xy_yyyyyzz, g_0_z_xy_yyyyzzz, g_0_z_xy_yyyzzzz, g_0_z_xy_yyzzzzz, g_0_z_xy_yzzzzzz, g_0_z_xy_zzzzzzz, g_0_z_y_xxxxxxx, g_0_z_y_xxxxxxxx, g_0_z_y_xxxxxxxy, g_0_z_y_xxxxxxxz, g_0_z_y_xxxxxxy, g_0_z_y_xxxxxxyy, g_0_z_y_xxxxxxyz, g_0_z_y_xxxxxxz, g_0_z_y_xxxxxxzz, g_0_z_y_xxxxxyy, g_0_z_y_xxxxxyyy, g_0_z_y_xxxxxyyz, g_0_z_y_xxxxxyz, g_0_z_y_xxxxxyzz, g_0_z_y_xxxxxzz, g_0_z_y_xxxxxzzz, g_0_z_y_xxxxyyy, g_0_z_y_xxxxyyyy, g_0_z_y_xxxxyyyz, g_0_z_y_xxxxyyz, g_0_z_y_xxxxyyzz, g_0_z_y_xxxxyzz, g_0_z_y_xxxxyzzz, g_0_z_y_xxxxzzz, g_0_z_y_xxxxzzzz, g_0_z_y_xxxyyyy, g_0_z_y_xxxyyyyy, g_0_z_y_xxxyyyyz, g_0_z_y_xxxyyyz, g_0_z_y_xxxyyyzz, g_0_z_y_xxxyyzz, g_0_z_y_xxxyyzzz, g_0_z_y_xxxyzzz, g_0_z_y_xxxyzzzz, g_0_z_y_xxxzzzz, g_0_z_y_xxxzzzzz, g_0_z_y_xxyyyyy, g_0_z_y_xxyyyyyy, g_0_z_y_xxyyyyyz, g_0_z_y_xxyyyyz, g_0_z_y_xxyyyyzz, g_0_z_y_xxyyyzz, g_0_z_y_xxyyyzzz, g_0_z_y_xxyyzzz, g_0_z_y_xxyyzzzz, g_0_z_y_xxyzzzz, g_0_z_y_xxyzzzzz, g_0_z_y_xxzzzzz, g_0_z_y_xxzzzzzz, g_0_z_y_xyyyyyy, g_0_z_y_xyyyyyyy, g_0_z_y_xyyyyyyz, g_0_z_y_xyyyyyz, g_0_z_y_xyyyyyzz, g_0_z_y_xyyyyzz, g_0_z_y_xyyyyzzz, g_0_z_y_xyyyzzz, g_0_z_y_xyyyzzzz, g_0_z_y_xyyzzzz, g_0_z_y_xyyzzzzz, g_0_z_y_xyzzzzz, g_0_z_y_xyzzzzzz, g_0_z_y_xzzzzzz, g_0_z_y_xzzzzzzz, g_0_z_y_yyyyyyy, g_0_z_y_yyyyyyz, g_0_z_y_yyyyyzz, g_0_z_y_yyyyzzz, g_0_z_y_yyyzzzz, g_0_z_y_yyzzzzz, g_0_z_y_yzzzzzz, g_0_z_y_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xy_xxxxxxx[k] = -g_0_z_y_xxxxxxx[k] * ab_x + g_0_z_y_xxxxxxxx[k];

                g_0_z_xy_xxxxxxy[k] = -g_0_z_y_xxxxxxy[k] * ab_x + g_0_z_y_xxxxxxxy[k];

                g_0_z_xy_xxxxxxz[k] = -g_0_z_y_xxxxxxz[k] * ab_x + g_0_z_y_xxxxxxxz[k];

                g_0_z_xy_xxxxxyy[k] = -g_0_z_y_xxxxxyy[k] * ab_x + g_0_z_y_xxxxxxyy[k];

                g_0_z_xy_xxxxxyz[k] = -g_0_z_y_xxxxxyz[k] * ab_x + g_0_z_y_xxxxxxyz[k];

                g_0_z_xy_xxxxxzz[k] = -g_0_z_y_xxxxxzz[k] * ab_x + g_0_z_y_xxxxxxzz[k];

                g_0_z_xy_xxxxyyy[k] = -g_0_z_y_xxxxyyy[k] * ab_x + g_0_z_y_xxxxxyyy[k];

                g_0_z_xy_xxxxyyz[k] = -g_0_z_y_xxxxyyz[k] * ab_x + g_0_z_y_xxxxxyyz[k];

                g_0_z_xy_xxxxyzz[k] = -g_0_z_y_xxxxyzz[k] * ab_x + g_0_z_y_xxxxxyzz[k];

                g_0_z_xy_xxxxzzz[k] = -g_0_z_y_xxxxzzz[k] * ab_x + g_0_z_y_xxxxxzzz[k];

                g_0_z_xy_xxxyyyy[k] = -g_0_z_y_xxxyyyy[k] * ab_x + g_0_z_y_xxxxyyyy[k];

                g_0_z_xy_xxxyyyz[k] = -g_0_z_y_xxxyyyz[k] * ab_x + g_0_z_y_xxxxyyyz[k];

                g_0_z_xy_xxxyyzz[k] = -g_0_z_y_xxxyyzz[k] * ab_x + g_0_z_y_xxxxyyzz[k];

                g_0_z_xy_xxxyzzz[k] = -g_0_z_y_xxxyzzz[k] * ab_x + g_0_z_y_xxxxyzzz[k];

                g_0_z_xy_xxxzzzz[k] = -g_0_z_y_xxxzzzz[k] * ab_x + g_0_z_y_xxxxzzzz[k];

                g_0_z_xy_xxyyyyy[k] = -g_0_z_y_xxyyyyy[k] * ab_x + g_0_z_y_xxxyyyyy[k];

                g_0_z_xy_xxyyyyz[k] = -g_0_z_y_xxyyyyz[k] * ab_x + g_0_z_y_xxxyyyyz[k];

                g_0_z_xy_xxyyyzz[k] = -g_0_z_y_xxyyyzz[k] * ab_x + g_0_z_y_xxxyyyzz[k];

                g_0_z_xy_xxyyzzz[k] = -g_0_z_y_xxyyzzz[k] * ab_x + g_0_z_y_xxxyyzzz[k];

                g_0_z_xy_xxyzzzz[k] = -g_0_z_y_xxyzzzz[k] * ab_x + g_0_z_y_xxxyzzzz[k];

                g_0_z_xy_xxzzzzz[k] = -g_0_z_y_xxzzzzz[k] * ab_x + g_0_z_y_xxxzzzzz[k];

                g_0_z_xy_xyyyyyy[k] = -g_0_z_y_xyyyyyy[k] * ab_x + g_0_z_y_xxyyyyyy[k];

                g_0_z_xy_xyyyyyz[k] = -g_0_z_y_xyyyyyz[k] * ab_x + g_0_z_y_xxyyyyyz[k];

                g_0_z_xy_xyyyyzz[k] = -g_0_z_y_xyyyyzz[k] * ab_x + g_0_z_y_xxyyyyzz[k];

                g_0_z_xy_xyyyzzz[k] = -g_0_z_y_xyyyzzz[k] * ab_x + g_0_z_y_xxyyyzzz[k];

                g_0_z_xy_xyyzzzz[k] = -g_0_z_y_xyyzzzz[k] * ab_x + g_0_z_y_xxyyzzzz[k];

                g_0_z_xy_xyzzzzz[k] = -g_0_z_y_xyzzzzz[k] * ab_x + g_0_z_y_xxyzzzzz[k];

                g_0_z_xy_xzzzzzz[k] = -g_0_z_y_xzzzzzz[k] * ab_x + g_0_z_y_xxzzzzzz[k];

                g_0_z_xy_yyyyyyy[k] = -g_0_z_y_yyyyyyy[k] * ab_x + g_0_z_y_xyyyyyyy[k];

                g_0_z_xy_yyyyyyz[k] = -g_0_z_y_yyyyyyz[k] * ab_x + g_0_z_y_xyyyyyyz[k];

                g_0_z_xy_yyyyyzz[k] = -g_0_z_y_yyyyyzz[k] * ab_x + g_0_z_y_xyyyyyzz[k];

                g_0_z_xy_yyyyzzz[k] = -g_0_z_y_yyyyzzz[k] * ab_x + g_0_z_y_xyyyyzzz[k];

                g_0_z_xy_yyyzzzz[k] = -g_0_z_y_yyyzzzz[k] * ab_x + g_0_z_y_xyyyzzzz[k];

                g_0_z_xy_yyzzzzz[k] = -g_0_z_y_yyzzzzz[k] * ab_x + g_0_z_y_xyyzzzzz[k];

                g_0_z_xy_yzzzzzz[k] = -g_0_z_y_yzzzzzz[k] * ab_x + g_0_z_y_xyzzzzzz[k];

                g_0_z_xy_zzzzzzz[k] = -g_0_z_y_zzzzzzz[k] * ab_x + g_0_z_y_xzzzzzzz[k];
            }

            /// Set up 504-540 components of targeted buffer : cbuffer.data(

            auto g_0_z_xz_xxxxxxx = cbuffer.data(dk_geom_01_off + 504 * ccomps * dcomps);

            auto g_0_z_xz_xxxxxxy = cbuffer.data(dk_geom_01_off + 505 * ccomps * dcomps);

            auto g_0_z_xz_xxxxxxz = cbuffer.data(dk_geom_01_off + 506 * ccomps * dcomps);

            auto g_0_z_xz_xxxxxyy = cbuffer.data(dk_geom_01_off + 507 * ccomps * dcomps);

            auto g_0_z_xz_xxxxxyz = cbuffer.data(dk_geom_01_off + 508 * ccomps * dcomps);

            auto g_0_z_xz_xxxxxzz = cbuffer.data(dk_geom_01_off + 509 * ccomps * dcomps);

            auto g_0_z_xz_xxxxyyy = cbuffer.data(dk_geom_01_off + 510 * ccomps * dcomps);

            auto g_0_z_xz_xxxxyyz = cbuffer.data(dk_geom_01_off + 511 * ccomps * dcomps);

            auto g_0_z_xz_xxxxyzz = cbuffer.data(dk_geom_01_off + 512 * ccomps * dcomps);

            auto g_0_z_xz_xxxxzzz = cbuffer.data(dk_geom_01_off + 513 * ccomps * dcomps);

            auto g_0_z_xz_xxxyyyy = cbuffer.data(dk_geom_01_off + 514 * ccomps * dcomps);

            auto g_0_z_xz_xxxyyyz = cbuffer.data(dk_geom_01_off + 515 * ccomps * dcomps);

            auto g_0_z_xz_xxxyyzz = cbuffer.data(dk_geom_01_off + 516 * ccomps * dcomps);

            auto g_0_z_xz_xxxyzzz = cbuffer.data(dk_geom_01_off + 517 * ccomps * dcomps);

            auto g_0_z_xz_xxxzzzz = cbuffer.data(dk_geom_01_off + 518 * ccomps * dcomps);

            auto g_0_z_xz_xxyyyyy = cbuffer.data(dk_geom_01_off + 519 * ccomps * dcomps);

            auto g_0_z_xz_xxyyyyz = cbuffer.data(dk_geom_01_off + 520 * ccomps * dcomps);

            auto g_0_z_xz_xxyyyzz = cbuffer.data(dk_geom_01_off + 521 * ccomps * dcomps);

            auto g_0_z_xz_xxyyzzz = cbuffer.data(dk_geom_01_off + 522 * ccomps * dcomps);

            auto g_0_z_xz_xxyzzzz = cbuffer.data(dk_geom_01_off + 523 * ccomps * dcomps);

            auto g_0_z_xz_xxzzzzz = cbuffer.data(dk_geom_01_off + 524 * ccomps * dcomps);

            auto g_0_z_xz_xyyyyyy = cbuffer.data(dk_geom_01_off + 525 * ccomps * dcomps);

            auto g_0_z_xz_xyyyyyz = cbuffer.data(dk_geom_01_off + 526 * ccomps * dcomps);

            auto g_0_z_xz_xyyyyzz = cbuffer.data(dk_geom_01_off + 527 * ccomps * dcomps);

            auto g_0_z_xz_xyyyzzz = cbuffer.data(dk_geom_01_off + 528 * ccomps * dcomps);

            auto g_0_z_xz_xyyzzzz = cbuffer.data(dk_geom_01_off + 529 * ccomps * dcomps);

            auto g_0_z_xz_xyzzzzz = cbuffer.data(dk_geom_01_off + 530 * ccomps * dcomps);

            auto g_0_z_xz_xzzzzzz = cbuffer.data(dk_geom_01_off + 531 * ccomps * dcomps);

            auto g_0_z_xz_yyyyyyy = cbuffer.data(dk_geom_01_off + 532 * ccomps * dcomps);

            auto g_0_z_xz_yyyyyyz = cbuffer.data(dk_geom_01_off + 533 * ccomps * dcomps);

            auto g_0_z_xz_yyyyyzz = cbuffer.data(dk_geom_01_off + 534 * ccomps * dcomps);

            auto g_0_z_xz_yyyyzzz = cbuffer.data(dk_geom_01_off + 535 * ccomps * dcomps);

            auto g_0_z_xz_yyyzzzz = cbuffer.data(dk_geom_01_off + 536 * ccomps * dcomps);

            auto g_0_z_xz_yyzzzzz = cbuffer.data(dk_geom_01_off + 537 * ccomps * dcomps);

            auto g_0_z_xz_yzzzzzz = cbuffer.data(dk_geom_01_off + 538 * ccomps * dcomps);

            auto g_0_z_xz_zzzzzzz = cbuffer.data(dk_geom_01_off + 539 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xz_xxxxxxx, g_0_z_xz_xxxxxxy, g_0_z_xz_xxxxxxz, g_0_z_xz_xxxxxyy, g_0_z_xz_xxxxxyz, g_0_z_xz_xxxxxzz, g_0_z_xz_xxxxyyy, g_0_z_xz_xxxxyyz, g_0_z_xz_xxxxyzz, g_0_z_xz_xxxxzzz, g_0_z_xz_xxxyyyy, g_0_z_xz_xxxyyyz, g_0_z_xz_xxxyyzz, g_0_z_xz_xxxyzzz, g_0_z_xz_xxxzzzz, g_0_z_xz_xxyyyyy, g_0_z_xz_xxyyyyz, g_0_z_xz_xxyyyzz, g_0_z_xz_xxyyzzz, g_0_z_xz_xxyzzzz, g_0_z_xz_xxzzzzz, g_0_z_xz_xyyyyyy, g_0_z_xz_xyyyyyz, g_0_z_xz_xyyyyzz, g_0_z_xz_xyyyzzz, g_0_z_xz_xyyzzzz, g_0_z_xz_xyzzzzz, g_0_z_xz_xzzzzzz, g_0_z_xz_yyyyyyy, g_0_z_xz_yyyyyyz, g_0_z_xz_yyyyyzz, g_0_z_xz_yyyyzzz, g_0_z_xz_yyyzzzz, g_0_z_xz_yyzzzzz, g_0_z_xz_yzzzzzz, g_0_z_xz_zzzzzzz, g_0_z_z_xxxxxxx, g_0_z_z_xxxxxxxx, g_0_z_z_xxxxxxxy, g_0_z_z_xxxxxxxz, g_0_z_z_xxxxxxy, g_0_z_z_xxxxxxyy, g_0_z_z_xxxxxxyz, g_0_z_z_xxxxxxz, g_0_z_z_xxxxxxzz, g_0_z_z_xxxxxyy, g_0_z_z_xxxxxyyy, g_0_z_z_xxxxxyyz, g_0_z_z_xxxxxyz, g_0_z_z_xxxxxyzz, g_0_z_z_xxxxxzz, g_0_z_z_xxxxxzzz, g_0_z_z_xxxxyyy, g_0_z_z_xxxxyyyy, g_0_z_z_xxxxyyyz, g_0_z_z_xxxxyyz, g_0_z_z_xxxxyyzz, g_0_z_z_xxxxyzz, g_0_z_z_xxxxyzzz, g_0_z_z_xxxxzzz, g_0_z_z_xxxxzzzz, g_0_z_z_xxxyyyy, g_0_z_z_xxxyyyyy, g_0_z_z_xxxyyyyz, g_0_z_z_xxxyyyz, g_0_z_z_xxxyyyzz, g_0_z_z_xxxyyzz, g_0_z_z_xxxyyzzz, g_0_z_z_xxxyzzz, g_0_z_z_xxxyzzzz, g_0_z_z_xxxzzzz, g_0_z_z_xxxzzzzz, g_0_z_z_xxyyyyy, g_0_z_z_xxyyyyyy, g_0_z_z_xxyyyyyz, g_0_z_z_xxyyyyz, g_0_z_z_xxyyyyzz, g_0_z_z_xxyyyzz, g_0_z_z_xxyyyzzz, g_0_z_z_xxyyzzz, g_0_z_z_xxyyzzzz, g_0_z_z_xxyzzzz, g_0_z_z_xxyzzzzz, g_0_z_z_xxzzzzz, g_0_z_z_xxzzzzzz, g_0_z_z_xyyyyyy, g_0_z_z_xyyyyyyy, g_0_z_z_xyyyyyyz, g_0_z_z_xyyyyyz, g_0_z_z_xyyyyyzz, g_0_z_z_xyyyyzz, g_0_z_z_xyyyyzzz, g_0_z_z_xyyyzzz, g_0_z_z_xyyyzzzz, g_0_z_z_xyyzzzz, g_0_z_z_xyyzzzzz, g_0_z_z_xyzzzzz, g_0_z_z_xyzzzzzz, g_0_z_z_xzzzzzz, g_0_z_z_xzzzzzzz, g_0_z_z_yyyyyyy, g_0_z_z_yyyyyyz, g_0_z_z_yyyyyzz, g_0_z_z_yyyyzzz, g_0_z_z_yyyzzzz, g_0_z_z_yyzzzzz, g_0_z_z_yzzzzzz, g_0_z_z_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xz_xxxxxxx[k] = -g_0_z_z_xxxxxxx[k] * ab_x + g_0_z_z_xxxxxxxx[k];

                g_0_z_xz_xxxxxxy[k] = -g_0_z_z_xxxxxxy[k] * ab_x + g_0_z_z_xxxxxxxy[k];

                g_0_z_xz_xxxxxxz[k] = -g_0_z_z_xxxxxxz[k] * ab_x + g_0_z_z_xxxxxxxz[k];

                g_0_z_xz_xxxxxyy[k] = -g_0_z_z_xxxxxyy[k] * ab_x + g_0_z_z_xxxxxxyy[k];

                g_0_z_xz_xxxxxyz[k] = -g_0_z_z_xxxxxyz[k] * ab_x + g_0_z_z_xxxxxxyz[k];

                g_0_z_xz_xxxxxzz[k] = -g_0_z_z_xxxxxzz[k] * ab_x + g_0_z_z_xxxxxxzz[k];

                g_0_z_xz_xxxxyyy[k] = -g_0_z_z_xxxxyyy[k] * ab_x + g_0_z_z_xxxxxyyy[k];

                g_0_z_xz_xxxxyyz[k] = -g_0_z_z_xxxxyyz[k] * ab_x + g_0_z_z_xxxxxyyz[k];

                g_0_z_xz_xxxxyzz[k] = -g_0_z_z_xxxxyzz[k] * ab_x + g_0_z_z_xxxxxyzz[k];

                g_0_z_xz_xxxxzzz[k] = -g_0_z_z_xxxxzzz[k] * ab_x + g_0_z_z_xxxxxzzz[k];

                g_0_z_xz_xxxyyyy[k] = -g_0_z_z_xxxyyyy[k] * ab_x + g_0_z_z_xxxxyyyy[k];

                g_0_z_xz_xxxyyyz[k] = -g_0_z_z_xxxyyyz[k] * ab_x + g_0_z_z_xxxxyyyz[k];

                g_0_z_xz_xxxyyzz[k] = -g_0_z_z_xxxyyzz[k] * ab_x + g_0_z_z_xxxxyyzz[k];

                g_0_z_xz_xxxyzzz[k] = -g_0_z_z_xxxyzzz[k] * ab_x + g_0_z_z_xxxxyzzz[k];

                g_0_z_xz_xxxzzzz[k] = -g_0_z_z_xxxzzzz[k] * ab_x + g_0_z_z_xxxxzzzz[k];

                g_0_z_xz_xxyyyyy[k] = -g_0_z_z_xxyyyyy[k] * ab_x + g_0_z_z_xxxyyyyy[k];

                g_0_z_xz_xxyyyyz[k] = -g_0_z_z_xxyyyyz[k] * ab_x + g_0_z_z_xxxyyyyz[k];

                g_0_z_xz_xxyyyzz[k] = -g_0_z_z_xxyyyzz[k] * ab_x + g_0_z_z_xxxyyyzz[k];

                g_0_z_xz_xxyyzzz[k] = -g_0_z_z_xxyyzzz[k] * ab_x + g_0_z_z_xxxyyzzz[k];

                g_0_z_xz_xxyzzzz[k] = -g_0_z_z_xxyzzzz[k] * ab_x + g_0_z_z_xxxyzzzz[k];

                g_0_z_xz_xxzzzzz[k] = -g_0_z_z_xxzzzzz[k] * ab_x + g_0_z_z_xxxzzzzz[k];

                g_0_z_xz_xyyyyyy[k] = -g_0_z_z_xyyyyyy[k] * ab_x + g_0_z_z_xxyyyyyy[k];

                g_0_z_xz_xyyyyyz[k] = -g_0_z_z_xyyyyyz[k] * ab_x + g_0_z_z_xxyyyyyz[k];

                g_0_z_xz_xyyyyzz[k] = -g_0_z_z_xyyyyzz[k] * ab_x + g_0_z_z_xxyyyyzz[k];

                g_0_z_xz_xyyyzzz[k] = -g_0_z_z_xyyyzzz[k] * ab_x + g_0_z_z_xxyyyzzz[k];

                g_0_z_xz_xyyzzzz[k] = -g_0_z_z_xyyzzzz[k] * ab_x + g_0_z_z_xxyyzzzz[k];

                g_0_z_xz_xyzzzzz[k] = -g_0_z_z_xyzzzzz[k] * ab_x + g_0_z_z_xxyzzzzz[k];

                g_0_z_xz_xzzzzzz[k] = -g_0_z_z_xzzzzzz[k] * ab_x + g_0_z_z_xxzzzzzz[k];

                g_0_z_xz_yyyyyyy[k] = -g_0_z_z_yyyyyyy[k] * ab_x + g_0_z_z_xyyyyyyy[k];

                g_0_z_xz_yyyyyyz[k] = -g_0_z_z_yyyyyyz[k] * ab_x + g_0_z_z_xyyyyyyz[k];

                g_0_z_xz_yyyyyzz[k] = -g_0_z_z_yyyyyzz[k] * ab_x + g_0_z_z_xyyyyyzz[k];

                g_0_z_xz_yyyyzzz[k] = -g_0_z_z_yyyyzzz[k] * ab_x + g_0_z_z_xyyyyzzz[k];

                g_0_z_xz_yyyzzzz[k] = -g_0_z_z_yyyzzzz[k] * ab_x + g_0_z_z_xyyyzzzz[k];

                g_0_z_xz_yyzzzzz[k] = -g_0_z_z_yyzzzzz[k] * ab_x + g_0_z_z_xyyzzzzz[k];

                g_0_z_xz_yzzzzzz[k] = -g_0_z_z_yzzzzzz[k] * ab_x + g_0_z_z_xyzzzzzz[k];

                g_0_z_xz_zzzzzzz[k] = -g_0_z_z_zzzzzzz[k] * ab_x + g_0_z_z_xzzzzzzz[k];
            }

            /// Set up 540-576 components of targeted buffer : cbuffer.data(

            auto g_0_z_yy_xxxxxxx = cbuffer.data(dk_geom_01_off + 540 * ccomps * dcomps);

            auto g_0_z_yy_xxxxxxy = cbuffer.data(dk_geom_01_off + 541 * ccomps * dcomps);

            auto g_0_z_yy_xxxxxxz = cbuffer.data(dk_geom_01_off + 542 * ccomps * dcomps);

            auto g_0_z_yy_xxxxxyy = cbuffer.data(dk_geom_01_off + 543 * ccomps * dcomps);

            auto g_0_z_yy_xxxxxyz = cbuffer.data(dk_geom_01_off + 544 * ccomps * dcomps);

            auto g_0_z_yy_xxxxxzz = cbuffer.data(dk_geom_01_off + 545 * ccomps * dcomps);

            auto g_0_z_yy_xxxxyyy = cbuffer.data(dk_geom_01_off + 546 * ccomps * dcomps);

            auto g_0_z_yy_xxxxyyz = cbuffer.data(dk_geom_01_off + 547 * ccomps * dcomps);

            auto g_0_z_yy_xxxxyzz = cbuffer.data(dk_geom_01_off + 548 * ccomps * dcomps);

            auto g_0_z_yy_xxxxzzz = cbuffer.data(dk_geom_01_off + 549 * ccomps * dcomps);

            auto g_0_z_yy_xxxyyyy = cbuffer.data(dk_geom_01_off + 550 * ccomps * dcomps);

            auto g_0_z_yy_xxxyyyz = cbuffer.data(dk_geom_01_off + 551 * ccomps * dcomps);

            auto g_0_z_yy_xxxyyzz = cbuffer.data(dk_geom_01_off + 552 * ccomps * dcomps);

            auto g_0_z_yy_xxxyzzz = cbuffer.data(dk_geom_01_off + 553 * ccomps * dcomps);

            auto g_0_z_yy_xxxzzzz = cbuffer.data(dk_geom_01_off + 554 * ccomps * dcomps);

            auto g_0_z_yy_xxyyyyy = cbuffer.data(dk_geom_01_off + 555 * ccomps * dcomps);

            auto g_0_z_yy_xxyyyyz = cbuffer.data(dk_geom_01_off + 556 * ccomps * dcomps);

            auto g_0_z_yy_xxyyyzz = cbuffer.data(dk_geom_01_off + 557 * ccomps * dcomps);

            auto g_0_z_yy_xxyyzzz = cbuffer.data(dk_geom_01_off + 558 * ccomps * dcomps);

            auto g_0_z_yy_xxyzzzz = cbuffer.data(dk_geom_01_off + 559 * ccomps * dcomps);

            auto g_0_z_yy_xxzzzzz = cbuffer.data(dk_geom_01_off + 560 * ccomps * dcomps);

            auto g_0_z_yy_xyyyyyy = cbuffer.data(dk_geom_01_off + 561 * ccomps * dcomps);

            auto g_0_z_yy_xyyyyyz = cbuffer.data(dk_geom_01_off + 562 * ccomps * dcomps);

            auto g_0_z_yy_xyyyyzz = cbuffer.data(dk_geom_01_off + 563 * ccomps * dcomps);

            auto g_0_z_yy_xyyyzzz = cbuffer.data(dk_geom_01_off + 564 * ccomps * dcomps);

            auto g_0_z_yy_xyyzzzz = cbuffer.data(dk_geom_01_off + 565 * ccomps * dcomps);

            auto g_0_z_yy_xyzzzzz = cbuffer.data(dk_geom_01_off + 566 * ccomps * dcomps);

            auto g_0_z_yy_xzzzzzz = cbuffer.data(dk_geom_01_off + 567 * ccomps * dcomps);

            auto g_0_z_yy_yyyyyyy = cbuffer.data(dk_geom_01_off + 568 * ccomps * dcomps);

            auto g_0_z_yy_yyyyyyz = cbuffer.data(dk_geom_01_off + 569 * ccomps * dcomps);

            auto g_0_z_yy_yyyyyzz = cbuffer.data(dk_geom_01_off + 570 * ccomps * dcomps);

            auto g_0_z_yy_yyyyzzz = cbuffer.data(dk_geom_01_off + 571 * ccomps * dcomps);

            auto g_0_z_yy_yyyzzzz = cbuffer.data(dk_geom_01_off + 572 * ccomps * dcomps);

            auto g_0_z_yy_yyzzzzz = cbuffer.data(dk_geom_01_off + 573 * ccomps * dcomps);

            auto g_0_z_yy_yzzzzzz = cbuffer.data(dk_geom_01_off + 574 * ccomps * dcomps);

            auto g_0_z_yy_zzzzzzz = cbuffer.data(dk_geom_01_off + 575 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_y_xxxxxxx, g_0_z_y_xxxxxxxy, g_0_z_y_xxxxxxy, g_0_z_y_xxxxxxyy, g_0_z_y_xxxxxxyz, g_0_z_y_xxxxxxz, g_0_z_y_xxxxxyy, g_0_z_y_xxxxxyyy, g_0_z_y_xxxxxyyz, g_0_z_y_xxxxxyz, g_0_z_y_xxxxxyzz, g_0_z_y_xxxxxzz, g_0_z_y_xxxxyyy, g_0_z_y_xxxxyyyy, g_0_z_y_xxxxyyyz, g_0_z_y_xxxxyyz, g_0_z_y_xxxxyyzz, g_0_z_y_xxxxyzz, g_0_z_y_xxxxyzzz, g_0_z_y_xxxxzzz, g_0_z_y_xxxyyyy, g_0_z_y_xxxyyyyy, g_0_z_y_xxxyyyyz, g_0_z_y_xxxyyyz, g_0_z_y_xxxyyyzz, g_0_z_y_xxxyyzz, g_0_z_y_xxxyyzzz, g_0_z_y_xxxyzzz, g_0_z_y_xxxyzzzz, g_0_z_y_xxxzzzz, g_0_z_y_xxyyyyy, g_0_z_y_xxyyyyyy, g_0_z_y_xxyyyyyz, g_0_z_y_xxyyyyz, g_0_z_y_xxyyyyzz, g_0_z_y_xxyyyzz, g_0_z_y_xxyyyzzz, g_0_z_y_xxyyzzz, g_0_z_y_xxyyzzzz, g_0_z_y_xxyzzzz, g_0_z_y_xxyzzzzz, g_0_z_y_xxzzzzz, g_0_z_y_xyyyyyy, g_0_z_y_xyyyyyyy, g_0_z_y_xyyyyyyz, g_0_z_y_xyyyyyz, g_0_z_y_xyyyyyzz, g_0_z_y_xyyyyzz, g_0_z_y_xyyyyzzz, g_0_z_y_xyyyzzz, g_0_z_y_xyyyzzzz, g_0_z_y_xyyzzzz, g_0_z_y_xyyzzzzz, g_0_z_y_xyzzzzz, g_0_z_y_xyzzzzzz, g_0_z_y_xzzzzzz, g_0_z_y_yyyyyyy, g_0_z_y_yyyyyyyy, g_0_z_y_yyyyyyyz, g_0_z_y_yyyyyyz, g_0_z_y_yyyyyyzz, g_0_z_y_yyyyyzz, g_0_z_y_yyyyyzzz, g_0_z_y_yyyyzzz, g_0_z_y_yyyyzzzz, g_0_z_y_yyyzzzz, g_0_z_y_yyyzzzzz, g_0_z_y_yyzzzzz, g_0_z_y_yyzzzzzz, g_0_z_y_yzzzzzz, g_0_z_y_yzzzzzzz, g_0_z_y_zzzzzzz, g_0_z_yy_xxxxxxx, g_0_z_yy_xxxxxxy, g_0_z_yy_xxxxxxz, g_0_z_yy_xxxxxyy, g_0_z_yy_xxxxxyz, g_0_z_yy_xxxxxzz, g_0_z_yy_xxxxyyy, g_0_z_yy_xxxxyyz, g_0_z_yy_xxxxyzz, g_0_z_yy_xxxxzzz, g_0_z_yy_xxxyyyy, g_0_z_yy_xxxyyyz, g_0_z_yy_xxxyyzz, g_0_z_yy_xxxyzzz, g_0_z_yy_xxxzzzz, g_0_z_yy_xxyyyyy, g_0_z_yy_xxyyyyz, g_0_z_yy_xxyyyzz, g_0_z_yy_xxyyzzz, g_0_z_yy_xxyzzzz, g_0_z_yy_xxzzzzz, g_0_z_yy_xyyyyyy, g_0_z_yy_xyyyyyz, g_0_z_yy_xyyyyzz, g_0_z_yy_xyyyzzz, g_0_z_yy_xyyzzzz, g_0_z_yy_xyzzzzz, g_0_z_yy_xzzzzzz, g_0_z_yy_yyyyyyy, g_0_z_yy_yyyyyyz, g_0_z_yy_yyyyyzz, g_0_z_yy_yyyyzzz, g_0_z_yy_yyyzzzz, g_0_z_yy_yyzzzzz, g_0_z_yy_yzzzzzz, g_0_z_yy_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yy_xxxxxxx[k] = -g_0_z_y_xxxxxxx[k] * ab_y + g_0_z_y_xxxxxxxy[k];

                g_0_z_yy_xxxxxxy[k] = -g_0_z_y_xxxxxxy[k] * ab_y + g_0_z_y_xxxxxxyy[k];

                g_0_z_yy_xxxxxxz[k] = -g_0_z_y_xxxxxxz[k] * ab_y + g_0_z_y_xxxxxxyz[k];

                g_0_z_yy_xxxxxyy[k] = -g_0_z_y_xxxxxyy[k] * ab_y + g_0_z_y_xxxxxyyy[k];

                g_0_z_yy_xxxxxyz[k] = -g_0_z_y_xxxxxyz[k] * ab_y + g_0_z_y_xxxxxyyz[k];

                g_0_z_yy_xxxxxzz[k] = -g_0_z_y_xxxxxzz[k] * ab_y + g_0_z_y_xxxxxyzz[k];

                g_0_z_yy_xxxxyyy[k] = -g_0_z_y_xxxxyyy[k] * ab_y + g_0_z_y_xxxxyyyy[k];

                g_0_z_yy_xxxxyyz[k] = -g_0_z_y_xxxxyyz[k] * ab_y + g_0_z_y_xxxxyyyz[k];

                g_0_z_yy_xxxxyzz[k] = -g_0_z_y_xxxxyzz[k] * ab_y + g_0_z_y_xxxxyyzz[k];

                g_0_z_yy_xxxxzzz[k] = -g_0_z_y_xxxxzzz[k] * ab_y + g_0_z_y_xxxxyzzz[k];

                g_0_z_yy_xxxyyyy[k] = -g_0_z_y_xxxyyyy[k] * ab_y + g_0_z_y_xxxyyyyy[k];

                g_0_z_yy_xxxyyyz[k] = -g_0_z_y_xxxyyyz[k] * ab_y + g_0_z_y_xxxyyyyz[k];

                g_0_z_yy_xxxyyzz[k] = -g_0_z_y_xxxyyzz[k] * ab_y + g_0_z_y_xxxyyyzz[k];

                g_0_z_yy_xxxyzzz[k] = -g_0_z_y_xxxyzzz[k] * ab_y + g_0_z_y_xxxyyzzz[k];

                g_0_z_yy_xxxzzzz[k] = -g_0_z_y_xxxzzzz[k] * ab_y + g_0_z_y_xxxyzzzz[k];

                g_0_z_yy_xxyyyyy[k] = -g_0_z_y_xxyyyyy[k] * ab_y + g_0_z_y_xxyyyyyy[k];

                g_0_z_yy_xxyyyyz[k] = -g_0_z_y_xxyyyyz[k] * ab_y + g_0_z_y_xxyyyyyz[k];

                g_0_z_yy_xxyyyzz[k] = -g_0_z_y_xxyyyzz[k] * ab_y + g_0_z_y_xxyyyyzz[k];

                g_0_z_yy_xxyyzzz[k] = -g_0_z_y_xxyyzzz[k] * ab_y + g_0_z_y_xxyyyzzz[k];

                g_0_z_yy_xxyzzzz[k] = -g_0_z_y_xxyzzzz[k] * ab_y + g_0_z_y_xxyyzzzz[k];

                g_0_z_yy_xxzzzzz[k] = -g_0_z_y_xxzzzzz[k] * ab_y + g_0_z_y_xxyzzzzz[k];

                g_0_z_yy_xyyyyyy[k] = -g_0_z_y_xyyyyyy[k] * ab_y + g_0_z_y_xyyyyyyy[k];

                g_0_z_yy_xyyyyyz[k] = -g_0_z_y_xyyyyyz[k] * ab_y + g_0_z_y_xyyyyyyz[k];

                g_0_z_yy_xyyyyzz[k] = -g_0_z_y_xyyyyzz[k] * ab_y + g_0_z_y_xyyyyyzz[k];

                g_0_z_yy_xyyyzzz[k] = -g_0_z_y_xyyyzzz[k] * ab_y + g_0_z_y_xyyyyzzz[k];

                g_0_z_yy_xyyzzzz[k] = -g_0_z_y_xyyzzzz[k] * ab_y + g_0_z_y_xyyyzzzz[k];

                g_0_z_yy_xyzzzzz[k] = -g_0_z_y_xyzzzzz[k] * ab_y + g_0_z_y_xyyzzzzz[k];

                g_0_z_yy_xzzzzzz[k] = -g_0_z_y_xzzzzzz[k] * ab_y + g_0_z_y_xyzzzzzz[k];

                g_0_z_yy_yyyyyyy[k] = -g_0_z_y_yyyyyyy[k] * ab_y + g_0_z_y_yyyyyyyy[k];

                g_0_z_yy_yyyyyyz[k] = -g_0_z_y_yyyyyyz[k] * ab_y + g_0_z_y_yyyyyyyz[k];

                g_0_z_yy_yyyyyzz[k] = -g_0_z_y_yyyyyzz[k] * ab_y + g_0_z_y_yyyyyyzz[k];

                g_0_z_yy_yyyyzzz[k] = -g_0_z_y_yyyyzzz[k] * ab_y + g_0_z_y_yyyyyzzz[k];

                g_0_z_yy_yyyzzzz[k] = -g_0_z_y_yyyzzzz[k] * ab_y + g_0_z_y_yyyyzzzz[k];

                g_0_z_yy_yyzzzzz[k] = -g_0_z_y_yyzzzzz[k] * ab_y + g_0_z_y_yyyzzzzz[k];

                g_0_z_yy_yzzzzzz[k] = -g_0_z_y_yzzzzzz[k] * ab_y + g_0_z_y_yyzzzzzz[k];

                g_0_z_yy_zzzzzzz[k] = -g_0_z_y_zzzzzzz[k] * ab_y + g_0_z_y_yzzzzzzz[k];
            }

            /// Set up 576-612 components of targeted buffer : cbuffer.data(

            auto g_0_z_yz_xxxxxxx = cbuffer.data(dk_geom_01_off + 576 * ccomps * dcomps);

            auto g_0_z_yz_xxxxxxy = cbuffer.data(dk_geom_01_off + 577 * ccomps * dcomps);

            auto g_0_z_yz_xxxxxxz = cbuffer.data(dk_geom_01_off + 578 * ccomps * dcomps);

            auto g_0_z_yz_xxxxxyy = cbuffer.data(dk_geom_01_off + 579 * ccomps * dcomps);

            auto g_0_z_yz_xxxxxyz = cbuffer.data(dk_geom_01_off + 580 * ccomps * dcomps);

            auto g_0_z_yz_xxxxxzz = cbuffer.data(dk_geom_01_off + 581 * ccomps * dcomps);

            auto g_0_z_yz_xxxxyyy = cbuffer.data(dk_geom_01_off + 582 * ccomps * dcomps);

            auto g_0_z_yz_xxxxyyz = cbuffer.data(dk_geom_01_off + 583 * ccomps * dcomps);

            auto g_0_z_yz_xxxxyzz = cbuffer.data(dk_geom_01_off + 584 * ccomps * dcomps);

            auto g_0_z_yz_xxxxzzz = cbuffer.data(dk_geom_01_off + 585 * ccomps * dcomps);

            auto g_0_z_yz_xxxyyyy = cbuffer.data(dk_geom_01_off + 586 * ccomps * dcomps);

            auto g_0_z_yz_xxxyyyz = cbuffer.data(dk_geom_01_off + 587 * ccomps * dcomps);

            auto g_0_z_yz_xxxyyzz = cbuffer.data(dk_geom_01_off + 588 * ccomps * dcomps);

            auto g_0_z_yz_xxxyzzz = cbuffer.data(dk_geom_01_off + 589 * ccomps * dcomps);

            auto g_0_z_yz_xxxzzzz = cbuffer.data(dk_geom_01_off + 590 * ccomps * dcomps);

            auto g_0_z_yz_xxyyyyy = cbuffer.data(dk_geom_01_off + 591 * ccomps * dcomps);

            auto g_0_z_yz_xxyyyyz = cbuffer.data(dk_geom_01_off + 592 * ccomps * dcomps);

            auto g_0_z_yz_xxyyyzz = cbuffer.data(dk_geom_01_off + 593 * ccomps * dcomps);

            auto g_0_z_yz_xxyyzzz = cbuffer.data(dk_geom_01_off + 594 * ccomps * dcomps);

            auto g_0_z_yz_xxyzzzz = cbuffer.data(dk_geom_01_off + 595 * ccomps * dcomps);

            auto g_0_z_yz_xxzzzzz = cbuffer.data(dk_geom_01_off + 596 * ccomps * dcomps);

            auto g_0_z_yz_xyyyyyy = cbuffer.data(dk_geom_01_off + 597 * ccomps * dcomps);

            auto g_0_z_yz_xyyyyyz = cbuffer.data(dk_geom_01_off + 598 * ccomps * dcomps);

            auto g_0_z_yz_xyyyyzz = cbuffer.data(dk_geom_01_off + 599 * ccomps * dcomps);

            auto g_0_z_yz_xyyyzzz = cbuffer.data(dk_geom_01_off + 600 * ccomps * dcomps);

            auto g_0_z_yz_xyyzzzz = cbuffer.data(dk_geom_01_off + 601 * ccomps * dcomps);

            auto g_0_z_yz_xyzzzzz = cbuffer.data(dk_geom_01_off + 602 * ccomps * dcomps);

            auto g_0_z_yz_xzzzzzz = cbuffer.data(dk_geom_01_off + 603 * ccomps * dcomps);

            auto g_0_z_yz_yyyyyyy = cbuffer.data(dk_geom_01_off + 604 * ccomps * dcomps);

            auto g_0_z_yz_yyyyyyz = cbuffer.data(dk_geom_01_off + 605 * ccomps * dcomps);

            auto g_0_z_yz_yyyyyzz = cbuffer.data(dk_geom_01_off + 606 * ccomps * dcomps);

            auto g_0_z_yz_yyyyzzz = cbuffer.data(dk_geom_01_off + 607 * ccomps * dcomps);

            auto g_0_z_yz_yyyzzzz = cbuffer.data(dk_geom_01_off + 608 * ccomps * dcomps);

            auto g_0_z_yz_yyzzzzz = cbuffer.data(dk_geom_01_off + 609 * ccomps * dcomps);

            auto g_0_z_yz_yzzzzzz = cbuffer.data(dk_geom_01_off + 610 * ccomps * dcomps);

            auto g_0_z_yz_zzzzzzz = cbuffer.data(dk_geom_01_off + 611 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yz_xxxxxxx, g_0_z_yz_xxxxxxy, g_0_z_yz_xxxxxxz, g_0_z_yz_xxxxxyy, g_0_z_yz_xxxxxyz, g_0_z_yz_xxxxxzz, g_0_z_yz_xxxxyyy, g_0_z_yz_xxxxyyz, g_0_z_yz_xxxxyzz, g_0_z_yz_xxxxzzz, g_0_z_yz_xxxyyyy, g_0_z_yz_xxxyyyz, g_0_z_yz_xxxyyzz, g_0_z_yz_xxxyzzz, g_0_z_yz_xxxzzzz, g_0_z_yz_xxyyyyy, g_0_z_yz_xxyyyyz, g_0_z_yz_xxyyyzz, g_0_z_yz_xxyyzzz, g_0_z_yz_xxyzzzz, g_0_z_yz_xxzzzzz, g_0_z_yz_xyyyyyy, g_0_z_yz_xyyyyyz, g_0_z_yz_xyyyyzz, g_0_z_yz_xyyyzzz, g_0_z_yz_xyyzzzz, g_0_z_yz_xyzzzzz, g_0_z_yz_xzzzzzz, g_0_z_yz_yyyyyyy, g_0_z_yz_yyyyyyz, g_0_z_yz_yyyyyzz, g_0_z_yz_yyyyzzz, g_0_z_yz_yyyzzzz, g_0_z_yz_yyzzzzz, g_0_z_yz_yzzzzzz, g_0_z_yz_zzzzzzz, g_0_z_z_xxxxxxx, g_0_z_z_xxxxxxxy, g_0_z_z_xxxxxxy, g_0_z_z_xxxxxxyy, g_0_z_z_xxxxxxyz, g_0_z_z_xxxxxxz, g_0_z_z_xxxxxyy, g_0_z_z_xxxxxyyy, g_0_z_z_xxxxxyyz, g_0_z_z_xxxxxyz, g_0_z_z_xxxxxyzz, g_0_z_z_xxxxxzz, g_0_z_z_xxxxyyy, g_0_z_z_xxxxyyyy, g_0_z_z_xxxxyyyz, g_0_z_z_xxxxyyz, g_0_z_z_xxxxyyzz, g_0_z_z_xxxxyzz, g_0_z_z_xxxxyzzz, g_0_z_z_xxxxzzz, g_0_z_z_xxxyyyy, g_0_z_z_xxxyyyyy, g_0_z_z_xxxyyyyz, g_0_z_z_xxxyyyz, g_0_z_z_xxxyyyzz, g_0_z_z_xxxyyzz, g_0_z_z_xxxyyzzz, g_0_z_z_xxxyzzz, g_0_z_z_xxxyzzzz, g_0_z_z_xxxzzzz, g_0_z_z_xxyyyyy, g_0_z_z_xxyyyyyy, g_0_z_z_xxyyyyyz, g_0_z_z_xxyyyyz, g_0_z_z_xxyyyyzz, g_0_z_z_xxyyyzz, g_0_z_z_xxyyyzzz, g_0_z_z_xxyyzzz, g_0_z_z_xxyyzzzz, g_0_z_z_xxyzzzz, g_0_z_z_xxyzzzzz, g_0_z_z_xxzzzzz, g_0_z_z_xyyyyyy, g_0_z_z_xyyyyyyy, g_0_z_z_xyyyyyyz, g_0_z_z_xyyyyyz, g_0_z_z_xyyyyyzz, g_0_z_z_xyyyyzz, g_0_z_z_xyyyyzzz, g_0_z_z_xyyyzzz, g_0_z_z_xyyyzzzz, g_0_z_z_xyyzzzz, g_0_z_z_xyyzzzzz, g_0_z_z_xyzzzzz, g_0_z_z_xyzzzzzz, g_0_z_z_xzzzzzz, g_0_z_z_yyyyyyy, g_0_z_z_yyyyyyyy, g_0_z_z_yyyyyyyz, g_0_z_z_yyyyyyz, g_0_z_z_yyyyyyzz, g_0_z_z_yyyyyzz, g_0_z_z_yyyyyzzz, g_0_z_z_yyyyzzz, g_0_z_z_yyyyzzzz, g_0_z_z_yyyzzzz, g_0_z_z_yyyzzzzz, g_0_z_z_yyzzzzz, g_0_z_z_yyzzzzzz, g_0_z_z_yzzzzzz, g_0_z_z_yzzzzzzz, g_0_z_z_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yz_xxxxxxx[k] = -g_0_z_z_xxxxxxx[k] * ab_y + g_0_z_z_xxxxxxxy[k];

                g_0_z_yz_xxxxxxy[k] = -g_0_z_z_xxxxxxy[k] * ab_y + g_0_z_z_xxxxxxyy[k];

                g_0_z_yz_xxxxxxz[k] = -g_0_z_z_xxxxxxz[k] * ab_y + g_0_z_z_xxxxxxyz[k];

                g_0_z_yz_xxxxxyy[k] = -g_0_z_z_xxxxxyy[k] * ab_y + g_0_z_z_xxxxxyyy[k];

                g_0_z_yz_xxxxxyz[k] = -g_0_z_z_xxxxxyz[k] * ab_y + g_0_z_z_xxxxxyyz[k];

                g_0_z_yz_xxxxxzz[k] = -g_0_z_z_xxxxxzz[k] * ab_y + g_0_z_z_xxxxxyzz[k];

                g_0_z_yz_xxxxyyy[k] = -g_0_z_z_xxxxyyy[k] * ab_y + g_0_z_z_xxxxyyyy[k];

                g_0_z_yz_xxxxyyz[k] = -g_0_z_z_xxxxyyz[k] * ab_y + g_0_z_z_xxxxyyyz[k];

                g_0_z_yz_xxxxyzz[k] = -g_0_z_z_xxxxyzz[k] * ab_y + g_0_z_z_xxxxyyzz[k];

                g_0_z_yz_xxxxzzz[k] = -g_0_z_z_xxxxzzz[k] * ab_y + g_0_z_z_xxxxyzzz[k];

                g_0_z_yz_xxxyyyy[k] = -g_0_z_z_xxxyyyy[k] * ab_y + g_0_z_z_xxxyyyyy[k];

                g_0_z_yz_xxxyyyz[k] = -g_0_z_z_xxxyyyz[k] * ab_y + g_0_z_z_xxxyyyyz[k];

                g_0_z_yz_xxxyyzz[k] = -g_0_z_z_xxxyyzz[k] * ab_y + g_0_z_z_xxxyyyzz[k];

                g_0_z_yz_xxxyzzz[k] = -g_0_z_z_xxxyzzz[k] * ab_y + g_0_z_z_xxxyyzzz[k];

                g_0_z_yz_xxxzzzz[k] = -g_0_z_z_xxxzzzz[k] * ab_y + g_0_z_z_xxxyzzzz[k];

                g_0_z_yz_xxyyyyy[k] = -g_0_z_z_xxyyyyy[k] * ab_y + g_0_z_z_xxyyyyyy[k];

                g_0_z_yz_xxyyyyz[k] = -g_0_z_z_xxyyyyz[k] * ab_y + g_0_z_z_xxyyyyyz[k];

                g_0_z_yz_xxyyyzz[k] = -g_0_z_z_xxyyyzz[k] * ab_y + g_0_z_z_xxyyyyzz[k];

                g_0_z_yz_xxyyzzz[k] = -g_0_z_z_xxyyzzz[k] * ab_y + g_0_z_z_xxyyyzzz[k];

                g_0_z_yz_xxyzzzz[k] = -g_0_z_z_xxyzzzz[k] * ab_y + g_0_z_z_xxyyzzzz[k];

                g_0_z_yz_xxzzzzz[k] = -g_0_z_z_xxzzzzz[k] * ab_y + g_0_z_z_xxyzzzzz[k];

                g_0_z_yz_xyyyyyy[k] = -g_0_z_z_xyyyyyy[k] * ab_y + g_0_z_z_xyyyyyyy[k];

                g_0_z_yz_xyyyyyz[k] = -g_0_z_z_xyyyyyz[k] * ab_y + g_0_z_z_xyyyyyyz[k];

                g_0_z_yz_xyyyyzz[k] = -g_0_z_z_xyyyyzz[k] * ab_y + g_0_z_z_xyyyyyzz[k];

                g_0_z_yz_xyyyzzz[k] = -g_0_z_z_xyyyzzz[k] * ab_y + g_0_z_z_xyyyyzzz[k];

                g_0_z_yz_xyyzzzz[k] = -g_0_z_z_xyyzzzz[k] * ab_y + g_0_z_z_xyyyzzzz[k];

                g_0_z_yz_xyzzzzz[k] = -g_0_z_z_xyzzzzz[k] * ab_y + g_0_z_z_xyyzzzzz[k];

                g_0_z_yz_xzzzzzz[k] = -g_0_z_z_xzzzzzz[k] * ab_y + g_0_z_z_xyzzzzzz[k];

                g_0_z_yz_yyyyyyy[k] = -g_0_z_z_yyyyyyy[k] * ab_y + g_0_z_z_yyyyyyyy[k];

                g_0_z_yz_yyyyyyz[k] = -g_0_z_z_yyyyyyz[k] * ab_y + g_0_z_z_yyyyyyyz[k];

                g_0_z_yz_yyyyyzz[k] = -g_0_z_z_yyyyyzz[k] * ab_y + g_0_z_z_yyyyyyzz[k];

                g_0_z_yz_yyyyzzz[k] = -g_0_z_z_yyyyzzz[k] * ab_y + g_0_z_z_yyyyyzzz[k];

                g_0_z_yz_yyyzzzz[k] = -g_0_z_z_yyyzzzz[k] * ab_y + g_0_z_z_yyyyzzzz[k];

                g_0_z_yz_yyzzzzz[k] = -g_0_z_z_yyzzzzz[k] * ab_y + g_0_z_z_yyyzzzzz[k];

                g_0_z_yz_yzzzzzz[k] = -g_0_z_z_yzzzzzz[k] * ab_y + g_0_z_z_yyzzzzzz[k];

                g_0_z_yz_zzzzzzz[k] = -g_0_z_z_zzzzzzz[k] * ab_y + g_0_z_z_yzzzzzzz[k];
            }

            /// Set up 612-648 components of targeted buffer : cbuffer.data(

            auto g_0_z_zz_xxxxxxx = cbuffer.data(dk_geom_01_off + 612 * ccomps * dcomps);

            auto g_0_z_zz_xxxxxxy = cbuffer.data(dk_geom_01_off + 613 * ccomps * dcomps);

            auto g_0_z_zz_xxxxxxz = cbuffer.data(dk_geom_01_off + 614 * ccomps * dcomps);

            auto g_0_z_zz_xxxxxyy = cbuffer.data(dk_geom_01_off + 615 * ccomps * dcomps);

            auto g_0_z_zz_xxxxxyz = cbuffer.data(dk_geom_01_off + 616 * ccomps * dcomps);

            auto g_0_z_zz_xxxxxzz = cbuffer.data(dk_geom_01_off + 617 * ccomps * dcomps);

            auto g_0_z_zz_xxxxyyy = cbuffer.data(dk_geom_01_off + 618 * ccomps * dcomps);

            auto g_0_z_zz_xxxxyyz = cbuffer.data(dk_geom_01_off + 619 * ccomps * dcomps);

            auto g_0_z_zz_xxxxyzz = cbuffer.data(dk_geom_01_off + 620 * ccomps * dcomps);

            auto g_0_z_zz_xxxxzzz = cbuffer.data(dk_geom_01_off + 621 * ccomps * dcomps);

            auto g_0_z_zz_xxxyyyy = cbuffer.data(dk_geom_01_off + 622 * ccomps * dcomps);

            auto g_0_z_zz_xxxyyyz = cbuffer.data(dk_geom_01_off + 623 * ccomps * dcomps);

            auto g_0_z_zz_xxxyyzz = cbuffer.data(dk_geom_01_off + 624 * ccomps * dcomps);

            auto g_0_z_zz_xxxyzzz = cbuffer.data(dk_geom_01_off + 625 * ccomps * dcomps);

            auto g_0_z_zz_xxxzzzz = cbuffer.data(dk_geom_01_off + 626 * ccomps * dcomps);

            auto g_0_z_zz_xxyyyyy = cbuffer.data(dk_geom_01_off + 627 * ccomps * dcomps);

            auto g_0_z_zz_xxyyyyz = cbuffer.data(dk_geom_01_off + 628 * ccomps * dcomps);

            auto g_0_z_zz_xxyyyzz = cbuffer.data(dk_geom_01_off + 629 * ccomps * dcomps);

            auto g_0_z_zz_xxyyzzz = cbuffer.data(dk_geom_01_off + 630 * ccomps * dcomps);

            auto g_0_z_zz_xxyzzzz = cbuffer.data(dk_geom_01_off + 631 * ccomps * dcomps);

            auto g_0_z_zz_xxzzzzz = cbuffer.data(dk_geom_01_off + 632 * ccomps * dcomps);

            auto g_0_z_zz_xyyyyyy = cbuffer.data(dk_geom_01_off + 633 * ccomps * dcomps);

            auto g_0_z_zz_xyyyyyz = cbuffer.data(dk_geom_01_off + 634 * ccomps * dcomps);

            auto g_0_z_zz_xyyyyzz = cbuffer.data(dk_geom_01_off + 635 * ccomps * dcomps);

            auto g_0_z_zz_xyyyzzz = cbuffer.data(dk_geom_01_off + 636 * ccomps * dcomps);

            auto g_0_z_zz_xyyzzzz = cbuffer.data(dk_geom_01_off + 637 * ccomps * dcomps);

            auto g_0_z_zz_xyzzzzz = cbuffer.data(dk_geom_01_off + 638 * ccomps * dcomps);

            auto g_0_z_zz_xzzzzzz = cbuffer.data(dk_geom_01_off + 639 * ccomps * dcomps);

            auto g_0_z_zz_yyyyyyy = cbuffer.data(dk_geom_01_off + 640 * ccomps * dcomps);

            auto g_0_z_zz_yyyyyyz = cbuffer.data(dk_geom_01_off + 641 * ccomps * dcomps);

            auto g_0_z_zz_yyyyyzz = cbuffer.data(dk_geom_01_off + 642 * ccomps * dcomps);

            auto g_0_z_zz_yyyyzzz = cbuffer.data(dk_geom_01_off + 643 * ccomps * dcomps);

            auto g_0_z_zz_yyyzzzz = cbuffer.data(dk_geom_01_off + 644 * ccomps * dcomps);

            auto g_0_z_zz_yyzzzzz = cbuffer.data(dk_geom_01_off + 645 * ccomps * dcomps);

            auto g_0_z_zz_yzzzzzz = cbuffer.data(dk_geom_01_off + 646 * ccomps * dcomps);

            auto g_0_z_zz_zzzzzzz = cbuffer.data(dk_geom_01_off + 647 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_z_xxxxxxx, g_0_z_z_xxxxxxxz, g_0_z_z_xxxxxxy, g_0_z_z_xxxxxxyz, g_0_z_z_xxxxxxz, g_0_z_z_xxxxxxzz, g_0_z_z_xxxxxyy, g_0_z_z_xxxxxyyz, g_0_z_z_xxxxxyz, g_0_z_z_xxxxxyzz, g_0_z_z_xxxxxzz, g_0_z_z_xxxxxzzz, g_0_z_z_xxxxyyy, g_0_z_z_xxxxyyyz, g_0_z_z_xxxxyyz, g_0_z_z_xxxxyyzz, g_0_z_z_xxxxyzz, g_0_z_z_xxxxyzzz, g_0_z_z_xxxxzzz, g_0_z_z_xxxxzzzz, g_0_z_z_xxxyyyy, g_0_z_z_xxxyyyyz, g_0_z_z_xxxyyyz, g_0_z_z_xxxyyyzz, g_0_z_z_xxxyyzz, g_0_z_z_xxxyyzzz, g_0_z_z_xxxyzzz, g_0_z_z_xxxyzzzz, g_0_z_z_xxxzzzz, g_0_z_z_xxxzzzzz, g_0_z_z_xxyyyyy, g_0_z_z_xxyyyyyz, g_0_z_z_xxyyyyz, g_0_z_z_xxyyyyzz, g_0_z_z_xxyyyzz, g_0_z_z_xxyyyzzz, g_0_z_z_xxyyzzz, g_0_z_z_xxyyzzzz, g_0_z_z_xxyzzzz, g_0_z_z_xxyzzzzz, g_0_z_z_xxzzzzz, g_0_z_z_xxzzzzzz, g_0_z_z_xyyyyyy, g_0_z_z_xyyyyyyz, g_0_z_z_xyyyyyz, g_0_z_z_xyyyyyzz, g_0_z_z_xyyyyzz, g_0_z_z_xyyyyzzz, g_0_z_z_xyyyzzz, g_0_z_z_xyyyzzzz, g_0_z_z_xyyzzzz, g_0_z_z_xyyzzzzz, g_0_z_z_xyzzzzz, g_0_z_z_xyzzzzzz, g_0_z_z_xzzzzzz, g_0_z_z_xzzzzzzz, g_0_z_z_yyyyyyy, g_0_z_z_yyyyyyyz, g_0_z_z_yyyyyyz, g_0_z_z_yyyyyyzz, g_0_z_z_yyyyyzz, g_0_z_z_yyyyyzzz, g_0_z_z_yyyyzzz, g_0_z_z_yyyyzzzz, g_0_z_z_yyyzzzz, g_0_z_z_yyyzzzzz, g_0_z_z_yyzzzzz, g_0_z_z_yyzzzzzz, g_0_z_z_yzzzzzz, g_0_z_z_yzzzzzzz, g_0_z_z_zzzzzzz, g_0_z_z_zzzzzzzz, g_0_z_zz_xxxxxxx, g_0_z_zz_xxxxxxy, g_0_z_zz_xxxxxxz, g_0_z_zz_xxxxxyy, g_0_z_zz_xxxxxyz, g_0_z_zz_xxxxxzz, g_0_z_zz_xxxxyyy, g_0_z_zz_xxxxyyz, g_0_z_zz_xxxxyzz, g_0_z_zz_xxxxzzz, g_0_z_zz_xxxyyyy, g_0_z_zz_xxxyyyz, g_0_z_zz_xxxyyzz, g_0_z_zz_xxxyzzz, g_0_z_zz_xxxzzzz, g_0_z_zz_xxyyyyy, g_0_z_zz_xxyyyyz, g_0_z_zz_xxyyyzz, g_0_z_zz_xxyyzzz, g_0_z_zz_xxyzzzz, g_0_z_zz_xxzzzzz, g_0_z_zz_xyyyyyy, g_0_z_zz_xyyyyyz, g_0_z_zz_xyyyyzz, g_0_z_zz_xyyyzzz, g_0_z_zz_xyyzzzz, g_0_z_zz_xyzzzzz, g_0_z_zz_xzzzzzz, g_0_z_zz_yyyyyyy, g_0_z_zz_yyyyyyz, g_0_z_zz_yyyyyzz, g_0_z_zz_yyyyzzz, g_0_z_zz_yyyzzzz, g_0_z_zz_yyzzzzz, g_0_z_zz_yzzzzzz, g_0_z_zz_zzzzzzz, g_z_xxxxxxx, g_z_xxxxxxy, g_z_xxxxxxz, g_z_xxxxxyy, g_z_xxxxxyz, g_z_xxxxxzz, g_z_xxxxyyy, g_z_xxxxyyz, g_z_xxxxyzz, g_z_xxxxzzz, g_z_xxxyyyy, g_z_xxxyyyz, g_z_xxxyyzz, g_z_xxxyzzz, g_z_xxxzzzz, g_z_xxyyyyy, g_z_xxyyyyz, g_z_xxyyyzz, g_z_xxyyzzz, g_z_xxyzzzz, g_z_xxzzzzz, g_z_xyyyyyy, g_z_xyyyyyz, g_z_xyyyyzz, g_z_xyyyzzz, g_z_xyyzzzz, g_z_xyzzzzz, g_z_xzzzzzz, g_z_yyyyyyy, g_z_yyyyyyz, g_z_yyyyyzz, g_z_yyyyzzz, g_z_yyyzzzz, g_z_yyzzzzz, g_z_yzzzzzz, g_z_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_zz_xxxxxxx[k] = g_z_xxxxxxx[k] - g_0_z_z_xxxxxxx[k] * ab_z + g_0_z_z_xxxxxxxz[k];

                g_0_z_zz_xxxxxxy[k] = g_z_xxxxxxy[k] - g_0_z_z_xxxxxxy[k] * ab_z + g_0_z_z_xxxxxxyz[k];

                g_0_z_zz_xxxxxxz[k] = g_z_xxxxxxz[k] - g_0_z_z_xxxxxxz[k] * ab_z + g_0_z_z_xxxxxxzz[k];

                g_0_z_zz_xxxxxyy[k] = g_z_xxxxxyy[k] - g_0_z_z_xxxxxyy[k] * ab_z + g_0_z_z_xxxxxyyz[k];

                g_0_z_zz_xxxxxyz[k] = g_z_xxxxxyz[k] - g_0_z_z_xxxxxyz[k] * ab_z + g_0_z_z_xxxxxyzz[k];

                g_0_z_zz_xxxxxzz[k] = g_z_xxxxxzz[k] - g_0_z_z_xxxxxzz[k] * ab_z + g_0_z_z_xxxxxzzz[k];

                g_0_z_zz_xxxxyyy[k] = g_z_xxxxyyy[k] - g_0_z_z_xxxxyyy[k] * ab_z + g_0_z_z_xxxxyyyz[k];

                g_0_z_zz_xxxxyyz[k] = g_z_xxxxyyz[k] - g_0_z_z_xxxxyyz[k] * ab_z + g_0_z_z_xxxxyyzz[k];

                g_0_z_zz_xxxxyzz[k] = g_z_xxxxyzz[k] - g_0_z_z_xxxxyzz[k] * ab_z + g_0_z_z_xxxxyzzz[k];

                g_0_z_zz_xxxxzzz[k] = g_z_xxxxzzz[k] - g_0_z_z_xxxxzzz[k] * ab_z + g_0_z_z_xxxxzzzz[k];

                g_0_z_zz_xxxyyyy[k] = g_z_xxxyyyy[k] - g_0_z_z_xxxyyyy[k] * ab_z + g_0_z_z_xxxyyyyz[k];

                g_0_z_zz_xxxyyyz[k] = g_z_xxxyyyz[k] - g_0_z_z_xxxyyyz[k] * ab_z + g_0_z_z_xxxyyyzz[k];

                g_0_z_zz_xxxyyzz[k] = g_z_xxxyyzz[k] - g_0_z_z_xxxyyzz[k] * ab_z + g_0_z_z_xxxyyzzz[k];

                g_0_z_zz_xxxyzzz[k] = g_z_xxxyzzz[k] - g_0_z_z_xxxyzzz[k] * ab_z + g_0_z_z_xxxyzzzz[k];

                g_0_z_zz_xxxzzzz[k] = g_z_xxxzzzz[k] - g_0_z_z_xxxzzzz[k] * ab_z + g_0_z_z_xxxzzzzz[k];

                g_0_z_zz_xxyyyyy[k] = g_z_xxyyyyy[k] - g_0_z_z_xxyyyyy[k] * ab_z + g_0_z_z_xxyyyyyz[k];

                g_0_z_zz_xxyyyyz[k] = g_z_xxyyyyz[k] - g_0_z_z_xxyyyyz[k] * ab_z + g_0_z_z_xxyyyyzz[k];

                g_0_z_zz_xxyyyzz[k] = g_z_xxyyyzz[k] - g_0_z_z_xxyyyzz[k] * ab_z + g_0_z_z_xxyyyzzz[k];

                g_0_z_zz_xxyyzzz[k] = g_z_xxyyzzz[k] - g_0_z_z_xxyyzzz[k] * ab_z + g_0_z_z_xxyyzzzz[k];

                g_0_z_zz_xxyzzzz[k] = g_z_xxyzzzz[k] - g_0_z_z_xxyzzzz[k] * ab_z + g_0_z_z_xxyzzzzz[k];

                g_0_z_zz_xxzzzzz[k] = g_z_xxzzzzz[k] - g_0_z_z_xxzzzzz[k] * ab_z + g_0_z_z_xxzzzzzz[k];

                g_0_z_zz_xyyyyyy[k] = g_z_xyyyyyy[k] - g_0_z_z_xyyyyyy[k] * ab_z + g_0_z_z_xyyyyyyz[k];

                g_0_z_zz_xyyyyyz[k] = g_z_xyyyyyz[k] - g_0_z_z_xyyyyyz[k] * ab_z + g_0_z_z_xyyyyyzz[k];

                g_0_z_zz_xyyyyzz[k] = g_z_xyyyyzz[k] - g_0_z_z_xyyyyzz[k] * ab_z + g_0_z_z_xyyyyzzz[k];

                g_0_z_zz_xyyyzzz[k] = g_z_xyyyzzz[k] - g_0_z_z_xyyyzzz[k] * ab_z + g_0_z_z_xyyyzzzz[k];

                g_0_z_zz_xyyzzzz[k] = g_z_xyyzzzz[k] - g_0_z_z_xyyzzzz[k] * ab_z + g_0_z_z_xyyzzzzz[k];

                g_0_z_zz_xyzzzzz[k] = g_z_xyzzzzz[k] - g_0_z_z_xyzzzzz[k] * ab_z + g_0_z_z_xyzzzzzz[k];

                g_0_z_zz_xzzzzzz[k] = g_z_xzzzzzz[k] - g_0_z_z_xzzzzzz[k] * ab_z + g_0_z_z_xzzzzzzz[k];

                g_0_z_zz_yyyyyyy[k] = g_z_yyyyyyy[k] - g_0_z_z_yyyyyyy[k] * ab_z + g_0_z_z_yyyyyyyz[k];

                g_0_z_zz_yyyyyyz[k] = g_z_yyyyyyz[k] - g_0_z_z_yyyyyyz[k] * ab_z + g_0_z_z_yyyyyyzz[k];

                g_0_z_zz_yyyyyzz[k] = g_z_yyyyyzz[k] - g_0_z_z_yyyyyzz[k] * ab_z + g_0_z_z_yyyyyzzz[k];

                g_0_z_zz_yyyyzzz[k] = g_z_yyyyzzz[k] - g_0_z_z_yyyyzzz[k] * ab_z + g_0_z_z_yyyyzzzz[k];

                g_0_z_zz_yyyzzzz[k] = g_z_yyyzzzz[k] - g_0_z_z_yyyzzzz[k] * ab_z + g_0_z_z_yyyzzzzz[k];

                g_0_z_zz_yyzzzzz[k] = g_z_yyzzzzz[k] - g_0_z_z_yyzzzzz[k] * ab_z + g_0_z_z_yyzzzzzz[k];

                g_0_z_zz_yzzzzzz[k] = g_z_yzzzzzz[k] - g_0_z_z_yzzzzzz[k] * ab_z + g_0_z_z_yzzzzzzz[k];

                g_0_z_zz_zzzzzzz[k] = g_z_zzzzzzz[k] - g_0_z_z_zzzzzzz[k] * ab_z + g_0_z_z_zzzzzzzz[k];
            }
        }
    }
}

} // erirec namespace

