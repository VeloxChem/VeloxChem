#include "ElectronRepulsionGeom0100ContrRecIDXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_idxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_idxx,
                                            const size_t idx_hdxx,
                                            const size_t idx_geom_01_hdxx,
                                            const size_t idx_geom_01_hfxx,
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
            /// Set up components of auxilary buffer : HDSS

            const auto hd_off = idx_hdxx + i * dcomps + j;

            auto g_xxxxx_xx = cbuffer.data(hd_off + 0 * ccomps * dcomps);

            auto g_xxxxx_xy = cbuffer.data(hd_off + 1 * ccomps * dcomps);

            auto g_xxxxx_xz = cbuffer.data(hd_off + 2 * ccomps * dcomps);

            auto g_xxxxx_yy = cbuffer.data(hd_off + 3 * ccomps * dcomps);

            auto g_xxxxx_yz = cbuffer.data(hd_off + 4 * ccomps * dcomps);

            auto g_xxxxx_zz = cbuffer.data(hd_off + 5 * ccomps * dcomps);

            auto g_xxxxy_xx = cbuffer.data(hd_off + 6 * ccomps * dcomps);

            auto g_xxxxy_xy = cbuffer.data(hd_off + 7 * ccomps * dcomps);

            auto g_xxxxy_xz = cbuffer.data(hd_off + 8 * ccomps * dcomps);

            auto g_xxxxy_yy = cbuffer.data(hd_off + 9 * ccomps * dcomps);

            auto g_xxxxy_yz = cbuffer.data(hd_off + 10 * ccomps * dcomps);

            auto g_xxxxy_zz = cbuffer.data(hd_off + 11 * ccomps * dcomps);

            auto g_xxxxz_xx = cbuffer.data(hd_off + 12 * ccomps * dcomps);

            auto g_xxxxz_xy = cbuffer.data(hd_off + 13 * ccomps * dcomps);

            auto g_xxxxz_xz = cbuffer.data(hd_off + 14 * ccomps * dcomps);

            auto g_xxxxz_yy = cbuffer.data(hd_off + 15 * ccomps * dcomps);

            auto g_xxxxz_yz = cbuffer.data(hd_off + 16 * ccomps * dcomps);

            auto g_xxxxz_zz = cbuffer.data(hd_off + 17 * ccomps * dcomps);

            auto g_xxxyy_xx = cbuffer.data(hd_off + 18 * ccomps * dcomps);

            auto g_xxxyy_xy = cbuffer.data(hd_off + 19 * ccomps * dcomps);

            auto g_xxxyy_xz = cbuffer.data(hd_off + 20 * ccomps * dcomps);

            auto g_xxxyy_yy = cbuffer.data(hd_off + 21 * ccomps * dcomps);

            auto g_xxxyy_yz = cbuffer.data(hd_off + 22 * ccomps * dcomps);

            auto g_xxxyy_zz = cbuffer.data(hd_off + 23 * ccomps * dcomps);

            auto g_xxxyz_xx = cbuffer.data(hd_off + 24 * ccomps * dcomps);

            auto g_xxxyz_xy = cbuffer.data(hd_off + 25 * ccomps * dcomps);

            auto g_xxxyz_xz = cbuffer.data(hd_off + 26 * ccomps * dcomps);

            auto g_xxxyz_yy = cbuffer.data(hd_off + 27 * ccomps * dcomps);

            auto g_xxxyz_yz = cbuffer.data(hd_off + 28 * ccomps * dcomps);

            auto g_xxxyz_zz = cbuffer.data(hd_off + 29 * ccomps * dcomps);

            auto g_xxxzz_xx = cbuffer.data(hd_off + 30 * ccomps * dcomps);

            auto g_xxxzz_xy = cbuffer.data(hd_off + 31 * ccomps * dcomps);

            auto g_xxxzz_xz = cbuffer.data(hd_off + 32 * ccomps * dcomps);

            auto g_xxxzz_yy = cbuffer.data(hd_off + 33 * ccomps * dcomps);

            auto g_xxxzz_yz = cbuffer.data(hd_off + 34 * ccomps * dcomps);

            auto g_xxxzz_zz = cbuffer.data(hd_off + 35 * ccomps * dcomps);

            auto g_xxyyy_xx = cbuffer.data(hd_off + 36 * ccomps * dcomps);

            auto g_xxyyy_xy = cbuffer.data(hd_off + 37 * ccomps * dcomps);

            auto g_xxyyy_xz = cbuffer.data(hd_off + 38 * ccomps * dcomps);

            auto g_xxyyy_yy = cbuffer.data(hd_off + 39 * ccomps * dcomps);

            auto g_xxyyy_yz = cbuffer.data(hd_off + 40 * ccomps * dcomps);

            auto g_xxyyy_zz = cbuffer.data(hd_off + 41 * ccomps * dcomps);

            auto g_xxyyz_xx = cbuffer.data(hd_off + 42 * ccomps * dcomps);

            auto g_xxyyz_xy = cbuffer.data(hd_off + 43 * ccomps * dcomps);

            auto g_xxyyz_xz = cbuffer.data(hd_off + 44 * ccomps * dcomps);

            auto g_xxyyz_yy = cbuffer.data(hd_off + 45 * ccomps * dcomps);

            auto g_xxyyz_yz = cbuffer.data(hd_off + 46 * ccomps * dcomps);

            auto g_xxyyz_zz = cbuffer.data(hd_off + 47 * ccomps * dcomps);

            auto g_xxyzz_xx = cbuffer.data(hd_off + 48 * ccomps * dcomps);

            auto g_xxyzz_xy = cbuffer.data(hd_off + 49 * ccomps * dcomps);

            auto g_xxyzz_xz = cbuffer.data(hd_off + 50 * ccomps * dcomps);

            auto g_xxyzz_yy = cbuffer.data(hd_off + 51 * ccomps * dcomps);

            auto g_xxyzz_yz = cbuffer.data(hd_off + 52 * ccomps * dcomps);

            auto g_xxyzz_zz = cbuffer.data(hd_off + 53 * ccomps * dcomps);

            auto g_xxzzz_xx = cbuffer.data(hd_off + 54 * ccomps * dcomps);

            auto g_xxzzz_xy = cbuffer.data(hd_off + 55 * ccomps * dcomps);

            auto g_xxzzz_xz = cbuffer.data(hd_off + 56 * ccomps * dcomps);

            auto g_xxzzz_yy = cbuffer.data(hd_off + 57 * ccomps * dcomps);

            auto g_xxzzz_yz = cbuffer.data(hd_off + 58 * ccomps * dcomps);

            auto g_xxzzz_zz = cbuffer.data(hd_off + 59 * ccomps * dcomps);

            auto g_xyyyy_xx = cbuffer.data(hd_off + 60 * ccomps * dcomps);

            auto g_xyyyy_xy = cbuffer.data(hd_off + 61 * ccomps * dcomps);

            auto g_xyyyy_xz = cbuffer.data(hd_off + 62 * ccomps * dcomps);

            auto g_xyyyy_yy = cbuffer.data(hd_off + 63 * ccomps * dcomps);

            auto g_xyyyy_yz = cbuffer.data(hd_off + 64 * ccomps * dcomps);

            auto g_xyyyy_zz = cbuffer.data(hd_off + 65 * ccomps * dcomps);

            auto g_xyyyz_xx = cbuffer.data(hd_off + 66 * ccomps * dcomps);

            auto g_xyyyz_xy = cbuffer.data(hd_off + 67 * ccomps * dcomps);

            auto g_xyyyz_xz = cbuffer.data(hd_off + 68 * ccomps * dcomps);

            auto g_xyyyz_yy = cbuffer.data(hd_off + 69 * ccomps * dcomps);

            auto g_xyyyz_yz = cbuffer.data(hd_off + 70 * ccomps * dcomps);

            auto g_xyyyz_zz = cbuffer.data(hd_off + 71 * ccomps * dcomps);

            auto g_xyyzz_xx = cbuffer.data(hd_off + 72 * ccomps * dcomps);

            auto g_xyyzz_xy = cbuffer.data(hd_off + 73 * ccomps * dcomps);

            auto g_xyyzz_xz = cbuffer.data(hd_off + 74 * ccomps * dcomps);

            auto g_xyyzz_yy = cbuffer.data(hd_off + 75 * ccomps * dcomps);

            auto g_xyyzz_yz = cbuffer.data(hd_off + 76 * ccomps * dcomps);

            auto g_xyyzz_zz = cbuffer.data(hd_off + 77 * ccomps * dcomps);

            auto g_xyzzz_xx = cbuffer.data(hd_off + 78 * ccomps * dcomps);

            auto g_xyzzz_xy = cbuffer.data(hd_off + 79 * ccomps * dcomps);

            auto g_xyzzz_xz = cbuffer.data(hd_off + 80 * ccomps * dcomps);

            auto g_xyzzz_yy = cbuffer.data(hd_off + 81 * ccomps * dcomps);

            auto g_xyzzz_yz = cbuffer.data(hd_off + 82 * ccomps * dcomps);

            auto g_xyzzz_zz = cbuffer.data(hd_off + 83 * ccomps * dcomps);

            auto g_xzzzz_xx = cbuffer.data(hd_off + 84 * ccomps * dcomps);

            auto g_xzzzz_xy = cbuffer.data(hd_off + 85 * ccomps * dcomps);

            auto g_xzzzz_xz = cbuffer.data(hd_off + 86 * ccomps * dcomps);

            auto g_xzzzz_yy = cbuffer.data(hd_off + 87 * ccomps * dcomps);

            auto g_xzzzz_yz = cbuffer.data(hd_off + 88 * ccomps * dcomps);

            auto g_xzzzz_zz = cbuffer.data(hd_off + 89 * ccomps * dcomps);

            auto g_yyyyy_xx = cbuffer.data(hd_off + 90 * ccomps * dcomps);

            auto g_yyyyy_xy = cbuffer.data(hd_off + 91 * ccomps * dcomps);

            auto g_yyyyy_xz = cbuffer.data(hd_off + 92 * ccomps * dcomps);

            auto g_yyyyy_yy = cbuffer.data(hd_off + 93 * ccomps * dcomps);

            auto g_yyyyy_yz = cbuffer.data(hd_off + 94 * ccomps * dcomps);

            auto g_yyyyy_zz = cbuffer.data(hd_off + 95 * ccomps * dcomps);

            auto g_yyyyz_xx = cbuffer.data(hd_off + 96 * ccomps * dcomps);

            auto g_yyyyz_xy = cbuffer.data(hd_off + 97 * ccomps * dcomps);

            auto g_yyyyz_xz = cbuffer.data(hd_off + 98 * ccomps * dcomps);

            auto g_yyyyz_yy = cbuffer.data(hd_off + 99 * ccomps * dcomps);

            auto g_yyyyz_yz = cbuffer.data(hd_off + 100 * ccomps * dcomps);

            auto g_yyyyz_zz = cbuffer.data(hd_off + 101 * ccomps * dcomps);

            auto g_yyyzz_xx = cbuffer.data(hd_off + 102 * ccomps * dcomps);

            auto g_yyyzz_xy = cbuffer.data(hd_off + 103 * ccomps * dcomps);

            auto g_yyyzz_xz = cbuffer.data(hd_off + 104 * ccomps * dcomps);

            auto g_yyyzz_yy = cbuffer.data(hd_off + 105 * ccomps * dcomps);

            auto g_yyyzz_yz = cbuffer.data(hd_off + 106 * ccomps * dcomps);

            auto g_yyyzz_zz = cbuffer.data(hd_off + 107 * ccomps * dcomps);

            auto g_yyzzz_xx = cbuffer.data(hd_off + 108 * ccomps * dcomps);

            auto g_yyzzz_xy = cbuffer.data(hd_off + 109 * ccomps * dcomps);

            auto g_yyzzz_xz = cbuffer.data(hd_off + 110 * ccomps * dcomps);

            auto g_yyzzz_yy = cbuffer.data(hd_off + 111 * ccomps * dcomps);

            auto g_yyzzz_yz = cbuffer.data(hd_off + 112 * ccomps * dcomps);

            auto g_yyzzz_zz = cbuffer.data(hd_off + 113 * ccomps * dcomps);

            auto g_yzzzz_xx = cbuffer.data(hd_off + 114 * ccomps * dcomps);

            auto g_yzzzz_xy = cbuffer.data(hd_off + 115 * ccomps * dcomps);

            auto g_yzzzz_xz = cbuffer.data(hd_off + 116 * ccomps * dcomps);

            auto g_yzzzz_yy = cbuffer.data(hd_off + 117 * ccomps * dcomps);

            auto g_yzzzz_yz = cbuffer.data(hd_off + 118 * ccomps * dcomps);

            auto g_yzzzz_zz = cbuffer.data(hd_off + 119 * ccomps * dcomps);

            auto g_zzzzz_xx = cbuffer.data(hd_off + 120 * ccomps * dcomps);

            auto g_zzzzz_xy = cbuffer.data(hd_off + 121 * ccomps * dcomps);

            auto g_zzzzz_xz = cbuffer.data(hd_off + 122 * ccomps * dcomps);

            auto g_zzzzz_yy = cbuffer.data(hd_off + 123 * ccomps * dcomps);

            auto g_zzzzz_yz = cbuffer.data(hd_off + 124 * ccomps * dcomps);

            auto g_zzzzz_zz = cbuffer.data(hd_off + 125 * ccomps * dcomps);

            /// Set up components of auxilary buffer : HDSS

            const auto hd_geom_01_off = idx_geom_01_hdxx + i * dcomps + j;

            auto g_0_x_xxxxx_xx = cbuffer.data(hd_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxxx_xy = cbuffer.data(hd_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxxx_xz = cbuffer.data(hd_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxxx_yy = cbuffer.data(hd_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxxx_yz = cbuffer.data(hd_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxxx_zz = cbuffer.data(hd_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxxxy_xx = cbuffer.data(hd_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxxy_xy = cbuffer.data(hd_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxxy_xz = cbuffer.data(hd_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxxxy_yy = cbuffer.data(hd_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxxxy_yz = cbuffer.data(hd_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxxy_zz = cbuffer.data(hd_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxxxz_xx = cbuffer.data(hd_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxxz_xy = cbuffer.data(hd_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxxz_xz = cbuffer.data(hd_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxxxz_yy = cbuffer.data(hd_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxxz_yz = cbuffer.data(hd_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxxz_zz = cbuffer.data(hd_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxxyy_xx = cbuffer.data(hd_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxxyy_xy = cbuffer.data(hd_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxxyy_xz = cbuffer.data(hd_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxxyy_yy = cbuffer.data(hd_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxxyy_yz = cbuffer.data(hd_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxxyy_zz = cbuffer.data(hd_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxxyz_xx = cbuffer.data(hd_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxxyz_xy = cbuffer.data(hd_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxxyz_xz = cbuffer.data(hd_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxxyz_yy = cbuffer.data(hd_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxxyz_yz = cbuffer.data(hd_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxxyz_zz = cbuffer.data(hd_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xxxzz_xx = cbuffer.data(hd_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxxzz_xy = cbuffer.data(hd_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxxzz_xz = cbuffer.data(hd_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxxzz_yy = cbuffer.data(hd_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxxzz_yz = cbuffer.data(hd_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxxzz_zz = cbuffer.data(hd_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xxyyy_xx = cbuffer.data(hd_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxyyy_xy = cbuffer.data(hd_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxyyy_xz = cbuffer.data(hd_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxyyy_yy = cbuffer.data(hd_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxyyy_yz = cbuffer.data(hd_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxyyy_zz = cbuffer.data(hd_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xxyyz_xx = cbuffer.data(hd_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxyyz_xy = cbuffer.data(hd_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxyyz_xz = cbuffer.data(hd_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xxyyz_yy = cbuffer.data(hd_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xxyyz_yz = cbuffer.data(hd_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xxyyz_zz = cbuffer.data(hd_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xxyzz_xx = cbuffer.data(hd_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xxyzz_xy = cbuffer.data(hd_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xxyzz_xz = cbuffer.data(hd_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xxyzz_yy = cbuffer.data(hd_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xxyzz_yz = cbuffer.data(hd_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xxyzz_zz = cbuffer.data(hd_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xxzzz_xx = cbuffer.data(hd_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xxzzz_xy = cbuffer.data(hd_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xxzzz_xz = cbuffer.data(hd_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xxzzz_yy = cbuffer.data(hd_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xxzzz_yz = cbuffer.data(hd_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xxzzz_zz = cbuffer.data(hd_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_xyyyy_xx = cbuffer.data(hd_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xyyyy_xy = cbuffer.data(hd_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xyyyy_xz = cbuffer.data(hd_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xyyyy_yy = cbuffer.data(hd_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xyyyy_yz = cbuffer.data(hd_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xyyyy_zz = cbuffer.data(hd_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xyyyz_xx = cbuffer.data(hd_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xyyyz_xy = cbuffer.data(hd_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xyyyz_xz = cbuffer.data(hd_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xyyyz_yy = cbuffer.data(hd_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xyyyz_yz = cbuffer.data(hd_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xyyyz_zz = cbuffer.data(hd_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_xyyzz_xx = cbuffer.data(hd_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xyyzz_xy = cbuffer.data(hd_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xyyzz_xz = cbuffer.data(hd_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xyyzz_yy = cbuffer.data(hd_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xyyzz_yz = cbuffer.data(hd_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xyyzz_zz = cbuffer.data(hd_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xyzzz_xx = cbuffer.data(hd_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xyzzz_xy = cbuffer.data(hd_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xyzzz_xz = cbuffer.data(hd_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xyzzz_yy = cbuffer.data(hd_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xyzzz_yz = cbuffer.data(hd_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xyzzz_zz = cbuffer.data(hd_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xzzzz_xx = cbuffer.data(hd_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xzzzz_xy = cbuffer.data(hd_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xzzzz_xz = cbuffer.data(hd_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xzzzz_yy = cbuffer.data(hd_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xzzzz_yz = cbuffer.data(hd_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xzzzz_zz = cbuffer.data(hd_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_yyyyy_xx = cbuffer.data(hd_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_yyyyy_xy = cbuffer.data(hd_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_yyyyy_xz = cbuffer.data(hd_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_yyyyy_yy = cbuffer.data(hd_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_yyyyy_yz = cbuffer.data(hd_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_yyyyy_zz = cbuffer.data(hd_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_yyyyz_xx = cbuffer.data(hd_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_yyyyz_xy = cbuffer.data(hd_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_yyyyz_xz = cbuffer.data(hd_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_yyyyz_yy = cbuffer.data(hd_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_yyyyz_yz = cbuffer.data(hd_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_yyyyz_zz = cbuffer.data(hd_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_yyyzz_xx = cbuffer.data(hd_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_yyyzz_xy = cbuffer.data(hd_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_yyyzz_xz = cbuffer.data(hd_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_yyyzz_yy = cbuffer.data(hd_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_yyyzz_yz = cbuffer.data(hd_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_yyyzz_zz = cbuffer.data(hd_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_yyzzz_xx = cbuffer.data(hd_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_yyzzz_xy = cbuffer.data(hd_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_yyzzz_xz = cbuffer.data(hd_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_yyzzz_yy = cbuffer.data(hd_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_yyzzz_yz = cbuffer.data(hd_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_yyzzz_zz = cbuffer.data(hd_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_yzzzz_xx = cbuffer.data(hd_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_yzzzz_xy = cbuffer.data(hd_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_yzzzz_xz = cbuffer.data(hd_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_yzzzz_yy = cbuffer.data(hd_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_yzzzz_yz = cbuffer.data(hd_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_yzzzz_zz = cbuffer.data(hd_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_x_zzzzz_xx = cbuffer.data(hd_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_zzzzz_xy = cbuffer.data(hd_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_zzzzz_xz = cbuffer.data(hd_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_zzzzz_yy = cbuffer.data(hd_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_zzzzz_yz = cbuffer.data(hd_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_zzzzz_zz = cbuffer.data(hd_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_y_xxxxx_xx = cbuffer.data(hd_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_y_xxxxx_xy = cbuffer.data(hd_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_y_xxxxx_xz = cbuffer.data(hd_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_y_xxxxx_yy = cbuffer.data(hd_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_y_xxxxx_yz = cbuffer.data(hd_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_y_xxxxx_zz = cbuffer.data(hd_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_y_xxxxy_xx = cbuffer.data(hd_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_y_xxxxy_xy = cbuffer.data(hd_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_y_xxxxy_xz = cbuffer.data(hd_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_y_xxxxy_yy = cbuffer.data(hd_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_y_xxxxy_yz = cbuffer.data(hd_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_y_xxxxy_zz = cbuffer.data(hd_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_y_xxxxz_xx = cbuffer.data(hd_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_y_xxxxz_xy = cbuffer.data(hd_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_y_xxxxz_xz = cbuffer.data(hd_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_y_xxxxz_yy = cbuffer.data(hd_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_y_xxxxz_yz = cbuffer.data(hd_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_y_xxxxz_zz = cbuffer.data(hd_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_y_xxxyy_xx = cbuffer.data(hd_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_y_xxxyy_xy = cbuffer.data(hd_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_y_xxxyy_xz = cbuffer.data(hd_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_y_xxxyy_yy = cbuffer.data(hd_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_y_xxxyy_yz = cbuffer.data(hd_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_y_xxxyy_zz = cbuffer.data(hd_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_y_xxxyz_xx = cbuffer.data(hd_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_y_xxxyz_xy = cbuffer.data(hd_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_y_xxxyz_xz = cbuffer.data(hd_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_y_xxxyz_yy = cbuffer.data(hd_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_y_xxxyz_yz = cbuffer.data(hd_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_y_xxxyz_zz = cbuffer.data(hd_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_y_xxxzz_xx = cbuffer.data(hd_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_y_xxxzz_xy = cbuffer.data(hd_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_y_xxxzz_xz = cbuffer.data(hd_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_y_xxxzz_yy = cbuffer.data(hd_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_y_xxxzz_yz = cbuffer.data(hd_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_y_xxxzz_zz = cbuffer.data(hd_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_y_xxyyy_xx = cbuffer.data(hd_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_y_xxyyy_xy = cbuffer.data(hd_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_y_xxyyy_xz = cbuffer.data(hd_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_y_xxyyy_yy = cbuffer.data(hd_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_y_xxyyy_yz = cbuffer.data(hd_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_y_xxyyy_zz = cbuffer.data(hd_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_y_xxyyz_xx = cbuffer.data(hd_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_y_xxyyz_xy = cbuffer.data(hd_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_y_xxyyz_xz = cbuffer.data(hd_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_y_xxyyz_yy = cbuffer.data(hd_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_y_xxyyz_yz = cbuffer.data(hd_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_y_xxyyz_zz = cbuffer.data(hd_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_y_xxyzz_xx = cbuffer.data(hd_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_y_xxyzz_xy = cbuffer.data(hd_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_y_xxyzz_xz = cbuffer.data(hd_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_y_xxyzz_yy = cbuffer.data(hd_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_y_xxyzz_yz = cbuffer.data(hd_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_y_xxyzz_zz = cbuffer.data(hd_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_y_xxzzz_xx = cbuffer.data(hd_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_y_xxzzz_xy = cbuffer.data(hd_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_y_xxzzz_xz = cbuffer.data(hd_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_y_xxzzz_yy = cbuffer.data(hd_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_y_xxzzz_yz = cbuffer.data(hd_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_y_xxzzz_zz = cbuffer.data(hd_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_y_xyyyy_xx = cbuffer.data(hd_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_y_xyyyy_xy = cbuffer.data(hd_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_y_xyyyy_xz = cbuffer.data(hd_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_y_xyyyy_yy = cbuffer.data(hd_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_y_xyyyy_yz = cbuffer.data(hd_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_y_xyyyy_zz = cbuffer.data(hd_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_y_xyyyz_xx = cbuffer.data(hd_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_y_xyyyz_xy = cbuffer.data(hd_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_y_xyyyz_xz = cbuffer.data(hd_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_y_xyyyz_yy = cbuffer.data(hd_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_y_xyyyz_yz = cbuffer.data(hd_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_y_xyyyz_zz = cbuffer.data(hd_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_y_xyyzz_xx = cbuffer.data(hd_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_y_xyyzz_xy = cbuffer.data(hd_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_y_xyyzz_xz = cbuffer.data(hd_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_y_xyyzz_yy = cbuffer.data(hd_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_y_xyyzz_yz = cbuffer.data(hd_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_y_xyyzz_zz = cbuffer.data(hd_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_y_xyzzz_xx = cbuffer.data(hd_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_y_xyzzz_xy = cbuffer.data(hd_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_y_xyzzz_xz = cbuffer.data(hd_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_y_xyzzz_yy = cbuffer.data(hd_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_y_xyzzz_yz = cbuffer.data(hd_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_y_xyzzz_zz = cbuffer.data(hd_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_y_xzzzz_xx = cbuffer.data(hd_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_y_xzzzz_xy = cbuffer.data(hd_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_y_xzzzz_xz = cbuffer.data(hd_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_y_xzzzz_yy = cbuffer.data(hd_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_y_xzzzz_yz = cbuffer.data(hd_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_y_xzzzz_zz = cbuffer.data(hd_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_y_yyyyy_xx = cbuffer.data(hd_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_y_yyyyy_xy = cbuffer.data(hd_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_y_yyyyy_xz = cbuffer.data(hd_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_y_yyyyy_yy = cbuffer.data(hd_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_y_yyyyy_yz = cbuffer.data(hd_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_y_yyyyy_zz = cbuffer.data(hd_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_y_yyyyz_xx = cbuffer.data(hd_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_y_yyyyz_xy = cbuffer.data(hd_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_y_yyyyz_xz = cbuffer.data(hd_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_y_yyyyz_yy = cbuffer.data(hd_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_y_yyyyz_yz = cbuffer.data(hd_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_y_yyyyz_zz = cbuffer.data(hd_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_y_yyyzz_xx = cbuffer.data(hd_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_y_yyyzz_xy = cbuffer.data(hd_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_y_yyyzz_xz = cbuffer.data(hd_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_y_yyyzz_yy = cbuffer.data(hd_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_y_yyyzz_yz = cbuffer.data(hd_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_y_yyyzz_zz = cbuffer.data(hd_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_y_yyzzz_xx = cbuffer.data(hd_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_y_yyzzz_xy = cbuffer.data(hd_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_y_yyzzz_xz = cbuffer.data(hd_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_y_yyzzz_yy = cbuffer.data(hd_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_y_yyzzz_yz = cbuffer.data(hd_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_y_yyzzz_zz = cbuffer.data(hd_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_y_yzzzz_xx = cbuffer.data(hd_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_y_yzzzz_xy = cbuffer.data(hd_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_y_yzzzz_xz = cbuffer.data(hd_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_y_yzzzz_yy = cbuffer.data(hd_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_y_yzzzz_yz = cbuffer.data(hd_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_y_yzzzz_zz = cbuffer.data(hd_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_y_zzzzz_xx = cbuffer.data(hd_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_y_zzzzz_xy = cbuffer.data(hd_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_y_zzzzz_xz = cbuffer.data(hd_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_y_zzzzz_yy = cbuffer.data(hd_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_y_zzzzz_yz = cbuffer.data(hd_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_y_zzzzz_zz = cbuffer.data(hd_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_z_xxxxx_xx = cbuffer.data(hd_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_z_xxxxx_xy = cbuffer.data(hd_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_z_xxxxx_xz = cbuffer.data(hd_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_z_xxxxx_yy = cbuffer.data(hd_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_z_xxxxx_yz = cbuffer.data(hd_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_z_xxxxx_zz = cbuffer.data(hd_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_z_xxxxy_xx = cbuffer.data(hd_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_z_xxxxy_xy = cbuffer.data(hd_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_z_xxxxy_xz = cbuffer.data(hd_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_z_xxxxy_yy = cbuffer.data(hd_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_z_xxxxy_yz = cbuffer.data(hd_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_z_xxxxy_zz = cbuffer.data(hd_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_z_xxxxz_xx = cbuffer.data(hd_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_z_xxxxz_xy = cbuffer.data(hd_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_z_xxxxz_xz = cbuffer.data(hd_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_z_xxxxz_yy = cbuffer.data(hd_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_z_xxxxz_yz = cbuffer.data(hd_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_z_xxxxz_zz = cbuffer.data(hd_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_z_xxxyy_xx = cbuffer.data(hd_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_z_xxxyy_xy = cbuffer.data(hd_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_z_xxxyy_xz = cbuffer.data(hd_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_z_xxxyy_yy = cbuffer.data(hd_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_z_xxxyy_yz = cbuffer.data(hd_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_z_xxxyy_zz = cbuffer.data(hd_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_z_xxxyz_xx = cbuffer.data(hd_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_z_xxxyz_xy = cbuffer.data(hd_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_z_xxxyz_xz = cbuffer.data(hd_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_z_xxxyz_yy = cbuffer.data(hd_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_z_xxxyz_yz = cbuffer.data(hd_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_z_xxxyz_zz = cbuffer.data(hd_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_z_xxxzz_xx = cbuffer.data(hd_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_z_xxxzz_xy = cbuffer.data(hd_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_z_xxxzz_xz = cbuffer.data(hd_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_z_xxxzz_yy = cbuffer.data(hd_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_z_xxxzz_yz = cbuffer.data(hd_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_z_xxxzz_zz = cbuffer.data(hd_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_z_xxyyy_xx = cbuffer.data(hd_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_z_xxyyy_xy = cbuffer.data(hd_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_z_xxyyy_xz = cbuffer.data(hd_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_z_xxyyy_yy = cbuffer.data(hd_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_z_xxyyy_yz = cbuffer.data(hd_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_z_xxyyy_zz = cbuffer.data(hd_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_z_xxyyz_xx = cbuffer.data(hd_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_z_xxyyz_xy = cbuffer.data(hd_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_z_xxyyz_xz = cbuffer.data(hd_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_z_xxyyz_yy = cbuffer.data(hd_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_z_xxyyz_yz = cbuffer.data(hd_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_z_xxyyz_zz = cbuffer.data(hd_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_z_xxyzz_xx = cbuffer.data(hd_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_z_xxyzz_xy = cbuffer.data(hd_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_z_xxyzz_xz = cbuffer.data(hd_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_z_xxyzz_yy = cbuffer.data(hd_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_z_xxyzz_yz = cbuffer.data(hd_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_z_xxyzz_zz = cbuffer.data(hd_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_z_xxzzz_xx = cbuffer.data(hd_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_z_xxzzz_xy = cbuffer.data(hd_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_z_xxzzz_xz = cbuffer.data(hd_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_z_xxzzz_yy = cbuffer.data(hd_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_z_xxzzz_yz = cbuffer.data(hd_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_z_xxzzz_zz = cbuffer.data(hd_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_z_xyyyy_xx = cbuffer.data(hd_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_z_xyyyy_xy = cbuffer.data(hd_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_z_xyyyy_xz = cbuffer.data(hd_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_z_xyyyy_yy = cbuffer.data(hd_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_z_xyyyy_yz = cbuffer.data(hd_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_z_xyyyy_zz = cbuffer.data(hd_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_z_xyyyz_xx = cbuffer.data(hd_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_z_xyyyz_xy = cbuffer.data(hd_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_z_xyyyz_xz = cbuffer.data(hd_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_z_xyyyz_yy = cbuffer.data(hd_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_z_xyyyz_yz = cbuffer.data(hd_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_z_xyyyz_zz = cbuffer.data(hd_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_z_xyyzz_xx = cbuffer.data(hd_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_z_xyyzz_xy = cbuffer.data(hd_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_z_xyyzz_xz = cbuffer.data(hd_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_z_xyyzz_yy = cbuffer.data(hd_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_z_xyyzz_yz = cbuffer.data(hd_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_z_xyyzz_zz = cbuffer.data(hd_geom_01_off + 329 * ccomps * dcomps);

            auto g_0_z_xyzzz_xx = cbuffer.data(hd_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_z_xyzzz_xy = cbuffer.data(hd_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_z_xyzzz_xz = cbuffer.data(hd_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_z_xyzzz_yy = cbuffer.data(hd_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_z_xyzzz_yz = cbuffer.data(hd_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_z_xyzzz_zz = cbuffer.data(hd_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_z_xzzzz_xx = cbuffer.data(hd_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_z_xzzzz_xy = cbuffer.data(hd_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_z_xzzzz_xz = cbuffer.data(hd_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_z_xzzzz_yy = cbuffer.data(hd_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_z_xzzzz_yz = cbuffer.data(hd_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_z_xzzzz_zz = cbuffer.data(hd_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_z_yyyyy_xx = cbuffer.data(hd_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_z_yyyyy_xy = cbuffer.data(hd_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_z_yyyyy_xz = cbuffer.data(hd_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_z_yyyyy_yy = cbuffer.data(hd_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_z_yyyyy_yz = cbuffer.data(hd_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_z_yyyyy_zz = cbuffer.data(hd_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_z_yyyyz_xx = cbuffer.data(hd_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_z_yyyyz_xy = cbuffer.data(hd_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_z_yyyyz_xz = cbuffer.data(hd_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_z_yyyyz_yy = cbuffer.data(hd_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_z_yyyyz_yz = cbuffer.data(hd_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_z_yyyyz_zz = cbuffer.data(hd_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_z_yyyzz_xx = cbuffer.data(hd_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_z_yyyzz_xy = cbuffer.data(hd_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_z_yyyzz_xz = cbuffer.data(hd_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_z_yyyzz_yy = cbuffer.data(hd_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_z_yyyzz_yz = cbuffer.data(hd_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_z_yyyzz_zz = cbuffer.data(hd_geom_01_off + 359 * ccomps * dcomps);

            auto g_0_z_yyzzz_xx = cbuffer.data(hd_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_z_yyzzz_xy = cbuffer.data(hd_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_z_yyzzz_xz = cbuffer.data(hd_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_z_yyzzz_yy = cbuffer.data(hd_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_z_yyzzz_yz = cbuffer.data(hd_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_z_yyzzz_zz = cbuffer.data(hd_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_z_yzzzz_xx = cbuffer.data(hd_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_z_yzzzz_xy = cbuffer.data(hd_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_z_yzzzz_xz = cbuffer.data(hd_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_z_yzzzz_yy = cbuffer.data(hd_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_z_yzzzz_yz = cbuffer.data(hd_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_z_yzzzz_zz = cbuffer.data(hd_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_z_zzzzz_xx = cbuffer.data(hd_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_z_zzzzz_xy = cbuffer.data(hd_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_z_zzzzz_xz = cbuffer.data(hd_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_z_zzzzz_yy = cbuffer.data(hd_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_z_zzzzz_yz = cbuffer.data(hd_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_z_zzzzz_zz = cbuffer.data(hd_geom_01_off + 377 * ccomps * dcomps);

            /// Set up components of auxilary buffer : HFSS

            const auto hf_geom_01_off = idx_geom_01_hfxx + i * dcomps + j;

            auto g_0_x_xxxxx_xxx = cbuffer.data(hf_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxy = cbuffer.data(hf_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxz = cbuffer.data(hf_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxxx_xyy = cbuffer.data(hf_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxxx_xyz = cbuffer.data(hf_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxxx_xzz = cbuffer.data(hf_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxxxx_yyy = cbuffer.data(hf_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxxx_yyz = cbuffer.data(hf_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxxx_yzz = cbuffer.data(hf_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxxxx_zzz = cbuffer.data(hf_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxx = cbuffer.data(hf_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxy = cbuffer.data(hf_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxz = cbuffer.data(hf_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxxy_xyy = cbuffer.data(hf_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxxy_xyz = cbuffer.data(hf_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxxxy_xzz = cbuffer.data(hf_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxxy_yyy = cbuffer.data(hf_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxxy_yyz = cbuffer.data(hf_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxxxy_yzz = cbuffer.data(hf_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxxxy_zzz = cbuffer.data(hf_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxx = cbuffer.data(hf_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxy = cbuffer.data(hf_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxz = cbuffer.data(hf_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxxxz_xyy = cbuffer.data(hf_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxxxz_xyz = cbuffer.data(hf_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxxxz_xzz = cbuffer.data(hf_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxxxz_yyy = cbuffer.data(hf_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxxxz_yyz = cbuffer.data(hf_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxxxz_yzz = cbuffer.data(hf_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxxxz_zzz = cbuffer.data(hf_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxx = cbuffer.data(hf_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxy = cbuffer.data(hf_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxz = cbuffer.data(hf_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxxyy_xyy = cbuffer.data(hf_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxxyy_xyz = cbuffer.data(hf_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxxyy_xzz = cbuffer.data(hf_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xxxyy_yyy = cbuffer.data(hf_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxxyy_yyz = cbuffer.data(hf_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxxyy_yzz = cbuffer.data(hf_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxxyy_zzz = cbuffer.data(hf_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxx = cbuffer.data(hf_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxy = cbuffer.data(hf_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxz = cbuffer.data(hf_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxxyz_xyy = cbuffer.data(hf_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxxyz_xyz = cbuffer.data(hf_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xxxyz_xzz = cbuffer.data(hf_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xxxyz_yyy = cbuffer.data(hf_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xxxyz_yyz = cbuffer.data(hf_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xxxyz_yzz = cbuffer.data(hf_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xxxyz_zzz = cbuffer.data(hf_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxx = cbuffer.data(hf_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxy = cbuffer.data(hf_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxz = cbuffer.data(hf_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xxxzz_xyy = cbuffer.data(hf_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xxxzz_xyz = cbuffer.data(hf_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xxxzz_xzz = cbuffer.data(hf_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xxxzz_yyy = cbuffer.data(hf_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xxxzz_yyz = cbuffer.data(hf_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xxxzz_yzz = cbuffer.data(hf_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xxxzz_zzz = cbuffer.data(hf_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxx = cbuffer.data(hf_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxy = cbuffer.data(hf_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxz = cbuffer.data(hf_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xxyyy_xyy = cbuffer.data(hf_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xxyyy_xyz = cbuffer.data(hf_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xxyyy_xzz = cbuffer.data(hf_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xxyyy_yyy = cbuffer.data(hf_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xxyyy_yyz = cbuffer.data(hf_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xxyyy_yzz = cbuffer.data(hf_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xxyyy_zzz = cbuffer.data(hf_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxx = cbuffer.data(hf_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxy = cbuffer.data(hf_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxz = cbuffer.data(hf_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xxyyz_xyy = cbuffer.data(hf_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xxyyz_xyz = cbuffer.data(hf_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xxyyz_xzz = cbuffer.data(hf_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xxyyz_yyy = cbuffer.data(hf_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xxyyz_yyz = cbuffer.data(hf_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xxyyz_yzz = cbuffer.data(hf_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xxyyz_zzz = cbuffer.data(hf_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxx = cbuffer.data(hf_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxy = cbuffer.data(hf_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxz = cbuffer.data(hf_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xxyzz_xyy = cbuffer.data(hf_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xxyzz_xyz = cbuffer.data(hf_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xxyzz_xzz = cbuffer.data(hf_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xxyzz_yyy = cbuffer.data(hf_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xxyzz_yyz = cbuffer.data(hf_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xxyzz_yzz = cbuffer.data(hf_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xxyzz_zzz = cbuffer.data(hf_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxx = cbuffer.data(hf_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxy = cbuffer.data(hf_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxz = cbuffer.data(hf_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_xxzzz_xyy = cbuffer.data(hf_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_xxzzz_xyz = cbuffer.data(hf_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_xxzzz_xzz = cbuffer.data(hf_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_xxzzz_yyy = cbuffer.data(hf_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_xxzzz_yyz = cbuffer.data(hf_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_xxzzz_yzz = cbuffer.data(hf_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_xxzzz_zzz = cbuffer.data(hf_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxx = cbuffer.data(hf_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxy = cbuffer.data(hf_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxz = cbuffer.data(hf_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_xyyyy_xyy = cbuffer.data(hf_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_xyyyy_xyz = cbuffer.data(hf_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_xyyyy_xzz = cbuffer.data(hf_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_xyyyy_yyy = cbuffer.data(hf_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_xyyyy_yyz = cbuffer.data(hf_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_xyyyy_yzz = cbuffer.data(hf_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_xyyyy_zzz = cbuffer.data(hf_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxx = cbuffer.data(hf_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxy = cbuffer.data(hf_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxz = cbuffer.data(hf_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_xyyyz_xyy = cbuffer.data(hf_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_xyyyz_xyz = cbuffer.data(hf_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_xyyyz_xzz = cbuffer.data(hf_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_xyyyz_yyy = cbuffer.data(hf_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_xyyyz_yyz = cbuffer.data(hf_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_xyyyz_yzz = cbuffer.data(hf_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_xyyyz_zzz = cbuffer.data(hf_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxx = cbuffer.data(hf_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxy = cbuffer.data(hf_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxz = cbuffer.data(hf_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_xyyzz_xyy = cbuffer.data(hf_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_xyyzz_xyz = cbuffer.data(hf_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_xyyzz_xzz = cbuffer.data(hf_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_x_xyyzz_yyy = cbuffer.data(hf_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_xyyzz_yyz = cbuffer.data(hf_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_xyyzz_yzz = cbuffer.data(hf_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_xyyzz_zzz = cbuffer.data(hf_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxx = cbuffer.data(hf_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxy = cbuffer.data(hf_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxz = cbuffer.data(hf_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_xyzzz_xyy = cbuffer.data(hf_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_xyzzz_xyz = cbuffer.data(hf_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_x_xyzzz_xzz = cbuffer.data(hf_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_xyzzz_yyy = cbuffer.data(hf_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_xyzzz_yyz = cbuffer.data(hf_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_x_xyzzz_yzz = cbuffer.data(hf_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_xyzzz_zzz = cbuffer.data(hf_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxx = cbuffer.data(hf_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxy = cbuffer.data(hf_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxz = cbuffer.data(hf_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_xzzzz_xyy = cbuffer.data(hf_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_x_xzzzz_xyz = cbuffer.data(hf_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_xzzzz_xzz = cbuffer.data(hf_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_xzzzz_yyy = cbuffer.data(hf_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_xzzzz_yyz = cbuffer.data(hf_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_xzzzz_yzz = cbuffer.data(hf_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_xzzzz_zzz = cbuffer.data(hf_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxx = cbuffer.data(hf_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxy = cbuffer.data(hf_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxz = cbuffer.data(hf_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_x_yyyyy_xyy = cbuffer.data(hf_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_x_yyyyy_xyz = cbuffer.data(hf_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_x_yyyyy_xzz = cbuffer.data(hf_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_x_yyyyy_yyy = cbuffer.data(hf_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_x_yyyyy_yyz = cbuffer.data(hf_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_x_yyyyy_yzz = cbuffer.data(hf_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_x_yyyyy_zzz = cbuffer.data(hf_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxx = cbuffer.data(hf_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxy = cbuffer.data(hf_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxz = cbuffer.data(hf_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_x_yyyyz_xyy = cbuffer.data(hf_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_x_yyyyz_xyz = cbuffer.data(hf_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_x_yyyyz_xzz = cbuffer.data(hf_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_x_yyyyz_yyy = cbuffer.data(hf_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_x_yyyyz_yyz = cbuffer.data(hf_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_x_yyyyz_yzz = cbuffer.data(hf_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_x_yyyyz_zzz = cbuffer.data(hf_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxx = cbuffer.data(hf_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxy = cbuffer.data(hf_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxz = cbuffer.data(hf_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_x_yyyzz_xyy = cbuffer.data(hf_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_x_yyyzz_xyz = cbuffer.data(hf_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_x_yyyzz_xzz = cbuffer.data(hf_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_x_yyyzz_yyy = cbuffer.data(hf_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_x_yyyzz_yyz = cbuffer.data(hf_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_x_yyyzz_yzz = cbuffer.data(hf_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_x_yyyzz_zzz = cbuffer.data(hf_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxx = cbuffer.data(hf_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxy = cbuffer.data(hf_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxz = cbuffer.data(hf_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_x_yyzzz_xyy = cbuffer.data(hf_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_x_yyzzz_xyz = cbuffer.data(hf_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_x_yyzzz_xzz = cbuffer.data(hf_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_x_yyzzz_yyy = cbuffer.data(hf_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_x_yyzzz_yyz = cbuffer.data(hf_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_x_yyzzz_yzz = cbuffer.data(hf_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_x_yyzzz_zzz = cbuffer.data(hf_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxx = cbuffer.data(hf_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxy = cbuffer.data(hf_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxz = cbuffer.data(hf_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_x_yzzzz_xyy = cbuffer.data(hf_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_x_yzzzz_xyz = cbuffer.data(hf_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_x_yzzzz_xzz = cbuffer.data(hf_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_x_yzzzz_yyy = cbuffer.data(hf_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_x_yzzzz_yyz = cbuffer.data(hf_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_x_yzzzz_yzz = cbuffer.data(hf_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_x_yzzzz_zzz = cbuffer.data(hf_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxx = cbuffer.data(hf_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxy = cbuffer.data(hf_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxz = cbuffer.data(hf_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_x_zzzzz_xyy = cbuffer.data(hf_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_x_zzzzz_xyz = cbuffer.data(hf_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_x_zzzzz_xzz = cbuffer.data(hf_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_x_zzzzz_yyy = cbuffer.data(hf_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_x_zzzzz_yyz = cbuffer.data(hf_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_x_zzzzz_yzz = cbuffer.data(hf_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_x_zzzzz_zzz = cbuffer.data(hf_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxx = cbuffer.data(hf_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxy = cbuffer.data(hf_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxz = cbuffer.data(hf_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_y_xxxxx_xyy = cbuffer.data(hf_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_y_xxxxx_xyz = cbuffer.data(hf_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_y_xxxxx_xzz = cbuffer.data(hf_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_y_xxxxx_yyy = cbuffer.data(hf_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_y_xxxxx_yyz = cbuffer.data(hf_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_y_xxxxx_yzz = cbuffer.data(hf_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_y_xxxxx_zzz = cbuffer.data(hf_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxx = cbuffer.data(hf_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxy = cbuffer.data(hf_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxz = cbuffer.data(hf_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_y_xxxxy_xyy = cbuffer.data(hf_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_y_xxxxy_xyz = cbuffer.data(hf_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_y_xxxxy_xzz = cbuffer.data(hf_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_y_xxxxy_yyy = cbuffer.data(hf_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_y_xxxxy_yyz = cbuffer.data(hf_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_y_xxxxy_yzz = cbuffer.data(hf_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_y_xxxxy_zzz = cbuffer.data(hf_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxx = cbuffer.data(hf_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxy = cbuffer.data(hf_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxz = cbuffer.data(hf_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_y_xxxxz_xyy = cbuffer.data(hf_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_y_xxxxz_xyz = cbuffer.data(hf_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_y_xxxxz_xzz = cbuffer.data(hf_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_y_xxxxz_yyy = cbuffer.data(hf_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_y_xxxxz_yyz = cbuffer.data(hf_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_y_xxxxz_yzz = cbuffer.data(hf_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_y_xxxxz_zzz = cbuffer.data(hf_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxx = cbuffer.data(hf_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxy = cbuffer.data(hf_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxz = cbuffer.data(hf_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_y_xxxyy_xyy = cbuffer.data(hf_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_y_xxxyy_xyz = cbuffer.data(hf_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_y_xxxyy_xzz = cbuffer.data(hf_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_y_xxxyy_yyy = cbuffer.data(hf_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_y_xxxyy_yyz = cbuffer.data(hf_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_y_xxxyy_yzz = cbuffer.data(hf_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_y_xxxyy_zzz = cbuffer.data(hf_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxx = cbuffer.data(hf_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxy = cbuffer.data(hf_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxz = cbuffer.data(hf_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_y_xxxyz_xyy = cbuffer.data(hf_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_y_xxxyz_xyz = cbuffer.data(hf_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_y_xxxyz_xzz = cbuffer.data(hf_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_y_xxxyz_yyy = cbuffer.data(hf_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_y_xxxyz_yyz = cbuffer.data(hf_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_y_xxxyz_yzz = cbuffer.data(hf_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_y_xxxyz_zzz = cbuffer.data(hf_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxx = cbuffer.data(hf_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxy = cbuffer.data(hf_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxz = cbuffer.data(hf_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_y_xxxzz_xyy = cbuffer.data(hf_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_y_xxxzz_xyz = cbuffer.data(hf_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_y_xxxzz_xzz = cbuffer.data(hf_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_y_xxxzz_yyy = cbuffer.data(hf_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_y_xxxzz_yyz = cbuffer.data(hf_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_y_xxxzz_yzz = cbuffer.data(hf_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_y_xxxzz_zzz = cbuffer.data(hf_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxx = cbuffer.data(hf_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxy = cbuffer.data(hf_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxz = cbuffer.data(hf_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_y_xxyyy_xyy = cbuffer.data(hf_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_y_xxyyy_xyz = cbuffer.data(hf_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_y_xxyyy_xzz = cbuffer.data(hf_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_y_xxyyy_yyy = cbuffer.data(hf_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_y_xxyyy_yyz = cbuffer.data(hf_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_y_xxyyy_yzz = cbuffer.data(hf_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_y_xxyyy_zzz = cbuffer.data(hf_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxx = cbuffer.data(hf_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxy = cbuffer.data(hf_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxz = cbuffer.data(hf_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_y_xxyyz_xyy = cbuffer.data(hf_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_y_xxyyz_xyz = cbuffer.data(hf_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_y_xxyyz_xzz = cbuffer.data(hf_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_y_xxyyz_yyy = cbuffer.data(hf_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_y_xxyyz_yyz = cbuffer.data(hf_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_y_xxyyz_yzz = cbuffer.data(hf_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_y_xxyyz_zzz = cbuffer.data(hf_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxx = cbuffer.data(hf_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxy = cbuffer.data(hf_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxz = cbuffer.data(hf_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_y_xxyzz_xyy = cbuffer.data(hf_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_y_xxyzz_xyz = cbuffer.data(hf_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_y_xxyzz_xzz = cbuffer.data(hf_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_y_xxyzz_yyy = cbuffer.data(hf_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_y_xxyzz_yyz = cbuffer.data(hf_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_y_xxyzz_yzz = cbuffer.data(hf_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_y_xxyzz_zzz = cbuffer.data(hf_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxx = cbuffer.data(hf_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxy = cbuffer.data(hf_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxz = cbuffer.data(hf_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_y_xxzzz_xyy = cbuffer.data(hf_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_y_xxzzz_xyz = cbuffer.data(hf_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_y_xxzzz_xzz = cbuffer.data(hf_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_y_xxzzz_yyy = cbuffer.data(hf_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_y_xxzzz_yyz = cbuffer.data(hf_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_y_xxzzz_yzz = cbuffer.data(hf_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_y_xxzzz_zzz = cbuffer.data(hf_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxx = cbuffer.data(hf_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxy = cbuffer.data(hf_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxz = cbuffer.data(hf_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_y_xyyyy_xyy = cbuffer.data(hf_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_y_xyyyy_xyz = cbuffer.data(hf_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_y_xyyyy_xzz = cbuffer.data(hf_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_y_xyyyy_yyy = cbuffer.data(hf_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_y_xyyyy_yyz = cbuffer.data(hf_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_y_xyyyy_yzz = cbuffer.data(hf_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_y_xyyyy_zzz = cbuffer.data(hf_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxx = cbuffer.data(hf_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxy = cbuffer.data(hf_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxz = cbuffer.data(hf_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_y_xyyyz_xyy = cbuffer.data(hf_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_y_xyyyz_xyz = cbuffer.data(hf_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_y_xyyyz_xzz = cbuffer.data(hf_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_y_xyyyz_yyy = cbuffer.data(hf_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_y_xyyyz_yyz = cbuffer.data(hf_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_y_xyyyz_yzz = cbuffer.data(hf_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_y_xyyyz_zzz = cbuffer.data(hf_geom_01_off + 329 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxx = cbuffer.data(hf_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxy = cbuffer.data(hf_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxz = cbuffer.data(hf_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_y_xyyzz_xyy = cbuffer.data(hf_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_y_xyyzz_xyz = cbuffer.data(hf_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_y_xyyzz_xzz = cbuffer.data(hf_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_y_xyyzz_yyy = cbuffer.data(hf_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_y_xyyzz_yyz = cbuffer.data(hf_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_y_xyyzz_yzz = cbuffer.data(hf_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_y_xyyzz_zzz = cbuffer.data(hf_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxx = cbuffer.data(hf_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxy = cbuffer.data(hf_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxz = cbuffer.data(hf_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_y_xyzzz_xyy = cbuffer.data(hf_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_y_xyzzz_xyz = cbuffer.data(hf_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_y_xyzzz_xzz = cbuffer.data(hf_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_y_xyzzz_yyy = cbuffer.data(hf_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_y_xyzzz_yyz = cbuffer.data(hf_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_y_xyzzz_yzz = cbuffer.data(hf_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_y_xyzzz_zzz = cbuffer.data(hf_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxx = cbuffer.data(hf_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxy = cbuffer.data(hf_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxz = cbuffer.data(hf_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_y_xzzzz_xyy = cbuffer.data(hf_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_y_xzzzz_xyz = cbuffer.data(hf_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_y_xzzzz_xzz = cbuffer.data(hf_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_y_xzzzz_yyy = cbuffer.data(hf_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_y_xzzzz_yyz = cbuffer.data(hf_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_y_xzzzz_yzz = cbuffer.data(hf_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_y_xzzzz_zzz = cbuffer.data(hf_geom_01_off + 359 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxx = cbuffer.data(hf_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxy = cbuffer.data(hf_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxz = cbuffer.data(hf_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_y_yyyyy_xyy = cbuffer.data(hf_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_y_yyyyy_xyz = cbuffer.data(hf_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_y_yyyyy_xzz = cbuffer.data(hf_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_y_yyyyy_yyy = cbuffer.data(hf_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_y_yyyyy_yyz = cbuffer.data(hf_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_y_yyyyy_yzz = cbuffer.data(hf_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_y_yyyyy_zzz = cbuffer.data(hf_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxx = cbuffer.data(hf_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxy = cbuffer.data(hf_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxz = cbuffer.data(hf_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_y_yyyyz_xyy = cbuffer.data(hf_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_y_yyyyz_xyz = cbuffer.data(hf_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_y_yyyyz_xzz = cbuffer.data(hf_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_y_yyyyz_yyy = cbuffer.data(hf_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_y_yyyyz_yyz = cbuffer.data(hf_geom_01_off + 377 * ccomps * dcomps);

            auto g_0_y_yyyyz_yzz = cbuffer.data(hf_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_y_yyyyz_zzz = cbuffer.data(hf_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxx = cbuffer.data(hf_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxy = cbuffer.data(hf_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxz = cbuffer.data(hf_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_y_yyyzz_xyy = cbuffer.data(hf_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_y_yyyzz_xyz = cbuffer.data(hf_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_y_yyyzz_xzz = cbuffer.data(hf_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_y_yyyzz_yyy = cbuffer.data(hf_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_y_yyyzz_yyz = cbuffer.data(hf_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_y_yyyzz_yzz = cbuffer.data(hf_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_y_yyyzz_zzz = cbuffer.data(hf_geom_01_off + 389 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxx = cbuffer.data(hf_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxy = cbuffer.data(hf_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxz = cbuffer.data(hf_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_y_yyzzz_xyy = cbuffer.data(hf_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_y_yyzzz_xyz = cbuffer.data(hf_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_y_yyzzz_xzz = cbuffer.data(hf_geom_01_off + 395 * ccomps * dcomps);

            auto g_0_y_yyzzz_yyy = cbuffer.data(hf_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_y_yyzzz_yyz = cbuffer.data(hf_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_y_yyzzz_yzz = cbuffer.data(hf_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_y_yyzzz_zzz = cbuffer.data(hf_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxx = cbuffer.data(hf_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxy = cbuffer.data(hf_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxz = cbuffer.data(hf_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_y_yzzzz_xyy = cbuffer.data(hf_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_y_yzzzz_xyz = cbuffer.data(hf_geom_01_off + 404 * ccomps * dcomps);

            auto g_0_y_yzzzz_xzz = cbuffer.data(hf_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_y_yzzzz_yyy = cbuffer.data(hf_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_y_yzzzz_yyz = cbuffer.data(hf_geom_01_off + 407 * ccomps * dcomps);

            auto g_0_y_yzzzz_yzz = cbuffer.data(hf_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_y_yzzzz_zzz = cbuffer.data(hf_geom_01_off + 409 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxx = cbuffer.data(hf_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxy = cbuffer.data(hf_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxz = cbuffer.data(hf_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_y_zzzzz_xyy = cbuffer.data(hf_geom_01_off + 413 * ccomps * dcomps);

            auto g_0_y_zzzzz_xyz = cbuffer.data(hf_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_y_zzzzz_xzz = cbuffer.data(hf_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_y_zzzzz_yyy = cbuffer.data(hf_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_y_zzzzz_yyz = cbuffer.data(hf_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_y_zzzzz_yzz = cbuffer.data(hf_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_y_zzzzz_zzz = cbuffer.data(hf_geom_01_off + 419 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxx = cbuffer.data(hf_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxy = cbuffer.data(hf_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxz = cbuffer.data(hf_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_z_xxxxx_xyy = cbuffer.data(hf_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_z_xxxxx_xyz = cbuffer.data(hf_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_z_xxxxx_xzz = cbuffer.data(hf_geom_01_off + 425 * ccomps * dcomps);

            auto g_0_z_xxxxx_yyy = cbuffer.data(hf_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_z_xxxxx_yyz = cbuffer.data(hf_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_z_xxxxx_yzz = cbuffer.data(hf_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_z_xxxxx_zzz = cbuffer.data(hf_geom_01_off + 429 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxx = cbuffer.data(hf_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxy = cbuffer.data(hf_geom_01_off + 431 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxz = cbuffer.data(hf_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_z_xxxxy_xyy = cbuffer.data(hf_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_z_xxxxy_xyz = cbuffer.data(hf_geom_01_off + 434 * ccomps * dcomps);

            auto g_0_z_xxxxy_xzz = cbuffer.data(hf_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_z_xxxxy_yyy = cbuffer.data(hf_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_z_xxxxy_yyz = cbuffer.data(hf_geom_01_off + 437 * ccomps * dcomps);

            auto g_0_z_xxxxy_yzz = cbuffer.data(hf_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_z_xxxxy_zzz = cbuffer.data(hf_geom_01_off + 439 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxx = cbuffer.data(hf_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxy = cbuffer.data(hf_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxz = cbuffer.data(hf_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_z_xxxxz_xyy = cbuffer.data(hf_geom_01_off + 443 * ccomps * dcomps);

            auto g_0_z_xxxxz_xyz = cbuffer.data(hf_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_z_xxxxz_xzz = cbuffer.data(hf_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_z_xxxxz_yyy = cbuffer.data(hf_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_z_xxxxz_yyz = cbuffer.data(hf_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_z_xxxxz_yzz = cbuffer.data(hf_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_z_xxxxz_zzz = cbuffer.data(hf_geom_01_off + 449 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxx = cbuffer.data(hf_geom_01_off + 450 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxy = cbuffer.data(hf_geom_01_off + 451 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxz = cbuffer.data(hf_geom_01_off + 452 * ccomps * dcomps);

            auto g_0_z_xxxyy_xyy = cbuffer.data(hf_geom_01_off + 453 * ccomps * dcomps);

            auto g_0_z_xxxyy_xyz = cbuffer.data(hf_geom_01_off + 454 * ccomps * dcomps);

            auto g_0_z_xxxyy_xzz = cbuffer.data(hf_geom_01_off + 455 * ccomps * dcomps);

            auto g_0_z_xxxyy_yyy = cbuffer.data(hf_geom_01_off + 456 * ccomps * dcomps);

            auto g_0_z_xxxyy_yyz = cbuffer.data(hf_geom_01_off + 457 * ccomps * dcomps);

            auto g_0_z_xxxyy_yzz = cbuffer.data(hf_geom_01_off + 458 * ccomps * dcomps);

            auto g_0_z_xxxyy_zzz = cbuffer.data(hf_geom_01_off + 459 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxx = cbuffer.data(hf_geom_01_off + 460 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxy = cbuffer.data(hf_geom_01_off + 461 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxz = cbuffer.data(hf_geom_01_off + 462 * ccomps * dcomps);

            auto g_0_z_xxxyz_xyy = cbuffer.data(hf_geom_01_off + 463 * ccomps * dcomps);

            auto g_0_z_xxxyz_xyz = cbuffer.data(hf_geom_01_off + 464 * ccomps * dcomps);

            auto g_0_z_xxxyz_xzz = cbuffer.data(hf_geom_01_off + 465 * ccomps * dcomps);

            auto g_0_z_xxxyz_yyy = cbuffer.data(hf_geom_01_off + 466 * ccomps * dcomps);

            auto g_0_z_xxxyz_yyz = cbuffer.data(hf_geom_01_off + 467 * ccomps * dcomps);

            auto g_0_z_xxxyz_yzz = cbuffer.data(hf_geom_01_off + 468 * ccomps * dcomps);

            auto g_0_z_xxxyz_zzz = cbuffer.data(hf_geom_01_off + 469 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxx = cbuffer.data(hf_geom_01_off + 470 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxy = cbuffer.data(hf_geom_01_off + 471 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxz = cbuffer.data(hf_geom_01_off + 472 * ccomps * dcomps);

            auto g_0_z_xxxzz_xyy = cbuffer.data(hf_geom_01_off + 473 * ccomps * dcomps);

            auto g_0_z_xxxzz_xyz = cbuffer.data(hf_geom_01_off + 474 * ccomps * dcomps);

            auto g_0_z_xxxzz_xzz = cbuffer.data(hf_geom_01_off + 475 * ccomps * dcomps);

            auto g_0_z_xxxzz_yyy = cbuffer.data(hf_geom_01_off + 476 * ccomps * dcomps);

            auto g_0_z_xxxzz_yyz = cbuffer.data(hf_geom_01_off + 477 * ccomps * dcomps);

            auto g_0_z_xxxzz_yzz = cbuffer.data(hf_geom_01_off + 478 * ccomps * dcomps);

            auto g_0_z_xxxzz_zzz = cbuffer.data(hf_geom_01_off + 479 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxx = cbuffer.data(hf_geom_01_off + 480 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxy = cbuffer.data(hf_geom_01_off + 481 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxz = cbuffer.data(hf_geom_01_off + 482 * ccomps * dcomps);

            auto g_0_z_xxyyy_xyy = cbuffer.data(hf_geom_01_off + 483 * ccomps * dcomps);

            auto g_0_z_xxyyy_xyz = cbuffer.data(hf_geom_01_off + 484 * ccomps * dcomps);

            auto g_0_z_xxyyy_xzz = cbuffer.data(hf_geom_01_off + 485 * ccomps * dcomps);

            auto g_0_z_xxyyy_yyy = cbuffer.data(hf_geom_01_off + 486 * ccomps * dcomps);

            auto g_0_z_xxyyy_yyz = cbuffer.data(hf_geom_01_off + 487 * ccomps * dcomps);

            auto g_0_z_xxyyy_yzz = cbuffer.data(hf_geom_01_off + 488 * ccomps * dcomps);

            auto g_0_z_xxyyy_zzz = cbuffer.data(hf_geom_01_off + 489 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxx = cbuffer.data(hf_geom_01_off + 490 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxy = cbuffer.data(hf_geom_01_off + 491 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxz = cbuffer.data(hf_geom_01_off + 492 * ccomps * dcomps);

            auto g_0_z_xxyyz_xyy = cbuffer.data(hf_geom_01_off + 493 * ccomps * dcomps);

            auto g_0_z_xxyyz_xyz = cbuffer.data(hf_geom_01_off + 494 * ccomps * dcomps);

            auto g_0_z_xxyyz_xzz = cbuffer.data(hf_geom_01_off + 495 * ccomps * dcomps);

            auto g_0_z_xxyyz_yyy = cbuffer.data(hf_geom_01_off + 496 * ccomps * dcomps);

            auto g_0_z_xxyyz_yyz = cbuffer.data(hf_geom_01_off + 497 * ccomps * dcomps);

            auto g_0_z_xxyyz_yzz = cbuffer.data(hf_geom_01_off + 498 * ccomps * dcomps);

            auto g_0_z_xxyyz_zzz = cbuffer.data(hf_geom_01_off + 499 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxx = cbuffer.data(hf_geom_01_off + 500 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxy = cbuffer.data(hf_geom_01_off + 501 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxz = cbuffer.data(hf_geom_01_off + 502 * ccomps * dcomps);

            auto g_0_z_xxyzz_xyy = cbuffer.data(hf_geom_01_off + 503 * ccomps * dcomps);

            auto g_0_z_xxyzz_xyz = cbuffer.data(hf_geom_01_off + 504 * ccomps * dcomps);

            auto g_0_z_xxyzz_xzz = cbuffer.data(hf_geom_01_off + 505 * ccomps * dcomps);

            auto g_0_z_xxyzz_yyy = cbuffer.data(hf_geom_01_off + 506 * ccomps * dcomps);

            auto g_0_z_xxyzz_yyz = cbuffer.data(hf_geom_01_off + 507 * ccomps * dcomps);

            auto g_0_z_xxyzz_yzz = cbuffer.data(hf_geom_01_off + 508 * ccomps * dcomps);

            auto g_0_z_xxyzz_zzz = cbuffer.data(hf_geom_01_off + 509 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxx = cbuffer.data(hf_geom_01_off + 510 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxy = cbuffer.data(hf_geom_01_off + 511 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxz = cbuffer.data(hf_geom_01_off + 512 * ccomps * dcomps);

            auto g_0_z_xxzzz_xyy = cbuffer.data(hf_geom_01_off + 513 * ccomps * dcomps);

            auto g_0_z_xxzzz_xyz = cbuffer.data(hf_geom_01_off + 514 * ccomps * dcomps);

            auto g_0_z_xxzzz_xzz = cbuffer.data(hf_geom_01_off + 515 * ccomps * dcomps);

            auto g_0_z_xxzzz_yyy = cbuffer.data(hf_geom_01_off + 516 * ccomps * dcomps);

            auto g_0_z_xxzzz_yyz = cbuffer.data(hf_geom_01_off + 517 * ccomps * dcomps);

            auto g_0_z_xxzzz_yzz = cbuffer.data(hf_geom_01_off + 518 * ccomps * dcomps);

            auto g_0_z_xxzzz_zzz = cbuffer.data(hf_geom_01_off + 519 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxx = cbuffer.data(hf_geom_01_off + 520 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxy = cbuffer.data(hf_geom_01_off + 521 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxz = cbuffer.data(hf_geom_01_off + 522 * ccomps * dcomps);

            auto g_0_z_xyyyy_xyy = cbuffer.data(hf_geom_01_off + 523 * ccomps * dcomps);

            auto g_0_z_xyyyy_xyz = cbuffer.data(hf_geom_01_off + 524 * ccomps * dcomps);

            auto g_0_z_xyyyy_xzz = cbuffer.data(hf_geom_01_off + 525 * ccomps * dcomps);

            auto g_0_z_xyyyy_yyy = cbuffer.data(hf_geom_01_off + 526 * ccomps * dcomps);

            auto g_0_z_xyyyy_yyz = cbuffer.data(hf_geom_01_off + 527 * ccomps * dcomps);

            auto g_0_z_xyyyy_yzz = cbuffer.data(hf_geom_01_off + 528 * ccomps * dcomps);

            auto g_0_z_xyyyy_zzz = cbuffer.data(hf_geom_01_off + 529 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxx = cbuffer.data(hf_geom_01_off + 530 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxy = cbuffer.data(hf_geom_01_off + 531 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxz = cbuffer.data(hf_geom_01_off + 532 * ccomps * dcomps);

            auto g_0_z_xyyyz_xyy = cbuffer.data(hf_geom_01_off + 533 * ccomps * dcomps);

            auto g_0_z_xyyyz_xyz = cbuffer.data(hf_geom_01_off + 534 * ccomps * dcomps);

            auto g_0_z_xyyyz_xzz = cbuffer.data(hf_geom_01_off + 535 * ccomps * dcomps);

            auto g_0_z_xyyyz_yyy = cbuffer.data(hf_geom_01_off + 536 * ccomps * dcomps);

            auto g_0_z_xyyyz_yyz = cbuffer.data(hf_geom_01_off + 537 * ccomps * dcomps);

            auto g_0_z_xyyyz_yzz = cbuffer.data(hf_geom_01_off + 538 * ccomps * dcomps);

            auto g_0_z_xyyyz_zzz = cbuffer.data(hf_geom_01_off + 539 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxx = cbuffer.data(hf_geom_01_off + 540 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxy = cbuffer.data(hf_geom_01_off + 541 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxz = cbuffer.data(hf_geom_01_off + 542 * ccomps * dcomps);

            auto g_0_z_xyyzz_xyy = cbuffer.data(hf_geom_01_off + 543 * ccomps * dcomps);

            auto g_0_z_xyyzz_xyz = cbuffer.data(hf_geom_01_off + 544 * ccomps * dcomps);

            auto g_0_z_xyyzz_xzz = cbuffer.data(hf_geom_01_off + 545 * ccomps * dcomps);

            auto g_0_z_xyyzz_yyy = cbuffer.data(hf_geom_01_off + 546 * ccomps * dcomps);

            auto g_0_z_xyyzz_yyz = cbuffer.data(hf_geom_01_off + 547 * ccomps * dcomps);

            auto g_0_z_xyyzz_yzz = cbuffer.data(hf_geom_01_off + 548 * ccomps * dcomps);

            auto g_0_z_xyyzz_zzz = cbuffer.data(hf_geom_01_off + 549 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxx = cbuffer.data(hf_geom_01_off + 550 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxy = cbuffer.data(hf_geom_01_off + 551 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxz = cbuffer.data(hf_geom_01_off + 552 * ccomps * dcomps);

            auto g_0_z_xyzzz_xyy = cbuffer.data(hf_geom_01_off + 553 * ccomps * dcomps);

            auto g_0_z_xyzzz_xyz = cbuffer.data(hf_geom_01_off + 554 * ccomps * dcomps);

            auto g_0_z_xyzzz_xzz = cbuffer.data(hf_geom_01_off + 555 * ccomps * dcomps);

            auto g_0_z_xyzzz_yyy = cbuffer.data(hf_geom_01_off + 556 * ccomps * dcomps);

            auto g_0_z_xyzzz_yyz = cbuffer.data(hf_geom_01_off + 557 * ccomps * dcomps);

            auto g_0_z_xyzzz_yzz = cbuffer.data(hf_geom_01_off + 558 * ccomps * dcomps);

            auto g_0_z_xyzzz_zzz = cbuffer.data(hf_geom_01_off + 559 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxx = cbuffer.data(hf_geom_01_off + 560 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxy = cbuffer.data(hf_geom_01_off + 561 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxz = cbuffer.data(hf_geom_01_off + 562 * ccomps * dcomps);

            auto g_0_z_xzzzz_xyy = cbuffer.data(hf_geom_01_off + 563 * ccomps * dcomps);

            auto g_0_z_xzzzz_xyz = cbuffer.data(hf_geom_01_off + 564 * ccomps * dcomps);

            auto g_0_z_xzzzz_xzz = cbuffer.data(hf_geom_01_off + 565 * ccomps * dcomps);

            auto g_0_z_xzzzz_yyy = cbuffer.data(hf_geom_01_off + 566 * ccomps * dcomps);

            auto g_0_z_xzzzz_yyz = cbuffer.data(hf_geom_01_off + 567 * ccomps * dcomps);

            auto g_0_z_xzzzz_yzz = cbuffer.data(hf_geom_01_off + 568 * ccomps * dcomps);

            auto g_0_z_xzzzz_zzz = cbuffer.data(hf_geom_01_off + 569 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxx = cbuffer.data(hf_geom_01_off + 570 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxy = cbuffer.data(hf_geom_01_off + 571 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxz = cbuffer.data(hf_geom_01_off + 572 * ccomps * dcomps);

            auto g_0_z_yyyyy_xyy = cbuffer.data(hf_geom_01_off + 573 * ccomps * dcomps);

            auto g_0_z_yyyyy_xyz = cbuffer.data(hf_geom_01_off + 574 * ccomps * dcomps);

            auto g_0_z_yyyyy_xzz = cbuffer.data(hf_geom_01_off + 575 * ccomps * dcomps);

            auto g_0_z_yyyyy_yyy = cbuffer.data(hf_geom_01_off + 576 * ccomps * dcomps);

            auto g_0_z_yyyyy_yyz = cbuffer.data(hf_geom_01_off + 577 * ccomps * dcomps);

            auto g_0_z_yyyyy_yzz = cbuffer.data(hf_geom_01_off + 578 * ccomps * dcomps);

            auto g_0_z_yyyyy_zzz = cbuffer.data(hf_geom_01_off + 579 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxx = cbuffer.data(hf_geom_01_off + 580 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxy = cbuffer.data(hf_geom_01_off + 581 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxz = cbuffer.data(hf_geom_01_off + 582 * ccomps * dcomps);

            auto g_0_z_yyyyz_xyy = cbuffer.data(hf_geom_01_off + 583 * ccomps * dcomps);

            auto g_0_z_yyyyz_xyz = cbuffer.data(hf_geom_01_off + 584 * ccomps * dcomps);

            auto g_0_z_yyyyz_xzz = cbuffer.data(hf_geom_01_off + 585 * ccomps * dcomps);

            auto g_0_z_yyyyz_yyy = cbuffer.data(hf_geom_01_off + 586 * ccomps * dcomps);

            auto g_0_z_yyyyz_yyz = cbuffer.data(hf_geom_01_off + 587 * ccomps * dcomps);

            auto g_0_z_yyyyz_yzz = cbuffer.data(hf_geom_01_off + 588 * ccomps * dcomps);

            auto g_0_z_yyyyz_zzz = cbuffer.data(hf_geom_01_off + 589 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxx = cbuffer.data(hf_geom_01_off + 590 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxy = cbuffer.data(hf_geom_01_off + 591 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxz = cbuffer.data(hf_geom_01_off + 592 * ccomps * dcomps);

            auto g_0_z_yyyzz_xyy = cbuffer.data(hf_geom_01_off + 593 * ccomps * dcomps);

            auto g_0_z_yyyzz_xyz = cbuffer.data(hf_geom_01_off + 594 * ccomps * dcomps);

            auto g_0_z_yyyzz_xzz = cbuffer.data(hf_geom_01_off + 595 * ccomps * dcomps);

            auto g_0_z_yyyzz_yyy = cbuffer.data(hf_geom_01_off + 596 * ccomps * dcomps);

            auto g_0_z_yyyzz_yyz = cbuffer.data(hf_geom_01_off + 597 * ccomps * dcomps);

            auto g_0_z_yyyzz_yzz = cbuffer.data(hf_geom_01_off + 598 * ccomps * dcomps);

            auto g_0_z_yyyzz_zzz = cbuffer.data(hf_geom_01_off + 599 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxx = cbuffer.data(hf_geom_01_off + 600 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxy = cbuffer.data(hf_geom_01_off + 601 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxz = cbuffer.data(hf_geom_01_off + 602 * ccomps * dcomps);

            auto g_0_z_yyzzz_xyy = cbuffer.data(hf_geom_01_off + 603 * ccomps * dcomps);

            auto g_0_z_yyzzz_xyz = cbuffer.data(hf_geom_01_off + 604 * ccomps * dcomps);

            auto g_0_z_yyzzz_xzz = cbuffer.data(hf_geom_01_off + 605 * ccomps * dcomps);

            auto g_0_z_yyzzz_yyy = cbuffer.data(hf_geom_01_off + 606 * ccomps * dcomps);

            auto g_0_z_yyzzz_yyz = cbuffer.data(hf_geom_01_off + 607 * ccomps * dcomps);

            auto g_0_z_yyzzz_yzz = cbuffer.data(hf_geom_01_off + 608 * ccomps * dcomps);

            auto g_0_z_yyzzz_zzz = cbuffer.data(hf_geom_01_off + 609 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxx = cbuffer.data(hf_geom_01_off + 610 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxy = cbuffer.data(hf_geom_01_off + 611 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxz = cbuffer.data(hf_geom_01_off + 612 * ccomps * dcomps);

            auto g_0_z_yzzzz_xyy = cbuffer.data(hf_geom_01_off + 613 * ccomps * dcomps);

            auto g_0_z_yzzzz_xyz = cbuffer.data(hf_geom_01_off + 614 * ccomps * dcomps);

            auto g_0_z_yzzzz_xzz = cbuffer.data(hf_geom_01_off + 615 * ccomps * dcomps);

            auto g_0_z_yzzzz_yyy = cbuffer.data(hf_geom_01_off + 616 * ccomps * dcomps);

            auto g_0_z_yzzzz_yyz = cbuffer.data(hf_geom_01_off + 617 * ccomps * dcomps);

            auto g_0_z_yzzzz_yzz = cbuffer.data(hf_geom_01_off + 618 * ccomps * dcomps);

            auto g_0_z_yzzzz_zzz = cbuffer.data(hf_geom_01_off + 619 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxx = cbuffer.data(hf_geom_01_off + 620 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxy = cbuffer.data(hf_geom_01_off + 621 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxz = cbuffer.data(hf_geom_01_off + 622 * ccomps * dcomps);

            auto g_0_z_zzzzz_xyy = cbuffer.data(hf_geom_01_off + 623 * ccomps * dcomps);

            auto g_0_z_zzzzz_xyz = cbuffer.data(hf_geom_01_off + 624 * ccomps * dcomps);

            auto g_0_z_zzzzz_xzz = cbuffer.data(hf_geom_01_off + 625 * ccomps * dcomps);

            auto g_0_z_zzzzz_yyy = cbuffer.data(hf_geom_01_off + 626 * ccomps * dcomps);

            auto g_0_z_zzzzz_yyz = cbuffer.data(hf_geom_01_off + 627 * ccomps * dcomps);

            auto g_0_z_zzzzz_yzz = cbuffer.data(hf_geom_01_off + 628 * ccomps * dcomps);

            auto g_0_z_zzzzz_zzz = cbuffer.data(hf_geom_01_off + 629 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_idxx

            const auto id_geom_01_off = idx_geom_01_idxx + i * dcomps + j;

            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxx_xx = cbuffer.data(id_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xy = cbuffer.data(id_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xz = cbuffer.data(id_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxxxx_yy = cbuffer.data(id_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxxxx_yz = cbuffer.data(id_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxxxx_zz = cbuffer.data(id_geom_01_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxx_xx, g_0_x_xxxxx_xxx, g_0_x_xxxxx_xxy, g_0_x_xxxxx_xxz, g_0_x_xxxxx_xy, g_0_x_xxxxx_xyy, g_0_x_xxxxx_xyz, g_0_x_xxxxx_xz, g_0_x_xxxxx_xzz, g_0_x_xxxxx_yy, g_0_x_xxxxx_yz, g_0_x_xxxxx_zz, g_0_x_xxxxxx_xx, g_0_x_xxxxxx_xy, g_0_x_xxxxxx_xz, g_0_x_xxxxxx_yy, g_0_x_xxxxxx_yz, g_0_x_xxxxxx_zz, g_xxxxx_xx, g_xxxxx_xy, g_xxxxx_xz, g_xxxxx_yy, g_xxxxx_yz, g_xxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxx_xx[k] = g_xxxxx_xx[k] - g_0_x_xxxxx_xx[k] * ab_x + g_0_x_xxxxx_xxx[k];

                g_0_x_xxxxxx_xy[k] = g_xxxxx_xy[k] - g_0_x_xxxxx_xy[k] * ab_x + g_0_x_xxxxx_xxy[k];

                g_0_x_xxxxxx_xz[k] = g_xxxxx_xz[k] - g_0_x_xxxxx_xz[k] * ab_x + g_0_x_xxxxx_xxz[k];

                g_0_x_xxxxxx_yy[k] = g_xxxxx_yy[k] - g_0_x_xxxxx_yy[k] * ab_x + g_0_x_xxxxx_xyy[k];

                g_0_x_xxxxxx_yz[k] = g_xxxxx_yz[k] - g_0_x_xxxxx_yz[k] * ab_x + g_0_x_xxxxx_xyz[k];

                g_0_x_xxxxxx_zz[k] = g_xxxxx_zz[k] - g_0_x_xxxxx_zz[k] * ab_x + g_0_x_xxxxx_xzz[k];
            }

            /// Set up 6-12 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxy_xx = cbuffer.data(id_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xy = cbuffer.data(id_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xz = cbuffer.data(id_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxxxxy_yy = cbuffer.data(id_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxxxxy_yz = cbuffer.data(id_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxxxy_zz = cbuffer.data(id_geom_01_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxx_xx, g_0_x_xxxxx_xxy, g_0_x_xxxxx_xy, g_0_x_xxxxx_xyy, g_0_x_xxxxx_xyz, g_0_x_xxxxx_xz, g_0_x_xxxxx_yy, g_0_x_xxxxx_yyy, g_0_x_xxxxx_yyz, g_0_x_xxxxx_yz, g_0_x_xxxxx_yzz, g_0_x_xxxxx_zz, g_0_x_xxxxxy_xx, g_0_x_xxxxxy_xy, g_0_x_xxxxxy_xz, g_0_x_xxxxxy_yy, g_0_x_xxxxxy_yz, g_0_x_xxxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxy_xx[k] = -g_0_x_xxxxx_xx[k] * ab_y + g_0_x_xxxxx_xxy[k];

                g_0_x_xxxxxy_xy[k] = -g_0_x_xxxxx_xy[k] * ab_y + g_0_x_xxxxx_xyy[k];

                g_0_x_xxxxxy_xz[k] = -g_0_x_xxxxx_xz[k] * ab_y + g_0_x_xxxxx_xyz[k];

                g_0_x_xxxxxy_yy[k] = -g_0_x_xxxxx_yy[k] * ab_y + g_0_x_xxxxx_yyy[k];

                g_0_x_xxxxxy_yz[k] = -g_0_x_xxxxx_yz[k] * ab_y + g_0_x_xxxxx_yyz[k];

                g_0_x_xxxxxy_zz[k] = -g_0_x_xxxxx_zz[k] * ab_y + g_0_x_xxxxx_yzz[k];
            }

            /// Set up 12-18 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxz_xx = cbuffer.data(id_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xy = cbuffer.data(id_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xz = cbuffer.data(id_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxxxxz_yy = cbuffer.data(id_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxxxz_yz = cbuffer.data(id_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxxxz_zz = cbuffer.data(id_geom_01_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxx_xx, g_0_x_xxxxx_xxz, g_0_x_xxxxx_xy, g_0_x_xxxxx_xyz, g_0_x_xxxxx_xz, g_0_x_xxxxx_xzz, g_0_x_xxxxx_yy, g_0_x_xxxxx_yyz, g_0_x_xxxxx_yz, g_0_x_xxxxx_yzz, g_0_x_xxxxx_zz, g_0_x_xxxxx_zzz, g_0_x_xxxxxz_xx, g_0_x_xxxxxz_xy, g_0_x_xxxxxz_xz, g_0_x_xxxxxz_yy, g_0_x_xxxxxz_yz, g_0_x_xxxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxz_xx[k] = -g_0_x_xxxxx_xx[k] * ab_z + g_0_x_xxxxx_xxz[k];

                g_0_x_xxxxxz_xy[k] = -g_0_x_xxxxx_xy[k] * ab_z + g_0_x_xxxxx_xyz[k];

                g_0_x_xxxxxz_xz[k] = -g_0_x_xxxxx_xz[k] * ab_z + g_0_x_xxxxx_xzz[k];

                g_0_x_xxxxxz_yy[k] = -g_0_x_xxxxx_yy[k] * ab_z + g_0_x_xxxxx_yyz[k];

                g_0_x_xxxxxz_yz[k] = -g_0_x_xxxxx_yz[k] * ab_z + g_0_x_xxxxx_yzz[k];

                g_0_x_xxxxxz_zz[k] = -g_0_x_xxxxx_zz[k] * ab_z + g_0_x_xxxxx_zzz[k];
            }

            /// Set up 18-24 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxyy_xx = cbuffer.data(id_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xy = cbuffer.data(id_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xz = cbuffer.data(id_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxxxyy_yy = cbuffer.data(id_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxxxyy_yz = cbuffer.data(id_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxxxyy_zz = cbuffer.data(id_geom_01_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxy_xx, g_0_x_xxxxy_xxy, g_0_x_xxxxy_xy, g_0_x_xxxxy_xyy, g_0_x_xxxxy_xyz, g_0_x_xxxxy_xz, g_0_x_xxxxy_yy, g_0_x_xxxxy_yyy, g_0_x_xxxxy_yyz, g_0_x_xxxxy_yz, g_0_x_xxxxy_yzz, g_0_x_xxxxy_zz, g_0_x_xxxxyy_xx, g_0_x_xxxxyy_xy, g_0_x_xxxxyy_xz, g_0_x_xxxxyy_yy, g_0_x_xxxxyy_yz, g_0_x_xxxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxyy_xx[k] = -g_0_x_xxxxy_xx[k] * ab_y + g_0_x_xxxxy_xxy[k];

                g_0_x_xxxxyy_xy[k] = -g_0_x_xxxxy_xy[k] * ab_y + g_0_x_xxxxy_xyy[k];

                g_0_x_xxxxyy_xz[k] = -g_0_x_xxxxy_xz[k] * ab_y + g_0_x_xxxxy_xyz[k];

                g_0_x_xxxxyy_yy[k] = -g_0_x_xxxxy_yy[k] * ab_y + g_0_x_xxxxy_yyy[k];

                g_0_x_xxxxyy_yz[k] = -g_0_x_xxxxy_yz[k] * ab_y + g_0_x_xxxxy_yyz[k];

                g_0_x_xxxxyy_zz[k] = -g_0_x_xxxxy_zz[k] * ab_y + g_0_x_xxxxy_yzz[k];
            }

            /// Set up 24-30 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxyz_xx = cbuffer.data(id_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xy = cbuffer.data(id_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xz = cbuffer.data(id_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxxxyz_yy = cbuffer.data(id_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxxxyz_yz = cbuffer.data(id_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxxxyz_zz = cbuffer.data(id_geom_01_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxyz_xx, g_0_x_xxxxyz_xy, g_0_x_xxxxyz_xz, g_0_x_xxxxyz_yy, g_0_x_xxxxyz_yz, g_0_x_xxxxyz_zz, g_0_x_xxxxz_xx, g_0_x_xxxxz_xxy, g_0_x_xxxxz_xy, g_0_x_xxxxz_xyy, g_0_x_xxxxz_xyz, g_0_x_xxxxz_xz, g_0_x_xxxxz_yy, g_0_x_xxxxz_yyy, g_0_x_xxxxz_yyz, g_0_x_xxxxz_yz, g_0_x_xxxxz_yzz, g_0_x_xxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxyz_xx[k] = -g_0_x_xxxxz_xx[k] * ab_y + g_0_x_xxxxz_xxy[k];

                g_0_x_xxxxyz_xy[k] = -g_0_x_xxxxz_xy[k] * ab_y + g_0_x_xxxxz_xyy[k];

                g_0_x_xxxxyz_xz[k] = -g_0_x_xxxxz_xz[k] * ab_y + g_0_x_xxxxz_xyz[k];

                g_0_x_xxxxyz_yy[k] = -g_0_x_xxxxz_yy[k] * ab_y + g_0_x_xxxxz_yyy[k];

                g_0_x_xxxxyz_yz[k] = -g_0_x_xxxxz_yz[k] * ab_y + g_0_x_xxxxz_yyz[k];

                g_0_x_xxxxyz_zz[k] = -g_0_x_xxxxz_zz[k] * ab_y + g_0_x_xxxxz_yzz[k];
            }

            /// Set up 30-36 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxzz_xx = cbuffer.data(id_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xy = cbuffer.data(id_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xz = cbuffer.data(id_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxxxzz_yy = cbuffer.data(id_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxxxzz_yz = cbuffer.data(id_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxxxzz_zz = cbuffer.data(id_geom_01_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxz_xx, g_0_x_xxxxz_xxz, g_0_x_xxxxz_xy, g_0_x_xxxxz_xyz, g_0_x_xxxxz_xz, g_0_x_xxxxz_xzz, g_0_x_xxxxz_yy, g_0_x_xxxxz_yyz, g_0_x_xxxxz_yz, g_0_x_xxxxz_yzz, g_0_x_xxxxz_zz, g_0_x_xxxxz_zzz, g_0_x_xxxxzz_xx, g_0_x_xxxxzz_xy, g_0_x_xxxxzz_xz, g_0_x_xxxxzz_yy, g_0_x_xxxxzz_yz, g_0_x_xxxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxzz_xx[k] = -g_0_x_xxxxz_xx[k] * ab_z + g_0_x_xxxxz_xxz[k];

                g_0_x_xxxxzz_xy[k] = -g_0_x_xxxxz_xy[k] * ab_z + g_0_x_xxxxz_xyz[k];

                g_0_x_xxxxzz_xz[k] = -g_0_x_xxxxz_xz[k] * ab_z + g_0_x_xxxxz_xzz[k];

                g_0_x_xxxxzz_yy[k] = -g_0_x_xxxxz_yy[k] * ab_z + g_0_x_xxxxz_yyz[k];

                g_0_x_xxxxzz_yz[k] = -g_0_x_xxxxz_yz[k] * ab_z + g_0_x_xxxxz_yzz[k];

                g_0_x_xxxxzz_zz[k] = -g_0_x_xxxxz_zz[k] * ab_z + g_0_x_xxxxz_zzz[k];
            }

            /// Set up 36-42 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyyy_xx = cbuffer.data(id_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xy = cbuffer.data(id_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xz = cbuffer.data(id_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxxyyy_yy = cbuffer.data(id_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxxyyy_yz = cbuffer.data(id_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxxyyy_zz = cbuffer.data(id_geom_01_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyy_xx, g_0_x_xxxyy_xxy, g_0_x_xxxyy_xy, g_0_x_xxxyy_xyy, g_0_x_xxxyy_xyz, g_0_x_xxxyy_xz, g_0_x_xxxyy_yy, g_0_x_xxxyy_yyy, g_0_x_xxxyy_yyz, g_0_x_xxxyy_yz, g_0_x_xxxyy_yzz, g_0_x_xxxyy_zz, g_0_x_xxxyyy_xx, g_0_x_xxxyyy_xy, g_0_x_xxxyyy_xz, g_0_x_xxxyyy_yy, g_0_x_xxxyyy_yz, g_0_x_xxxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyyy_xx[k] = -g_0_x_xxxyy_xx[k] * ab_y + g_0_x_xxxyy_xxy[k];

                g_0_x_xxxyyy_xy[k] = -g_0_x_xxxyy_xy[k] * ab_y + g_0_x_xxxyy_xyy[k];

                g_0_x_xxxyyy_xz[k] = -g_0_x_xxxyy_xz[k] * ab_y + g_0_x_xxxyy_xyz[k];

                g_0_x_xxxyyy_yy[k] = -g_0_x_xxxyy_yy[k] * ab_y + g_0_x_xxxyy_yyy[k];

                g_0_x_xxxyyy_yz[k] = -g_0_x_xxxyy_yz[k] * ab_y + g_0_x_xxxyy_yyz[k];

                g_0_x_xxxyyy_zz[k] = -g_0_x_xxxyy_zz[k] * ab_y + g_0_x_xxxyy_yzz[k];
            }

            /// Set up 42-48 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyyz_xx = cbuffer.data(id_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xy = cbuffer.data(id_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xz = cbuffer.data(id_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xxxyyz_yy = cbuffer.data(id_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xxxyyz_yz = cbuffer.data(id_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xxxyyz_zz = cbuffer.data(id_geom_01_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyyz_xx, g_0_x_xxxyyz_xy, g_0_x_xxxyyz_xz, g_0_x_xxxyyz_yy, g_0_x_xxxyyz_yz, g_0_x_xxxyyz_zz, g_0_x_xxxyz_xx, g_0_x_xxxyz_xxy, g_0_x_xxxyz_xy, g_0_x_xxxyz_xyy, g_0_x_xxxyz_xyz, g_0_x_xxxyz_xz, g_0_x_xxxyz_yy, g_0_x_xxxyz_yyy, g_0_x_xxxyz_yyz, g_0_x_xxxyz_yz, g_0_x_xxxyz_yzz, g_0_x_xxxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyyz_xx[k] = -g_0_x_xxxyz_xx[k] * ab_y + g_0_x_xxxyz_xxy[k];

                g_0_x_xxxyyz_xy[k] = -g_0_x_xxxyz_xy[k] * ab_y + g_0_x_xxxyz_xyy[k];

                g_0_x_xxxyyz_xz[k] = -g_0_x_xxxyz_xz[k] * ab_y + g_0_x_xxxyz_xyz[k];

                g_0_x_xxxyyz_yy[k] = -g_0_x_xxxyz_yy[k] * ab_y + g_0_x_xxxyz_yyy[k];

                g_0_x_xxxyyz_yz[k] = -g_0_x_xxxyz_yz[k] * ab_y + g_0_x_xxxyz_yyz[k];

                g_0_x_xxxyyz_zz[k] = -g_0_x_xxxyz_zz[k] * ab_y + g_0_x_xxxyz_yzz[k];
            }

            /// Set up 48-54 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyzz_xx = cbuffer.data(id_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xy = cbuffer.data(id_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xz = cbuffer.data(id_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xxxyzz_yy = cbuffer.data(id_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xxxyzz_yz = cbuffer.data(id_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xxxyzz_zz = cbuffer.data(id_geom_01_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyzz_xx, g_0_x_xxxyzz_xy, g_0_x_xxxyzz_xz, g_0_x_xxxyzz_yy, g_0_x_xxxyzz_yz, g_0_x_xxxyzz_zz, g_0_x_xxxzz_xx, g_0_x_xxxzz_xxy, g_0_x_xxxzz_xy, g_0_x_xxxzz_xyy, g_0_x_xxxzz_xyz, g_0_x_xxxzz_xz, g_0_x_xxxzz_yy, g_0_x_xxxzz_yyy, g_0_x_xxxzz_yyz, g_0_x_xxxzz_yz, g_0_x_xxxzz_yzz, g_0_x_xxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyzz_xx[k] = -g_0_x_xxxzz_xx[k] * ab_y + g_0_x_xxxzz_xxy[k];

                g_0_x_xxxyzz_xy[k] = -g_0_x_xxxzz_xy[k] * ab_y + g_0_x_xxxzz_xyy[k];

                g_0_x_xxxyzz_xz[k] = -g_0_x_xxxzz_xz[k] * ab_y + g_0_x_xxxzz_xyz[k];

                g_0_x_xxxyzz_yy[k] = -g_0_x_xxxzz_yy[k] * ab_y + g_0_x_xxxzz_yyy[k];

                g_0_x_xxxyzz_yz[k] = -g_0_x_xxxzz_yz[k] * ab_y + g_0_x_xxxzz_yyz[k];

                g_0_x_xxxyzz_zz[k] = -g_0_x_xxxzz_zz[k] * ab_y + g_0_x_xxxzz_yzz[k];
            }

            /// Set up 54-60 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxzzz_xx = cbuffer.data(id_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xy = cbuffer.data(id_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xz = cbuffer.data(id_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xxxzzz_yy = cbuffer.data(id_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xxxzzz_yz = cbuffer.data(id_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xxxzzz_zz = cbuffer.data(id_geom_01_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxzz_xx, g_0_x_xxxzz_xxz, g_0_x_xxxzz_xy, g_0_x_xxxzz_xyz, g_0_x_xxxzz_xz, g_0_x_xxxzz_xzz, g_0_x_xxxzz_yy, g_0_x_xxxzz_yyz, g_0_x_xxxzz_yz, g_0_x_xxxzz_yzz, g_0_x_xxxzz_zz, g_0_x_xxxzz_zzz, g_0_x_xxxzzz_xx, g_0_x_xxxzzz_xy, g_0_x_xxxzzz_xz, g_0_x_xxxzzz_yy, g_0_x_xxxzzz_yz, g_0_x_xxxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxzzz_xx[k] = -g_0_x_xxxzz_xx[k] * ab_z + g_0_x_xxxzz_xxz[k];

                g_0_x_xxxzzz_xy[k] = -g_0_x_xxxzz_xy[k] * ab_z + g_0_x_xxxzz_xyz[k];

                g_0_x_xxxzzz_xz[k] = -g_0_x_xxxzz_xz[k] * ab_z + g_0_x_xxxzz_xzz[k];

                g_0_x_xxxzzz_yy[k] = -g_0_x_xxxzz_yy[k] * ab_z + g_0_x_xxxzz_yyz[k];

                g_0_x_xxxzzz_yz[k] = -g_0_x_xxxzz_yz[k] * ab_z + g_0_x_xxxzz_yzz[k];

                g_0_x_xxxzzz_zz[k] = -g_0_x_xxxzz_zz[k] * ab_z + g_0_x_xxxzz_zzz[k];
            }

            /// Set up 60-66 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyyy_xx = cbuffer.data(id_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xy = cbuffer.data(id_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xz = cbuffer.data(id_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xxyyyy_yy = cbuffer.data(id_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xxyyyy_yz = cbuffer.data(id_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xxyyyy_zz = cbuffer.data(id_geom_01_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyy_xx, g_0_x_xxyyy_xxy, g_0_x_xxyyy_xy, g_0_x_xxyyy_xyy, g_0_x_xxyyy_xyz, g_0_x_xxyyy_xz, g_0_x_xxyyy_yy, g_0_x_xxyyy_yyy, g_0_x_xxyyy_yyz, g_0_x_xxyyy_yz, g_0_x_xxyyy_yzz, g_0_x_xxyyy_zz, g_0_x_xxyyyy_xx, g_0_x_xxyyyy_xy, g_0_x_xxyyyy_xz, g_0_x_xxyyyy_yy, g_0_x_xxyyyy_yz, g_0_x_xxyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyyy_xx[k] = -g_0_x_xxyyy_xx[k] * ab_y + g_0_x_xxyyy_xxy[k];

                g_0_x_xxyyyy_xy[k] = -g_0_x_xxyyy_xy[k] * ab_y + g_0_x_xxyyy_xyy[k];

                g_0_x_xxyyyy_xz[k] = -g_0_x_xxyyy_xz[k] * ab_y + g_0_x_xxyyy_xyz[k];

                g_0_x_xxyyyy_yy[k] = -g_0_x_xxyyy_yy[k] * ab_y + g_0_x_xxyyy_yyy[k];

                g_0_x_xxyyyy_yz[k] = -g_0_x_xxyyy_yz[k] * ab_y + g_0_x_xxyyy_yyz[k];

                g_0_x_xxyyyy_zz[k] = -g_0_x_xxyyy_zz[k] * ab_y + g_0_x_xxyyy_yzz[k];
            }

            /// Set up 66-72 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyyz_xx = cbuffer.data(id_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xy = cbuffer.data(id_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xz = cbuffer.data(id_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xxyyyz_yy = cbuffer.data(id_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xxyyyz_yz = cbuffer.data(id_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xxyyyz_zz = cbuffer.data(id_geom_01_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyyz_xx, g_0_x_xxyyyz_xy, g_0_x_xxyyyz_xz, g_0_x_xxyyyz_yy, g_0_x_xxyyyz_yz, g_0_x_xxyyyz_zz, g_0_x_xxyyz_xx, g_0_x_xxyyz_xxy, g_0_x_xxyyz_xy, g_0_x_xxyyz_xyy, g_0_x_xxyyz_xyz, g_0_x_xxyyz_xz, g_0_x_xxyyz_yy, g_0_x_xxyyz_yyy, g_0_x_xxyyz_yyz, g_0_x_xxyyz_yz, g_0_x_xxyyz_yzz, g_0_x_xxyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyyz_xx[k] = -g_0_x_xxyyz_xx[k] * ab_y + g_0_x_xxyyz_xxy[k];

                g_0_x_xxyyyz_xy[k] = -g_0_x_xxyyz_xy[k] * ab_y + g_0_x_xxyyz_xyy[k];

                g_0_x_xxyyyz_xz[k] = -g_0_x_xxyyz_xz[k] * ab_y + g_0_x_xxyyz_xyz[k];

                g_0_x_xxyyyz_yy[k] = -g_0_x_xxyyz_yy[k] * ab_y + g_0_x_xxyyz_yyy[k];

                g_0_x_xxyyyz_yz[k] = -g_0_x_xxyyz_yz[k] * ab_y + g_0_x_xxyyz_yyz[k];

                g_0_x_xxyyyz_zz[k] = -g_0_x_xxyyz_zz[k] * ab_y + g_0_x_xxyyz_yzz[k];
            }

            /// Set up 72-78 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyzz_xx = cbuffer.data(id_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xy = cbuffer.data(id_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xz = cbuffer.data(id_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xxyyzz_yy = cbuffer.data(id_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xxyyzz_yz = cbuffer.data(id_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xxyyzz_zz = cbuffer.data(id_geom_01_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyzz_xx, g_0_x_xxyyzz_xy, g_0_x_xxyyzz_xz, g_0_x_xxyyzz_yy, g_0_x_xxyyzz_yz, g_0_x_xxyyzz_zz, g_0_x_xxyzz_xx, g_0_x_xxyzz_xxy, g_0_x_xxyzz_xy, g_0_x_xxyzz_xyy, g_0_x_xxyzz_xyz, g_0_x_xxyzz_xz, g_0_x_xxyzz_yy, g_0_x_xxyzz_yyy, g_0_x_xxyzz_yyz, g_0_x_xxyzz_yz, g_0_x_xxyzz_yzz, g_0_x_xxyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyzz_xx[k] = -g_0_x_xxyzz_xx[k] * ab_y + g_0_x_xxyzz_xxy[k];

                g_0_x_xxyyzz_xy[k] = -g_0_x_xxyzz_xy[k] * ab_y + g_0_x_xxyzz_xyy[k];

                g_0_x_xxyyzz_xz[k] = -g_0_x_xxyzz_xz[k] * ab_y + g_0_x_xxyzz_xyz[k];

                g_0_x_xxyyzz_yy[k] = -g_0_x_xxyzz_yy[k] * ab_y + g_0_x_xxyzz_yyy[k];

                g_0_x_xxyyzz_yz[k] = -g_0_x_xxyzz_yz[k] * ab_y + g_0_x_xxyzz_yyz[k];

                g_0_x_xxyyzz_zz[k] = -g_0_x_xxyzz_zz[k] * ab_y + g_0_x_xxyzz_yzz[k];
            }

            /// Set up 78-84 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyzzz_xx = cbuffer.data(id_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xy = cbuffer.data(id_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xz = cbuffer.data(id_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xxyzzz_yy = cbuffer.data(id_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xxyzzz_yz = cbuffer.data(id_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xxyzzz_zz = cbuffer.data(id_geom_01_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyzzz_xx, g_0_x_xxyzzz_xy, g_0_x_xxyzzz_xz, g_0_x_xxyzzz_yy, g_0_x_xxyzzz_yz, g_0_x_xxyzzz_zz, g_0_x_xxzzz_xx, g_0_x_xxzzz_xxy, g_0_x_xxzzz_xy, g_0_x_xxzzz_xyy, g_0_x_xxzzz_xyz, g_0_x_xxzzz_xz, g_0_x_xxzzz_yy, g_0_x_xxzzz_yyy, g_0_x_xxzzz_yyz, g_0_x_xxzzz_yz, g_0_x_xxzzz_yzz, g_0_x_xxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyzzz_xx[k] = -g_0_x_xxzzz_xx[k] * ab_y + g_0_x_xxzzz_xxy[k];

                g_0_x_xxyzzz_xy[k] = -g_0_x_xxzzz_xy[k] * ab_y + g_0_x_xxzzz_xyy[k];

                g_0_x_xxyzzz_xz[k] = -g_0_x_xxzzz_xz[k] * ab_y + g_0_x_xxzzz_xyz[k];

                g_0_x_xxyzzz_yy[k] = -g_0_x_xxzzz_yy[k] * ab_y + g_0_x_xxzzz_yyy[k];

                g_0_x_xxyzzz_yz[k] = -g_0_x_xxzzz_yz[k] * ab_y + g_0_x_xxzzz_yyz[k];

                g_0_x_xxyzzz_zz[k] = -g_0_x_xxzzz_zz[k] * ab_y + g_0_x_xxzzz_yzz[k];
            }

            /// Set up 84-90 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxzzzz_xx = cbuffer.data(id_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xy = cbuffer.data(id_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xz = cbuffer.data(id_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xxzzzz_yy = cbuffer.data(id_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xxzzzz_yz = cbuffer.data(id_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xxzzzz_zz = cbuffer.data(id_geom_01_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxzzz_xx, g_0_x_xxzzz_xxz, g_0_x_xxzzz_xy, g_0_x_xxzzz_xyz, g_0_x_xxzzz_xz, g_0_x_xxzzz_xzz, g_0_x_xxzzz_yy, g_0_x_xxzzz_yyz, g_0_x_xxzzz_yz, g_0_x_xxzzz_yzz, g_0_x_xxzzz_zz, g_0_x_xxzzz_zzz, g_0_x_xxzzzz_xx, g_0_x_xxzzzz_xy, g_0_x_xxzzzz_xz, g_0_x_xxzzzz_yy, g_0_x_xxzzzz_yz, g_0_x_xxzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxzzzz_xx[k] = -g_0_x_xxzzz_xx[k] * ab_z + g_0_x_xxzzz_xxz[k];

                g_0_x_xxzzzz_xy[k] = -g_0_x_xxzzz_xy[k] * ab_z + g_0_x_xxzzz_xyz[k];

                g_0_x_xxzzzz_xz[k] = -g_0_x_xxzzz_xz[k] * ab_z + g_0_x_xxzzz_xzz[k];

                g_0_x_xxzzzz_yy[k] = -g_0_x_xxzzz_yy[k] * ab_z + g_0_x_xxzzz_yyz[k];

                g_0_x_xxzzzz_yz[k] = -g_0_x_xxzzz_yz[k] * ab_z + g_0_x_xxzzz_yzz[k];

                g_0_x_xxzzzz_zz[k] = -g_0_x_xxzzz_zz[k] * ab_z + g_0_x_xxzzz_zzz[k];
            }

            /// Set up 90-96 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyyy_xx = cbuffer.data(id_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xy = cbuffer.data(id_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xz = cbuffer.data(id_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_xyyyyy_yy = cbuffer.data(id_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_xyyyyy_yz = cbuffer.data(id_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_xyyyyy_zz = cbuffer.data(id_geom_01_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyy_xx, g_0_x_xyyyy_xxy, g_0_x_xyyyy_xy, g_0_x_xyyyy_xyy, g_0_x_xyyyy_xyz, g_0_x_xyyyy_xz, g_0_x_xyyyy_yy, g_0_x_xyyyy_yyy, g_0_x_xyyyy_yyz, g_0_x_xyyyy_yz, g_0_x_xyyyy_yzz, g_0_x_xyyyy_zz, g_0_x_xyyyyy_xx, g_0_x_xyyyyy_xy, g_0_x_xyyyyy_xz, g_0_x_xyyyyy_yy, g_0_x_xyyyyy_yz, g_0_x_xyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyyy_xx[k] = -g_0_x_xyyyy_xx[k] * ab_y + g_0_x_xyyyy_xxy[k];

                g_0_x_xyyyyy_xy[k] = -g_0_x_xyyyy_xy[k] * ab_y + g_0_x_xyyyy_xyy[k];

                g_0_x_xyyyyy_xz[k] = -g_0_x_xyyyy_xz[k] * ab_y + g_0_x_xyyyy_xyz[k];

                g_0_x_xyyyyy_yy[k] = -g_0_x_xyyyy_yy[k] * ab_y + g_0_x_xyyyy_yyy[k];

                g_0_x_xyyyyy_yz[k] = -g_0_x_xyyyy_yz[k] * ab_y + g_0_x_xyyyy_yyz[k];

                g_0_x_xyyyyy_zz[k] = -g_0_x_xyyyy_zz[k] * ab_y + g_0_x_xyyyy_yzz[k];
            }

            /// Set up 96-102 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyyz_xx = cbuffer.data(id_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xy = cbuffer.data(id_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xz = cbuffer.data(id_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_xyyyyz_yy = cbuffer.data(id_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_xyyyyz_yz = cbuffer.data(id_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_xyyyyz_zz = cbuffer.data(id_geom_01_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyyz_xx, g_0_x_xyyyyz_xy, g_0_x_xyyyyz_xz, g_0_x_xyyyyz_yy, g_0_x_xyyyyz_yz, g_0_x_xyyyyz_zz, g_0_x_xyyyz_xx, g_0_x_xyyyz_xxy, g_0_x_xyyyz_xy, g_0_x_xyyyz_xyy, g_0_x_xyyyz_xyz, g_0_x_xyyyz_xz, g_0_x_xyyyz_yy, g_0_x_xyyyz_yyy, g_0_x_xyyyz_yyz, g_0_x_xyyyz_yz, g_0_x_xyyyz_yzz, g_0_x_xyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyyz_xx[k] = -g_0_x_xyyyz_xx[k] * ab_y + g_0_x_xyyyz_xxy[k];

                g_0_x_xyyyyz_xy[k] = -g_0_x_xyyyz_xy[k] * ab_y + g_0_x_xyyyz_xyy[k];

                g_0_x_xyyyyz_xz[k] = -g_0_x_xyyyz_xz[k] * ab_y + g_0_x_xyyyz_xyz[k];

                g_0_x_xyyyyz_yy[k] = -g_0_x_xyyyz_yy[k] * ab_y + g_0_x_xyyyz_yyy[k];

                g_0_x_xyyyyz_yz[k] = -g_0_x_xyyyz_yz[k] * ab_y + g_0_x_xyyyz_yyz[k];

                g_0_x_xyyyyz_zz[k] = -g_0_x_xyyyz_zz[k] * ab_y + g_0_x_xyyyz_yzz[k];
            }

            /// Set up 102-108 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyzz_xx = cbuffer.data(id_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xy = cbuffer.data(id_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xz = cbuffer.data(id_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_xyyyzz_yy = cbuffer.data(id_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_xyyyzz_yz = cbuffer.data(id_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_xyyyzz_zz = cbuffer.data(id_geom_01_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyzz_xx, g_0_x_xyyyzz_xy, g_0_x_xyyyzz_xz, g_0_x_xyyyzz_yy, g_0_x_xyyyzz_yz, g_0_x_xyyyzz_zz, g_0_x_xyyzz_xx, g_0_x_xyyzz_xxy, g_0_x_xyyzz_xy, g_0_x_xyyzz_xyy, g_0_x_xyyzz_xyz, g_0_x_xyyzz_xz, g_0_x_xyyzz_yy, g_0_x_xyyzz_yyy, g_0_x_xyyzz_yyz, g_0_x_xyyzz_yz, g_0_x_xyyzz_yzz, g_0_x_xyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyzz_xx[k] = -g_0_x_xyyzz_xx[k] * ab_y + g_0_x_xyyzz_xxy[k];

                g_0_x_xyyyzz_xy[k] = -g_0_x_xyyzz_xy[k] * ab_y + g_0_x_xyyzz_xyy[k];

                g_0_x_xyyyzz_xz[k] = -g_0_x_xyyzz_xz[k] * ab_y + g_0_x_xyyzz_xyz[k];

                g_0_x_xyyyzz_yy[k] = -g_0_x_xyyzz_yy[k] * ab_y + g_0_x_xyyzz_yyy[k];

                g_0_x_xyyyzz_yz[k] = -g_0_x_xyyzz_yz[k] * ab_y + g_0_x_xyyzz_yyz[k];

                g_0_x_xyyyzz_zz[k] = -g_0_x_xyyzz_zz[k] * ab_y + g_0_x_xyyzz_yzz[k];
            }

            /// Set up 108-114 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyzzz_xx = cbuffer.data(id_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xy = cbuffer.data(id_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xz = cbuffer.data(id_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_xyyzzz_yy = cbuffer.data(id_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_xyyzzz_yz = cbuffer.data(id_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_xyyzzz_zz = cbuffer.data(id_geom_01_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyzzz_xx, g_0_x_xyyzzz_xy, g_0_x_xyyzzz_xz, g_0_x_xyyzzz_yy, g_0_x_xyyzzz_yz, g_0_x_xyyzzz_zz, g_0_x_xyzzz_xx, g_0_x_xyzzz_xxy, g_0_x_xyzzz_xy, g_0_x_xyzzz_xyy, g_0_x_xyzzz_xyz, g_0_x_xyzzz_xz, g_0_x_xyzzz_yy, g_0_x_xyzzz_yyy, g_0_x_xyzzz_yyz, g_0_x_xyzzz_yz, g_0_x_xyzzz_yzz, g_0_x_xyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyzzz_xx[k] = -g_0_x_xyzzz_xx[k] * ab_y + g_0_x_xyzzz_xxy[k];

                g_0_x_xyyzzz_xy[k] = -g_0_x_xyzzz_xy[k] * ab_y + g_0_x_xyzzz_xyy[k];

                g_0_x_xyyzzz_xz[k] = -g_0_x_xyzzz_xz[k] * ab_y + g_0_x_xyzzz_xyz[k];

                g_0_x_xyyzzz_yy[k] = -g_0_x_xyzzz_yy[k] * ab_y + g_0_x_xyzzz_yyy[k];

                g_0_x_xyyzzz_yz[k] = -g_0_x_xyzzz_yz[k] * ab_y + g_0_x_xyzzz_yyz[k];

                g_0_x_xyyzzz_zz[k] = -g_0_x_xyzzz_zz[k] * ab_y + g_0_x_xyzzz_yzz[k];
            }

            /// Set up 114-120 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyzzzz_xx = cbuffer.data(id_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xy = cbuffer.data(id_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xz = cbuffer.data(id_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_xyzzzz_yy = cbuffer.data(id_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_xyzzzz_yz = cbuffer.data(id_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_xyzzzz_zz = cbuffer.data(id_geom_01_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyzzzz_xx, g_0_x_xyzzzz_xy, g_0_x_xyzzzz_xz, g_0_x_xyzzzz_yy, g_0_x_xyzzzz_yz, g_0_x_xyzzzz_zz, g_0_x_xzzzz_xx, g_0_x_xzzzz_xxy, g_0_x_xzzzz_xy, g_0_x_xzzzz_xyy, g_0_x_xzzzz_xyz, g_0_x_xzzzz_xz, g_0_x_xzzzz_yy, g_0_x_xzzzz_yyy, g_0_x_xzzzz_yyz, g_0_x_xzzzz_yz, g_0_x_xzzzz_yzz, g_0_x_xzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyzzzz_xx[k] = -g_0_x_xzzzz_xx[k] * ab_y + g_0_x_xzzzz_xxy[k];

                g_0_x_xyzzzz_xy[k] = -g_0_x_xzzzz_xy[k] * ab_y + g_0_x_xzzzz_xyy[k];

                g_0_x_xyzzzz_xz[k] = -g_0_x_xzzzz_xz[k] * ab_y + g_0_x_xzzzz_xyz[k];

                g_0_x_xyzzzz_yy[k] = -g_0_x_xzzzz_yy[k] * ab_y + g_0_x_xzzzz_yyy[k];

                g_0_x_xyzzzz_yz[k] = -g_0_x_xzzzz_yz[k] * ab_y + g_0_x_xzzzz_yyz[k];

                g_0_x_xyzzzz_zz[k] = -g_0_x_xzzzz_zz[k] * ab_y + g_0_x_xzzzz_yzz[k];
            }

            /// Set up 120-126 components of targeted buffer : cbuffer.data(

            auto g_0_x_xzzzzz_xx = cbuffer.data(id_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xy = cbuffer.data(id_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xz = cbuffer.data(id_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_xzzzzz_yy = cbuffer.data(id_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_xzzzzz_yz = cbuffer.data(id_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_xzzzzz_zz = cbuffer.data(id_geom_01_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xzzzz_xx, g_0_x_xzzzz_xxz, g_0_x_xzzzz_xy, g_0_x_xzzzz_xyz, g_0_x_xzzzz_xz, g_0_x_xzzzz_xzz, g_0_x_xzzzz_yy, g_0_x_xzzzz_yyz, g_0_x_xzzzz_yz, g_0_x_xzzzz_yzz, g_0_x_xzzzz_zz, g_0_x_xzzzz_zzz, g_0_x_xzzzzz_xx, g_0_x_xzzzzz_xy, g_0_x_xzzzzz_xz, g_0_x_xzzzzz_yy, g_0_x_xzzzzz_yz, g_0_x_xzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xzzzzz_xx[k] = -g_0_x_xzzzz_xx[k] * ab_z + g_0_x_xzzzz_xxz[k];

                g_0_x_xzzzzz_xy[k] = -g_0_x_xzzzz_xy[k] * ab_z + g_0_x_xzzzz_xyz[k];

                g_0_x_xzzzzz_xz[k] = -g_0_x_xzzzz_xz[k] * ab_z + g_0_x_xzzzz_xzz[k];

                g_0_x_xzzzzz_yy[k] = -g_0_x_xzzzz_yy[k] * ab_z + g_0_x_xzzzz_yyz[k];

                g_0_x_xzzzzz_yz[k] = -g_0_x_xzzzz_yz[k] * ab_z + g_0_x_xzzzz_yzz[k];

                g_0_x_xzzzzz_zz[k] = -g_0_x_xzzzz_zz[k] * ab_z + g_0_x_xzzzz_zzz[k];
            }

            /// Set up 126-132 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyyy_xx = cbuffer.data(id_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xy = cbuffer.data(id_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xz = cbuffer.data(id_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_yyyyyy_yy = cbuffer.data(id_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_yyyyyy_yz = cbuffer.data(id_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_yyyyyy_zz = cbuffer.data(id_geom_01_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyy_xx, g_0_x_yyyyy_xxy, g_0_x_yyyyy_xy, g_0_x_yyyyy_xyy, g_0_x_yyyyy_xyz, g_0_x_yyyyy_xz, g_0_x_yyyyy_yy, g_0_x_yyyyy_yyy, g_0_x_yyyyy_yyz, g_0_x_yyyyy_yz, g_0_x_yyyyy_yzz, g_0_x_yyyyy_zz, g_0_x_yyyyyy_xx, g_0_x_yyyyyy_xy, g_0_x_yyyyyy_xz, g_0_x_yyyyyy_yy, g_0_x_yyyyyy_yz, g_0_x_yyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyyy_xx[k] = -g_0_x_yyyyy_xx[k] * ab_y + g_0_x_yyyyy_xxy[k];

                g_0_x_yyyyyy_xy[k] = -g_0_x_yyyyy_xy[k] * ab_y + g_0_x_yyyyy_xyy[k];

                g_0_x_yyyyyy_xz[k] = -g_0_x_yyyyy_xz[k] * ab_y + g_0_x_yyyyy_xyz[k];

                g_0_x_yyyyyy_yy[k] = -g_0_x_yyyyy_yy[k] * ab_y + g_0_x_yyyyy_yyy[k];

                g_0_x_yyyyyy_yz[k] = -g_0_x_yyyyy_yz[k] * ab_y + g_0_x_yyyyy_yyz[k];

                g_0_x_yyyyyy_zz[k] = -g_0_x_yyyyy_zz[k] * ab_y + g_0_x_yyyyy_yzz[k];
            }

            /// Set up 132-138 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyyz_xx = cbuffer.data(id_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xy = cbuffer.data(id_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xz = cbuffer.data(id_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_x_yyyyyz_yy = cbuffer.data(id_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_yyyyyz_yz = cbuffer.data(id_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_yyyyyz_zz = cbuffer.data(id_geom_01_off + 137 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyyz_xx, g_0_x_yyyyyz_xy, g_0_x_yyyyyz_xz, g_0_x_yyyyyz_yy, g_0_x_yyyyyz_yz, g_0_x_yyyyyz_zz, g_0_x_yyyyz_xx, g_0_x_yyyyz_xxy, g_0_x_yyyyz_xy, g_0_x_yyyyz_xyy, g_0_x_yyyyz_xyz, g_0_x_yyyyz_xz, g_0_x_yyyyz_yy, g_0_x_yyyyz_yyy, g_0_x_yyyyz_yyz, g_0_x_yyyyz_yz, g_0_x_yyyyz_yzz, g_0_x_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyyz_xx[k] = -g_0_x_yyyyz_xx[k] * ab_y + g_0_x_yyyyz_xxy[k];

                g_0_x_yyyyyz_xy[k] = -g_0_x_yyyyz_xy[k] * ab_y + g_0_x_yyyyz_xyy[k];

                g_0_x_yyyyyz_xz[k] = -g_0_x_yyyyz_xz[k] * ab_y + g_0_x_yyyyz_xyz[k];

                g_0_x_yyyyyz_yy[k] = -g_0_x_yyyyz_yy[k] * ab_y + g_0_x_yyyyz_yyy[k];

                g_0_x_yyyyyz_yz[k] = -g_0_x_yyyyz_yz[k] * ab_y + g_0_x_yyyyz_yyz[k];

                g_0_x_yyyyyz_zz[k] = -g_0_x_yyyyz_zz[k] * ab_y + g_0_x_yyyyz_yzz[k];
            }

            /// Set up 138-144 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyzz_xx = cbuffer.data(id_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xy = cbuffer.data(id_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xz = cbuffer.data(id_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_yyyyzz_yy = cbuffer.data(id_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_yyyyzz_yz = cbuffer.data(id_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_yyyyzz_zz = cbuffer.data(id_geom_01_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyzz_xx, g_0_x_yyyyzz_xy, g_0_x_yyyyzz_xz, g_0_x_yyyyzz_yy, g_0_x_yyyyzz_yz, g_0_x_yyyyzz_zz, g_0_x_yyyzz_xx, g_0_x_yyyzz_xxy, g_0_x_yyyzz_xy, g_0_x_yyyzz_xyy, g_0_x_yyyzz_xyz, g_0_x_yyyzz_xz, g_0_x_yyyzz_yy, g_0_x_yyyzz_yyy, g_0_x_yyyzz_yyz, g_0_x_yyyzz_yz, g_0_x_yyyzz_yzz, g_0_x_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyzz_xx[k] = -g_0_x_yyyzz_xx[k] * ab_y + g_0_x_yyyzz_xxy[k];

                g_0_x_yyyyzz_xy[k] = -g_0_x_yyyzz_xy[k] * ab_y + g_0_x_yyyzz_xyy[k];

                g_0_x_yyyyzz_xz[k] = -g_0_x_yyyzz_xz[k] * ab_y + g_0_x_yyyzz_xyz[k];

                g_0_x_yyyyzz_yy[k] = -g_0_x_yyyzz_yy[k] * ab_y + g_0_x_yyyzz_yyy[k];

                g_0_x_yyyyzz_yz[k] = -g_0_x_yyyzz_yz[k] * ab_y + g_0_x_yyyzz_yyz[k];

                g_0_x_yyyyzz_zz[k] = -g_0_x_yyyzz_zz[k] * ab_y + g_0_x_yyyzz_yzz[k];
            }

            /// Set up 144-150 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyzzz_xx = cbuffer.data(id_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xy = cbuffer.data(id_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xz = cbuffer.data(id_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_yyyzzz_yy = cbuffer.data(id_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_yyyzzz_yz = cbuffer.data(id_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_yyyzzz_zz = cbuffer.data(id_geom_01_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyzzz_xx, g_0_x_yyyzzz_xy, g_0_x_yyyzzz_xz, g_0_x_yyyzzz_yy, g_0_x_yyyzzz_yz, g_0_x_yyyzzz_zz, g_0_x_yyzzz_xx, g_0_x_yyzzz_xxy, g_0_x_yyzzz_xy, g_0_x_yyzzz_xyy, g_0_x_yyzzz_xyz, g_0_x_yyzzz_xz, g_0_x_yyzzz_yy, g_0_x_yyzzz_yyy, g_0_x_yyzzz_yyz, g_0_x_yyzzz_yz, g_0_x_yyzzz_yzz, g_0_x_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyzzz_xx[k] = -g_0_x_yyzzz_xx[k] * ab_y + g_0_x_yyzzz_xxy[k];

                g_0_x_yyyzzz_xy[k] = -g_0_x_yyzzz_xy[k] * ab_y + g_0_x_yyzzz_xyy[k];

                g_0_x_yyyzzz_xz[k] = -g_0_x_yyzzz_xz[k] * ab_y + g_0_x_yyzzz_xyz[k];

                g_0_x_yyyzzz_yy[k] = -g_0_x_yyzzz_yy[k] * ab_y + g_0_x_yyzzz_yyy[k];

                g_0_x_yyyzzz_yz[k] = -g_0_x_yyzzz_yz[k] * ab_y + g_0_x_yyzzz_yyz[k];

                g_0_x_yyyzzz_zz[k] = -g_0_x_yyzzz_zz[k] * ab_y + g_0_x_yyzzz_yzz[k];
            }

            /// Set up 150-156 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyzzzz_xx = cbuffer.data(id_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xy = cbuffer.data(id_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xz = cbuffer.data(id_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_x_yyzzzz_yy = cbuffer.data(id_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_x_yyzzzz_yz = cbuffer.data(id_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_x_yyzzzz_zz = cbuffer.data(id_geom_01_off + 155 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyzzzz_xx, g_0_x_yyzzzz_xy, g_0_x_yyzzzz_xz, g_0_x_yyzzzz_yy, g_0_x_yyzzzz_yz, g_0_x_yyzzzz_zz, g_0_x_yzzzz_xx, g_0_x_yzzzz_xxy, g_0_x_yzzzz_xy, g_0_x_yzzzz_xyy, g_0_x_yzzzz_xyz, g_0_x_yzzzz_xz, g_0_x_yzzzz_yy, g_0_x_yzzzz_yyy, g_0_x_yzzzz_yyz, g_0_x_yzzzz_yz, g_0_x_yzzzz_yzz, g_0_x_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyzzzz_xx[k] = -g_0_x_yzzzz_xx[k] * ab_y + g_0_x_yzzzz_xxy[k];

                g_0_x_yyzzzz_xy[k] = -g_0_x_yzzzz_xy[k] * ab_y + g_0_x_yzzzz_xyy[k];

                g_0_x_yyzzzz_xz[k] = -g_0_x_yzzzz_xz[k] * ab_y + g_0_x_yzzzz_xyz[k];

                g_0_x_yyzzzz_yy[k] = -g_0_x_yzzzz_yy[k] * ab_y + g_0_x_yzzzz_yyy[k];

                g_0_x_yyzzzz_yz[k] = -g_0_x_yzzzz_yz[k] * ab_y + g_0_x_yzzzz_yyz[k];

                g_0_x_yyzzzz_zz[k] = -g_0_x_yzzzz_zz[k] * ab_y + g_0_x_yzzzz_yzz[k];
            }

            /// Set up 156-162 components of targeted buffer : cbuffer.data(

            auto g_0_x_yzzzzz_xx = cbuffer.data(id_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xy = cbuffer.data(id_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xz = cbuffer.data(id_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_x_yzzzzz_yy = cbuffer.data(id_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_x_yzzzzz_yz = cbuffer.data(id_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_x_yzzzzz_zz = cbuffer.data(id_geom_01_off + 161 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yzzzzz_xx, g_0_x_yzzzzz_xy, g_0_x_yzzzzz_xz, g_0_x_yzzzzz_yy, g_0_x_yzzzzz_yz, g_0_x_yzzzzz_zz, g_0_x_zzzzz_xx, g_0_x_zzzzz_xxy, g_0_x_zzzzz_xy, g_0_x_zzzzz_xyy, g_0_x_zzzzz_xyz, g_0_x_zzzzz_xz, g_0_x_zzzzz_yy, g_0_x_zzzzz_yyy, g_0_x_zzzzz_yyz, g_0_x_zzzzz_yz, g_0_x_zzzzz_yzz, g_0_x_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yzzzzz_xx[k] = -g_0_x_zzzzz_xx[k] * ab_y + g_0_x_zzzzz_xxy[k];

                g_0_x_yzzzzz_xy[k] = -g_0_x_zzzzz_xy[k] * ab_y + g_0_x_zzzzz_xyy[k];

                g_0_x_yzzzzz_xz[k] = -g_0_x_zzzzz_xz[k] * ab_y + g_0_x_zzzzz_xyz[k];

                g_0_x_yzzzzz_yy[k] = -g_0_x_zzzzz_yy[k] * ab_y + g_0_x_zzzzz_yyy[k];

                g_0_x_yzzzzz_yz[k] = -g_0_x_zzzzz_yz[k] * ab_y + g_0_x_zzzzz_yyz[k];

                g_0_x_yzzzzz_zz[k] = -g_0_x_zzzzz_zz[k] * ab_y + g_0_x_zzzzz_yzz[k];
            }

            /// Set up 162-168 components of targeted buffer : cbuffer.data(

            auto g_0_x_zzzzzz_xx = cbuffer.data(id_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xy = cbuffer.data(id_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xz = cbuffer.data(id_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_x_zzzzzz_yy = cbuffer.data(id_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_x_zzzzzz_yz = cbuffer.data(id_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_x_zzzzzz_zz = cbuffer.data(id_geom_01_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_zzzzz_xx, g_0_x_zzzzz_xxz, g_0_x_zzzzz_xy, g_0_x_zzzzz_xyz, g_0_x_zzzzz_xz, g_0_x_zzzzz_xzz, g_0_x_zzzzz_yy, g_0_x_zzzzz_yyz, g_0_x_zzzzz_yz, g_0_x_zzzzz_yzz, g_0_x_zzzzz_zz, g_0_x_zzzzz_zzz, g_0_x_zzzzzz_xx, g_0_x_zzzzzz_xy, g_0_x_zzzzzz_xz, g_0_x_zzzzzz_yy, g_0_x_zzzzzz_yz, g_0_x_zzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_zzzzzz_xx[k] = -g_0_x_zzzzz_xx[k] * ab_z + g_0_x_zzzzz_xxz[k];

                g_0_x_zzzzzz_xy[k] = -g_0_x_zzzzz_xy[k] * ab_z + g_0_x_zzzzz_xyz[k];

                g_0_x_zzzzzz_xz[k] = -g_0_x_zzzzz_xz[k] * ab_z + g_0_x_zzzzz_xzz[k];

                g_0_x_zzzzzz_yy[k] = -g_0_x_zzzzz_yy[k] * ab_z + g_0_x_zzzzz_yyz[k];

                g_0_x_zzzzzz_yz[k] = -g_0_x_zzzzz_yz[k] * ab_z + g_0_x_zzzzz_yzz[k];

                g_0_x_zzzzzz_zz[k] = -g_0_x_zzzzz_zz[k] * ab_z + g_0_x_zzzzz_zzz[k];
            }

            /// Set up 168-174 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxx_xx = cbuffer.data(id_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xy = cbuffer.data(id_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xz = cbuffer.data(id_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_y_xxxxxx_yy = cbuffer.data(id_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_y_xxxxxx_yz = cbuffer.data(id_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_y_xxxxxx_zz = cbuffer.data(id_geom_01_off + 173 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxx_xx, g_0_y_xxxxx_xxx, g_0_y_xxxxx_xxy, g_0_y_xxxxx_xxz, g_0_y_xxxxx_xy, g_0_y_xxxxx_xyy, g_0_y_xxxxx_xyz, g_0_y_xxxxx_xz, g_0_y_xxxxx_xzz, g_0_y_xxxxx_yy, g_0_y_xxxxx_yz, g_0_y_xxxxx_zz, g_0_y_xxxxxx_xx, g_0_y_xxxxxx_xy, g_0_y_xxxxxx_xz, g_0_y_xxxxxx_yy, g_0_y_xxxxxx_yz, g_0_y_xxxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxx_xx[k] = -g_0_y_xxxxx_xx[k] * ab_x + g_0_y_xxxxx_xxx[k];

                g_0_y_xxxxxx_xy[k] = -g_0_y_xxxxx_xy[k] * ab_x + g_0_y_xxxxx_xxy[k];

                g_0_y_xxxxxx_xz[k] = -g_0_y_xxxxx_xz[k] * ab_x + g_0_y_xxxxx_xxz[k];

                g_0_y_xxxxxx_yy[k] = -g_0_y_xxxxx_yy[k] * ab_x + g_0_y_xxxxx_xyy[k];

                g_0_y_xxxxxx_yz[k] = -g_0_y_xxxxx_yz[k] * ab_x + g_0_y_xxxxx_xyz[k];

                g_0_y_xxxxxx_zz[k] = -g_0_y_xxxxx_zz[k] * ab_x + g_0_y_xxxxx_xzz[k];
            }

            /// Set up 174-180 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxy_xx = cbuffer.data(id_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xy = cbuffer.data(id_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xz = cbuffer.data(id_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_y_xxxxxy_yy = cbuffer.data(id_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_y_xxxxxy_yz = cbuffer.data(id_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_y_xxxxxy_zz = cbuffer.data(id_geom_01_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxy_xx, g_0_y_xxxxxy_xy, g_0_y_xxxxxy_xz, g_0_y_xxxxxy_yy, g_0_y_xxxxxy_yz, g_0_y_xxxxxy_zz, g_0_y_xxxxy_xx, g_0_y_xxxxy_xxx, g_0_y_xxxxy_xxy, g_0_y_xxxxy_xxz, g_0_y_xxxxy_xy, g_0_y_xxxxy_xyy, g_0_y_xxxxy_xyz, g_0_y_xxxxy_xz, g_0_y_xxxxy_xzz, g_0_y_xxxxy_yy, g_0_y_xxxxy_yz, g_0_y_xxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxy_xx[k] = -g_0_y_xxxxy_xx[k] * ab_x + g_0_y_xxxxy_xxx[k];

                g_0_y_xxxxxy_xy[k] = -g_0_y_xxxxy_xy[k] * ab_x + g_0_y_xxxxy_xxy[k];

                g_0_y_xxxxxy_xz[k] = -g_0_y_xxxxy_xz[k] * ab_x + g_0_y_xxxxy_xxz[k];

                g_0_y_xxxxxy_yy[k] = -g_0_y_xxxxy_yy[k] * ab_x + g_0_y_xxxxy_xyy[k];

                g_0_y_xxxxxy_yz[k] = -g_0_y_xxxxy_yz[k] * ab_x + g_0_y_xxxxy_xyz[k];

                g_0_y_xxxxxy_zz[k] = -g_0_y_xxxxy_zz[k] * ab_x + g_0_y_xxxxy_xzz[k];
            }

            /// Set up 180-186 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxz_xx = cbuffer.data(id_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xy = cbuffer.data(id_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xz = cbuffer.data(id_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_y_xxxxxz_yy = cbuffer.data(id_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_y_xxxxxz_yz = cbuffer.data(id_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_y_xxxxxz_zz = cbuffer.data(id_geom_01_off + 185 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxz_xx, g_0_y_xxxxxz_xy, g_0_y_xxxxxz_xz, g_0_y_xxxxxz_yy, g_0_y_xxxxxz_yz, g_0_y_xxxxxz_zz, g_0_y_xxxxz_xx, g_0_y_xxxxz_xxx, g_0_y_xxxxz_xxy, g_0_y_xxxxz_xxz, g_0_y_xxxxz_xy, g_0_y_xxxxz_xyy, g_0_y_xxxxz_xyz, g_0_y_xxxxz_xz, g_0_y_xxxxz_xzz, g_0_y_xxxxz_yy, g_0_y_xxxxz_yz, g_0_y_xxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxz_xx[k] = -g_0_y_xxxxz_xx[k] * ab_x + g_0_y_xxxxz_xxx[k];

                g_0_y_xxxxxz_xy[k] = -g_0_y_xxxxz_xy[k] * ab_x + g_0_y_xxxxz_xxy[k];

                g_0_y_xxxxxz_xz[k] = -g_0_y_xxxxz_xz[k] * ab_x + g_0_y_xxxxz_xxz[k];

                g_0_y_xxxxxz_yy[k] = -g_0_y_xxxxz_yy[k] * ab_x + g_0_y_xxxxz_xyy[k];

                g_0_y_xxxxxz_yz[k] = -g_0_y_xxxxz_yz[k] * ab_x + g_0_y_xxxxz_xyz[k];

                g_0_y_xxxxxz_zz[k] = -g_0_y_xxxxz_zz[k] * ab_x + g_0_y_xxxxz_xzz[k];
            }

            /// Set up 186-192 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxyy_xx = cbuffer.data(id_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xy = cbuffer.data(id_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xz = cbuffer.data(id_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_y_xxxxyy_yy = cbuffer.data(id_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_y_xxxxyy_yz = cbuffer.data(id_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_y_xxxxyy_zz = cbuffer.data(id_geom_01_off + 191 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxyy_xx, g_0_y_xxxxyy_xy, g_0_y_xxxxyy_xz, g_0_y_xxxxyy_yy, g_0_y_xxxxyy_yz, g_0_y_xxxxyy_zz, g_0_y_xxxyy_xx, g_0_y_xxxyy_xxx, g_0_y_xxxyy_xxy, g_0_y_xxxyy_xxz, g_0_y_xxxyy_xy, g_0_y_xxxyy_xyy, g_0_y_xxxyy_xyz, g_0_y_xxxyy_xz, g_0_y_xxxyy_xzz, g_0_y_xxxyy_yy, g_0_y_xxxyy_yz, g_0_y_xxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxyy_xx[k] = -g_0_y_xxxyy_xx[k] * ab_x + g_0_y_xxxyy_xxx[k];

                g_0_y_xxxxyy_xy[k] = -g_0_y_xxxyy_xy[k] * ab_x + g_0_y_xxxyy_xxy[k];

                g_0_y_xxxxyy_xz[k] = -g_0_y_xxxyy_xz[k] * ab_x + g_0_y_xxxyy_xxz[k];

                g_0_y_xxxxyy_yy[k] = -g_0_y_xxxyy_yy[k] * ab_x + g_0_y_xxxyy_xyy[k];

                g_0_y_xxxxyy_yz[k] = -g_0_y_xxxyy_yz[k] * ab_x + g_0_y_xxxyy_xyz[k];

                g_0_y_xxxxyy_zz[k] = -g_0_y_xxxyy_zz[k] * ab_x + g_0_y_xxxyy_xzz[k];
            }

            /// Set up 192-198 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxyz_xx = cbuffer.data(id_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xy = cbuffer.data(id_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xz = cbuffer.data(id_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_y_xxxxyz_yy = cbuffer.data(id_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_y_xxxxyz_yz = cbuffer.data(id_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_y_xxxxyz_zz = cbuffer.data(id_geom_01_off + 197 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxyz_xx, g_0_y_xxxxyz_xy, g_0_y_xxxxyz_xz, g_0_y_xxxxyz_yy, g_0_y_xxxxyz_yz, g_0_y_xxxxyz_zz, g_0_y_xxxyz_xx, g_0_y_xxxyz_xxx, g_0_y_xxxyz_xxy, g_0_y_xxxyz_xxz, g_0_y_xxxyz_xy, g_0_y_xxxyz_xyy, g_0_y_xxxyz_xyz, g_0_y_xxxyz_xz, g_0_y_xxxyz_xzz, g_0_y_xxxyz_yy, g_0_y_xxxyz_yz, g_0_y_xxxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxyz_xx[k] = -g_0_y_xxxyz_xx[k] * ab_x + g_0_y_xxxyz_xxx[k];

                g_0_y_xxxxyz_xy[k] = -g_0_y_xxxyz_xy[k] * ab_x + g_0_y_xxxyz_xxy[k];

                g_0_y_xxxxyz_xz[k] = -g_0_y_xxxyz_xz[k] * ab_x + g_0_y_xxxyz_xxz[k];

                g_0_y_xxxxyz_yy[k] = -g_0_y_xxxyz_yy[k] * ab_x + g_0_y_xxxyz_xyy[k];

                g_0_y_xxxxyz_yz[k] = -g_0_y_xxxyz_yz[k] * ab_x + g_0_y_xxxyz_xyz[k];

                g_0_y_xxxxyz_zz[k] = -g_0_y_xxxyz_zz[k] * ab_x + g_0_y_xxxyz_xzz[k];
            }

            /// Set up 198-204 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxzz_xx = cbuffer.data(id_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xy = cbuffer.data(id_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xz = cbuffer.data(id_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_y_xxxxzz_yy = cbuffer.data(id_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_y_xxxxzz_yz = cbuffer.data(id_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_y_xxxxzz_zz = cbuffer.data(id_geom_01_off + 203 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxzz_xx, g_0_y_xxxxzz_xy, g_0_y_xxxxzz_xz, g_0_y_xxxxzz_yy, g_0_y_xxxxzz_yz, g_0_y_xxxxzz_zz, g_0_y_xxxzz_xx, g_0_y_xxxzz_xxx, g_0_y_xxxzz_xxy, g_0_y_xxxzz_xxz, g_0_y_xxxzz_xy, g_0_y_xxxzz_xyy, g_0_y_xxxzz_xyz, g_0_y_xxxzz_xz, g_0_y_xxxzz_xzz, g_0_y_xxxzz_yy, g_0_y_xxxzz_yz, g_0_y_xxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxzz_xx[k] = -g_0_y_xxxzz_xx[k] * ab_x + g_0_y_xxxzz_xxx[k];

                g_0_y_xxxxzz_xy[k] = -g_0_y_xxxzz_xy[k] * ab_x + g_0_y_xxxzz_xxy[k];

                g_0_y_xxxxzz_xz[k] = -g_0_y_xxxzz_xz[k] * ab_x + g_0_y_xxxzz_xxz[k];

                g_0_y_xxxxzz_yy[k] = -g_0_y_xxxzz_yy[k] * ab_x + g_0_y_xxxzz_xyy[k];

                g_0_y_xxxxzz_yz[k] = -g_0_y_xxxzz_yz[k] * ab_x + g_0_y_xxxzz_xyz[k];

                g_0_y_xxxxzz_zz[k] = -g_0_y_xxxzz_zz[k] * ab_x + g_0_y_xxxzz_xzz[k];
            }

            /// Set up 204-210 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyyy_xx = cbuffer.data(id_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xy = cbuffer.data(id_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xz = cbuffer.data(id_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_y_xxxyyy_yy = cbuffer.data(id_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_y_xxxyyy_yz = cbuffer.data(id_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_y_xxxyyy_zz = cbuffer.data(id_geom_01_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyyy_xx, g_0_y_xxxyyy_xy, g_0_y_xxxyyy_xz, g_0_y_xxxyyy_yy, g_0_y_xxxyyy_yz, g_0_y_xxxyyy_zz, g_0_y_xxyyy_xx, g_0_y_xxyyy_xxx, g_0_y_xxyyy_xxy, g_0_y_xxyyy_xxz, g_0_y_xxyyy_xy, g_0_y_xxyyy_xyy, g_0_y_xxyyy_xyz, g_0_y_xxyyy_xz, g_0_y_xxyyy_xzz, g_0_y_xxyyy_yy, g_0_y_xxyyy_yz, g_0_y_xxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyyy_xx[k] = -g_0_y_xxyyy_xx[k] * ab_x + g_0_y_xxyyy_xxx[k];

                g_0_y_xxxyyy_xy[k] = -g_0_y_xxyyy_xy[k] * ab_x + g_0_y_xxyyy_xxy[k];

                g_0_y_xxxyyy_xz[k] = -g_0_y_xxyyy_xz[k] * ab_x + g_0_y_xxyyy_xxz[k];

                g_0_y_xxxyyy_yy[k] = -g_0_y_xxyyy_yy[k] * ab_x + g_0_y_xxyyy_xyy[k];

                g_0_y_xxxyyy_yz[k] = -g_0_y_xxyyy_yz[k] * ab_x + g_0_y_xxyyy_xyz[k];

                g_0_y_xxxyyy_zz[k] = -g_0_y_xxyyy_zz[k] * ab_x + g_0_y_xxyyy_xzz[k];
            }

            /// Set up 210-216 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyyz_xx = cbuffer.data(id_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xy = cbuffer.data(id_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xz = cbuffer.data(id_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_y_xxxyyz_yy = cbuffer.data(id_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_y_xxxyyz_yz = cbuffer.data(id_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_y_xxxyyz_zz = cbuffer.data(id_geom_01_off + 215 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyyz_xx, g_0_y_xxxyyz_xy, g_0_y_xxxyyz_xz, g_0_y_xxxyyz_yy, g_0_y_xxxyyz_yz, g_0_y_xxxyyz_zz, g_0_y_xxyyz_xx, g_0_y_xxyyz_xxx, g_0_y_xxyyz_xxy, g_0_y_xxyyz_xxz, g_0_y_xxyyz_xy, g_0_y_xxyyz_xyy, g_0_y_xxyyz_xyz, g_0_y_xxyyz_xz, g_0_y_xxyyz_xzz, g_0_y_xxyyz_yy, g_0_y_xxyyz_yz, g_0_y_xxyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyyz_xx[k] = -g_0_y_xxyyz_xx[k] * ab_x + g_0_y_xxyyz_xxx[k];

                g_0_y_xxxyyz_xy[k] = -g_0_y_xxyyz_xy[k] * ab_x + g_0_y_xxyyz_xxy[k];

                g_0_y_xxxyyz_xz[k] = -g_0_y_xxyyz_xz[k] * ab_x + g_0_y_xxyyz_xxz[k];

                g_0_y_xxxyyz_yy[k] = -g_0_y_xxyyz_yy[k] * ab_x + g_0_y_xxyyz_xyy[k];

                g_0_y_xxxyyz_yz[k] = -g_0_y_xxyyz_yz[k] * ab_x + g_0_y_xxyyz_xyz[k];

                g_0_y_xxxyyz_zz[k] = -g_0_y_xxyyz_zz[k] * ab_x + g_0_y_xxyyz_xzz[k];
            }

            /// Set up 216-222 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyzz_xx = cbuffer.data(id_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xy = cbuffer.data(id_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xz = cbuffer.data(id_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_y_xxxyzz_yy = cbuffer.data(id_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_y_xxxyzz_yz = cbuffer.data(id_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_y_xxxyzz_zz = cbuffer.data(id_geom_01_off + 221 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyzz_xx, g_0_y_xxxyzz_xy, g_0_y_xxxyzz_xz, g_0_y_xxxyzz_yy, g_0_y_xxxyzz_yz, g_0_y_xxxyzz_zz, g_0_y_xxyzz_xx, g_0_y_xxyzz_xxx, g_0_y_xxyzz_xxy, g_0_y_xxyzz_xxz, g_0_y_xxyzz_xy, g_0_y_xxyzz_xyy, g_0_y_xxyzz_xyz, g_0_y_xxyzz_xz, g_0_y_xxyzz_xzz, g_0_y_xxyzz_yy, g_0_y_xxyzz_yz, g_0_y_xxyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyzz_xx[k] = -g_0_y_xxyzz_xx[k] * ab_x + g_0_y_xxyzz_xxx[k];

                g_0_y_xxxyzz_xy[k] = -g_0_y_xxyzz_xy[k] * ab_x + g_0_y_xxyzz_xxy[k];

                g_0_y_xxxyzz_xz[k] = -g_0_y_xxyzz_xz[k] * ab_x + g_0_y_xxyzz_xxz[k];

                g_0_y_xxxyzz_yy[k] = -g_0_y_xxyzz_yy[k] * ab_x + g_0_y_xxyzz_xyy[k];

                g_0_y_xxxyzz_yz[k] = -g_0_y_xxyzz_yz[k] * ab_x + g_0_y_xxyzz_xyz[k];

                g_0_y_xxxyzz_zz[k] = -g_0_y_xxyzz_zz[k] * ab_x + g_0_y_xxyzz_xzz[k];
            }

            /// Set up 222-228 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxzzz_xx = cbuffer.data(id_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xy = cbuffer.data(id_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xz = cbuffer.data(id_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_y_xxxzzz_yy = cbuffer.data(id_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_y_xxxzzz_yz = cbuffer.data(id_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_y_xxxzzz_zz = cbuffer.data(id_geom_01_off + 227 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxzzz_xx, g_0_y_xxxzzz_xy, g_0_y_xxxzzz_xz, g_0_y_xxxzzz_yy, g_0_y_xxxzzz_yz, g_0_y_xxxzzz_zz, g_0_y_xxzzz_xx, g_0_y_xxzzz_xxx, g_0_y_xxzzz_xxy, g_0_y_xxzzz_xxz, g_0_y_xxzzz_xy, g_0_y_xxzzz_xyy, g_0_y_xxzzz_xyz, g_0_y_xxzzz_xz, g_0_y_xxzzz_xzz, g_0_y_xxzzz_yy, g_0_y_xxzzz_yz, g_0_y_xxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxzzz_xx[k] = -g_0_y_xxzzz_xx[k] * ab_x + g_0_y_xxzzz_xxx[k];

                g_0_y_xxxzzz_xy[k] = -g_0_y_xxzzz_xy[k] * ab_x + g_0_y_xxzzz_xxy[k];

                g_0_y_xxxzzz_xz[k] = -g_0_y_xxzzz_xz[k] * ab_x + g_0_y_xxzzz_xxz[k];

                g_0_y_xxxzzz_yy[k] = -g_0_y_xxzzz_yy[k] * ab_x + g_0_y_xxzzz_xyy[k];

                g_0_y_xxxzzz_yz[k] = -g_0_y_xxzzz_yz[k] * ab_x + g_0_y_xxzzz_xyz[k];

                g_0_y_xxxzzz_zz[k] = -g_0_y_xxzzz_zz[k] * ab_x + g_0_y_xxzzz_xzz[k];
            }

            /// Set up 228-234 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyyy_xx = cbuffer.data(id_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xy = cbuffer.data(id_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xz = cbuffer.data(id_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_y_xxyyyy_yy = cbuffer.data(id_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_y_xxyyyy_yz = cbuffer.data(id_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_y_xxyyyy_zz = cbuffer.data(id_geom_01_off + 233 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyyy_xx, g_0_y_xxyyyy_xy, g_0_y_xxyyyy_xz, g_0_y_xxyyyy_yy, g_0_y_xxyyyy_yz, g_0_y_xxyyyy_zz, g_0_y_xyyyy_xx, g_0_y_xyyyy_xxx, g_0_y_xyyyy_xxy, g_0_y_xyyyy_xxz, g_0_y_xyyyy_xy, g_0_y_xyyyy_xyy, g_0_y_xyyyy_xyz, g_0_y_xyyyy_xz, g_0_y_xyyyy_xzz, g_0_y_xyyyy_yy, g_0_y_xyyyy_yz, g_0_y_xyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyyy_xx[k] = -g_0_y_xyyyy_xx[k] * ab_x + g_0_y_xyyyy_xxx[k];

                g_0_y_xxyyyy_xy[k] = -g_0_y_xyyyy_xy[k] * ab_x + g_0_y_xyyyy_xxy[k];

                g_0_y_xxyyyy_xz[k] = -g_0_y_xyyyy_xz[k] * ab_x + g_0_y_xyyyy_xxz[k];

                g_0_y_xxyyyy_yy[k] = -g_0_y_xyyyy_yy[k] * ab_x + g_0_y_xyyyy_xyy[k];

                g_0_y_xxyyyy_yz[k] = -g_0_y_xyyyy_yz[k] * ab_x + g_0_y_xyyyy_xyz[k];

                g_0_y_xxyyyy_zz[k] = -g_0_y_xyyyy_zz[k] * ab_x + g_0_y_xyyyy_xzz[k];
            }

            /// Set up 234-240 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyyz_xx = cbuffer.data(id_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xy = cbuffer.data(id_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xz = cbuffer.data(id_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_y_xxyyyz_yy = cbuffer.data(id_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_y_xxyyyz_yz = cbuffer.data(id_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_y_xxyyyz_zz = cbuffer.data(id_geom_01_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyyz_xx, g_0_y_xxyyyz_xy, g_0_y_xxyyyz_xz, g_0_y_xxyyyz_yy, g_0_y_xxyyyz_yz, g_0_y_xxyyyz_zz, g_0_y_xyyyz_xx, g_0_y_xyyyz_xxx, g_0_y_xyyyz_xxy, g_0_y_xyyyz_xxz, g_0_y_xyyyz_xy, g_0_y_xyyyz_xyy, g_0_y_xyyyz_xyz, g_0_y_xyyyz_xz, g_0_y_xyyyz_xzz, g_0_y_xyyyz_yy, g_0_y_xyyyz_yz, g_0_y_xyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyyz_xx[k] = -g_0_y_xyyyz_xx[k] * ab_x + g_0_y_xyyyz_xxx[k];

                g_0_y_xxyyyz_xy[k] = -g_0_y_xyyyz_xy[k] * ab_x + g_0_y_xyyyz_xxy[k];

                g_0_y_xxyyyz_xz[k] = -g_0_y_xyyyz_xz[k] * ab_x + g_0_y_xyyyz_xxz[k];

                g_0_y_xxyyyz_yy[k] = -g_0_y_xyyyz_yy[k] * ab_x + g_0_y_xyyyz_xyy[k];

                g_0_y_xxyyyz_yz[k] = -g_0_y_xyyyz_yz[k] * ab_x + g_0_y_xyyyz_xyz[k];

                g_0_y_xxyyyz_zz[k] = -g_0_y_xyyyz_zz[k] * ab_x + g_0_y_xyyyz_xzz[k];
            }

            /// Set up 240-246 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyzz_xx = cbuffer.data(id_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xy = cbuffer.data(id_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xz = cbuffer.data(id_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_y_xxyyzz_yy = cbuffer.data(id_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_y_xxyyzz_yz = cbuffer.data(id_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_y_xxyyzz_zz = cbuffer.data(id_geom_01_off + 245 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyzz_xx, g_0_y_xxyyzz_xy, g_0_y_xxyyzz_xz, g_0_y_xxyyzz_yy, g_0_y_xxyyzz_yz, g_0_y_xxyyzz_zz, g_0_y_xyyzz_xx, g_0_y_xyyzz_xxx, g_0_y_xyyzz_xxy, g_0_y_xyyzz_xxz, g_0_y_xyyzz_xy, g_0_y_xyyzz_xyy, g_0_y_xyyzz_xyz, g_0_y_xyyzz_xz, g_0_y_xyyzz_xzz, g_0_y_xyyzz_yy, g_0_y_xyyzz_yz, g_0_y_xyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyzz_xx[k] = -g_0_y_xyyzz_xx[k] * ab_x + g_0_y_xyyzz_xxx[k];

                g_0_y_xxyyzz_xy[k] = -g_0_y_xyyzz_xy[k] * ab_x + g_0_y_xyyzz_xxy[k];

                g_0_y_xxyyzz_xz[k] = -g_0_y_xyyzz_xz[k] * ab_x + g_0_y_xyyzz_xxz[k];

                g_0_y_xxyyzz_yy[k] = -g_0_y_xyyzz_yy[k] * ab_x + g_0_y_xyyzz_xyy[k];

                g_0_y_xxyyzz_yz[k] = -g_0_y_xyyzz_yz[k] * ab_x + g_0_y_xyyzz_xyz[k];

                g_0_y_xxyyzz_zz[k] = -g_0_y_xyyzz_zz[k] * ab_x + g_0_y_xyyzz_xzz[k];
            }

            /// Set up 246-252 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyzzz_xx = cbuffer.data(id_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xy = cbuffer.data(id_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xz = cbuffer.data(id_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_y_xxyzzz_yy = cbuffer.data(id_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_y_xxyzzz_yz = cbuffer.data(id_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_y_xxyzzz_zz = cbuffer.data(id_geom_01_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyzzz_xx, g_0_y_xxyzzz_xy, g_0_y_xxyzzz_xz, g_0_y_xxyzzz_yy, g_0_y_xxyzzz_yz, g_0_y_xxyzzz_zz, g_0_y_xyzzz_xx, g_0_y_xyzzz_xxx, g_0_y_xyzzz_xxy, g_0_y_xyzzz_xxz, g_0_y_xyzzz_xy, g_0_y_xyzzz_xyy, g_0_y_xyzzz_xyz, g_0_y_xyzzz_xz, g_0_y_xyzzz_xzz, g_0_y_xyzzz_yy, g_0_y_xyzzz_yz, g_0_y_xyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyzzz_xx[k] = -g_0_y_xyzzz_xx[k] * ab_x + g_0_y_xyzzz_xxx[k];

                g_0_y_xxyzzz_xy[k] = -g_0_y_xyzzz_xy[k] * ab_x + g_0_y_xyzzz_xxy[k];

                g_0_y_xxyzzz_xz[k] = -g_0_y_xyzzz_xz[k] * ab_x + g_0_y_xyzzz_xxz[k];

                g_0_y_xxyzzz_yy[k] = -g_0_y_xyzzz_yy[k] * ab_x + g_0_y_xyzzz_xyy[k];

                g_0_y_xxyzzz_yz[k] = -g_0_y_xyzzz_yz[k] * ab_x + g_0_y_xyzzz_xyz[k];

                g_0_y_xxyzzz_zz[k] = -g_0_y_xyzzz_zz[k] * ab_x + g_0_y_xyzzz_xzz[k];
            }

            /// Set up 252-258 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxzzzz_xx = cbuffer.data(id_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xy = cbuffer.data(id_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xz = cbuffer.data(id_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_y_xxzzzz_yy = cbuffer.data(id_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_y_xxzzzz_yz = cbuffer.data(id_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_y_xxzzzz_zz = cbuffer.data(id_geom_01_off + 257 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxzzzz_xx, g_0_y_xxzzzz_xy, g_0_y_xxzzzz_xz, g_0_y_xxzzzz_yy, g_0_y_xxzzzz_yz, g_0_y_xxzzzz_zz, g_0_y_xzzzz_xx, g_0_y_xzzzz_xxx, g_0_y_xzzzz_xxy, g_0_y_xzzzz_xxz, g_0_y_xzzzz_xy, g_0_y_xzzzz_xyy, g_0_y_xzzzz_xyz, g_0_y_xzzzz_xz, g_0_y_xzzzz_xzz, g_0_y_xzzzz_yy, g_0_y_xzzzz_yz, g_0_y_xzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxzzzz_xx[k] = -g_0_y_xzzzz_xx[k] * ab_x + g_0_y_xzzzz_xxx[k];

                g_0_y_xxzzzz_xy[k] = -g_0_y_xzzzz_xy[k] * ab_x + g_0_y_xzzzz_xxy[k];

                g_0_y_xxzzzz_xz[k] = -g_0_y_xzzzz_xz[k] * ab_x + g_0_y_xzzzz_xxz[k];

                g_0_y_xxzzzz_yy[k] = -g_0_y_xzzzz_yy[k] * ab_x + g_0_y_xzzzz_xyy[k];

                g_0_y_xxzzzz_yz[k] = -g_0_y_xzzzz_yz[k] * ab_x + g_0_y_xzzzz_xyz[k];

                g_0_y_xxzzzz_zz[k] = -g_0_y_xzzzz_zz[k] * ab_x + g_0_y_xzzzz_xzz[k];
            }

            /// Set up 258-264 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyyy_xx = cbuffer.data(id_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xy = cbuffer.data(id_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xz = cbuffer.data(id_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_y_xyyyyy_yy = cbuffer.data(id_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_y_xyyyyy_yz = cbuffer.data(id_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_y_xyyyyy_zz = cbuffer.data(id_geom_01_off + 263 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyyy_xx, g_0_y_xyyyyy_xy, g_0_y_xyyyyy_xz, g_0_y_xyyyyy_yy, g_0_y_xyyyyy_yz, g_0_y_xyyyyy_zz, g_0_y_yyyyy_xx, g_0_y_yyyyy_xxx, g_0_y_yyyyy_xxy, g_0_y_yyyyy_xxz, g_0_y_yyyyy_xy, g_0_y_yyyyy_xyy, g_0_y_yyyyy_xyz, g_0_y_yyyyy_xz, g_0_y_yyyyy_xzz, g_0_y_yyyyy_yy, g_0_y_yyyyy_yz, g_0_y_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyyy_xx[k] = -g_0_y_yyyyy_xx[k] * ab_x + g_0_y_yyyyy_xxx[k];

                g_0_y_xyyyyy_xy[k] = -g_0_y_yyyyy_xy[k] * ab_x + g_0_y_yyyyy_xxy[k];

                g_0_y_xyyyyy_xz[k] = -g_0_y_yyyyy_xz[k] * ab_x + g_0_y_yyyyy_xxz[k];

                g_0_y_xyyyyy_yy[k] = -g_0_y_yyyyy_yy[k] * ab_x + g_0_y_yyyyy_xyy[k];

                g_0_y_xyyyyy_yz[k] = -g_0_y_yyyyy_yz[k] * ab_x + g_0_y_yyyyy_xyz[k];

                g_0_y_xyyyyy_zz[k] = -g_0_y_yyyyy_zz[k] * ab_x + g_0_y_yyyyy_xzz[k];
            }

            /// Set up 264-270 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyyz_xx = cbuffer.data(id_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xy = cbuffer.data(id_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xz = cbuffer.data(id_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_y_xyyyyz_yy = cbuffer.data(id_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_y_xyyyyz_yz = cbuffer.data(id_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_y_xyyyyz_zz = cbuffer.data(id_geom_01_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyyz_xx, g_0_y_xyyyyz_xy, g_0_y_xyyyyz_xz, g_0_y_xyyyyz_yy, g_0_y_xyyyyz_yz, g_0_y_xyyyyz_zz, g_0_y_yyyyz_xx, g_0_y_yyyyz_xxx, g_0_y_yyyyz_xxy, g_0_y_yyyyz_xxz, g_0_y_yyyyz_xy, g_0_y_yyyyz_xyy, g_0_y_yyyyz_xyz, g_0_y_yyyyz_xz, g_0_y_yyyyz_xzz, g_0_y_yyyyz_yy, g_0_y_yyyyz_yz, g_0_y_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyyz_xx[k] = -g_0_y_yyyyz_xx[k] * ab_x + g_0_y_yyyyz_xxx[k];

                g_0_y_xyyyyz_xy[k] = -g_0_y_yyyyz_xy[k] * ab_x + g_0_y_yyyyz_xxy[k];

                g_0_y_xyyyyz_xz[k] = -g_0_y_yyyyz_xz[k] * ab_x + g_0_y_yyyyz_xxz[k];

                g_0_y_xyyyyz_yy[k] = -g_0_y_yyyyz_yy[k] * ab_x + g_0_y_yyyyz_xyy[k];

                g_0_y_xyyyyz_yz[k] = -g_0_y_yyyyz_yz[k] * ab_x + g_0_y_yyyyz_xyz[k];

                g_0_y_xyyyyz_zz[k] = -g_0_y_yyyyz_zz[k] * ab_x + g_0_y_yyyyz_xzz[k];
            }

            /// Set up 270-276 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyzz_xx = cbuffer.data(id_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xy = cbuffer.data(id_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xz = cbuffer.data(id_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_y_xyyyzz_yy = cbuffer.data(id_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_y_xyyyzz_yz = cbuffer.data(id_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_y_xyyyzz_zz = cbuffer.data(id_geom_01_off + 275 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyzz_xx, g_0_y_xyyyzz_xy, g_0_y_xyyyzz_xz, g_0_y_xyyyzz_yy, g_0_y_xyyyzz_yz, g_0_y_xyyyzz_zz, g_0_y_yyyzz_xx, g_0_y_yyyzz_xxx, g_0_y_yyyzz_xxy, g_0_y_yyyzz_xxz, g_0_y_yyyzz_xy, g_0_y_yyyzz_xyy, g_0_y_yyyzz_xyz, g_0_y_yyyzz_xz, g_0_y_yyyzz_xzz, g_0_y_yyyzz_yy, g_0_y_yyyzz_yz, g_0_y_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyzz_xx[k] = -g_0_y_yyyzz_xx[k] * ab_x + g_0_y_yyyzz_xxx[k];

                g_0_y_xyyyzz_xy[k] = -g_0_y_yyyzz_xy[k] * ab_x + g_0_y_yyyzz_xxy[k];

                g_0_y_xyyyzz_xz[k] = -g_0_y_yyyzz_xz[k] * ab_x + g_0_y_yyyzz_xxz[k];

                g_0_y_xyyyzz_yy[k] = -g_0_y_yyyzz_yy[k] * ab_x + g_0_y_yyyzz_xyy[k];

                g_0_y_xyyyzz_yz[k] = -g_0_y_yyyzz_yz[k] * ab_x + g_0_y_yyyzz_xyz[k];

                g_0_y_xyyyzz_zz[k] = -g_0_y_yyyzz_zz[k] * ab_x + g_0_y_yyyzz_xzz[k];
            }

            /// Set up 276-282 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyzzz_xx = cbuffer.data(id_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xy = cbuffer.data(id_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xz = cbuffer.data(id_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_y_xyyzzz_yy = cbuffer.data(id_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_y_xyyzzz_yz = cbuffer.data(id_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_y_xyyzzz_zz = cbuffer.data(id_geom_01_off + 281 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyzzz_xx, g_0_y_xyyzzz_xy, g_0_y_xyyzzz_xz, g_0_y_xyyzzz_yy, g_0_y_xyyzzz_yz, g_0_y_xyyzzz_zz, g_0_y_yyzzz_xx, g_0_y_yyzzz_xxx, g_0_y_yyzzz_xxy, g_0_y_yyzzz_xxz, g_0_y_yyzzz_xy, g_0_y_yyzzz_xyy, g_0_y_yyzzz_xyz, g_0_y_yyzzz_xz, g_0_y_yyzzz_xzz, g_0_y_yyzzz_yy, g_0_y_yyzzz_yz, g_0_y_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyzzz_xx[k] = -g_0_y_yyzzz_xx[k] * ab_x + g_0_y_yyzzz_xxx[k];

                g_0_y_xyyzzz_xy[k] = -g_0_y_yyzzz_xy[k] * ab_x + g_0_y_yyzzz_xxy[k];

                g_0_y_xyyzzz_xz[k] = -g_0_y_yyzzz_xz[k] * ab_x + g_0_y_yyzzz_xxz[k];

                g_0_y_xyyzzz_yy[k] = -g_0_y_yyzzz_yy[k] * ab_x + g_0_y_yyzzz_xyy[k];

                g_0_y_xyyzzz_yz[k] = -g_0_y_yyzzz_yz[k] * ab_x + g_0_y_yyzzz_xyz[k];

                g_0_y_xyyzzz_zz[k] = -g_0_y_yyzzz_zz[k] * ab_x + g_0_y_yyzzz_xzz[k];
            }

            /// Set up 282-288 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyzzzz_xx = cbuffer.data(id_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xy = cbuffer.data(id_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xz = cbuffer.data(id_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_y_xyzzzz_yy = cbuffer.data(id_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_y_xyzzzz_yz = cbuffer.data(id_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_y_xyzzzz_zz = cbuffer.data(id_geom_01_off + 287 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyzzzz_xx, g_0_y_xyzzzz_xy, g_0_y_xyzzzz_xz, g_0_y_xyzzzz_yy, g_0_y_xyzzzz_yz, g_0_y_xyzzzz_zz, g_0_y_yzzzz_xx, g_0_y_yzzzz_xxx, g_0_y_yzzzz_xxy, g_0_y_yzzzz_xxz, g_0_y_yzzzz_xy, g_0_y_yzzzz_xyy, g_0_y_yzzzz_xyz, g_0_y_yzzzz_xz, g_0_y_yzzzz_xzz, g_0_y_yzzzz_yy, g_0_y_yzzzz_yz, g_0_y_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyzzzz_xx[k] = -g_0_y_yzzzz_xx[k] * ab_x + g_0_y_yzzzz_xxx[k];

                g_0_y_xyzzzz_xy[k] = -g_0_y_yzzzz_xy[k] * ab_x + g_0_y_yzzzz_xxy[k];

                g_0_y_xyzzzz_xz[k] = -g_0_y_yzzzz_xz[k] * ab_x + g_0_y_yzzzz_xxz[k];

                g_0_y_xyzzzz_yy[k] = -g_0_y_yzzzz_yy[k] * ab_x + g_0_y_yzzzz_xyy[k];

                g_0_y_xyzzzz_yz[k] = -g_0_y_yzzzz_yz[k] * ab_x + g_0_y_yzzzz_xyz[k];

                g_0_y_xyzzzz_zz[k] = -g_0_y_yzzzz_zz[k] * ab_x + g_0_y_yzzzz_xzz[k];
            }

            /// Set up 288-294 components of targeted buffer : cbuffer.data(

            auto g_0_y_xzzzzz_xx = cbuffer.data(id_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xy = cbuffer.data(id_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xz = cbuffer.data(id_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_y_xzzzzz_yy = cbuffer.data(id_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_y_xzzzzz_yz = cbuffer.data(id_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_y_xzzzzz_zz = cbuffer.data(id_geom_01_off + 293 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xzzzzz_xx, g_0_y_xzzzzz_xy, g_0_y_xzzzzz_xz, g_0_y_xzzzzz_yy, g_0_y_xzzzzz_yz, g_0_y_xzzzzz_zz, g_0_y_zzzzz_xx, g_0_y_zzzzz_xxx, g_0_y_zzzzz_xxy, g_0_y_zzzzz_xxz, g_0_y_zzzzz_xy, g_0_y_zzzzz_xyy, g_0_y_zzzzz_xyz, g_0_y_zzzzz_xz, g_0_y_zzzzz_xzz, g_0_y_zzzzz_yy, g_0_y_zzzzz_yz, g_0_y_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xzzzzz_xx[k] = -g_0_y_zzzzz_xx[k] * ab_x + g_0_y_zzzzz_xxx[k];

                g_0_y_xzzzzz_xy[k] = -g_0_y_zzzzz_xy[k] * ab_x + g_0_y_zzzzz_xxy[k];

                g_0_y_xzzzzz_xz[k] = -g_0_y_zzzzz_xz[k] * ab_x + g_0_y_zzzzz_xxz[k];

                g_0_y_xzzzzz_yy[k] = -g_0_y_zzzzz_yy[k] * ab_x + g_0_y_zzzzz_xyy[k];

                g_0_y_xzzzzz_yz[k] = -g_0_y_zzzzz_yz[k] * ab_x + g_0_y_zzzzz_xyz[k];

                g_0_y_xzzzzz_zz[k] = -g_0_y_zzzzz_zz[k] * ab_x + g_0_y_zzzzz_xzz[k];
            }

            /// Set up 294-300 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyyy_xx = cbuffer.data(id_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xy = cbuffer.data(id_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xz = cbuffer.data(id_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_y_yyyyyy_yy = cbuffer.data(id_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_y_yyyyyy_yz = cbuffer.data(id_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_y_yyyyyy_zz = cbuffer.data(id_geom_01_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyyy_xx, g_0_y_yyyyy_xxy, g_0_y_yyyyy_xy, g_0_y_yyyyy_xyy, g_0_y_yyyyy_xyz, g_0_y_yyyyy_xz, g_0_y_yyyyy_yy, g_0_y_yyyyy_yyy, g_0_y_yyyyy_yyz, g_0_y_yyyyy_yz, g_0_y_yyyyy_yzz, g_0_y_yyyyy_zz, g_0_y_yyyyyy_xx, g_0_y_yyyyyy_xy, g_0_y_yyyyyy_xz, g_0_y_yyyyyy_yy, g_0_y_yyyyyy_yz, g_0_y_yyyyyy_zz, g_yyyyy_xx, g_yyyyy_xy, g_yyyyy_xz, g_yyyyy_yy, g_yyyyy_yz, g_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyyy_xx[k] = g_yyyyy_xx[k] - g_0_y_yyyyy_xx[k] * ab_y + g_0_y_yyyyy_xxy[k];

                g_0_y_yyyyyy_xy[k] = g_yyyyy_xy[k] - g_0_y_yyyyy_xy[k] * ab_y + g_0_y_yyyyy_xyy[k];

                g_0_y_yyyyyy_xz[k] = g_yyyyy_xz[k] - g_0_y_yyyyy_xz[k] * ab_y + g_0_y_yyyyy_xyz[k];

                g_0_y_yyyyyy_yy[k] = g_yyyyy_yy[k] - g_0_y_yyyyy_yy[k] * ab_y + g_0_y_yyyyy_yyy[k];

                g_0_y_yyyyyy_yz[k] = g_yyyyy_yz[k] - g_0_y_yyyyy_yz[k] * ab_y + g_0_y_yyyyy_yyz[k];

                g_0_y_yyyyyy_zz[k] = g_yyyyy_zz[k] - g_0_y_yyyyy_zz[k] * ab_y + g_0_y_yyyyy_yzz[k];
            }

            /// Set up 300-306 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyyz_xx = cbuffer.data(id_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xy = cbuffer.data(id_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xz = cbuffer.data(id_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_y_yyyyyz_yy = cbuffer.data(id_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_y_yyyyyz_yz = cbuffer.data(id_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_y_yyyyyz_zz = cbuffer.data(id_geom_01_off + 305 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyyy_xx, g_0_y_yyyyy_xxz, g_0_y_yyyyy_xy, g_0_y_yyyyy_xyz, g_0_y_yyyyy_xz, g_0_y_yyyyy_xzz, g_0_y_yyyyy_yy, g_0_y_yyyyy_yyz, g_0_y_yyyyy_yz, g_0_y_yyyyy_yzz, g_0_y_yyyyy_zz, g_0_y_yyyyy_zzz, g_0_y_yyyyyz_xx, g_0_y_yyyyyz_xy, g_0_y_yyyyyz_xz, g_0_y_yyyyyz_yy, g_0_y_yyyyyz_yz, g_0_y_yyyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyyz_xx[k] = -g_0_y_yyyyy_xx[k] * ab_z + g_0_y_yyyyy_xxz[k];

                g_0_y_yyyyyz_xy[k] = -g_0_y_yyyyy_xy[k] * ab_z + g_0_y_yyyyy_xyz[k];

                g_0_y_yyyyyz_xz[k] = -g_0_y_yyyyy_xz[k] * ab_z + g_0_y_yyyyy_xzz[k];

                g_0_y_yyyyyz_yy[k] = -g_0_y_yyyyy_yy[k] * ab_z + g_0_y_yyyyy_yyz[k];

                g_0_y_yyyyyz_yz[k] = -g_0_y_yyyyy_yz[k] * ab_z + g_0_y_yyyyy_yzz[k];

                g_0_y_yyyyyz_zz[k] = -g_0_y_yyyyy_zz[k] * ab_z + g_0_y_yyyyy_zzz[k];
            }

            /// Set up 306-312 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyzz_xx = cbuffer.data(id_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xy = cbuffer.data(id_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xz = cbuffer.data(id_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_y_yyyyzz_yy = cbuffer.data(id_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_y_yyyyzz_yz = cbuffer.data(id_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_y_yyyyzz_zz = cbuffer.data(id_geom_01_off + 311 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyyz_xx, g_0_y_yyyyz_xxz, g_0_y_yyyyz_xy, g_0_y_yyyyz_xyz, g_0_y_yyyyz_xz, g_0_y_yyyyz_xzz, g_0_y_yyyyz_yy, g_0_y_yyyyz_yyz, g_0_y_yyyyz_yz, g_0_y_yyyyz_yzz, g_0_y_yyyyz_zz, g_0_y_yyyyz_zzz, g_0_y_yyyyzz_xx, g_0_y_yyyyzz_xy, g_0_y_yyyyzz_xz, g_0_y_yyyyzz_yy, g_0_y_yyyyzz_yz, g_0_y_yyyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyzz_xx[k] = -g_0_y_yyyyz_xx[k] * ab_z + g_0_y_yyyyz_xxz[k];

                g_0_y_yyyyzz_xy[k] = -g_0_y_yyyyz_xy[k] * ab_z + g_0_y_yyyyz_xyz[k];

                g_0_y_yyyyzz_xz[k] = -g_0_y_yyyyz_xz[k] * ab_z + g_0_y_yyyyz_xzz[k];

                g_0_y_yyyyzz_yy[k] = -g_0_y_yyyyz_yy[k] * ab_z + g_0_y_yyyyz_yyz[k];

                g_0_y_yyyyzz_yz[k] = -g_0_y_yyyyz_yz[k] * ab_z + g_0_y_yyyyz_yzz[k];

                g_0_y_yyyyzz_zz[k] = -g_0_y_yyyyz_zz[k] * ab_z + g_0_y_yyyyz_zzz[k];
            }

            /// Set up 312-318 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyzzz_xx = cbuffer.data(id_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xy = cbuffer.data(id_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xz = cbuffer.data(id_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_y_yyyzzz_yy = cbuffer.data(id_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_y_yyyzzz_yz = cbuffer.data(id_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_y_yyyzzz_zz = cbuffer.data(id_geom_01_off + 317 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyzz_xx, g_0_y_yyyzz_xxz, g_0_y_yyyzz_xy, g_0_y_yyyzz_xyz, g_0_y_yyyzz_xz, g_0_y_yyyzz_xzz, g_0_y_yyyzz_yy, g_0_y_yyyzz_yyz, g_0_y_yyyzz_yz, g_0_y_yyyzz_yzz, g_0_y_yyyzz_zz, g_0_y_yyyzz_zzz, g_0_y_yyyzzz_xx, g_0_y_yyyzzz_xy, g_0_y_yyyzzz_xz, g_0_y_yyyzzz_yy, g_0_y_yyyzzz_yz, g_0_y_yyyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyzzz_xx[k] = -g_0_y_yyyzz_xx[k] * ab_z + g_0_y_yyyzz_xxz[k];

                g_0_y_yyyzzz_xy[k] = -g_0_y_yyyzz_xy[k] * ab_z + g_0_y_yyyzz_xyz[k];

                g_0_y_yyyzzz_xz[k] = -g_0_y_yyyzz_xz[k] * ab_z + g_0_y_yyyzz_xzz[k];

                g_0_y_yyyzzz_yy[k] = -g_0_y_yyyzz_yy[k] * ab_z + g_0_y_yyyzz_yyz[k];

                g_0_y_yyyzzz_yz[k] = -g_0_y_yyyzz_yz[k] * ab_z + g_0_y_yyyzz_yzz[k];

                g_0_y_yyyzzz_zz[k] = -g_0_y_yyyzz_zz[k] * ab_z + g_0_y_yyyzz_zzz[k];
            }

            /// Set up 318-324 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyzzzz_xx = cbuffer.data(id_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xy = cbuffer.data(id_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xz = cbuffer.data(id_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_y_yyzzzz_yy = cbuffer.data(id_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_y_yyzzzz_yz = cbuffer.data(id_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_y_yyzzzz_zz = cbuffer.data(id_geom_01_off + 323 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyzzz_xx, g_0_y_yyzzz_xxz, g_0_y_yyzzz_xy, g_0_y_yyzzz_xyz, g_0_y_yyzzz_xz, g_0_y_yyzzz_xzz, g_0_y_yyzzz_yy, g_0_y_yyzzz_yyz, g_0_y_yyzzz_yz, g_0_y_yyzzz_yzz, g_0_y_yyzzz_zz, g_0_y_yyzzz_zzz, g_0_y_yyzzzz_xx, g_0_y_yyzzzz_xy, g_0_y_yyzzzz_xz, g_0_y_yyzzzz_yy, g_0_y_yyzzzz_yz, g_0_y_yyzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyzzzz_xx[k] = -g_0_y_yyzzz_xx[k] * ab_z + g_0_y_yyzzz_xxz[k];

                g_0_y_yyzzzz_xy[k] = -g_0_y_yyzzz_xy[k] * ab_z + g_0_y_yyzzz_xyz[k];

                g_0_y_yyzzzz_xz[k] = -g_0_y_yyzzz_xz[k] * ab_z + g_0_y_yyzzz_xzz[k];

                g_0_y_yyzzzz_yy[k] = -g_0_y_yyzzz_yy[k] * ab_z + g_0_y_yyzzz_yyz[k];

                g_0_y_yyzzzz_yz[k] = -g_0_y_yyzzz_yz[k] * ab_z + g_0_y_yyzzz_yzz[k];

                g_0_y_yyzzzz_zz[k] = -g_0_y_yyzzz_zz[k] * ab_z + g_0_y_yyzzz_zzz[k];
            }

            /// Set up 324-330 components of targeted buffer : cbuffer.data(

            auto g_0_y_yzzzzz_xx = cbuffer.data(id_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xy = cbuffer.data(id_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xz = cbuffer.data(id_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_y_yzzzzz_yy = cbuffer.data(id_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_y_yzzzzz_yz = cbuffer.data(id_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_y_yzzzzz_zz = cbuffer.data(id_geom_01_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yzzzz_xx, g_0_y_yzzzz_xxz, g_0_y_yzzzz_xy, g_0_y_yzzzz_xyz, g_0_y_yzzzz_xz, g_0_y_yzzzz_xzz, g_0_y_yzzzz_yy, g_0_y_yzzzz_yyz, g_0_y_yzzzz_yz, g_0_y_yzzzz_yzz, g_0_y_yzzzz_zz, g_0_y_yzzzz_zzz, g_0_y_yzzzzz_xx, g_0_y_yzzzzz_xy, g_0_y_yzzzzz_xz, g_0_y_yzzzzz_yy, g_0_y_yzzzzz_yz, g_0_y_yzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yzzzzz_xx[k] = -g_0_y_yzzzz_xx[k] * ab_z + g_0_y_yzzzz_xxz[k];

                g_0_y_yzzzzz_xy[k] = -g_0_y_yzzzz_xy[k] * ab_z + g_0_y_yzzzz_xyz[k];

                g_0_y_yzzzzz_xz[k] = -g_0_y_yzzzz_xz[k] * ab_z + g_0_y_yzzzz_xzz[k];

                g_0_y_yzzzzz_yy[k] = -g_0_y_yzzzz_yy[k] * ab_z + g_0_y_yzzzz_yyz[k];

                g_0_y_yzzzzz_yz[k] = -g_0_y_yzzzz_yz[k] * ab_z + g_0_y_yzzzz_yzz[k];

                g_0_y_yzzzzz_zz[k] = -g_0_y_yzzzz_zz[k] * ab_z + g_0_y_yzzzz_zzz[k];
            }

            /// Set up 330-336 components of targeted buffer : cbuffer.data(

            auto g_0_y_zzzzzz_xx = cbuffer.data(id_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xy = cbuffer.data(id_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xz = cbuffer.data(id_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_y_zzzzzz_yy = cbuffer.data(id_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_y_zzzzzz_yz = cbuffer.data(id_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_y_zzzzzz_zz = cbuffer.data(id_geom_01_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_zzzzz_xx, g_0_y_zzzzz_xxz, g_0_y_zzzzz_xy, g_0_y_zzzzz_xyz, g_0_y_zzzzz_xz, g_0_y_zzzzz_xzz, g_0_y_zzzzz_yy, g_0_y_zzzzz_yyz, g_0_y_zzzzz_yz, g_0_y_zzzzz_yzz, g_0_y_zzzzz_zz, g_0_y_zzzzz_zzz, g_0_y_zzzzzz_xx, g_0_y_zzzzzz_xy, g_0_y_zzzzzz_xz, g_0_y_zzzzzz_yy, g_0_y_zzzzzz_yz, g_0_y_zzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_zzzzzz_xx[k] = -g_0_y_zzzzz_xx[k] * ab_z + g_0_y_zzzzz_xxz[k];

                g_0_y_zzzzzz_xy[k] = -g_0_y_zzzzz_xy[k] * ab_z + g_0_y_zzzzz_xyz[k];

                g_0_y_zzzzzz_xz[k] = -g_0_y_zzzzz_xz[k] * ab_z + g_0_y_zzzzz_xzz[k];

                g_0_y_zzzzzz_yy[k] = -g_0_y_zzzzz_yy[k] * ab_z + g_0_y_zzzzz_yyz[k];

                g_0_y_zzzzzz_yz[k] = -g_0_y_zzzzz_yz[k] * ab_z + g_0_y_zzzzz_yzz[k];

                g_0_y_zzzzzz_zz[k] = -g_0_y_zzzzz_zz[k] * ab_z + g_0_y_zzzzz_zzz[k];
            }

            /// Set up 336-342 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxx_xx = cbuffer.data(id_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xy = cbuffer.data(id_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xz = cbuffer.data(id_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_z_xxxxxx_yy = cbuffer.data(id_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_z_xxxxxx_yz = cbuffer.data(id_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_z_xxxxxx_zz = cbuffer.data(id_geom_01_off + 341 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxx_xx, g_0_z_xxxxx_xxx, g_0_z_xxxxx_xxy, g_0_z_xxxxx_xxz, g_0_z_xxxxx_xy, g_0_z_xxxxx_xyy, g_0_z_xxxxx_xyz, g_0_z_xxxxx_xz, g_0_z_xxxxx_xzz, g_0_z_xxxxx_yy, g_0_z_xxxxx_yz, g_0_z_xxxxx_zz, g_0_z_xxxxxx_xx, g_0_z_xxxxxx_xy, g_0_z_xxxxxx_xz, g_0_z_xxxxxx_yy, g_0_z_xxxxxx_yz, g_0_z_xxxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxx_xx[k] = -g_0_z_xxxxx_xx[k] * ab_x + g_0_z_xxxxx_xxx[k];

                g_0_z_xxxxxx_xy[k] = -g_0_z_xxxxx_xy[k] * ab_x + g_0_z_xxxxx_xxy[k];

                g_0_z_xxxxxx_xz[k] = -g_0_z_xxxxx_xz[k] * ab_x + g_0_z_xxxxx_xxz[k];

                g_0_z_xxxxxx_yy[k] = -g_0_z_xxxxx_yy[k] * ab_x + g_0_z_xxxxx_xyy[k];

                g_0_z_xxxxxx_yz[k] = -g_0_z_xxxxx_yz[k] * ab_x + g_0_z_xxxxx_xyz[k];

                g_0_z_xxxxxx_zz[k] = -g_0_z_xxxxx_zz[k] * ab_x + g_0_z_xxxxx_xzz[k];
            }

            /// Set up 342-348 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxy_xx = cbuffer.data(id_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xy = cbuffer.data(id_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xz = cbuffer.data(id_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_z_xxxxxy_yy = cbuffer.data(id_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_z_xxxxxy_yz = cbuffer.data(id_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_z_xxxxxy_zz = cbuffer.data(id_geom_01_off + 347 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxy_xx, g_0_z_xxxxxy_xy, g_0_z_xxxxxy_xz, g_0_z_xxxxxy_yy, g_0_z_xxxxxy_yz, g_0_z_xxxxxy_zz, g_0_z_xxxxy_xx, g_0_z_xxxxy_xxx, g_0_z_xxxxy_xxy, g_0_z_xxxxy_xxz, g_0_z_xxxxy_xy, g_0_z_xxxxy_xyy, g_0_z_xxxxy_xyz, g_0_z_xxxxy_xz, g_0_z_xxxxy_xzz, g_0_z_xxxxy_yy, g_0_z_xxxxy_yz, g_0_z_xxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxy_xx[k] = -g_0_z_xxxxy_xx[k] * ab_x + g_0_z_xxxxy_xxx[k];

                g_0_z_xxxxxy_xy[k] = -g_0_z_xxxxy_xy[k] * ab_x + g_0_z_xxxxy_xxy[k];

                g_0_z_xxxxxy_xz[k] = -g_0_z_xxxxy_xz[k] * ab_x + g_0_z_xxxxy_xxz[k];

                g_0_z_xxxxxy_yy[k] = -g_0_z_xxxxy_yy[k] * ab_x + g_0_z_xxxxy_xyy[k];

                g_0_z_xxxxxy_yz[k] = -g_0_z_xxxxy_yz[k] * ab_x + g_0_z_xxxxy_xyz[k];

                g_0_z_xxxxxy_zz[k] = -g_0_z_xxxxy_zz[k] * ab_x + g_0_z_xxxxy_xzz[k];
            }

            /// Set up 348-354 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxz_xx = cbuffer.data(id_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xy = cbuffer.data(id_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xz = cbuffer.data(id_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_z_xxxxxz_yy = cbuffer.data(id_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_z_xxxxxz_yz = cbuffer.data(id_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_z_xxxxxz_zz = cbuffer.data(id_geom_01_off + 353 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxz_xx, g_0_z_xxxxxz_xy, g_0_z_xxxxxz_xz, g_0_z_xxxxxz_yy, g_0_z_xxxxxz_yz, g_0_z_xxxxxz_zz, g_0_z_xxxxz_xx, g_0_z_xxxxz_xxx, g_0_z_xxxxz_xxy, g_0_z_xxxxz_xxz, g_0_z_xxxxz_xy, g_0_z_xxxxz_xyy, g_0_z_xxxxz_xyz, g_0_z_xxxxz_xz, g_0_z_xxxxz_xzz, g_0_z_xxxxz_yy, g_0_z_xxxxz_yz, g_0_z_xxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxz_xx[k] = -g_0_z_xxxxz_xx[k] * ab_x + g_0_z_xxxxz_xxx[k];

                g_0_z_xxxxxz_xy[k] = -g_0_z_xxxxz_xy[k] * ab_x + g_0_z_xxxxz_xxy[k];

                g_0_z_xxxxxz_xz[k] = -g_0_z_xxxxz_xz[k] * ab_x + g_0_z_xxxxz_xxz[k];

                g_0_z_xxxxxz_yy[k] = -g_0_z_xxxxz_yy[k] * ab_x + g_0_z_xxxxz_xyy[k];

                g_0_z_xxxxxz_yz[k] = -g_0_z_xxxxz_yz[k] * ab_x + g_0_z_xxxxz_xyz[k];

                g_0_z_xxxxxz_zz[k] = -g_0_z_xxxxz_zz[k] * ab_x + g_0_z_xxxxz_xzz[k];
            }

            /// Set up 354-360 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxyy_xx = cbuffer.data(id_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xy = cbuffer.data(id_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xz = cbuffer.data(id_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_z_xxxxyy_yy = cbuffer.data(id_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_z_xxxxyy_yz = cbuffer.data(id_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_z_xxxxyy_zz = cbuffer.data(id_geom_01_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxyy_xx, g_0_z_xxxxyy_xy, g_0_z_xxxxyy_xz, g_0_z_xxxxyy_yy, g_0_z_xxxxyy_yz, g_0_z_xxxxyy_zz, g_0_z_xxxyy_xx, g_0_z_xxxyy_xxx, g_0_z_xxxyy_xxy, g_0_z_xxxyy_xxz, g_0_z_xxxyy_xy, g_0_z_xxxyy_xyy, g_0_z_xxxyy_xyz, g_0_z_xxxyy_xz, g_0_z_xxxyy_xzz, g_0_z_xxxyy_yy, g_0_z_xxxyy_yz, g_0_z_xxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxyy_xx[k] = -g_0_z_xxxyy_xx[k] * ab_x + g_0_z_xxxyy_xxx[k];

                g_0_z_xxxxyy_xy[k] = -g_0_z_xxxyy_xy[k] * ab_x + g_0_z_xxxyy_xxy[k];

                g_0_z_xxxxyy_xz[k] = -g_0_z_xxxyy_xz[k] * ab_x + g_0_z_xxxyy_xxz[k];

                g_0_z_xxxxyy_yy[k] = -g_0_z_xxxyy_yy[k] * ab_x + g_0_z_xxxyy_xyy[k];

                g_0_z_xxxxyy_yz[k] = -g_0_z_xxxyy_yz[k] * ab_x + g_0_z_xxxyy_xyz[k];

                g_0_z_xxxxyy_zz[k] = -g_0_z_xxxyy_zz[k] * ab_x + g_0_z_xxxyy_xzz[k];
            }

            /// Set up 360-366 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxyz_xx = cbuffer.data(id_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xy = cbuffer.data(id_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xz = cbuffer.data(id_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_z_xxxxyz_yy = cbuffer.data(id_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_z_xxxxyz_yz = cbuffer.data(id_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_z_xxxxyz_zz = cbuffer.data(id_geom_01_off + 365 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxyz_xx, g_0_z_xxxxyz_xy, g_0_z_xxxxyz_xz, g_0_z_xxxxyz_yy, g_0_z_xxxxyz_yz, g_0_z_xxxxyz_zz, g_0_z_xxxyz_xx, g_0_z_xxxyz_xxx, g_0_z_xxxyz_xxy, g_0_z_xxxyz_xxz, g_0_z_xxxyz_xy, g_0_z_xxxyz_xyy, g_0_z_xxxyz_xyz, g_0_z_xxxyz_xz, g_0_z_xxxyz_xzz, g_0_z_xxxyz_yy, g_0_z_xxxyz_yz, g_0_z_xxxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxyz_xx[k] = -g_0_z_xxxyz_xx[k] * ab_x + g_0_z_xxxyz_xxx[k];

                g_0_z_xxxxyz_xy[k] = -g_0_z_xxxyz_xy[k] * ab_x + g_0_z_xxxyz_xxy[k];

                g_0_z_xxxxyz_xz[k] = -g_0_z_xxxyz_xz[k] * ab_x + g_0_z_xxxyz_xxz[k];

                g_0_z_xxxxyz_yy[k] = -g_0_z_xxxyz_yy[k] * ab_x + g_0_z_xxxyz_xyy[k];

                g_0_z_xxxxyz_yz[k] = -g_0_z_xxxyz_yz[k] * ab_x + g_0_z_xxxyz_xyz[k];

                g_0_z_xxxxyz_zz[k] = -g_0_z_xxxyz_zz[k] * ab_x + g_0_z_xxxyz_xzz[k];
            }

            /// Set up 366-372 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxzz_xx = cbuffer.data(id_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xy = cbuffer.data(id_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xz = cbuffer.data(id_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_z_xxxxzz_yy = cbuffer.data(id_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_z_xxxxzz_yz = cbuffer.data(id_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_z_xxxxzz_zz = cbuffer.data(id_geom_01_off + 371 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxzz_xx, g_0_z_xxxxzz_xy, g_0_z_xxxxzz_xz, g_0_z_xxxxzz_yy, g_0_z_xxxxzz_yz, g_0_z_xxxxzz_zz, g_0_z_xxxzz_xx, g_0_z_xxxzz_xxx, g_0_z_xxxzz_xxy, g_0_z_xxxzz_xxz, g_0_z_xxxzz_xy, g_0_z_xxxzz_xyy, g_0_z_xxxzz_xyz, g_0_z_xxxzz_xz, g_0_z_xxxzz_xzz, g_0_z_xxxzz_yy, g_0_z_xxxzz_yz, g_0_z_xxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxzz_xx[k] = -g_0_z_xxxzz_xx[k] * ab_x + g_0_z_xxxzz_xxx[k];

                g_0_z_xxxxzz_xy[k] = -g_0_z_xxxzz_xy[k] * ab_x + g_0_z_xxxzz_xxy[k];

                g_0_z_xxxxzz_xz[k] = -g_0_z_xxxzz_xz[k] * ab_x + g_0_z_xxxzz_xxz[k];

                g_0_z_xxxxzz_yy[k] = -g_0_z_xxxzz_yy[k] * ab_x + g_0_z_xxxzz_xyy[k];

                g_0_z_xxxxzz_yz[k] = -g_0_z_xxxzz_yz[k] * ab_x + g_0_z_xxxzz_xyz[k];

                g_0_z_xxxxzz_zz[k] = -g_0_z_xxxzz_zz[k] * ab_x + g_0_z_xxxzz_xzz[k];
            }

            /// Set up 372-378 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyyy_xx = cbuffer.data(id_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xy = cbuffer.data(id_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xz = cbuffer.data(id_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_z_xxxyyy_yy = cbuffer.data(id_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_z_xxxyyy_yz = cbuffer.data(id_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_z_xxxyyy_zz = cbuffer.data(id_geom_01_off + 377 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyyy_xx, g_0_z_xxxyyy_xy, g_0_z_xxxyyy_xz, g_0_z_xxxyyy_yy, g_0_z_xxxyyy_yz, g_0_z_xxxyyy_zz, g_0_z_xxyyy_xx, g_0_z_xxyyy_xxx, g_0_z_xxyyy_xxy, g_0_z_xxyyy_xxz, g_0_z_xxyyy_xy, g_0_z_xxyyy_xyy, g_0_z_xxyyy_xyz, g_0_z_xxyyy_xz, g_0_z_xxyyy_xzz, g_0_z_xxyyy_yy, g_0_z_xxyyy_yz, g_0_z_xxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyyy_xx[k] = -g_0_z_xxyyy_xx[k] * ab_x + g_0_z_xxyyy_xxx[k];

                g_0_z_xxxyyy_xy[k] = -g_0_z_xxyyy_xy[k] * ab_x + g_0_z_xxyyy_xxy[k];

                g_0_z_xxxyyy_xz[k] = -g_0_z_xxyyy_xz[k] * ab_x + g_0_z_xxyyy_xxz[k];

                g_0_z_xxxyyy_yy[k] = -g_0_z_xxyyy_yy[k] * ab_x + g_0_z_xxyyy_xyy[k];

                g_0_z_xxxyyy_yz[k] = -g_0_z_xxyyy_yz[k] * ab_x + g_0_z_xxyyy_xyz[k];

                g_0_z_xxxyyy_zz[k] = -g_0_z_xxyyy_zz[k] * ab_x + g_0_z_xxyyy_xzz[k];
            }

            /// Set up 378-384 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyyz_xx = cbuffer.data(id_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xy = cbuffer.data(id_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xz = cbuffer.data(id_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_z_xxxyyz_yy = cbuffer.data(id_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_z_xxxyyz_yz = cbuffer.data(id_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_z_xxxyyz_zz = cbuffer.data(id_geom_01_off + 383 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyyz_xx, g_0_z_xxxyyz_xy, g_0_z_xxxyyz_xz, g_0_z_xxxyyz_yy, g_0_z_xxxyyz_yz, g_0_z_xxxyyz_zz, g_0_z_xxyyz_xx, g_0_z_xxyyz_xxx, g_0_z_xxyyz_xxy, g_0_z_xxyyz_xxz, g_0_z_xxyyz_xy, g_0_z_xxyyz_xyy, g_0_z_xxyyz_xyz, g_0_z_xxyyz_xz, g_0_z_xxyyz_xzz, g_0_z_xxyyz_yy, g_0_z_xxyyz_yz, g_0_z_xxyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyyz_xx[k] = -g_0_z_xxyyz_xx[k] * ab_x + g_0_z_xxyyz_xxx[k];

                g_0_z_xxxyyz_xy[k] = -g_0_z_xxyyz_xy[k] * ab_x + g_0_z_xxyyz_xxy[k];

                g_0_z_xxxyyz_xz[k] = -g_0_z_xxyyz_xz[k] * ab_x + g_0_z_xxyyz_xxz[k];

                g_0_z_xxxyyz_yy[k] = -g_0_z_xxyyz_yy[k] * ab_x + g_0_z_xxyyz_xyy[k];

                g_0_z_xxxyyz_yz[k] = -g_0_z_xxyyz_yz[k] * ab_x + g_0_z_xxyyz_xyz[k];

                g_0_z_xxxyyz_zz[k] = -g_0_z_xxyyz_zz[k] * ab_x + g_0_z_xxyyz_xzz[k];
            }

            /// Set up 384-390 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyzz_xx = cbuffer.data(id_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xy = cbuffer.data(id_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xz = cbuffer.data(id_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_z_xxxyzz_yy = cbuffer.data(id_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_z_xxxyzz_yz = cbuffer.data(id_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_z_xxxyzz_zz = cbuffer.data(id_geom_01_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyzz_xx, g_0_z_xxxyzz_xy, g_0_z_xxxyzz_xz, g_0_z_xxxyzz_yy, g_0_z_xxxyzz_yz, g_0_z_xxxyzz_zz, g_0_z_xxyzz_xx, g_0_z_xxyzz_xxx, g_0_z_xxyzz_xxy, g_0_z_xxyzz_xxz, g_0_z_xxyzz_xy, g_0_z_xxyzz_xyy, g_0_z_xxyzz_xyz, g_0_z_xxyzz_xz, g_0_z_xxyzz_xzz, g_0_z_xxyzz_yy, g_0_z_xxyzz_yz, g_0_z_xxyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyzz_xx[k] = -g_0_z_xxyzz_xx[k] * ab_x + g_0_z_xxyzz_xxx[k];

                g_0_z_xxxyzz_xy[k] = -g_0_z_xxyzz_xy[k] * ab_x + g_0_z_xxyzz_xxy[k];

                g_0_z_xxxyzz_xz[k] = -g_0_z_xxyzz_xz[k] * ab_x + g_0_z_xxyzz_xxz[k];

                g_0_z_xxxyzz_yy[k] = -g_0_z_xxyzz_yy[k] * ab_x + g_0_z_xxyzz_xyy[k];

                g_0_z_xxxyzz_yz[k] = -g_0_z_xxyzz_yz[k] * ab_x + g_0_z_xxyzz_xyz[k];

                g_0_z_xxxyzz_zz[k] = -g_0_z_xxyzz_zz[k] * ab_x + g_0_z_xxyzz_xzz[k];
            }

            /// Set up 390-396 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxzzz_xx = cbuffer.data(id_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xy = cbuffer.data(id_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xz = cbuffer.data(id_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_z_xxxzzz_yy = cbuffer.data(id_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_z_xxxzzz_yz = cbuffer.data(id_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_z_xxxzzz_zz = cbuffer.data(id_geom_01_off + 395 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxzzz_xx, g_0_z_xxxzzz_xy, g_0_z_xxxzzz_xz, g_0_z_xxxzzz_yy, g_0_z_xxxzzz_yz, g_0_z_xxxzzz_zz, g_0_z_xxzzz_xx, g_0_z_xxzzz_xxx, g_0_z_xxzzz_xxy, g_0_z_xxzzz_xxz, g_0_z_xxzzz_xy, g_0_z_xxzzz_xyy, g_0_z_xxzzz_xyz, g_0_z_xxzzz_xz, g_0_z_xxzzz_xzz, g_0_z_xxzzz_yy, g_0_z_xxzzz_yz, g_0_z_xxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxzzz_xx[k] = -g_0_z_xxzzz_xx[k] * ab_x + g_0_z_xxzzz_xxx[k];

                g_0_z_xxxzzz_xy[k] = -g_0_z_xxzzz_xy[k] * ab_x + g_0_z_xxzzz_xxy[k];

                g_0_z_xxxzzz_xz[k] = -g_0_z_xxzzz_xz[k] * ab_x + g_0_z_xxzzz_xxz[k];

                g_0_z_xxxzzz_yy[k] = -g_0_z_xxzzz_yy[k] * ab_x + g_0_z_xxzzz_xyy[k];

                g_0_z_xxxzzz_yz[k] = -g_0_z_xxzzz_yz[k] * ab_x + g_0_z_xxzzz_xyz[k];

                g_0_z_xxxzzz_zz[k] = -g_0_z_xxzzz_zz[k] * ab_x + g_0_z_xxzzz_xzz[k];
            }

            /// Set up 396-402 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyyy_xx = cbuffer.data(id_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xy = cbuffer.data(id_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xz = cbuffer.data(id_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_z_xxyyyy_yy = cbuffer.data(id_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_z_xxyyyy_yz = cbuffer.data(id_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_z_xxyyyy_zz = cbuffer.data(id_geom_01_off + 401 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyyy_xx, g_0_z_xxyyyy_xy, g_0_z_xxyyyy_xz, g_0_z_xxyyyy_yy, g_0_z_xxyyyy_yz, g_0_z_xxyyyy_zz, g_0_z_xyyyy_xx, g_0_z_xyyyy_xxx, g_0_z_xyyyy_xxy, g_0_z_xyyyy_xxz, g_0_z_xyyyy_xy, g_0_z_xyyyy_xyy, g_0_z_xyyyy_xyz, g_0_z_xyyyy_xz, g_0_z_xyyyy_xzz, g_0_z_xyyyy_yy, g_0_z_xyyyy_yz, g_0_z_xyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyyy_xx[k] = -g_0_z_xyyyy_xx[k] * ab_x + g_0_z_xyyyy_xxx[k];

                g_0_z_xxyyyy_xy[k] = -g_0_z_xyyyy_xy[k] * ab_x + g_0_z_xyyyy_xxy[k];

                g_0_z_xxyyyy_xz[k] = -g_0_z_xyyyy_xz[k] * ab_x + g_0_z_xyyyy_xxz[k];

                g_0_z_xxyyyy_yy[k] = -g_0_z_xyyyy_yy[k] * ab_x + g_0_z_xyyyy_xyy[k];

                g_0_z_xxyyyy_yz[k] = -g_0_z_xyyyy_yz[k] * ab_x + g_0_z_xyyyy_xyz[k];

                g_0_z_xxyyyy_zz[k] = -g_0_z_xyyyy_zz[k] * ab_x + g_0_z_xyyyy_xzz[k];
            }

            /// Set up 402-408 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyyz_xx = cbuffer.data(id_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xy = cbuffer.data(id_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xz = cbuffer.data(id_geom_01_off + 404 * ccomps * dcomps);

            auto g_0_z_xxyyyz_yy = cbuffer.data(id_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_z_xxyyyz_yz = cbuffer.data(id_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_z_xxyyyz_zz = cbuffer.data(id_geom_01_off + 407 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyyz_xx, g_0_z_xxyyyz_xy, g_0_z_xxyyyz_xz, g_0_z_xxyyyz_yy, g_0_z_xxyyyz_yz, g_0_z_xxyyyz_zz, g_0_z_xyyyz_xx, g_0_z_xyyyz_xxx, g_0_z_xyyyz_xxy, g_0_z_xyyyz_xxz, g_0_z_xyyyz_xy, g_0_z_xyyyz_xyy, g_0_z_xyyyz_xyz, g_0_z_xyyyz_xz, g_0_z_xyyyz_xzz, g_0_z_xyyyz_yy, g_0_z_xyyyz_yz, g_0_z_xyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyyz_xx[k] = -g_0_z_xyyyz_xx[k] * ab_x + g_0_z_xyyyz_xxx[k];

                g_0_z_xxyyyz_xy[k] = -g_0_z_xyyyz_xy[k] * ab_x + g_0_z_xyyyz_xxy[k];

                g_0_z_xxyyyz_xz[k] = -g_0_z_xyyyz_xz[k] * ab_x + g_0_z_xyyyz_xxz[k];

                g_0_z_xxyyyz_yy[k] = -g_0_z_xyyyz_yy[k] * ab_x + g_0_z_xyyyz_xyy[k];

                g_0_z_xxyyyz_yz[k] = -g_0_z_xyyyz_yz[k] * ab_x + g_0_z_xyyyz_xyz[k];

                g_0_z_xxyyyz_zz[k] = -g_0_z_xyyyz_zz[k] * ab_x + g_0_z_xyyyz_xzz[k];
            }

            /// Set up 408-414 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyzz_xx = cbuffer.data(id_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xy = cbuffer.data(id_geom_01_off + 409 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xz = cbuffer.data(id_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_z_xxyyzz_yy = cbuffer.data(id_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_z_xxyyzz_yz = cbuffer.data(id_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_z_xxyyzz_zz = cbuffer.data(id_geom_01_off + 413 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyzz_xx, g_0_z_xxyyzz_xy, g_0_z_xxyyzz_xz, g_0_z_xxyyzz_yy, g_0_z_xxyyzz_yz, g_0_z_xxyyzz_zz, g_0_z_xyyzz_xx, g_0_z_xyyzz_xxx, g_0_z_xyyzz_xxy, g_0_z_xyyzz_xxz, g_0_z_xyyzz_xy, g_0_z_xyyzz_xyy, g_0_z_xyyzz_xyz, g_0_z_xyyzz_xz, g_0_z_xyyzz_xzz, g_0_z_xyyzz_yy, g_0_z_xyyzz_yz, g_0_z_xyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyzz_xx[k] = -g_0_z_xyyzz_xx[k] * ab_x + g_0_z_xyyzz_xxx[k];

                g_0_z_xxyyzz_xy[k] = -g_0_z_xyyzz_xy[k] * ab_x + g_0_z_xyyzz_xxy[k];

                g_0_z_xxyyzz_xz[k] = -g_0_z_xyyzz_xz[k] * ab_x + g_0_z_xyyzz_xxz[k];

                g_0_z_xxyyzz_yy[k] = -g_0_z_xyyzz_yy[k] * ab_x + g_0_z_xyyzz_xyy[k];

                g_0_z_xxyyzz_yz[k] = -g_0_z_xyyzz_yz[k] * ab_x + g_0_z_xyyzz_xyz[k];

                g_0_z_xxyyzz_zz[k] = -g_0_z_xyyzz_zz[k] * ab_x + g_0_z_xyyzz_xzz[k];
            }

            /// Set up 414-420 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyzzz_xx = cbuffer.data(id_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xy = cbuffer.data(id_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xz = cbuffer.data(id_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_z_xxyzzz_yy = cbuffer.data(id_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_z_xxyzzz_yz = cbuffer.data(id_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_z_xxyzzz_zz = cbuffer.data(id_geom_01_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyzzz_xx, g_0_z_xxyzzz_xy, g_0_z_xxyzzz_xz, g_0_z_xxyzzz_yy, g_0_z_xxyzzz_yz, g_0_z_xxyzzz_zz, g_0_z_xyzzz_xx, g_0_z_xyzzz_xxx, g_0_z_xyzzz_xxy, g_0_z_xyzzz_xxz, g_0_z_xyzzz_xy, g_0_z_xyzzz_xyy, g_0_z_xyzzz_xyz, g_0_z_xyzzz_xz, g_0_z_xyzzz_xzz, g_0_z_xyzzz_yy, g_0_z_xyzzz_yz, g_0_z_xyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyzzz_xx[k] = -g_0_z_xyzzz_xx[k] * ab_x + g_0_z_xyzzz_xxx[k];

                g_0_z_xxyzzz_xy[k] = -g_0_z_xyzzz_xy[k] * ab_x + g_0_z_xyzzz_xxy[k];

                g_0_z_xxyzzz_xz[k] = -g_0_z_xyzzz_xz[k] * ab_x + g_0_z_xyzzz_xxz[k];

                g_0_z_xxyzzz_yy[k] = -g_0_z_xyzzz_yy[k] * ab_x + g_0_z_xyzzz_xyy[k];

                g_0_z_xxyzzz_yz[k] = -g_0_z_xyzzz_yz[k] * ab_x + g_0_z_xyzzz_xyz[k];

                g_0_z_xxyzzz_zz[k] = -g_0_z_xyzzz_zz[k] * ab_x + g_0_z_xyzzz_xzz[k];
            }

            /// Set up 420-426 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxzzzz_xx = cbuffer.data(id_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xy = cbuffer.data(id_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xz = cbuffer.data(id_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_z_xxzzzz_yy = cbuffer.data(id_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_z_xxzzzz_yz = cbuffer.data(id_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_z_xxzzzz_zz = cbuffer.data(id_geom_01_off + 425 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxzzzz_xx, g_0_z_xxzzzz_xy, g_0_z_xxzzzz_xz, g_0_z_xxzzzz_yy, g_0_z_xxzzzz_yz, g_0_z_xxzzzz_zz, g_0_z_xzzzz_xx, g_0_z_xzzzz_xxx, g_0_z_xzzzz_xxy, g_0_z_xzzzz_xxz, g_0_z_xzzzz_xy, g_0_z_xzzzz_xyy, g_0_z_xzzzz_xyz, g_0_z_xzzzz_xz, g_0_z_xzzzz_xzz, g_0_z_xzzzz_yy, g_0_z_xzzzz_yz, g_0_z_xzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxzzzz_xx[k] = -g_0_z_xzzzz_xx[k] * ab_x + g_0_z_xzzzz_xxx[k];

                g_0_z_xxzzzz_xy[k] = -g_0_z_xzzzz_xy[k] * ab_x + g_0_z_xzzzz_xxy[k];

                g_0_z_xxzzzz_xz[k] = -g_0_z_xzzzz_xz[k] * ab_x + g_0_z_xzzzz_xxz[k];

                g_0_z_xxzzzz_yy[k] = -g_0_z_xzzzz_yy[k] * ab_x + g_0_z_xzzzz_xyy[k];

                g_0_z_xxzzzz_yz[k] = -g_0_z_xzzzz_yz[k] * ab_x + g_0_z_xzzzz_xyz[k];

                g_0_z_xxzzzz_zz[k] = -g_0_z_xzzzz_zz[k] * ab_x + g_0_z_xzzzz_xzz[k];
            }

            /// Set up 426-432 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyyy_xx = cbuffer.data(id_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xy = cbuffer.data(id_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xz = cbuffer.data(id_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_z_xyyyyy_yy = cbuffer.data(id_geom_01_off + 429 * ccomps * dcomps);

            auto g_0_z_xyyyyy_yz = cbuffer.data(id_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_z_xyyyyy_zz = cbuffer.data(id_geom_01_off + 431 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyyy_xx, g_0_z_xyyyyy_xy, g_0_z_xyyyyy_xz, g_0_z_xyyyyy_yy, g_0_z_xyyyyy_yz, g_0_z_xyyyyy_zz, g_0_z_yyyyy_xx, g_0_z_yyyyy_xxx, g_0_z_yyyyy_xxy, g_0_z_yyyyy_xxz, g_0_z_yyyyy_xy, g_0_z_yyyyy_xyy, g_0_z_yyyyy_xyz, g_0_z_yyyyy_xz, g_0_z_yyyyy_xzz, g_0_z_yyyyy_yy, g_0_z_yyyyy_yz, g_0_z_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyyy_xx[k] = -g_0_z_yyyyy_xx[k] * ab_x + g_0_z_yyyyy_xxx[k];

                g_0_z_xyyyyy_xy[k] = -g_0_z_yyyyy_xy[k] * ab_x + g_0_z_yyyyy_xxy[k];

                g_0_z_xyyyyy_xz[k] = -g_0_z_yyyyy_xz[k] * ab_x + g_0_z_yyyyy_xxz[k];

                g_0_z_xyyyyy_yy[k] = -g_0_z_yyyyy_yy[k] * ab_x + g_0_z_yyyyy_xyy[k];

                g_0_z_xyyyyy_yz[k] = -g_0_z_yyyyy_yz[k] * ab_x + g_0_z_yyyyy_xyz[k];

                g_0_z_xyyyyy_zz[k] = -g_0_z_yyyyy_zz[k] * ab_x + g_0_z_yyyyy_xzz[k];
            }

            /// Set up 432-438 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyyz_xx = cbuffer.data(id_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xy = cbuffer.data(id_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xz = cbuffer.data(id_geom_01_off + 434 * ccomps * dcomps);

            auto g_0_z_xyyyyz_yy = cbuffer.data(id_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_z_xyyyyz_yz = cbuffer.data(id_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_z_xyyyyz_zz = cbuffer.data(id_geom_01_off + 437 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyyz_xx, g_0_z_xyyyyz_xy, g_0_z_xyyyyz_xz, g_0_z_xyyyyz_yy, g_0_z_xyyyyz_yz, g_0_z_xyyyyz_zz, g_0_z_yyyyz_xx, g_0_z_yyyyz_xxx, g_0_z_yyyyz_xxy, g_0_z_yyyyz_xxz, g_0_z_yyyyz_xy, g_0_z_yyyyz_xyy, g_0_z_yyyyz_xyz, g_0_z_yyyyz_xz, g_0_z_yyyyz_xzz, g_0_z_yyyyz_yy, g_0_z_yyyyz_yz, g_0_z_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyyz_xx[k] = -g_0_z_yyyyz_xx[k] * ab_x + g_0_z_yyyyz_xxx[k];

                g_0_z_xyyyyz_xy[k] = -g_0_z_yyyyz_xy[k] * ab_x + g_0_z_yyyyz_xxy[k];

                g_0_z_xyyyyz_xz[k] = -g_0_z_yyyyz_xz[k] * ab_x + g_0_z_yyyyz_xxz[k];

                g_0_z_xyyyyz_yy[k] = -g_0_z_yyyyz_yy[k] * ab_x + g_0_z_yyyyz_xyy[k];

                g_0_z_xyyyyz_yz[k] = -g_0_z_yyyyz_yz[k] * ab_x + g_0_z_yyyyz_xyz[k];

                g_0_z_xyyyyz_zz[k] = -g_0_z_yyyyz_zz[k] * ab_x + g_0_z_yyyyz_xzz[k];
            }

            /// Set up 438-444 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyzz_xx = cbuffer.data(id_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xy = cbuffer.data(id_geom_01_off + 439 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xz = cbuffer.data(id_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_z_xyyyzz_yy = cbuffer.data(id_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_z_xyyyzz_yz = cbuffer.data(id_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_z_xyyyzz_zz = cbuffer.data(id_geom_01_off + 443 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyzz_xx, g_0_z_xyyyzz_xy, g_0_z_xyyyzz_xz, g_0_z_xyyyzz_yy, g_0_z_xyyyzz_yz, g_0_z_xyyyzz_zz, g_0_z_yyyzz_xx, g_0_z_yyyzz_xxx, g_0_z_yyyzz_xxy, g_0_z_yyyzz_xxz, g_0_z_yyyzz_xy, g_0_z_yyyzz_xyy, g_0_z_yyyzz_xyz, g_0_z_yyyzz_xz, g_0_z_yyyzz_xzz, g_0_z_yyyzz_yy, g_0_z_yyyzz_yz, g_0_z_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyzz_xx[k] = -g_0_z_yyyzz_xx[k] * ab_x + g_0_z_yyyzz_xxx[k];

                g_0_z_xyyyzz_xy[k] = -g_0_z_yyyzz_xy[k] * ab_x + g_0_z_yyyzz_xxy[k];

                g_0_z_xyyyzz_xz[k] = -g_0_z_yyyzz_xz[k] * ab_x + g_0_z_yyyzz_xxz[k];

                g_0_z_xyyyzz_yy[k] = -g_0_z_yyyzz_yy[k] * ab_x + g_0_z_yyyzz_xyy[k];

                g_0_z_xyyyzz_yz[k] = -g_0_z_yyyzz_yz[k] * ab_x + g_0_z_yyyzz_xyz[k];

                g_0_z_xyyyzz_zz[k] = -g_0_z_yyyzz_zz[k] * ab_x + g_0_z_yyyzz_xzz[k];
            }

            /// Set up 444-450 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyzzz_xx = cbuffer.data(id_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xy = cbuffer.data(id_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xz = cbuffer.data(id_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_z_xyyzzz_yy = cbuffer.data(id_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_z_xyyzzz_yz = cbuffer.data(id_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_z_xyyzzz_zz = cbuffer.data(id_geom_01_off + 449 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyzzz_xx, g_0_z_xyyzzz_xy, g_0_z_xyyzzz_xz, g_0_z_xyyzzz_yy, g_0_z_xyyzzz_yz, g_0_z_xyyzzz_zz, g_0_z_yyzzz_xx, g_0_z_yyzzz_xxx, g_0_z_yyzzz_xxy, g_0_z_yyzzz_xxz, g_0_z_yyzzz_xy, g_0_z_yyzzz_xyy, g_0_z_yyzzz_xyz, g_0_z_yyzzz_xz, g_0_z_yyzzz_xzz, g_0_z_yyzzz_yy, g_0_z_yyzzz_yz, g_0_z_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyzzz_xx[k] = -g_0_z_yyzzz_xx[k] * ab_x + g_0_z_yyzzz_xxx[k];

                g_0_z_xyyzzz_xy[k] = -g_0_z_yyzzz_xy[k] * ab_x + g_0_z_yyzzz_xxy[k];

                g_0_z_xyyzzz_xz[k] = -g_0_z_yyzzz_xz[k] * ab_x + g_0_z_yyzzz_xxz[k];

                g_0_z_xyyzzz_yy[k] = -g_0_z_yyzzz_yy[k] * ab_x + g_0_z_yyzzz_xyy[k];

                g_0_z_xyyzzz_yz[k] = -g_0_z_yyzzz_yz[k] * ab_x + g_0_z_yyzzz_xyz[k];

                g_0_z_xyyzzz_zz[k] = -g_0_z_yyzzz_zz[k] * ab_x + g_0_z_yyzzz_xzz[k];
            }

            /// Set up 450-456 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyzzzz_xx = cbuffer.data(id_geom_01_off + 450 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xy = cbuffer.data(id_geom_01_off + 451 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xz = cbuffer.data(id_geom_01_off + 452 * ccomps * dcomps);

            auto g_0_z_xyzzzz_yy = cbuffer.data(id_geom_01_off + 453 * ccomps * dcomps);

            auto g_0_z_xyzzzz_yz = cbuffer.data(id_geom_01_off + 454 * ccomps * dcomps);

            auto g_0_z_xyzzzz_zz = cbuffer.data(id_geom_01_off + 455 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyzzzz_xx, g_0_z_xyzzzz_xy, g_0_z_xyzzzz_xz, g_0_z_xyzzzz_yy, g_0_z_xyzzzz_yz, g_0_z_xyzzzz_zz, g_0_z_yzzzz_xx, g_0_z_yzzzz_xxx, g_0_z_yzzzz_xxy, g_0_z_yzzzz_xxz, g_0_z_yzzzz_xy, g_0_z_yzzzz_xyy, g_0_z_yzzzz_xyz, g_0_z_yzzzz_xz, g_0_z_yzzzz_xzz, g_0_z_yzzzz_yy, g_0_z_yzzzz_yz, g_0_z_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyzzzz_xx[k] = -g_0_z_yzzzz_xx[k] * ab_x + g_0_z_yzzzz_xxx[k];

                g_0_z_xyzzzz_xy[k] = -g_0_z_yzzzz_xy[k] * ab_x + g_0_z_yzzzz_xxy[k];

                g_0_z_xyzzzz_xz[k] = -g_0_z_yzzzz_xz[k] * ab_x + g_0_z_yzzzz_xxz[k];

                g_0_z_xyzzzz_yy[k] = -g_0_z_yzzzz_yy[k] * ab_x + g_0_z_yzzzz_xyy[k];

                g_0_z_xyzzzz_yz[k] = -g_0_z_yzzzz_yz[k] * ab_x + g_0_z_yzzzz_xyz[k];

                g_0_z_xyzzzz_zz[k] = -g_0_z_yzzzz_zz[k] * ab_x + g_0_z_yzzzz_xzz[k];
            }

            /// Set up 456-462 components of targeted buffer : cbuffer.data(

            auto g_0_z_xzzzzz_xx = cbuffer.data(id_geom_01_off + 456 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xy = cbuffer.data(id_geom_01_off + 457 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xz = cbuffer.data(id_geom_01_off + 458 * ccomps * dcomps);

            auto g_0_z_xzzzzz_yy = cbuffer.data(id_geom_01_off + 459 * ccomps * dcomps);

            auto g_0_z_xzzzzz_yz = cbuffer.data(id_geom_01_off + 460 * ccomps * dcomps);

            auto g_0_z_xzzzzz_zz = cbuffer.data(id_geom_01_off + 461 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xzzzzz_xx, g_0_z_xzzzzz_xy, g_0_z_xzzzzz_xz, g_0_z_xzzzzz_yy, g_0_z_xzzzzz_yz, g_0_z_xzzzzz_zz, g_0_z_zzzzz_xx, g_0_z_zzzzz_xxx, g_0_z_zzzzz_xxy, g_0_z_zzzzz_xxz, g_0_z_zzzzz_xy, g_0_z_zzzzz_xyy, g_0_z_zzzzz_xyz, g_0_z_zzzzz_xz, g_0_z_zzzzz_xzz, g_0_z_zzzzz_yy, g_0_z_zzzzz_yz, g_0_z_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xzzzzz_xx[k] = -g_0_z_zzzzz_xx[k] * ab_x + g_0_z_zzzzz_xxx[k];

                g_0_z_xzzzzz_xy[k] = -g_0_z_zzzzz_xy[k] * ab_x + g_0_z_zzzzz_xxy[k];

                g_0_z_xzzzzz_xz[k] = -g_0_z_zzzzz_xz[k] * ab_x + g_0_z_zzzzz_xxz[k];

                g_0_z_xzzzzz_yy[k] = -g_0_z_zzzzz_yy[k] * ab_x + g_0_z_zzzzz_xyy[k];

                g_0_z_xzzzzz_yz[k] = -g_0_z_zzzzz_yz[k] * ab_x + g_0_z_zzzzz_xyz[k];

                g_0_z_xzzzzz_zz[k] = -g_0_z_zzzzz_zz[k] * ab_x + g_0_z_zzzzz_xzz[k];
            }

            /// Set up 462-468 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyyy_xx = cbuffer.data(id_geom_01_off + 462 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xy = cbuffer.data(id_geom_01_off + 463 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xz = cbuffer.data(id_geom_01_off + 464 * ccomps * dcomps);

            auto g_0_z_yyyyyy_yy = cbuffer.data(id_geom_01_off + 465 * ccomps * dcomps);

            auto g_0_z_yyyyyy_yz = cbuffer.data(id_geom_01_off + 466 * ccomps * dcomps);

            auto g_0_z_yyyyyy_zz = cbuffer.data(id_geom_01_off + 467 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyy_xx, g_0_z_yyyyy_xxy, g_0_z_yyyyy_xy, g_0_z_yyyyy_xyy, g_0_z_yyyyy_xyz, g_0_z_yyyyy_xz, g_0_z_yyyyy_yy, g_0_z_yyyyy_yyy, g_0_z_yyyyy_yyz, g_0_z_yyyyy_yz, g_0_z_yyyyy_yzz, g_0_z_yyyyy_zz, g_0_z_yyyyyy_xx, g_0_z_yyyyyy_xy, g_0_z_yyyyyy_xz, g_0_z_yyyyyy_yy, g_0_z_yyyyyy_yz, g_0_z_yyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyyy_xx[k] = -g_0_z_yyyyy_xx[k] * ab_y + g_0_z_yyyyy_xxy[k];

                g_0_z_yyyyyy_xy[k] = -g_0_z_yyyyy_xy[k] * ab_y + g_0_z_yyyyy_xyy[k];

                g_0_z_yyyyyy_xz[k] = -g_0_z_yyyyy_xz[k] * ab_y + g_0_z_yyyyy_xyz[k];

                g_0_z_yyyyyy_yy[k] = -g_0_z_yyyyy_yy[k] * ab_y + g_0_z_yyyyy_yyy[k];

                g_0_z_yyyyyy_yz[k] = -g_0_z_yyyyy_yz[k] * ab_y + g_0_z_yyyyy_yyz[k];

                g_0_z_yyyyyy_zz[k] = -g_0_z_yyyyy_zz[k] * ab_y + g_0_z_yyyyy_yzz[k];
            }

            /// Set up 468-474 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyyz_xx = cbuffer.data(id_geom_01_off + 468 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xy = cbuffer.data(id_geom_01_off + 469 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xz = cbuffer.data(id_geom_01_off + 470 * ccomps * dcomps);

            auto g_0_z_yyyyyz_yy = cbuffer.data(id_geom_01_off + 471 * ccomps * dcomps);

            auto g_0_z_yyyyyz_yz = cbuffer.data(id_geom_01_off + 472 * ccomps * dcomps);

            auto g_0_z_yyyyyz_zz = cbuffer.data(id_geom_01_off + 473 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyyz_xx, g_0_z_yyyyyz_xy, g_0_z_yyyyyz_xz, g_0_z_yyyyyz_yy, g_0_z_yyyyyz_yz, g_0_z_yyyyyz_zz, g_0_z_yyyyz_xx, g_0_z_yyyyz_xxy, g_0_z_yyyyz_xy, g_0_z_yyyyz_xyy, g_0_z_yyyyz_xyz, g_0_z_yyyyz_xz, g_0_z_yyyyz_yy, g_0_z_yyyyz_yyy, g_0_z_yyyyz_yyz, g_0_z_yyyyz_yz, g_0_z_yyyyz_yzz, g_0_z_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyyz_xx[k] = -g_0_z_yyyyz_xx[k] * ab_y + g_0_z_yyyyz_xxy[k];

                g_0_z_yyyyyz_xy[k] = -g_0_z_yyyyz_xy[k] * ab_y + g_0_z_yyyyz_xyy[k];

                g_0_z_yyyyyz_xz[k] = -g_0_z_yyyyz_xz[k] * ab_y + g_0_z_yyyyz_xyz[k];

                g_0_z_yyyyyz_yy[k] = -g_0_z_yyyyz_yy[k] * ab_y + g_0_z_yyyyz_yyy[k];

                g_0_z_yyyyyz_yz[k] = -g_0_z_yyyyz_yz[k] * ab_y + g_0_z_yyyyz_yyz[k];

                g_0_z_yyyyyz_zz[k] = -g_0_z_yyyyz_zz[k] * ab_y + g_0_z_yyyyz_yzz[k];
            }

            /// Set up 474-480 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyzz_xx = cbuffer.data(id_geom_01_off + 474 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xy = cbuffer.data(id_geom_01_off + 475 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xz = cbuffer.data(id_geom_01_off + 476 * ccomps * dcomps);

            auto g_0_z_yyyyzz_yy = cbuffer.data(id_geom_01_off + 477 * ccomps * dcomps);

            auto g_0_z_yyyyzz_yz = cbuffer.data(id_geom_01_off + 478 * ccomps * dcomps);

            auto g_0_z_yyyyzz_zz = cbuffer.data(id_geom_01_off + 479 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyzz_xx, g_0_z_yyyyzz_xy, g_0_z_yyyyzz_xz, g_0_z_yyyyzz_yy, g_0_z_yyyyzz_yz, g_0_z_yyyyzz_zz, g_0_z_yyyzz_xx, g_0_z_yyyzz_xxy, g_0_z_yyyzz_xy, g_0_z_yyyzz_xyy, g_0_z_yyyzz_xyz, g_0_z_yyyzz_xz, g_0_z_yyyzz_yy, g_0_z_yyyzz_yyy, g_0_z_yyyzz_yyz, g_0_z_yyyzz_yz, g_0_z_yyyzz_yzz, g_0_z_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyzz_xx[k] = -g_0_z_yyyzz_xx[k] * ab_y + g_0_z_yyyzz_xxy[k];

                g_0_z_yyyyzz_xy[k] = -g_0_z_yyyzz_xy[k] * ab_y + g_0_z_yyyzz_xyy[k];

                g_0_z_yyyyzz_xz[k] = -g_0_z_yyyzz_xz[k] * ab_y + g_0_z_yyyzz_xyz[k];

                g_0_z_yyyyzz_yy[k] = -g_0_z_yyyzz_yy[k] * ab_y + g_0_z_yyyzz_yyy[k];

                g_0_z_yyyyzz_yz[k] = -g_0_z_yyyzz_yz[k] * ab_y + g_0_z_yyyzz_yyz[k];

                g_0_z_yyyyzz_zz[k] = -g_0_z_yyyzz_zz[k] * ab_y + g_0_z_yyyzz_yzz[k];
            }

            /// Set up 480-486 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyzzz_xx = cbuffer.data(id_geom_01_off + 480 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xy = cbuffer.data(id_geom_01_off + 481 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xz = cbuffer.data(id_geom_01_off + 482 * ccomps * dcomps);

            auto g_0_z_yyyzzz_yy = cbuffer.data(id_geom_01_off + 483 * ccomps * dcomps);

            auto g_0_z_yyyzzz_yz = cbuffer.data(id_geom_01_off + 484 * ccomps * dcomps);

            auto g_0_z_yyyzzz_zz = cbuffer.data(id_geom_01_off + 485 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyzzz_xx, g_0_z_yyyzzz_xy, g_0_z_yyyzzz_xz, g_0_z_yyyzzz_yy, g_0_z_yyyzzz_yz, g_0_z_yyyzzz_zz, g_0_z_yyzzz_xx, g_0_z_yyzzz_xxy, g_0_z_yyzzz_xy, g_0_z_yyzzz_xyy, g_0_z_yyzzz_xyz, g_0_z_yyzzz_xz, g_0_z_yyzzz_yy, g_0_z_yyzzz_yyy, g_0_z_yyzzz_yyz, g_0_z_yyzzz_yz, g_0_z_yyzzz_yzz, g_0_z_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyzzz_xx[k] = -g_0_z_yyzzz_xx[k] * ab_y + g_0_z_yyzzz_xxy[k];

                g_0_z_yyyzzz_xy[k] = -g_0_z_yyzzz_xy[k] * ab_y + g_0_z_yyzzz_xyy[k];

                g_0_z_yyyzzz_xz[k] = -g_0_z_yyzzz_xz[k] * ab_y + g_0_z_yyzzz_xyz[k];

                g_0_z_yyyzzz_yy[k] = -g_0_z_yyzzz_yy[k] * ab_y + g_0_z_yyzzz_yyy[k];

                g_0_z_yyyzzz_yz[k] = -g_0_z_yyzzz_yz[k] * ab_y + g_0_z_yyzzz_yyz[k];

                g_0_z_yyyzzz_zz[k] = -g_0_z_yyzzz_zz[k] * ab_y + g_0_z_yyzzz_yzz[k];
            }

            /// Set up 486-492 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyzzzz_xx = cbuffer.data(id_geom_01_off + 486 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xy = cbuffer.data(id_geom_01_off + 487 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xz = cbuffer.data(id_geom_01_off + 488 * ccomps * dcomps);

            auto g_0_z_yyzzzz_yy = cbuffer.data(id_geom_01_off + 489 * ccomps * dcomps);

            auto g_0_z_yyzzzz_yz = cbuffer.data(id_geom_01_off + 490 * ccomps * dcomps);

            auto g_0_z_yyzzzz_zz = cbuffer.data(id_geom_01_off + 491 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyzzzz_xx, g_0_z_yyzzzz_xy, g_0_z_yyzzzz_xz, g_0_z_yyzzzz_yy, g_0_z_yyzzzz_yz, g_0_z_yyzzzz_zz, g_0_z_yzzzz_xx, g_0_z_yzzzz_xxy, g_0_z_yzzzz_xy, g_0_z_yzzzz_xyy, g_0_z_yzzzz_xyz, g_0_z_yzzzz_xz, g_0_z_yzzzz_yy, g_0_z_yzzzz_yyy, g_0_z_yzzzz_yyz, g_0_z_yzzzz_yz, g_0_z_yzzzz_yzz, g_0_z_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyzzzz_xx[k] = -g_0_z_yzzzz_xx[k] * ab_y + g_0_z_yzzzz_xxy[k];

                g_0_z_yyzzzz_xy[k] = -g_0_z_yzzzz_xy[k] * ab_y + g_0_z_yzzzz_xyy[k];

                g_0_z_yyzzzz_xz[k] = -g_0_z_yzzzz_xz[k] * ab_y + g_0_z_yzzzz_xyz[k];

                g_0_z_yyzzzz_yy[k] = -g_0_z_yzzzz_yy[k] * ab_y + g_0_z_yzzzz_yyy[k];

                g_0_z_yyzzzz_yz[k] = -g_0_z_yzzzz_yz[k] * ab_y + g_0_z_yzzzz_yyz[k];

                g_0_z_yyzzzz_zz[k] = -g_0_z_yzzzz_zz[k] * ab_y + g_0_z_yzzzz_yzz[k];
            }

            /// Set up 492-498 components of targeted buffer : cbuffer.data(

            auto g_0_z_yzzzzz_xx = cbuffer.data(id_geom_01_off + 492 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xy = cbuffer.data(id_geom_01_off + 493 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xz = cbuffer.data(id_geom_01_off + 494 * ccomps * dcomps);

            auto g_0_z_yzzzzz_yy = cbuffer.data(id_geom_01_off + 495 * ccomps * dcomps);

            auto g_0_z_yzzzzz_yz = cbuffer.data(id_geom_01_off + 496 * ccomps * dcomps);

            auto g_0_z_yzzzzz_zz = cbuffer.data(id_geom_01_off + 497 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yzzzzz_xx, g_0_z_yzzzzz_xy, g_0_z_yzzzzz_xz, g_0_z_yzzzzz_yy, g_0_z_yzzzzz_yz, g_0_z_yzzzzz_zz, g_0_z_zzzzz_xx, g_0_z_zzzzz_xxy, g_0_z_zzzzz_xy, g_0_z_zzzzz_xyy, g_0_z_zzzzz_xyz, g_0_z_zzzzz_xz, g_0_z_zzzzz_yy, g_0_z_zzzzz_yyy, g_0_z_zzzzz_yyz, g_0_z_zzzzz_yz, g_0_z_zzzzz_yzz, g_0_z_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yzzzzz_xx[k] = -g_0_z_zzzzz_xx[k] * ab_y + g_0_z_zzzzz_xxy[k];

                g_0_z_yzzzzz_xy[k] = -g_0_z_zzzzz_xy[k] * ab_y + g_0_z_zzzzz_xyy[k];

                g_0_z_yzzzzz_xz[k] = -g_0_z_zzzzz_xz[k] * ab_y + g_0_z_zzzzz_xyz[k];

                g_0_z_yzzzzz_yy[k] = -g_0_z_zzzzz_yy[k] * ab_y + g_0_z_zzzzz_yyy[k];

                g_0_z_yzzzzz_yz[k] = -g_0_z_zzzzz_yz[k] * ab_y + g_0_z_zzzzz_yyz[k];

                g_0_z_yzzzzz_zz[k] = -g_0_z_zzzzz_zz[k] * ab_y + g_0_z_zzzzz_yzz[k];
            }

            /// Set up 498-504 components of targeted buffer : cbuffer.data(

            auto g_0_z_zzzzzz_xx = cbuffer.data(id_geom_01_off + 498 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xy = cbuffer.data(id_geom_01_off + 499 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xz = cbuffer.data(id_geom_01_off + 500 * ccomps * dcomps);

            auto g_0_z_zzzzzz_yy = cbuffer.data(id_geom_01_off + 501 * ccomps * dcomps);

            auto g_0_z_zzzzzz_yz = cbuffer.data(id_geom_01_off + 502 * ccomps * dcomps);

            auto g_0_z_zzzzzz_zz = cbuffer.data(id_geom_01_off + 503 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zzzzz_xx, g_0_z_zzzzz_xxz, g_0_z_zzzzz_xy, g_0_z_zzzzz_xyz, g_0_z_zzzzz_xz, g_0_z_zzzzz_xzz, g_0_z_zzzzz_yy, g_0_z_zzzzz_yyz, g_0_z_zzzzz_yz, g_0_z_zzzzz_yzz, g_0_z_zzzzz_zz, g_0_z_zzzzz_zzz, g_0_z_zzzzzz_xx, g_0_z_zzzzzz_xy, g_0_z_zzzzzz_xz, g_0_z_zzzzzz_yy, g_0_z_zzzzzz_yz, g_0_z_zzzzzz_zz, g_zzzzz_xx, g_zzzzz_xy, g_zzzzz_xz, g_zzzzz_yy, g_zzzzz_yz, g_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_zzzzzz_xx[k] = g_zzzzz_xx[k] - g_0_z_zzzzz_xx[k] * ab_z + g_0_z_zzzzz_xxz[k];

                g_0_z_zzzzzz_xy[k] = g_zzzzz_xy[k] - g_0_z_zzzzz_xy[k] * ab_z + g_0_z_zzzzz_xyz[k];

                g_0_z_zzzzzz_xz[k] = g_zzzzz_xz[k] - g_0_z_zzzzz_xz[k] * ab_z + g_0_z_zzzzz_xzz[k];

                g_0_z_zzzzzz_yy[k] = g_zzzzz_yy[k] - g_0_z_zzzzz_yy[k] * ab_z + g_0_z_zzzzz_yyz[k];

                g_0_z_zzzzzz_yz[k] = g_zzzzz_yz[k] - g_0_z_zzzzz_yz[k] * ab_z + g_0_z_zzzzz_yzz[k];

                g_0_z_zzzzzz_zz[k] = g_zzzzz_zz[k] - g_0_z_zzzzz_zz[k] * ab_z + g_0_z_zzzzz_zzz[k];
            }
        }
    }
}

} // erirec namespace

