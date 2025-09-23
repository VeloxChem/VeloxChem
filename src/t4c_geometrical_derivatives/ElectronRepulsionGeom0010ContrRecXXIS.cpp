#include "ElectronRepulsionGeom0010ContrRecXXIS.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxis(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxis,
                                            const size_t idx_xxhs,
                                            const size_t idx_geom_10_xxhs,
                                            const size_t idx_geom_10_xxhp,
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
            /// Set up components of auxilary buffer : SSHS

            const auto hs_off = idx_xxhs + (i * bcomps + j) * 21;

            auto g_xxxxx_0 = cbuffer.data(hs_off + 0);

            auto g_xxxxy_0 = cbuffer.data(hs_off + 1);

            auto g_xxxxz_0 = cbuffer.data(hs_off + 2);

            auto g_xxxyy_0 = cbuffer.data(hs_off + 3);

            auto g_xxxyz_0 = cbuffer.data(hs_off + 4);

            auto g_xxxzz_0 = cbuffer.data(hs_off + 5);

            auto g_xxyyy_0 = cbuffer.data(hs_off + 6);

            auto g_xxyyz_0 = cbuffer.data(hs_off + 7);

            auto g_xxyzz_0 = cbuffer.data(hs_off + 8);

            auto g_xxzzz_0 = cbuffer.data(hs_off + 9);

            auto g_xyyyy_0 = cbuffer.data(hs_off + 10);

            auto g_xyyyz_0 = cbuffer.data(hs_off + 11);

            auto g_xyyzz_0 = cbuffer.data(hs_off + 12);

            auto g_xyzzz_0 = cbuffer.data(hs_off + 13);

            auto g_xzzzz_0 = cbuffer.data(hs_off + 14);

            auto g_yyyyy_0 = cbuffer.data(hs_off + 15);

            auto g_yyyyz_0 = cbuffer.data(hs_off + 16);

            auto g_yyyzz_0 = cbuffer.data(hs_off + 17);

            auto g_yyzzz_0 = cbuffer.data(hs_off + 18);

            auto g_yzzzz_0 = cbuffer.data(hs_off + 19);

            auto g_zzzzz_0 = cbuffer.data(hs_off + 20);

            /// Set up components of auxilary buffer : SSHS

            const auto hs_geom_10_off = idx_geom_10_xxhs + (i * bcomps + j) * 21;

            auto g_x_0_xxxxx_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxxy_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxxz_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxyy_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxyz_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxzz_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxyyy_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxyyz_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxyzz_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxzzz_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xyyyy_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xyyyz_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xyyzz_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xyzzz_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xzzzz_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_yyyyy_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_yyyyz_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_yyyzz_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_yyzzz_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_yzzzz_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_zzzzz_0 = cbuffer.data(hs_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_y_0_xxxxx_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 0);

            auto g_y_0_xxxxy_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 1);

            auto g_y_0_xxxxz_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 2);

            auto g_y_0_xxxyy_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 3);

            auto g_y_0_xxxyz_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 4);

            auto g_y_0_xxxzz_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 5);

            auto g_y_0_xxyyy_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 6);

            auto g_y_0_xxyyz_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 7);

            auto g_y_0_xxyzz_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 8);

            auto g_y_0_xxzzz_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 9);

            auto g_y_0_xyyyy_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 10);

            auto g_y_0_xyyyz_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 11);

            auto g_y_0_xyyzz_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 12);

            auto g_y_0_xyzzz_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 13);

            auto g_y_0_xzzzz_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 14);

            auto g_y_0_yyyyy_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 15);

            auto g_y_0_yyyyz_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 16);

            auto g_y_0_yyyzz_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 17);

            auto g_y_0_yyzzz_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 18);

            auto g_y_0_yzzzz_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 19);

            auto g_y_0_zzzzz_0 = cbuffer.data(hs_geom_10_off + 21 * acomps * bcomps + 20);

            auto g_z_0_xxxxx_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 0);

            auto g_z_0_xxxxy_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 1);

            auto g_z_0_xxxxz_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 2);

            auto g_z_0_xxxyy_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 3);

            auto g_z_0_xxxyz_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 4);

            auto g_z_0_xxxzz_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 5);

            auto g_z_0_xxyyy_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 6);

            auto g_z_0_xxyyz_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 7);

            auto g_z_0_xxyzz_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 8);

            auto g_z_0_xxzzz_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 9);

            auto g_z_0_xyyyy_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 10);

            auto g_z_0_xyyyz_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 11);

            auto g_z_0_xyyzz_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 12);

            auto g_z_0_xyzzz_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 13);

            auto g_z_0_xzzzz_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 14);

            auto g_z_0_yyyyy_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 15);

            auto g_z_0_yyyyz_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 16);

            auto g_z_0_yyyzz_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 17);

            auto g_z_0_yyzzz_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 18);

            auto g_z_0_yzzzz_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 19);

            auto g_z_0_zzzzz_0 = cbuffer.data(hs_geom_10_off + 42 * acomps * bcomps + 20);

            /// Set up components of auxilary buffer : SSHP

            const auto hp_geom_10_off = idx_geom_10_xxhp + (i * bcomps + j) * 63;

            auto g_x_0_xxxxx_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxxx_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxxx_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxxy_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxxy_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxxy_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxxxz_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxxz_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxxz_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxxyy_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxxyy_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxxyy_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxxyz_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxxyz_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxxyz_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxxzz_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxxzz_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxzz_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xxyyy_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxyyy_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxyyy_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xxyyz_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxyyz_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxyyz_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xxyzz_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxyzz_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxyzz_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxzzz_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxzzz_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxzzz_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_xyyyy_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xyyyy_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xyyyy_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xyyyz_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xyyyz_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xyyyz_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xyyzz_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xyyzz_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xyyzz_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xyzzz_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xyzzz_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xyzzz_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_xzzzz_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xzzzz_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xzzzz_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_yyyyy_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_yyyyy_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_yyyyy_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_yyyyz_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_yyyyz_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_yyyyz_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_yyyzz_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_yyyzz_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_yyyzz_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_yyzzz_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_yyzzz_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_yyzzz_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_yzzzz_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_yzzzz_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_yzzzz_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 59);

            auto g_x_0_zzzzz_x = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_zzzzz_y = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_zzzzz_z = cbuffer.data(hp_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_y_0_xxxxx_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 0);

            auto g_y_0_xxxxx_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 1);

            auto g_y_0_xxxxx_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 2);

            auto g_y_0_xxxxy_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 3);

            auto g_y_0_xxxxy_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 4);

            auto g_y_0_xxxxy_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 5);

            auto g_y_0_xxxxz_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 6);

            auto g_y_0_xxxxz_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 7);

            auto g_y_0_xxxxz_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 8);

            auto g_y_0_xxxyy_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 9);

            auto g_y_0_xxxyy_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 10);

            auto g_y_0_xxxyy_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 11);

            auto g_y_0_xxxyz_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 12);

            auto g_y_0_xxxyz_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 13);

            auto g_y_0_xxxyz_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 14);

            auto g_y_0_xxxzz_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 15);

            auto g_y_0_xxxzz_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 16);

            auto g_y_0_xxxzz_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 17);

            auto g_y_0_xxyyy_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 18);

            auto g_y_0_xxyyy_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 19);

            auto g_y_0_xxyyy_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 20);

            auto g_y_0_xxyyz_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 21);

            auto g_y_0_xxyyz_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 22);

            auto g_y_0_xxyyz_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 23);

            auto g_y_0_xxyzz_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 24);

            auto g_y_0_xxyzz_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 25);

            auto g_y_0_xxyzz_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 26);

            auto g_y_0_xxzzz_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 27);

            auto g_y_0_xxzzz_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 28);

            auto g_y_0_xxzzz_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 29);

            auto g_y_0_xyyyy_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 30);

            auto g_y_0_xyyyy_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 31);

            auto g_y_0_xyyyy_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 32);

            auto g_y_0_xyyyz_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 33);

            auto g_y_0_xyyyz_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 34);

            auto g_y_0_xyyyz_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 35);

            auto g_y_0_xyyzz_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 36);

            auto g_y_0_xyyzz_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 37);

            auto g_y_0_xyyzz_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 38);

            auto g_y_0_xyzzz_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 39);

            auto g_y_0_xyzzz_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 40);

            auto g_y_0_xyzzz_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 41);

            auto g_y_0_xzzzz_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 42);

            auto g_y_0_xzzzz_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 43);

            auto g_y_0_xzzzz_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 44);

            auto g_y_0_yyyyy_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 45);

            auto g_y_0_yyyyy_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 46);

            auto g_y_0_yyyyy_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 47);

            auto g_y_0_yyyyz_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 48);

            auto g_y_0_yyyyz_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 49);

            auto g_y_0_yyyyz_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 50);

            auto g_y_0_yyyzz_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 51);

            auto g_y_0_yyyzz_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 52);

            auto g_y_0_yyyzz_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 53);

            auto g_y_0_yyzzz_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 54);

            auto g_y_0_yyzzz_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 55);

            auto g_y_0_yyzzz_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 56);

            auto g_y_0_yzzzz_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 57);

            auto g_y_0_yzzzz_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 58);

            auto g_y_0_yzzzz_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 59);

            auto g_y_0_zzzzz_x = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 60);

            auto g_y_0_zzzzz_y = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 61);

            auto g_y_0_zzzzz_z = cbuffer.data(hp_geom_10_off + 63 * acomps * bcomps + 62);

            auto g_z_0_xxxxx_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 0);

            auto g_z_0_xxxxx_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 1);

            auto g_z_0_xxxxx_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 2);

            auto g_z_0_xxxxy_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 3);

            auto g_z_0_xxxxy_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 4);

            auto g_z_0_xxxxy_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 5);

            auto g_z_0_xxxxz_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 6);

            auto g_z_0_xxxxz_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 7);

            auto g_z_0_xxxxz_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 8);

            auto g_z_0_xxxyy_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 9);

            auto g_z_0_xxxyy_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 10);

            auto g_z_0_xxxyy_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 11);

            auto g_z_0_xxxyz_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 12);

            auto g_z_0_xxxyz_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 13);

            auto g_z_0_xxxyz_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 14);

            auto g_z_0_xxxzz_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 15);

            auto g_z_0_xxxzz_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 16);

            auto g_z_0_xxxzz_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 17);

            auto g_z_0_xxyyy_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 18);

            auto g_z_0_xxyyy_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 19);

            auto g_z_0_xxyyy_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 20);

            auto g_z_0_xxyyz_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 21);

            auto g_z_0_xxyyz_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 22);

            auto g_z_0_xxyyz_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 23);

            auto g_z_0_xxyzz_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 24);

            auto g_z_0_xxyzz_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 25);

            auto g_z_0_xxyzz_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 26);

            auto g_z_0_xxzzz_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 27);

            auto g_z_0_xxzzz_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 28);

            auto g_z_0_xxzzz_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 29);

            auto g_z_0_xyyyy_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 30);

            auto g_z_0_xyyyy_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 31);

            auto g_z_0_xyyyy_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 32);

            auto g_z_0_xyyyz_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 33);

            auto g_z_0_xyyyz_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 34);

            auto g_z_0_xyyyz_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 35);

            auto g_z_0_xyyzz_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 36);

            auto g_z_0_xyyzz_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 37);

            auto g_z_0_xyyzz_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 38);

            auto g_z_0_xyzzz_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 39);

            auto g_z_0_xyzzz_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 40);

            auto g_z_0_xyzzz_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 41);

            auto g_z_0_xzzzz_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 42);

            auto g_z_0_xzzzz_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 43);

            auto g_z_0_xzzzz_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 44);

            auto g_z_0_yyyyy_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 45);

            auto g_z_0_yyyyy_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 46);

            auto g_z_0_yyyyy_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 47);

            auto g_z_0_yyyyz_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 48);

            auto g_z_0_yyyyz_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 49);

            auto g_z_0_yyyyz_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 50);

            auto g_z_0_yyyzz_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 51);

            auto g_z_0_yyyzz_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 52);

            auto g_z_0_yyyzz_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 53);

            auto g_z_0_yyzzz_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 54);

            auto g_z_0_yyzzz_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 55);

            auto g_z_0_yyzzz_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 56);

            auto g_z_0_yzzzz_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 57);

            auto g_z_0_yzzzz_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 58);

            auto g_z_0_yzzzz_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 59);

            auto g_z_0_zzzzz_x = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 60);

            auto g_z_0_zzzzz_y = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 61);

            auto g_z_0_zzzzz_z = cbuffer.data(hp_geom_10_off + 126 * acomps * bcomps + 62);

            /// set up bra offset for contr_buffer_xxis

            const auto is_geom_10_off = idx_geom_10_xxis + (i * bcomps + j) * 28;

            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxx_0 = cbuffer.data(is_geom_10_off + 0 * acomps * bcomps + 0);

            #pragma omp simd aligned(cd_x, g_x_0_xxxxx_0, g_x_0_xxxxx_x, g_x_0_xxxxxx_0, g_xxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxx_0[k] = -g_xxxxx_0[k] - g_x_0_xxxxx_0[k] * cd_x[k] + g_x_0_xxxxx_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxy_0 = cbuffer.data(is_geom_10_off + 0 * acomps * bcomps + 1);

            #pragma omp simd aligned(cd_y, g_x_0_xxxxx_0, g_x_0_xxxxx_y, g_x_0_xxxxxy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxy_0[k] = -g_x_0_xxxxx_0[k] * cd_y[k] + g_x_0_xxxxx_y[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxz_0 = cbuffer.data(is_geom_10_off + 0 * acomps * bcomps + 2);

            #pragma omp simd aligned(cd_z, g_x_0_xxxxx_0, g_x_0_xxxxx_z, g_x_0_xxxxxz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxz_0[k] = -g_x_0_xxxxx_0[k] * cd_z[k] + g_x_0_xxxxx_z[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxyy_0 = cbuffer.data(is_geom_10_off + 0 * acomps * bcomps + 3);

            #pragma omp simd aligned(cd_y, g_x_0_xxxxy_0, g_x_0_xxxxy_y, g_x_0_xxxxyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxyy_0[k] = -g_x_0_xxxxy_0[k] * cd_y[k] + g_x_0_xxxxy_y[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxyz_0 = cbuffer.data(is_geom_10_off + 0 * acomps * bcomps + 4);

            #pragma omp simd aligned(cd_y, g_x_0_xxxxyz_0, g_x_0_xxxxz_0, g_x_0_xxxxz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxyz_0[k] = -g_x_0_xxxxz_0[k] * cd_y[k] + g_x_0_xxxxz_y[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxzz_0 = cbuffer.data(is_geom_10_off + 0 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_z, g_x_0_xxxxz_0, g_x_0_xxxxz_z, g_x_0_xxxxzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxzz_0[k] = -g_x_0_xxxxz_0[k] * cd_z[k] + g_x_0_xxxxz_z[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyyy_0 = cbuffer.data(is_geom_10_off + 0 * acomps * bcomps + 6);

            #pragma omp simd aligned(cd_y, g_x_0_xxxyy_0, g_x_0_xxxyy_y, g_x_0_xxxyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyyy_0[k] = -g_x_0_xxxyy_0[k] * cd_y[k] + g_x_0_xxxyy_y[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyyz_0 = cbuffer.data(is_geom_10_off + 0 * acomps * bcomps + 7);

            #pragma omp simd aligned(cd_y, g_x_0_xxxyyz_0, g_x_0_xxxyz_0, g_x_0_xxxyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyyz_0[k] = -g_x_0_xxxyz_0[k] * cd_y[k] + g_x_0_xxxyz_y[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyzz_0 = cbuffer.data(is_geom_10_off + 0 * acomps * bcomps + 8);

            #pragma omp simd aligned(cd_y, g_x_0_xxxyzz_0, g_x_0_xxxzz_0, g_x_0_xxxzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyzz_0[k] = -g_x_0_xxxzz_0[k] * cd_y[k] + g_x_0_xxxzz_y[k];
            }

            /// Set up 9-10 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxzzz_0 = cbuffer.data(is_geom_10_off + 0 * acomps * bcomps + 9);

            #pragma omp simd aligned(cd_z, g_x_0_xxxzz_0, g_x_0_xxxzz_z, g_x_0_xxxzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxzzz_0[k] = -g_x_0_xxxzz_0[k] * cd_z[k] + g_x_0_xxxzz_z[k];
            }

            /// Set up 10-11 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyyy_0 = cbuffer.data(is_geom_10_off + 0 * acomps * bcomps + 10);

            #pragma omp simd aligned(cd_y, g_x_0_xxyyy_0, g_x_0_xxyyy_y, g_x_0_xxyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyyy_0[k] = -g_x_0_xxyyy_0[k] * cd_y[k] + g_x_0_xxyyy_y[k];
            }

            /// Set up 11-12 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyyz_0 = cbuffer.data(is_geom_10_off + 0 * acomps * bcomps + 11);

            #pragma omp simd aligned(cd_y, g_x_0_xxyyyz_0, g_x_0_xxyyz_0, g_x_0_xxyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyyz_0[k] = -g_x_0_xxyyz_0[k] * cd_y[k] + g_x_0_xxyyz_y[k];
            }

            /// Set up 12-13 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyzz_0 = cbuffer.data(is_geom_10_off + 0 * acomps * bcomps + 12);

            #pragma omp simd aligned(cd_y, g_x_0_xxyyzz_0, g_x_0_xxyzz_0, g_x_0_xxyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyzz_0[k] = -g_x_0_xxyzz_0[k] * cd_y[k] + g_x_0_xxyzz_y[k];
            }

            /// Set up 13-14 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyzzz_0 = cbuffer.data(is_geom_10_off + 0 * acomps * bcomps + 13);

            #pragma omp simd aligned(cd_y, g_x_0_xxyzzz_0, g_x_0_xxzzz_0, g_x_0_xxzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyzzz_0[k] = -g_x_0_xxzzz_0[k] * cd_y[k] + g_x_0_xxzzz_y[k];
            }

            /// Set up 14-15 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxzzzz_0 = cbuffer.data(is_geom_10_off + 0 * acomps * bcomps + 14);

            #pragma omp simd aligned(cd_z, g_x_0_xxzzz_0, g_x_0_xxzzz_z, g_x_0_xxzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzzzz_0[k] = -g_x_0_xxzzz_0[k] * cd_z[k] + g_x_0_xxzzz_z[k];
            }

            /// Set up 15-16 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyyy_0 = cbuffer.data(is_geom_10_off + 0 * acomps * bcomps + 15);

            #pragma omp simd aligned(cd_y, g_x_0_xyyyy_0, g_x_0_xyyyy_y, g_x_0_xyyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyyy_0[k] = -g_x_0_xyyyy_0[k] * cd_y[k] + g_x_0_xyyyy_y[k];
            }

            /// Set up 16-17 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyyz_0 = cbuffer.data(is_geom_10_off + 0 * acomps * bcomps + 16);

            #pragma omp simd aligned(cd_y, g_x_0_xyyyyz_0, g_x_0_xyyyz_0, g_x_0_xyyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyyz_0[k] = -g_x_0_xyyyz_0[k] * cd_y[k] + g_x_0_xyyyz_y[k];
            }

            /// Set up 17-18 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyzz_0 = cbuffer.data(is_geom_10_off + 0 * acomps * bcomps + 17);

            #pragma omp simd aligned(cd_y, g_x_0_xyyyzz_0, g_x_0_xyyzz_0, g_x_0_xyyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyzz_0[k] = -g_x_0_xyyzz_0[k] * cd_y[k] + g_x_0_xyyzz_y[k];
            }

            /// Set up 18-19 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyzzz_0 = cbuffer.data(is_geom_10_off + 0 * acomps * bcomps + 18);

            #pragma omp simd aligned(cd_y, g_x_0_xyyzzz_0, g_x_0_xyzzz_0, g_x_0_xyzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyzzz_0[k] = -g_x_0_xyzzz_0[k] * cd_y[k] + g_x_0_xyzzz_y[k];
            }

            /// Set up 19-20 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyzzzz_0 = cbuffer.data(is_geom_10_off + 0 * acomps * bcomps + 19);

            #pragma omp simd aligned(cd_y, g_x_0_xyzzzz_0, g_x_0_xzzzz_0, g_x_0_xzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzzzz_0[k] = -g_x_0_xzzzz_0[k] * cd_y[k] + g_x_0_xzzzz_y[k];
            }

            /// Set up 20-21 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzzzzz_0 = cbuffer.data(is_geom_10_off + 0 * acomps * bcomps + 20);

            #pragma omp simd aligned(cd_z, g_x_0_xzzzz_0, g_x_0_xzzzz_z, g_x_0_xzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzzzz_0[k] = -g_x_0_xzzzz_0[k] * cd_z[k] + g_x_0_xzzzz_z[k];
            }

            /// Set up 21-22 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyyy_0 = cbuffer.data(is_geom_10_off + 0 * acomps * bcomps + 21);

            #pragma omp simd aligned(cd_y, g_x_0_yyyyy_0, g_x_0_yyyyy_y, g_x_0_yyyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyyy_0[k] = -g_x_0_yyyyy_0[k] * cd_y[k] + g_x_0_yyyyy_y[k];
            }

            /// Set up 22-23 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyyz_0 = cbuffer.data(is_geom_10_off + 0 * acomps * bcomps + 22);

            #pragma omp simd aligned(cd_y, g_x_0_yyyyyz_0, g_x_0_yyyyz_0, g_x_0_yyyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyyz_0[k] = -g_x_0_yyyyz_0[k] * cd_y[k] + g_x_0_yyyyz_y[k];
            }

            /// Set up 23-24 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyzz_0 = cbuffer.data(is_geom_10_off + 0 * acomps * bcomps + 23);

            #pragma omp simd aligned(cd_y, g_x_0_yyyyzz_0, g_x_0_yyyzz_0, g_x_0_yyyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyzz_0[k] = -g_x_0_yyyzz_0[k] * cd_y[k] + g_x_0_yyyzz_y[k];
            }

            /// Set up 24-25 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyzzz_0 = cbuffer.data(is_geom_10_off + 0 * acomps * bcomps + 24);

            #pragma omp simd aligned(cd_y, g_x_0_yyyzzz_0, g_x_0_yyzzz_0, g_x_0_yyzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyzzz_0[k] = -g_x_0_yyzzz_0[k] * cd_y[k] + g_x_0_yyzzz_y[k];
            }

            /// Set up 25-26 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyzzzz_0 = cbuffer.data(is_geom_10_off + 0 * acomps * bcomps + 25);

            #pragma omp simd aligned(cd_y, g_x_0_yyzzzz_0, g_x_0_yzzzz_0, g_x_0_yzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzzzz_0[k] = -g_x_0_yzzzz_0[k] * cd_y[k] + g_x_0_yzzzz_y[k];
            }

            /// Set up 26-27 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzzzzz_0 = cbuffer.data(is_geom_10_off + 0 * acomps * bcomps + 26);

            #pragma omp simd aligned(cd_y, g_x_0_yzzzzz_0, g_x_0_zzzzz_0, g_x_0_zzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzzzz_0[k] = -g_x_0_zzzzz_0[k] * cd_y[k] + g_x_0_zzzzz_y[k];
            }

            /// Set up 27-28 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzzzzz_0 = cbuffer.data(is_geom_10_off + 0 * acomps * bcomps + 27);

            #pragma omp simd aligned(cd_z, g_x_0_zzzzz_0, g_x_0_zzzzz_z, g_x_0_zzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzzzz_0[k] = -g_x_0_zzzzz_0[k] * cd_z[k] + g_x_0_zzzzz_z[k];
            }
            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxx_0 = cbuffer.data(is_geom_10_off + 28 * acomps * bcomps + 0);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxx_0, g_y_0_xxxxx_x, g_y_0_xxxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxx_0[k] = -g_y_0_xxxxx_0[k] * cd_x[k] + g_y_0_xxxxx_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxy_0 = cbuffer.data(is_geom_10_off + 28 * acomps * bcomps + 1);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxxy_0, g_y_0_xxxxy_0, g_y_0_xxxxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxy_0[k] = -g_y_0_xxxxy_0[k] * cd_x[k] + g_y_0_xxxxy_x[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxz_0 = cbuffer.data(is_geom_10_off + 28 * acomps * bcomps + 2);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxxz_0, g_y_0_xxxxz_0, g_y_0_xxxxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxz_0[k] = -g_y_0_xxxxz_0[k] * cd_x[k] + g_y_0_xxxxz_x[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxyy_0 = cbuffer.data(is_geom_10_off + 28 * acomps * bcomps + 3);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxyy_0, g_y_0_xxxyy_0, g_y_0_xxxyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxyy_0[k] = -g_y_0_xxxyy_0[k] * cd_x[k] + g_y_0_xxxyy_x[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxyz_0 = cbuffer.data(is_geom_10_off + 28 * acomps * bcomps + 4);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxyz_0, g_y_0_xxxyz_0, g_y_0_xxxyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxyz_0[k] = -g_y_0_xxxyz_0[k] * cd_x[k] + g_y_0_xxxyz_x[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxzz_0 = cbuffer.data(is_geom_10_off + 28 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxzz_0, g_y_0_xxxzz_0, g_y_0_xxxzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxzz_0[k] = -g_y_0_xxxzz_0[k] * cd_x[k] + g_y_0_xxxzz_x[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyyy_0 = cbuffer.data(is_geom_10_off + 28 * acomps * bcomps + 6);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyyy_0, g_y_0_xxyyy_0, g_y_0_xxyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyyy_0[k] = -g_y_0_xxyyy_0[k] * cd_x[k] + g_y_0_xxyyy_x[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyyz_0 = cbuffer.data(is_geom_10_off + 28 * acomps * bcomps + 7);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyyz_0, g_y_0_xxyyz_0, g_y_0_xxyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyyz_0[k] = -g_y_0_xxyyz_0[k] * cd_x[k] + g_y_0_xxyyz_x[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyzz_0 = cbuffer.data(is_geom_10_off + 28 * acomps * bcomps + 8);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyzz_0, g_y_0_xxyzz_0, g_y_0_xxyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyzz_0[k] = -g_y_0_xxyzz_0[k] * cd_x[k] + g_y_0_xxyzz_x[k];
            }

            /// Set up 9-10 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxzzz_0 = cbuffer.data(is_geom_10_off + 28 * acomps * bcomps + 9);

            #pragma omp simd aligned(cd_x, g_y_0_xxxzzz_0, g_y_0_xxzzz_0, g_y_0_xxzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxzzz_0[k] = -g_y_0_xxzzz_0[k] * cd_x[k] + g_y_0_xxzzz_x[k];
            }

            /// Set up 10-11 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyyy_0 = cbuffer.data(is_geom_10_off + 28 * acomps * bcomps + 10);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyyy_0, g_y_0_xyyyy_0, g_y_0_xyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyyy_0[k] = -g_y_0_xyyyy_0[k] * cd_x[k] + g_y_0_xyyyy_x[k];
            }

            /// Set up 11-12 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyyz_0 = cbuffer.data(is_geom_10_off + 28 * acomps * bcomps + 11);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyyz_0, g_y_0_xyyyz_0, g_y_0_xyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyyz_0[k] = -g_y_0_xyyyz_0[k] * cd_x[k] + g_y_0_xyyyz_x[k];
            }

            /// Set up 12-13 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyzz_0 = cbuffer.data(is_geom_10_off + 28 * acomps * bcomps + 12);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyzz_0, g_y_0_xyyzz_0, g_y_0_xyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyzz_0[k] = -g_y_0_xyyzz_0[k] * cd_x[k] + g_y_0_xyyzz_x[k];
            }

            /// Set up 13-14 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyzzz_0 = cbuffer.data(is_geom_10_off + 28 * acomps * bcomps + 13);

            #pragma omp simd aligned(cd_x, g_y_0_xxyzzz_0, g_y_0_xyzzz_0, g_y_0_xyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyzzz_0[k] = -g_y_0_xyzzz_0[k] * cd_x[k] + g_y_0_xyzzz_x[k];
            }

            /// Set up 14-15 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxzzzz_0 = cbuffer.data(is_geom_10_off + 28 * acomps * bcomps + 14);

            #pragma omp simd aligned(cd_x, g_y_0_xxzzzz_0, g_y_0_xzzzz_0, g_y_0_xzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzzzz_0[k] = -g_y_0_xzzzz_0[k] * cd_x[k] + g_y_0_xzzzz_x[k];
            }

            /// Set up 15-16 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyyy_0 = cbuffer.data(is_geom_10_off + 28 * acomps * bcomps + 15);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyyy_0, g_y_0_yyyyy_0, g_y_0_yyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyyy_0[k] = -g_y_0_yyyyy_0[k] * cd_x[k] + g_y_0_yyyyy_x[k];
            }

            /// Set up 16-17 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyyz_0 = cbuffer.data(is_geom_10_off + 28 * acomps * bcomps + 16);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyyz_0, g_y_0_yyyyz_0, g_y_0_yyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyyz_0[k] = -g_y_0_yyyyz_0[k] * cd_x[k] + g_y_0_yyyyz_x[k];
            }

            /// Set up 17-18 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyzz_0 = cbuffer.data(is_geom_10_off + 28 * acomps * bcomps + 17);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyzz_0, g_y_0_yyyzz_0, g_y_0_yyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyzz_0[k] = -g_y_0_yyyzz_0[k] * cd_x[k] + g_y_0_yyyzz_x[k];
            }

            /// Set up 18-19 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyzzz_0 = cbuffer.data(is_geom_10_off + 28 * acomps * bcomps + 18);

            #pragma omp simd aligned(cd_x, g_y_0_xyyzzz_0, g_y_0_yyzzz_0, g_y_0_yyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyzzz_0[k] = -g_y_0_yyzzz_0[k] * cd_x[k] + g_y_0_yyzzz_x[k];
            }

            /// Set up 19-20 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyzzzz_0 = cbuffer.data(is_geom_10_off + 28 * acomps * bcomps + 19);

            #pragma omp simd aligned(cd_x, g_y_0_xyzzzz_0, g_y_0_yzzzz_0, g_y_0_yzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzzzz_0[k] = -g_y_0_yzzzz_0[k] * cd_x[k] + g_y_0_yzzzz_x[k];
            }

            /// Set up 20-21 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzzzzz_0 = cbuffer.data(is_geom_10_off + 28 * acomps * bcomps + 20);

            #pragma omp simd aligned(cd_x, g_y_0_xzzzzz_0, g_y_0_zzzzz_0, g_y_0_zzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzzzz_0[k] = -g_y_0_zzzzz_0[k] * cd_x[k] + g_y_0_zzzzz_x[k];
            }

            /// Set up 21-22 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyyy_0 = cbuffer.data(is_geom_10_off + 28 * acomps * bcomps + 21);

            #pragma omp simd aligned(cd_y, g_y_0_yyyyy_0, g_y_0_yyyyy_y, g_y_0_yyyyyy_0, g_yyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyyy_0[k] = -g_yyyyy_0[k] - g_y_0_yyyyy_0[k] * cd_y[k] + g_y_0_yyyyy_y[k];
            }

            /// Set up 22-23 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyyz_0 = cbuffer.data(is_geom_10_off + 28 * acomps * bcomps + 22);

            #pragma omp simd aligned(cd_z, g_y_0_yyyyy_0, g_y_0_yyyyy_z, g_y_0_yyyyyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyyz_0[k] = -g_y_0_yyyyy_0[k] * cd_z[k] + g_y_0_yyyyy_z[k];
            }

            /// Set up 23-24 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyzz_0 = cbuffer.data(is_geom_10_off + 28 * acomps * bcomps + 23);

            #pragma omp simd aligned(cd_z, g_y_0_yyyyz_0, g_y_0_yyyyz_z, g_y_0_yyyyzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyzz_0[k] = -g_y_0_yyyyz_0[k] * cd_z[k] + g_y_0_yyyyz_z[k];
            }

            /// Set up 24-25 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyzzz_0 = cbuffer.data(is_geom_10_off + 28 * acomps * bcomps + 24);

            #pragma omp simd aligned(cd_z, g_y_0_yyyzz_0, g_y_0_yyyzz_z, g_y_0_yyyzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyzzz_0[k] = -g_y_0_yyyzz_0[k] * cd_z[k] + g_y_0_yyyzz_z[k];
            }

            /// Set up 25-26 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyzzzz_0 = cbuffer.data(is_geom_10_off + 28 * acomps * bcomps + 25);

            #pragma omp simd aligned(cd_z, g_y_0_yyzzz_0, g_y_0_yyzzz_z, g_y_0_yyzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzzzz_0[k] = -g_y_0_yyzzz_0[k] * cd_z[k] + g_y_0_yyzzz_z[k];
            }

            /// Set up 26-27 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzzzzz_0 = cbuffer.data(is_geom_10_off + 28 * acomps * bcomps + 26);

            #pragma omp simd aligned(cd_z, g_y_0_yzzzz_0, g_y_0_yzzzz_z, g_y_0_yzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzzzz_0[k] = -g_y_0_yzzzz_0[k] * cd_z[k] + g_y_0_yzzzz_z[k];
            }

            /// Set up 27-28 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzzzzz_0 = cbuffer.data(is_geom_10_off + 28 * acomps * bcomps + 27);

            #pragma omp simd aligned(cd_z, g_y_0_zzzzz_0, g_y_0_zzzzz_z, g_y_0_zzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzzzz_0[k] = -g_y_0_zzzzz_0[k] * cd_z[k] + g_y_0_zzzzz_z[k];
            }
            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxx_0 = cbuffer.data(is_geom_10_off + 56 * acomps * bcomps + 0);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxx_0, g_z_0_xxxxx_x, g_z_0_xxxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxx_0[k] = -g_z_0_xxxxx_0[k] * cd_x[k] + g_z_0_xxxxx_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxy_0 = cbuffer.data(is_geom_10_off + 56 * acomps * bcomps + 1);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxxy_0, g_z_0_xxxxy_0, g_z_0_xxxxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxy_0[k] = -g_z_0_xxxxy_0[k] * cd_x[k] + g_z_0_xxxxy_x[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxz_0 = cbuffer.data(is_geom_10_off + 56 * acomps * bcomps + 2);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxxz_0, g_z_0_xxxxz_0, g_z_0_xxxxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxz_0[k] = -g_z_0_xxxxz_0[k] * cd_x[k] + g_z_0_xxxxz_x[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxyy_0 = cbuffer.data(is_geom_10_off + 56 * acomps * bcomps + 3);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxyy_0, g_z_0_xxxyy_0, g_z_0_xxxyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxyy_0[k] = -g_z_0_xxxyy_0[k] * cd_x[k] + g_z_0_xxxyy_x[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxyz_0 = cbuffer.data(is_geom_10_off + 56 * acomps * bcomps + 4);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxyz_0, g_z_0_xxxyz_0, g_z_0_xxxyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxyz_0[k] = -g_z_0_xxxyz_0[k] * cd_x[k] + g_z_0_xxxyz_x[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxzz_0 = cbuffer.data(is_geom_10_off + 56 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxzz_0, g_z_0_xxxzz_0, g_z_0_xxxzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxzz_0[k] = -g_z_0_xxxzz_0[k] * cd_x[k] + g_z_0_xxxzz_x[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyyy_0 = cbuffer.data(is_geom_10_off + 56 * acomps * bcomps + 6);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyyy_0, g_z_0_xxyyy_0, g_z_0_xxyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyyy_0[k] = -g_z_0_xxyyy_0[k] * cd_x[k] + g_z_0_xxyyy_x[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyyz_0 = cbuffer.data(is_geom_10_off + 56 * acomps * bcomps + 7);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyyz_0, g_z_0_xxyyz_0, g_z_0_xxyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyyz_0[k] = -g_z_0_xxyyz_0[k] * cd_x[k] + g_z_0_xxyyz_x[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyzz_0 = cbuffer.data(is_geom_10_off + 56 * acomps * bcomps + 8);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyzz_0, g_z_0_xxyzz_0, g_z_0_xxyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyzz_0[k] = -g_z_0_xxyzz_0[k] * cd_x[k] + g_z_0_xxyzz_x[k];
            }

            /// Set up 9-10 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxzzz_0 = cbuffer.data(is_geom_10_off + 56 * acomps * bcomps + 9);

            #pragma omp simd aligned(cd_x, g_z_0_xxxzzz_0, g_z_0_xxzzz_0, g_z_0_xxzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxzzz_0[k] = -g_z_0_xxzzz_0[k] * cd_x[k] + g_z_0_xxzzz_x[k];
            }

            /// Set up 10-11 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyyy_0 = cbuffer.data(is_geom_10_off + 56 * acomps * bcomps + 10);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyyy_0, g_z_0_xyyyy_0, g_z_0_xyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyyy_0[k] = -g_z_0_xyyyy_0[k] * cd_x[k] + g_z_0_xyyyy_x[k];
            }

            /// Set up 11-12 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyyz_0 = cbuffer.data(is_geom_10_off + 56 * acomps * bcomps + 11);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyyz_0, g_z_0_xyyyz_0, g_z_0_xyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyyz_0[k] = -g_z_0_xyyyz_0[k] * cd_x[k] + g_z_0_xyyyz_x[k];
            }

            /// Set up 12-13 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyzz_0 = cbuffer.data(is_geom_10_off + 56 * acomps * bcomps + 12);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyzz_0, g_z_0_xyyzz_0, g_z_0_xyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyzz_0[k] = -g_z_0_xyyzz_0[k] * cd_x[k] + g_z_0_xyyzz_x[k];
            }

            /// Set up 13-14 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyzzz_0 = cbuffer.data(is_geom_10_off + 56 * acomps * bcomps + 13);

            #pragma omp simd aligned(cd_x, g_z_0_xxyzzz_0, g_z_0_xyzzz_0, g_z_0_xyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyzzz_0[k] = -g_z_0_xyzzz_0[k] * cd_x[k] + g_z_0_xyzzz_x[k];
            }

            /// Set up 14-15 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxzzzz_0 = cbuffer.data(is_geom_10_off + 56 * acomps * bcomps + 14);

            #pragma omp simd aligned(cd_x, g_z_0_xxzzzz_0, g_z_0_xzzzz_0, g_z_0_xzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzzzz_0[k] = -g_z_0_xzzzz_0[k] * cd_x[k] + g_z_0_xzzzz_x[k];
            }

            /// Set up 15-16 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyyy_0 = cbuffer.data(is_geom_10_off + 56 * acomps * bcomps + 15);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyyy_0, g_z_0_yyyyy_0, g_z_0_yyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyyy_0[k] = -g_z_0_yyyyy_0[k] * cd_x[k] + g_z_0_yyyyy_x[k];
            }

            /// Set up 16-17 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyyz_0 = cbuffer.data(is_geom_10_off + 56 * acomps * bcomps + 16);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyyz_0, g_z_0_yyyyz_0, g_z_0_yyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyyz_0[k] = -g_z_0_yyyyz_0[k] * cd_x[k] + g_z_0_yyyyz_x[k];
            }

            /// Set up 17-18 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyzz_0 = cbuffer.data(is_geom_10_off + 56 * acomps * bcomps + 17);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyzz_0, g_z_0_yyyzz_0, g_z_0_yyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyzz_0[k] = -g_z_0_yyyzz_0[k] * cd_x[k] + g_z_0_yyyzz_x[k];
            }

            /// Set up 18-19 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyzzz_0 = cbuffer.data(is_geom_10_off + 56 * acomps * bcomps + 18);

            #pragma omp simd aligned(cd_x, g_z_0_xyyzzz_0, g_z_0_yyzzz_0, g_z_0_yyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyzzz_0[k] = -g_z_0_yyzzz_0[k] * cd_x[k] + g_z_0_yyzzz_x[k];
            }

            /// Set up 19-20 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyzzzz_0 = cbuffer.data(is_geom_10_off + 56 * acomps * bcomps + 19);

            #pragma omp simd aligned(cd_x, g_z_0_xyzzzz_0, g_z_0_yzzzz_0, g_z_0_yzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzzzz_0[k] = -g_z_0_yzzzz_0[k] * cd_x[k] + g_z_0_yzzzz_x[k];
            }

            /// Set up 20-21 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzzzzz_0 = cbuffer.data(is_geom_10_off + 56 * acomps * bcomps + 20);

            #pragma omp simd aligned(cd_x, g_z_0_xzzzzz_0, g_z_0_zzzzz_0, g_z_0_zzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzzzz_0[k] = -g_z_0_zzzzz_0[k] * cd_x[k] + g_z_0_zzzzz_x[k];
            }

            /// Set up 21-22 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyyy_0 = cbuffer.data(is_geom_10_off + 56 * acomps * bcomps + 21);

            #pragma omp simd aligned(cd_y, g_z_0_yyyyy_0, g_z_0_yyyyy_y, g_z_0_yyyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyyy_0[k] = -g_z_0_yyyyy_0[k] * cd_y[k] + g_z_0_yyyyy_y[k];
            }

            /// Set up 22-23 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyyz_0 = cbuffer.data(is_geom_10_off + 56 * acomps * bcomps + 22);

            #pragma omp simd aligned(cd_y, g_z_0_yyyyyz_0, g_z_0_yyyyz_0, g_z_0_yyyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyyz_0[k] = -g_z_0_yyyyz_0[k] * cd_y[k] + g_z_0_yyyyz_y[k];
            }

            /// Set up 23-24 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyzz_0 = cbuffer.data(is_geom_10_off + 56 * acomps * bcomps + 23);

            #pragma omp simd aligned(cd_y, g_z_0_yyyyzz_0, g_z_0_yyyzz_0, g_z_0_yyyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyzz_0[k] = -g_z_0_yyyzz_0[k] * cd_y[k] + g_z_0_yyyzz_y[k];
            }

            /// Set up 24-25 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyzzz_0 = cbuffer.data(is_geom_10_off + 56 * acomps * bcomps + 24);

            #pragma omp simd aligned(cd_y, g_z_0_yyyzzz_0, g_z_0_yyzzz_0, g_z_0_yyzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyzzz_0[k] = -g_z_0_yyzzz_0[k] * cd_y[k] + g_z_0_yyzzz_y[k];
            }

            /// Set up 25-26 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyzzzz_0 = cbuffer.data(is_geom_10_off + 56 * acomps * bcomps + 25);

            #pragma omp simd aligned(cd_y, g_z_0_yyzzzz_0, g_z_0_yzzzz_0, g_z_0_yzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzzzz_0[k] = -g_z_0_yzzzz_0[k] * cd_y[k] + g_z_0_yzzzz_y[k];
            }

            /// Set up 26-27 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzzzzz_0 = cbuffer.data(is_geom_10_off + 56 * acomps * bcomps + 26);

            #pragma omp simd aligned(cd_y, g_z_0_yzzzzz_0, g_z_0_zzzzz_0, g_z_0_zzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzzzz_0[k] = -g_z_0_zzzzz_0[k] * cd_y[k] + g_z_0_zzzzz_y[k];
            }

            /// Set up 27-28 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzzzzz_0 = cbuffer.data(is_geom_10_off + 56 * acomps * bcomps + 27);

            #pragma omp simd aligned(cd_z, g_z_0_zzzzz_0, g_z_0_zzzzz_z, g_z_0_zzzzzz_0, g_zzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzzzz_0[k] = -g_zzzzz_0[k] - g_z_0_zzzzz_0[k] * cd_z[k] + g_z_0_zzzzz_z[k];
            }
        }
    }
}

} // erirec namespace

