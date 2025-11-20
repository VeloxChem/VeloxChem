#include "ElectronRepulsionGeom0010ContrRecXXIP.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxip(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxip,
                                            const size_t idx_xxhp,
                                            const size_t idx_geom_10_xxhp,
                                            const size_t idx_geom_10_xxhd,
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
            /// Set up components of auxilary buffer : SSHP

            const auto hp_off = idx_xxhp + (i * bcomps + j) * 63;

            auto g_xxxxx_x = cbuffer.data(hp_off + 0);

            auto g_xxxxx_y = cbuffer.data(hp_off + 1);

            auto g_xxxxx_z = cbuffer.data(hp_off + 2);

            auto g_xxxxy_x = cbuffer.data(hp_off + 3);

            auto g_xxxxy_y = cbuffer.data(hp_off + 4);

            auto g_xxxxy_z = cbuffer.data(hp_off + 5);

            auto g_xxxxz_x = cbuffer.data(hp_off + 6);

            auto g_xxxxz_y = cbuffer.data(hp_off + 7);

            auto g_xxxxz_z = cbuffer.data(hp_off + 8);

            auto g_xxxyy_x = cbuffer.data(hp_off + 9);

            auto g_xxxyy_y = cbuffer.data(hp_off + 10);

            auto g_xxxyy_z = cbuffer.data(hp_off + 11);

            auto g_xxxyz_x = cbuffer.data(hp_off + 12);

            auto g_xxxyz_y = cbuffer.data(hp_off + 13);

            auto g_xxxyz_z = cbuffer.data(hp_off + 14);

            auto g_xxxzz_x = cbuffer.data(hp_off + 15);

            auto g_xxxzz_y = cbuffer.data(hp_off + 16);

            auto g_xxxzz_z = cbuffer.data(hp_off + 17);

            auto g_xxyyy_x = cbuffer.data(hp_off + 18);

            auto g_xxyyy_y = cbuffer.data(hp_off + 19);

            auto g_xxyyy_z = cbuffer.data(hp_off + 20);

            auto g_xxyyz_x = cbuffer.data(hp_off + 21);

            auto g_xxyyz_y = cbuffer.data(hp_off + 22);

            auto g_xxyyz_z = cbuffer.data(hp_off + 23);

            auto g_xxyzz_x = cbuffer.data(hp_off + 24);

            auto g_xxyzz_y = cbuffer.data(hp_off + 25);

            auto g_xxyzz_z = cbuffer.data(hp_off + 26);

            auto g_xxzzz_x = cbuffer.data(hp_off + 27);

            auto g_xxzzz_y = cbuffer.data(hp_off + 28);

            auto g_xxzzz_z = cbuffer.data(hp_off + 29);

            auto g_xyyyy_x = cbuffer.data(hp_off + 30);

            auto g_xyyyy_y = cbuffer.data(hp_off + 31);

            auto g_xyyyy_z = cbuffer.data(hp_off + 32);

            auto g_xyyyz_x = cbuffer.data(hp_off + 33);

            auto g_xyyyz_y = cbuffer.data(hp_off + 34);

            auto g_xyyyz_z = cbuffer.data(hp_off + 35);

            auto g_xyyzz_x = cbuffer.data(hp_off + 36);

            auto g_xyyzz_y = cbuffer.data(hp_off + 37);

            auto g_xyyzz_z = cbuffer.data(hp_off + 38);

            auto g_xyzzz_x = cbuffer.data(hp_off + 39);

            auto g_xyzzz_y = cbuffer.data(hp_off + 40);

            auto g_xyzzz_z = cbuffer.data(hp_off + 41);

            auto g_xzzzz_x = cbuffer.data(hp_off + 42);

            auto g_xzzzz_y = cbuffer.data(hp_off + 43);

            auto g_xzzzz_z = cbuffer.data(hp_off + 44);

            auto g_yyyyy_x = cbuffer.data(hp_off + 45);

            auto g_yyyyy_y = cbuffer.data(hp_off + 46);

            auto g_yyyyy_z = cbuffer.data(hp_off + 47);

            auto g_yyyyz_x = cbuffer.data(hp_off + 48);

            auto g_yyyyz_y = cbuffer.data(hp_off + 49);

            auto g_yyyyz_z = cbuffer.data(hp_off + 50);

            auto g_yyyzz_x = cbuffer.data(hp_off + 51);

            auto g_yyyzz_y = cbuffer.data(hp_off + 52);

            auto g_yyyzz_z = cbuffer.data(hp_off + 53);

            auto g_yyzzz_x = cbuffer.data(hp_off + 54);

            auto g_yyzzz_y = cbuffer.data(hp_off + 55);

            auto g_yyzzz_z = cbuffer.data(hp_off + 56);

            auto g_yzzzz_x = cbuffer.data(hp_off + 57);

            auto g_yzzzz_y = cbuffer.data(hp_off + 58);

            auto g_yzzzz_z = cbuffer.data(hp_off + 59);

            auto g_zzzzz_x = cbuffer.data(hp_off + 60);

            auto g_zzzzz_y = cbuffer.data(hp_off + 61);

            auto g_zzzzz_z = cbuffer.data(hp_off + 62);

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

            /// Set up components of auxilary buffer : SSHD

            const auto hd_geom_10_off = idx_geom_10_xxhd + (i * bcomps + j) * 126;

            auto g_x_0_xxxxx_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxxx_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxxx_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxxx_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxxx_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxxx_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxxxy_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxxy_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxxy_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxxxy_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxxxy_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxxxy_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxxxz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxxxz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxxxz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxxxz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxxxz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxxz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xxxyy_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxxyy_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxxyy_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xxxyy_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxxyy_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxxyy_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xxxyz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxxyz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxxyz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxxyz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxxyz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxxyz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_xxxzz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xxxzz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxxzz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxxzz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxxzz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxxzz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xxyyy_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xxyyy_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xxyyy_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xxyyy_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xxyyy_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxyyy_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_xxyyz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xxyyz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxyyz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xxyyz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xxyyz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xxyyz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xxyzz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xxyzz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xxyzz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xxyzz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xxyzz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xxyzz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xxzzz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xxzzz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xxzzz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xxzzz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xxzzz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xxzzz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 59);

            auto g_x_0_xyyyy_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xyyyy_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xyyyy_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_x_0_xyyyy_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xyyyy_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xyyyy_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 65);

            auto g_x_0_xyyyz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xyyyz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xyyyz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xyyyz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_xyyyz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xyyyz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_xyyzz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xyyzz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xyyzz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_xyyzz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_xyyzz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xyyzz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_xyzzz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xyzzz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_xyzzz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xyzzz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xyzzz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xyzzz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 83);

            auto g_x_0_xzzzz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_xzzzz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xzzzz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_xzzzz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xzzzz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xzzzz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 89);

            auto g_x_0_yyyyy_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_yyyyy_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_yyyyy_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_yyyyy_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_yyyyy_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_yyyyy_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 95);

            auto g_x_0_yyyyz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_yyyyz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_yyyyz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_yyyyz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_x_0_yyyyz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_yyyyz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_yyyzz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_yyyzz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_yyyzz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 104);

            auto g_x_0_yyyzz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_yyyzz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_yyyzz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 107);

            auto g_x_0_yyzzz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_yyzzz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_yyzzz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_yyzzz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 111);

            auto g_x_0_yyzzz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_yyzzz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_yzzzz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_yzzzz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_yzzzz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_yzzzz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_yzzzz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_yzzzz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 119);

            auto g_x_0_zzzzz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_zzzzz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_zzzzz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_zzzzz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_zzzzz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_zzzzz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 125);

            auto g_y_0_xxxxx_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 0);

            auto g_y_0_xxxxx_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 1);

            auto g_y_0_xxxxx_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 2);

            auto g_y_0_xxxxx_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 3);

            auto g_y_0_xxxxx_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 4);

            auto g_y_0_xxxxx_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 5);

            auto g_y_0_xxxxy_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 6);

            auto g_y_0_xxxxy_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 7);

            auto g_y_0_xxxxy_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 8);

            auto g_y_0_xxxxy_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 9);

            auto g_y_0_xxxxy_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 10);

            auto g_y_0_xxxxy_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 11);

            auto g_y_0_xxxxz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 12);

            auto g_y_0_xxxxz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 13);

            auto g_y_0_xxxxz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 14);

            auto g_y_0_xxxxz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 15);

            auto g_y_0_xxxxz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 16);

            auto g_y_0_xxxxz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 17);

            auto g_y_0_xxxyy_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 18);

            auto g_y_0_xxxyy_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 19);

            auto g_y_0_xxxyy_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 20);

            auto g_y_0_xxxyy_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 21);

            auto g_y_0_xxxyy_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 22);

            auto g_y_0_xxxyy_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 23);

            auto g_y_0_xxxyz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 24);

            auto g_y_0_xxxyz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 25);

            auto g_y_0_xxxyz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 26);

            auto g_y_0_xxxyz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 27);

            auto g_y_0_xxxyz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 28);

            auto g_y_0_xxxyz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 29);

            auto g_y_0_xxxzz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 30);

            auto g_y_0_xxxzz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 31);

            auto g_y_0_xxxzz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 32);

            auto g_y_0_xxxzz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 33);

            auto g_y_0_xxxzz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 34);

            auto g_y_0_xxxzz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 35);

            auto g_y_0_xxyyy_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 36);

            auto g_y_0_xxyyy_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 37);

            auto g_y_0_xxyyy_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 38);

            auto g_y_0_xxyyy_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 39);

            auto g_y_0_xxyyy_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 40);

            auto g_y_0_xxyyy_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 41);

            auto g_y_0_xxyyz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 42);

            auto g_y_0_xxyyz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 43);

            auto g_y_0_xxyyz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 44);

            auto g_y_0_xxyyz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 45);

            auto g_y_0_xxyyz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 46);

            auto g_y_0_xxyyz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 47);

            auto g_y_0_xxyzz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 48);

            auto g_y_0_xxyzz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 49);

            auto g_y_0_xxyzz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 50);

            auto g_y_0_xxyzz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 51);

            auto g_y_0_xxyzz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 52);

            auto g_y_0_xxyzz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 53);

            auto g_y_0_xxzzz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 54);

            auto g_y_0_xxzzz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 55);

            auto g_y_0_xxzzz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 56);

            auto g_y_0_xxzzz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 57);

            auto g_y_0_xxzzz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 58);

            auto g_y_0_xxzzz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 59);

            auto g_y_0_xyyyy_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 60);

            auto g_y_0_xyyyy_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 61);

            auto g_y_0_xyyyy_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 62);

            auto g_y_0_xyyyy_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 63);

            auto g_y_0_xyyyy_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 64);

            auto g_y_0_xyyyy_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 65);

            auto g_y_0_xyyyz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 66);

            auto g_y_0_xyyyz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 67);

            auto g_y_0_xyyyz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 68);

            auto g_y_0_xyyyz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 69);

            auto g_y_0_xyyyz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 70);

            auto g_y_0_xyyyz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 71);

            auto g_y_0_xyyzz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 72);

            auto g_y_0_xyyzz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 73);

            auto g_y_0_xyyzz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 74);

            auto g_y_0_xyyzz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 75);

            auto g_y_0_xyyzz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 76);

            auto g_y_0_xyyzz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 77);

            auto g_y_0_xyzzz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 78);

            auto g_y_0_xyzzz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 79);

            auto g_y_0_xyzzz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 80);

            auto g_y_0_xyzzz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 81);

            auto g_y_0_xyzzz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 82);

            auto g_y_0_xyzzz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 83);

            auto g_y_0_xzzzz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 84);

            auto g_y_0_xzzzz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 85);

            auto g_y_0_xzzzz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 86);

            auto g_y_0_xzzzz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 87);

            auto g_y_0_xzzzz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 88);

            auto g_y_0_xzzzz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 89);

            auto g_y_0_yyyyy_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 90);

            auto g_y_0_yyyyy_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 91);

            auto g_y_0_yyyyy_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 92);

            auto g_y_0_yyyyy_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 93);

            auto g_y_0_yyyyy_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 94);

            auto g_y_0_yyyyy_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 95);

            auto g_y_0_yyyyz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 96);

            auto g_y_0_yyyyz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 97);

            auto g_y_0_yyyyz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 98);

            auto g_y_0_yyyyz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 99);

            auto g_y_0_yyyyz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 100);

            auto g_y_0_yyyyz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 101);

            auto g_y_0_yyyzz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 102);

            auto g_y_0_yyyzz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 103);

            auto g_y_0_yyyzz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 104);

            auto g_y_0_yyyzz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 105);

            auto g_y_0_yyyzz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 106);

            auto g_y_0_yyyzz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 107);

            auto g_y_0_yyzzz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 108);

            auto g_y_0_yyzzz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 109);

            auto g_y_0_yyzzz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 110);

            auto g_y_0_yyzzz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 111);

            auto g_y_0_yyzzz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 112);

            auto g_y_0_yyzzz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 113);

            auto g_y_0_yzzzz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 114);

            auto g_y_0_yzzzz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 115);

            auto g_y_0_yzzzz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 116);

            auto g_y_0_yzzzz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 117);

            auto g_y_0_yzzzz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 118);

            auto g_y_0_yzzzz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 119);

            auto g_y_0_zzzzz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 120);

            auto g_y_0_zzzzz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 121);

            auto g_y_0_zzzzz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 122);

            auto g_y_0_zzzzz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 123);

            auto g_y_0_zzzzz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 124);

            auto g_y_0_zzzzz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 125);

            auto g_z_0_xxxxx_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 0);

            auto g_z_0_xxxxx_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 1);

            auto g_z_0_xxxxx_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 2);

            auto g_z_0_xxxxx_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 3);

            auto g_z_0_xxxxx_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 4);

            auto g_z_0_xxxxx_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 5);

            auto g_z_0_xxxxy_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 6);

            auto g_z_0_xxxxy_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 7);

            auto g_z_0_xxxxy_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 8);

            auto g_z_0_xxxxy_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 9);

            auto g_z_0_xxxxy_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 10);

            auto g_z_0_xxxxy_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 11);

            auto g_z_0_xxxxz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 12);

            auto g_z_0_xxxxz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 13);

            auto g_z_0_xxxxz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 14);

            auto g_z_0_xxxxz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 15);

            auto g_z_0_xxxxz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 16);

            auto g_z_0_xxxxz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 17);

            auto g_z_0_xxxyy_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 18);

            auto g_z_0_xxxyy_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 19);

            auto g_z_0_xxxyy_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 20);

            auto g_z_0_xxxyy_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 21);

            auto g_z_0_xxxyy_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 22);

            auto g_z_0_xxxyy_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 23);

            auto g_z_0_xxxyz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 24);

            auto g_z_0_xxxyz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 25);

            auto g_z_0_xxxyz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 26);

            auto g_z_0_xxxyz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 27);

            auto g_z_0_xxxyz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 28);

            auto g_z_0_xxxyz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 29);

            auto g_z_0_xxxzz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 30);

            auto g_z_0_xxxzz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 31);

            auto g_z_0_xxxzz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 32);

            auto g_z_0_xxxzz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 33);

            auto g_z_0_xxxzz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 34);

            auto g_z_0_xxxzz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 35);

            auto g_z_0_xxyyy_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 36);

            auto g_z_0_xxyyy_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 37);

            auto g_z_0_xxyyy_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 38);

            auto g_z_0_xxyyy_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 39);

            auto g_z_0_xxyyy_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 40);

            auto g_z_0_xxyyy_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 41);

            auto g_z_0_xxyyz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 42);

            auto g_z_0_xxyyz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 43);

            auto g_z_0_xxyyz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 44);

            auto g_z_0_xxyyz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 45);

            auto g_z_0_xxyyz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 46);

            auto g_z_0_xxyyz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 47);

            auto g_z_0_xxyzz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 48);

            auto g_z_0_xxyzz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 49);

            auto g_z_0_xxyzz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 50);

            auto g_z_0_xxyzz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 51);

            auto g_z_0_xxyzz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 52);

            auto g_z_0_xxyzz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 53);

            auto g_z_0_xxzzz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 54);

            auto g_z_0_xxzzz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 55);

            auto g_z_0_xxzzz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 56);

            auto g_z_0_xxzzz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 57);

            auto g_z_0_xxzzz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 58);

            auto g_z_0_xxzzz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 59);

            auto g_z_0_xyyyy_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 60);

            auto g_z_0_xyyyy_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 61);

            auto g_z_0_xyyyy_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 62);

            auto g_z_0_xyyyy_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 63);

            auto g_z_0_xyyyy_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 64);

            auto g_z_0_xyyyy_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 65);

            auto g_z_0_xyyyz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 66);

            auto g_z_0_xyyyz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 67);

            auto g_z_0_xyyyz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 68);

            auto g_z_0_xyyyz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 69);

            auto g_z_0_xyyyz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 70);

            auto g_z_0_xyyyz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 71);

            auto g_z_0_xyyzz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 72);

            auto g_z_0_xyyzz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 73);

            auto g_z_0_xyyzz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 74);

            auto g_z_0_xyyzz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 75);

            auto g_z_0_xyyzz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 76);

            auto g_z_0_xyyzz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 77);

            auto g_z_0_xyzzz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 78);

            auto g_z_0_xyzzz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 79);

            auto g_z_0_xyzzz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 80);

            auto g_z_0_xyzzz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 81);

            auto g_z_0_xyzzz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 82);

            auto g_z_0_xyzzz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 83);

            auto g_z_0_xzzzz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 84);

            auto g_z_0_xzzzz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 85);

            auto g_z_0_xzzzz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 86);

            auto g_z_0_xzzzz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 87);

            auto g_z_0_xzzzz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 88);

            auto g_z_0_xzzzz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 89);

            auto g_z_0_yyyyy_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 90);

            auto g_z_0_yyyyy_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 91);

            auto g_z_0_yyyyy_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 92);

            auto g_z_0_yyyyy_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 93);

            auto g_z_0_yyyyy_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 94);

            auto g_z_0_yyyyy_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 95);

            auto g_z_0_yyyyz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 96);

            auto g_z_0_yyyyz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 97);

            auto g_z_0_yyyyz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 98);

            auto g_z_0_yyyyz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 99);

            auto g_z_0_yyyyz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 100);

            auto g_z_0_yyyyz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 101);

            auto g_z_0_yyyzz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 102);

            auto g_z_0_yyyzz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 103);

            auto g_z_0_yyyzz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 104);

            auto g_z_0_yyyzz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 105);

            auto g_z_0_yyyzz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 106);

            auto g_z_0_yyyzz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 107);

            auto g_z_0_yyzzz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 108);

            auto g_z_0_yyzzz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 109);

            auto g_z_0_yyzzz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 110);

            auto g_z_0_yyzzz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 111);

            auto g_z_0_yyzzz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 112);

            auto g_z_0_yyzzz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 113);

            auto g_z_0_yzzzz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 114);

            auto g_z_0_yzzzz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 115);

            auto g_z_0_yzzzz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 116);

            auto g_z_0_yzzzz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 117);

            auto g_z_0_yzzzz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 118);

            auto g_z_0_yzzzz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 119);

            auto g_z_0_zzzzz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 120);

            auto g_z_0_zzzzz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 121);

            auto g_z_0_zzzzz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 122);

            auto g_z_0_zzzzz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 123);

            auto g_z_0_zzzzz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 124);

            auto g_z_0_zzzzz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 125);

            /// set up bra offset for contr_buffer_xxip

            const auto ip_geom_10_off = idx_geom_10_xxip + (i * bcomps + j) * 84;

            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxx_x = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxxxx_y = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxxxx_z = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 2);

            #pragma omp simd aligned(cd_x, g_x_0_xxxxx_x, g_x_0_xxxxx_xx, g_x_0_xxxxx_xy, g_x_0_xxxxx_xz, g_x_0_xxxxx_y, g_x_0_xxxxx_z, g_x_0_xxxxxx_x, g_x_0_xxxxxx_y, g_x_0_xxxxxx_z, g_xxxxx_x, g_xxxxx_y, g_xxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxx_x[k] = -g_xxxxx_x[k] - g_x_0_xxxxx_x[k] * cd_x[k] + g_x_0_xxxxx_xx[k];

                g_x_0_xxxxxx_y[k] = -g_xxxxx_y[k] - g_x_0_xxxxx_y[k] * cd_x[k] + g_x_0_xxxxx_xy[k];

                g_x_0_xxxxxx_z[k] = -g_xxxxx_z[k] - g_x_0_xxxxx_z[k] * cd_x[k] + g_x_0_xxxxx_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxy_x = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxxxy_y = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxxxy_z = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_y, g_x_0_xxxxx_x, g_x_0_xxxxx_xy, g_x_0_xxxxx_y, g_x_0_xxxxx_yy, g_x_0_xxxxx_yz, g_x_0_xxxxx_z, g_x_0_xxxxxy_x, g_x_0_xxxxxy_y, g_x_0_xxxxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxy_x[k] = -g_x_0_xxxxx_x[k] * cd_y[k] + g_x_0_xxxxx_xy[k];

                g_x_0_xxxxxy_y[k] = -g_x_0_xxxxx_y[k] * cd_y[k] + g_x_0_xxxxx_yy[k];

                g_x_0_xxxxxy_z[k] = -g_x_0_xxxxx_z[k] * cd_y[k] + g_x_0_xxxxx_yz[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxz_x = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxxxz_y = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxxxz_z = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 8);

            #pragma omp simd aligned(cd_z, g_x_0_xxxxx_x, g_x_0_xxxxx_xz, g_x_0_xxxxx_y, g_x_0_xxxxx_yz, g_x_0_xxxxx_z, g_x_0_xxxxx_zz, g_x_0_xxxxxz_x, g_x_0_xxxxxz_y, g_x_0_xxxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxz_x[k] = -g_x_0_xxxxx_x[k] * cd_z[k] + g_x_0_xxxxx_xz[k];

                g_x_0_xxxxxz_y[k] = -g_x_0_xxxxx_y[k] * cd_z[k] + g_x_0_xxxxx_yz[k];

                g_x_0_xxxxxz_z[k] = -g_x_0_xxxxx_z[k] * cd_z[k] + g_x_0_xxxxx_zz[k];
            }

            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxyy_x = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxxxyy_y = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxxxyy_z = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 11);

            #pragma omp simd aligned(cd_y, g_x_0_xxxxy_x, g_x_0_xxxxy_xy, g_x_0_xxxxy_y, g_x_0_xxxxy_yy, g_x_0_xxxxy_yz, g_x_0_xxxxy_z, g_x_0_xxxxyy_x, g_x_0_xxxxyy_y, g_x_0_xxxxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxyy_x[k] = -g_x_0_xxxxy_x[k] * cd_y[k] + g_x_0_xxxxy_xy[k];

                g_x_0_xxxxyy_y[k] = -g_x_0_xxxxy_y[k] * cd_y[k] + g_x_0_xxxxy_yy[k];

                g_x_0_xxxxyy_z[k] = -g_x_0_xxxxy_z[k] * cd_y[k] + g_x_0_xxxxy_yz[k];
            }

            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxyz_x = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxxxyz_y = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxxxyz_z = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 14);

            #pragma omp simd aligned(cd_y, g_x_0_xxxxyz_x, g_x_0_xxxxyz_y, g_x_0_xxxxyz_z, g_x_0_xxxxz_x, g_x_0_xxxxz_xy, g_x_0_xxxxz_y, g_x_0_xxxxz_yy, g_x_0_xxxxz_yz, g_x_0_xxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxyz_x[k] = -g_x_0_xxxxz_x[k] * cd_y[k] + g_x_0_xxxxz_xy[k];

                g_x_0_xxxxyz_y[k] = -g_x_0_xxxxz_y[k] * cd_y[k] + g_x_0_xxxxz_yy[k];

                g_x_0_xxxxyz_z[k] = -g_x_0_xxxxz_z[k] * cd_y[k] + g_x_0_xxxxz_yz[k];
            }

            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxzz_x = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxxxzz_y = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxxzz_z = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 17);

            #pragma omp simd aligned(cd_z, g_x_0_xxxxz_x, g_x_0_xxxxz_xz, g_x_0_xxxxz_y, g_x_0_xxxxz_yz, g_x_0_xxxxz_z, g_x_0_xxxxz_zz, g_x_0_xxxxzz_x, g_x_0_xxxxzz_y, g_x_0_xxxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxzz_x[k] = -g_x_0_xxxxz_x[k] * cd_z[k] + g_x_0_xxxxz_xz[k];

                g_x_0_xxxxzz_y[k] = -g_x_0_xxxxz_y[k] * cd_z[k] + g_x_0_xxxxz_yz[k];

                g_x_0_xxxxzz_z[k] = -g_x_0_xxxxz_z[k] * cd_z[k] + g_x_0_xxxxz_zz[k];
            }

            /// Set up 18-21 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyyy_x = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxxyyy_y = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxxyyy_z = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 20);

            #pragma omp simd aligned(cd_y, g_x_0_xxxyy_x, g_x_0_xxxyy_xy, g_x_0_xxxyy_y, g_x_0_xxxyy_yy, g_x_0_xxxyy_yz, g_x_0_xxxyy_z, g_x_0_xxxyyy_x, g_x_0_xxxyyy_y, g_x_0_xxxyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyyy_x[k] = -g_x_0_xxxyy_x[k] * cd_y[k] + g_x_0_xxxyy_xy[k];

                g_x_0_xxxyyy_y[k] = -g_x_0_xxxyy_y[k] * cd_y[k] + g_x_0_xxxyy_yy[k];

                g_x_0_xxxyyy_z[k] = -g_x_0_xxxyy_z[k] * cd_y[k] + g_x_0_xxxyy_yz[k];
            }

            /// Set up 21-24 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyyz_x = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxxyyz_y = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxxyyz_z = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 23);

            #pragma omp simd aligned(cd_y, g_x_0_xxxyyz_x, g_x_0_xxxyyz_y, g_x_0_xxxyyz_z, g_x_0_xxxyz_x, g_x_0_xxxyz_xy, g_x_0_xxxyz_y, g_x_0_xxxyz_yy, g_x_0_xxxyz_yz, g_x_0_xxxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyyz_x[k] = -g_x_0_xxxyz_x[k] * cd_y[k] + g_x_0_xxxyz_xy[k];

                g_x_0_xxxyyz_y[k] = -g_x_0_xxxyz_y[k] * cd_y[k] + g_x_0_xxxyz_yy[k];

                g_x_0_xxxyyz_z[k] = -g_x_0_xxxyz_z[k] * cd_y[k] + g_x_0_xxxyz_yz[k];
            }

            /// Set up 24-27 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyzz_x = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxxyzz_y = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxxyzz_z = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 26);

            #pragma omp simd aligned(cd_y, g_x_0_xxxyzz_x, g_x_0_xxxyzz_y, g_x_0_xxxyzz_z, g_x_0_xxxzz_x, g_x_0_xxxzz_xy, g_x_0_xxxzz_y, g_x_0_xxxzz_yy, g_x_0_xxxzz_yz, g_x_0_xxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyzz_x[k] = -g_x_0_xxxzz_x[k] * cd_y[k] + g_x_0_xxxzz_xy[k];

                g_x_0_xxxyzz_y[k] = -g_x_0_xxxzz_y[k] * cd_y[k] + g_x_0_xxxzz_yy[k];

                g_x_0_xxxyzz_z[k] = -g_x_0_xxxzz_z[k] * cd_y[k] + g_x_0_xxxzz_yz[k];
            }

            /// Set up 27-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxzzz_x = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxxzzz_y = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxxzzz_z = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_z, g_x_0_xxxzz_x, g_x_0_xxxzz_xz, g_x_0_xxxzz_y, g_x_0_xxxzz_yz, g_x_0_xxxzz_z, g_x_0_xxxzz_zz, g_x_0_xxxzzz_x, g_x_0_xxxzzz_y, g_x_0_xxxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxzzz_x[k] = -g_x_0_xxxzz_x[k] * cd_z[k] + g_x_0_xxxzz_xz[k];

                g_x_0_xxxzzz_y[k] = -g_x_0_xxxzz_y[k] * cd_z[k] + g_x_0_xxxzz_yz[k];

                g_x_0_xxxzzz_z[k] = -g_x_0_xxxzz_z[k] * cd_z[k] + g_x_0_xxxzz_zz[k];
            }

            /// Set up 30-33 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyyy_x = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xxyyyy_y = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxyyyy_z = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 32);

            #pragma omp simd aligned(cd_y, g_x_0_xxyyy_x, g_x_0_xxyyy_xy, g_x_0_xxyyy_y, g_x_0_xxyyy_yy, g_x_0_xxyyy_yz, g_x_0_xxyyy_z, g_x_0_xxyyyy_x, g_x_0_xxyyyy_y, g_x_0_xxyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyyy_x[k] = -g_x_0_xxyyy_x[k] * cd_y[k] + g_x_0_xxyyy_xy[k];

                g_x_0_xxyyyy_y[k] = -g_x_0_xxyyy_y[k] * cd_y[k] + g_x_0_xxyyy_yy[k];

                g_x_0_xxyyyy_z[k] = -g_x_0_xxyyy_z[k] * cd_y[k] + g_x_0_xxyyy_yz[k];
            }

            /// Set up 33-36 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyyz_x = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxyyyz_y = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxyyyz_z = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 35);

            #pragma omp simd aligned(cd_y, g_x_0_xxyyyz_x, g_x_0_xxyyyz_y, g_x_0_xxyyyz_z, g_x_0_xxyyz_x, g_x_0_xxyyz_xy, g_x_0_xxyyz_y, g_x_0_xxyyz_yy, g_x_0_xxyyz_yz, g_x_0_xxyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyyz_x[k] = -g_x_0_xxyyz_x[k] * cd_y[k] + g_x_0_xxyyz_xy[k];

                g_x_0_xxyyyz_y[k] = -g_x_0_xxyyz_y[k] * cd_y[k] + g_x_0_xxyyz_yy[k];

                g_x_0_xxyyyz_z[k] = -g_x_0_xxyyz_z[k] * cd_y[k] + g_x_0_xxyyz_yz[k];
            }

            /// Set up 36-39 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyzz_x = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xxyyzz_y = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xxyyzz_z = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 38);

            #pragma omp simd aligned(cd_y, g_x_0_xxyyzz_x, g_x_0_xxyyzz_y, g_x_0_xxyyzz_z, g_x_0_xxyzz_x, g_x_0_xxyzz_xy, g_x_0_xxyzz_y, g_x_0_xxyzz_yy, g_x_0_xxyzz_yz, g_x_0_xxyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyzz_x[k] = -g_x_0_xxyzz_x[k] * cd_y[k] + g_x_0_xxyzz_xy[k];

                g_x_0_xxyyzz_y[k] = -g_x_0_xxyzz_y[k] * cd_y[k] + g_x_0_xxyzz_yy[k];

                g_x_0_xxyyzz_z[k] = -g_x_0_xxyzz_z[k] * cd_y[k] + g_x_0_xxyzz_yz[k];
            }

            /// Set up 39-42 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyzzz_x = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xxyzzz_y = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxyzzz_z = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 41);

            #pragma omp simd aligned(cd_y, g_x_0_xxyzzz_x, g_x_0_xxyzzz_y, g_x_0_xxyzzz_z, g_x_0_xxzzz_x, g_x_0_xxzzz_xy, g_x_0_xxzzz_y, g_x_0_xxzzz_yy, g_x_0_xxzzz_yz, g_x_0_xxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyzzz_x[k] = -g_x_0_xxzzz_x[k] * cd_y[k] + g_x_0_xxzzz_xy[k];

                g_x_0_xxyzzz_y[k] = -g_x_0_xxzzz_y[k] * cd_y[k] + g_x_0_xxzzz_yy[k];

                g_x_0_xxyzzz_z[k] = -g_x_0_xxzzz_z[k] * cd_y[k] + g_x_0_xxzzz_yz[k];
            }

            /// Set up 42-45 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxzzzz_x = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xxzzzz_y = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxzzzz_z = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 44);

            #pragma omp simd aligned(cd_z, g_x_0_xxzzz_x, g_x_0_xxzzz_xz, g_x_0_xxzzz_y, g_x_0_xxzzz_yz, g_x_0_xxzzz_z, g_x_0_xxzzz_zz, g_x_0_xxzzzz_x, g_x_0_xxzzzz_y, g_x_0_xxzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzzzz_x[k] = -g_x_0_xxzzz_x[k] * cd_z[k] + g_x_0_xxzzz_xz[k];

                g_x_0_xxzzzz_y[k] = -g_x_0_xxzzz_y[k] * cd_z[k] + g_x_0_xxzzz_yz[k];

                g_x_0_xxzzzz_z[k] = -g_x_0_xxzzz_z[k] * cd_z[k] + g_x_0_xxzzz_zz[k];
            }

            /// Set up 45-48 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyyy_x = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xyyyyy_y = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xyyyyy_z = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 47);

            #pragma omp simd aligned(cd_y, g_x_0_xyyyy_x, g_x_0_xyyyy_xy, g_x_0_xyyyy_y, g_x_0_xyyyy_yy, g_x_0_xyyyy_yz, g_x_0_xyyyy_z, g_x_0_xyyyyy_x, g_x_0_xyyyyy_y, g_x_0_xyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyyy_x[k] = -g_x_0_xyyyy_x[k] * cd_y[k] + g_x_0_xyyyy_xy[k];

                g_x_0_xyyyyy_y[k] = -g_x_0_xyyyy_y[k] * cd_y[k] + g_x_0_xyyyy_yy[k];

                g_x_0_xyyyyy_z[k] = -g_x_0_xyyyy_z[k] * cd_y[k] + g_x_0_xyyyy_yz[k];
            }

            /// Set up 48-51 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyyz_x = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xyyyyz_y = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xyyyyz_z = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 50);

            #pragma omp simd aligned(cd_y, g_x_0_xyyyyz_x, g_x_0_xyyyyz_y, g_x_0_xyyyyz_z, g_x_0_xyyyz_x, g_x_0_xyyyz_xy, g_x_0_xyyyz_y, g_x_0_xyyyz_yy, g_x_0_xyyyz_yz, g_x_0_xyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyyz_x[k] = -g_x_0_xyyyz_x[k] * cd_y[k] + g_x_0_xyyyz_xy[k];

                g_x_0_xyyyyz_y[k] = -g_x_0_xyyyz_y[k] * cd_y[k] + g_x_0_xyyyz_yy[k];

                g_x_0_xyyyyz_z[k] = -g_x_0_xyyyz_z[k] * cd_y[k] + g_x_0_xyyyz_yz[k];
            }

            /// Set up 51-54 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyzz_x = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xyyyzz_y = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xyyyzz_z = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 53);

            #pragma omp simd aligned(cd_y, g_x_0_xyyyzz_x, g_x_0_xyyyzz_y, g_x_0_xyyyzz_z, g_x_0_xyyzz_x, g_x_0_xyyzz_xy, g_x_0_xyyzz_y, g_x_0_xyyzz_yy, g_x_0_xyyzz_yz, g_x_0_xyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyzz_x[k] = -g_x_0_xyyzz_x[k] * cd_y[k] + g_x_0_xyyzz_xy[k];

                g_x_0_xyyyzz_y[k] = -g_x_0_xyyzz_y[k] * cd_y[k] + g_x_0_xyyzz_yy[k];

                g_x_0_xyyyzz_z[k] = -g_x_0_xyyzz_z[k] * cd_y[k] + g_x_0_xyyzz_yz[k];
            }

            /// Set up 54-57 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyzzz_x = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xyyzzz_y = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xyyzzz_z = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 56);

            #pragma omp simd aligned(cd_y, g_x_0_xyyzzz_x, g_x_0_xyyzzz_y, g_x_0_xyyzzz_z, g_x_0_xyzzz_x, g_x_0_xyzzz_xy, g_x_0_xyzzz_y, g_x_0_xyzzz_yy, g_x_0_xyzzz_yz, g_x_0_xyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyzzz_x[k] = -g_x_0_xyzzz_x[k] * cd_y[k] + g_x_0_xyzzz_xy[k];

                g_x_0_xyyzzz_y[k] = -g_x_0_xyzzz_y[k] * cd_y[k] + g_x_0_xyzzz_yy[k];

                g_x_0_xyyzzz_z[k] = -g_x_0_xyzzz_z[k] * cd_y[k] + g_x_0_xyzzz_yz[k];
            }

            /// Set up 57-60 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyzzzz_x = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xyzzzz_y = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xyzzzz_z = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 59);

            #pragma omp simd aligned(cd_y, g_x_0_xyzzzz_x, g_x_0_xyzzzz_y, g_x_0_xyzzzz_z, g_x_0_xzzzz_x, g_x_0_xzzzz_xy, g_x_0_xzzzz_y, g_x_0_xzzzz_yy, g_x_0_xzzzz_yz, g_x_0_xzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzzzz_x[k] = -g_x_0_xzzzz_x[k] * cd_y[k] + g_x_0_xzzzz_xy[k];

                g_x_0_xyzzzz_y[k] = -g_x_0_xzzzz_y[k] * cd_y[k] + g_x_0_xzzzz_yy[k];

                g_x_0_xyzzzz_z[k] = -g_x_0_xzzzz_z[k] * cd_y[k] + g_x_0_xzzzz_yz[k];
            }

            /// Set up 60-63 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzzzzz_x = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xzzzzz_y = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xzzzzz_z = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 62);

            #pragma omp simd aligned(cd_z, g_x_0_xzzzz_x, g_x_0_xzzzz_xz, g_x_0_xzzzz_y, g_x_0_xzzzz_yz, g_x_0_xzzzz_z, g_x_0_xzzzz_zz, g_x_0_xzzzzz_x, g_x_0_xzzzzz_y, g_x_0_xzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzzzz_x[k] = -g_x_0_xzzzz_x[k] * cd_z[k] + g_x_0_xzzzz_xz[k];

                g_x_0_xzzzzz_y[k] = -g_x_0_xzzzz_y[k] * cd_z[k] + g_x_0_xzzzz_yz[k];

                g_x_0_xzzzzz_z[k] = -g_x_0_xzzzz_z[k] * cd_z[k] + g_x_0_xzzzz_zz[k];
            }

            /// Set up 63-66 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyyy_x = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_yyyyyy_y = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_yyyyyy_z = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 65);

            #pragma omp simd aligned(cd_y, g_x_0_yyyyy_x, g_x_0_yyyyy_xy, g_x_0_yyyyy_y, g_x_0_yyyyy_yy, g_x_0_yyyyy_yz, g_x_0_yyyyy_z, g_x_0_yyyyyy_x, g_x_0_yyyyyy_y, g_x_0_yyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyyy_x[k] = -g_x_0_yyyyy_x[k] * cd_y[k] + g_x_0_yyyyy_xy[k];

                g_x_0_yyyyyy_y[k] = -g_x_0_yyyyy_y[k] * cd_y[k] + g_x_0_yyyyy_yy[k];

                g_x_0_yyyyyy_z[k] = -g_x_0_yyyyy_z[k] * cd_y[k] + g_x_0_yyyyy_yz[k];
            }

            /// Set up 66-69 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyyz_x = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_yyyyyz_y = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_yyyyyz_z = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 68);

            #pragma omp simd aligned(cd_y, g_x_0_yyyyyz_x, g_x_0_yyyyyz_y, g_x_0_yyyyyz_z, g_x_0_yyyyz_x, g_x_0_yyyyz_xy, g_x_0_yyyyz_y, g_x_0_yyyyz_yy, g_x_0_yyyyz_yz, g_x_0_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyyz_x[k] = -g_x_0_yyyyz_x[k] * cd_y[k] + g_x_0_yyyyz_xy[k];

                g_x_0_yyyyyz_y[k] = -g_x_0_yyyyz_y[k] * cd_y[k] + g_x_0_yyyyz_yy[k];

                g_x_0_yyyyyz_z[k] = -g_x_0_yyyyz_z[k] * cd_y[k] + g_x_0_yyyyz_yz[k];
            }

            /// Set up 69-72 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyzz_x = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_yyyyzz_y = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_yyyyzz_z = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 71);

            #pragma omp simd aligned(cd_y, g_x_0_yyyyzz_x, g_x_0_yyyyzz_y, g_x_0_yyyyzz_z, g_x_0_yyyzz_x, g_x_0_yyyzz_xy, g_x_0_yyyzz_y, g_x_0_yyyzz_yy, g_x_0_yyyzz_yz, g_x_0_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyzz_x[k] = -g_x_0_yyyzz_x[k] * cd_y[k] + g_x_0_yyyzz_xy[k];

                g_x_0_yyyyzz_y[k] = -g_x_0_yyyzz_y[k] * cd_y[k] + g_x_0_yyyzz_yy[k];

                g_x_0_yyyyzz_z[k] = -g_x_0_yyyzz_z[k] * cd_y[k] + g_x_0_yyyzz_yz[k];
            }

            /// Set up 72-75 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyzzz_x = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_yyyzzz_y = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_yyyzzz_z = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 74);

            #pragma omp simd aligned(cd_y, g_x_0_yyyzzz_x, g_x_0_yyyzzz_y, g_x_0_yyyzzz_z, g_x_0_yyzzz_x, g_x_0_yyzzz_xy, g_x_0_yyzzz_y, g_x_0_yyzzz_yy, g_x_0_yyzzz_yz, g_x_0_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyzzz_x[k] = -g_x_0_yyzzz_x[k] * cd_y[k] + g_x_0_yyzzz_xy[k];

                g_x_0_yyyzzz_y[k] = -g_x_0_yyzzz_y[k] * cd_y[k] + g_x_0_yyzzz_yy[k];

                g_x_0_yyyzzz_z[k] = -g_x_0_yyzzz_z[k] * cd_y[k] + g_x_0_yyzzz_yz[k];
            }

            /// Set up 75-78 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyzzzz_x = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_yyzzzz_y = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_yyzzzz_z = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 77);

            #pragma omp simd aligned(cd_y, g_x_0_yyzzzz_x, g_x_0_yyzzzz_y, g_x_0_yyzzzz_z, g_x_0_yzzzz_x, g_x_0_yzzzz_xy, g_x_0_yzzzz_y, g_x_0_yzzzz_yy, g_x_0_yzzzz_yz, g_x_0_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzzzz_x[k] = -g_x_0_yzzzz_x[k] * cd_y[k] + g_x_0_yzzzz_xy[k];

                g_x_0_yyzzzz_y[k] = -g_x_0_yzzzz_y[k] * cd_y[k] + g_x_0_yzzzz_yy[k];

                g_x_0_yyzzzz_z[k] = -g_x_0_yzzzz_z[k] * cd_y[k] + g_x_0_yzzzz_yz[k];
            }

            /// Set up 78-81 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzzzzz_x = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_yzzzzz_y = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_yzzzzz_z = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 80);

            #pragma omp simd aligned(cd_y, g_x_0_yzzzzz_x, g_x_0_yzzzzz_y, g_x_0_yzzzzz_z, g_x_0_zzzzz_x, g_x_0_zzzzz_xy, g_x_0_zzzzz_y, g_x_0_zzzzz_yy, g_x_0_zzzzz_yz, g_x_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzzzz_x[k] = -g_x_0_zzzzz_x[k] * cd_y[k] + g_x_0_zzzzz_xy[k];

                g_x_0_yzzzzz_y[k] = -g_x_0_zzzzz_y[k] * cd_y[k] + g_x_0_zzzzz_yy[k];

                g_x_0_yzzzzz_z[k] = -g_x_0_zzzzz_z[k] * cd_y[k] + g_x_0_zzzzz_yz[k];
            }

            /// Set up 81-84 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzzzzz_x = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_zzzzzz_y = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_zzzzzz_z = cbuffer.data(ip_geom_10_off + 0 * acomps * bcomps + 83);

            #pragma omp simd aligned(cd_z, g_x_0_zzzzz_x, g_x_0_zzzzz_xz, g_x_0_zzzzz_y, g_x_0_zzzzz_yz, g_x_0_zzzzz_z, g_x_0_zzzzz_zz, g_x_0_zzzzzz_x, g_x_0_zzzzzz_y, g_x_0_zzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzzzz_x[k] = -g_x_0_zzzzz_x[k] * cd_z[k] + g_x_0_zzzzz_xz[k];

                g_x_0_zzzzzz_y[k] = -g_x_0_zzzzz_y[k] * cd_z[k] + g_x_0_zzzzz_yz[k];

                g_x_0_zzzzzz_z[k] = -g_x_0_zzzzz_z[k] * cd_z[k] + g_x_0_zzzzz_zz[k];
            }
            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxx_x = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 0);

            auto g_y_0_xxxxxx_y = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 1);

            auto g_y_0_xxxxxx_z = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 2);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxx_x, g_y_0_xxxxx_xx, g_y_0_xxxxx_xy, g_y_0_xxxxx_xz, g_y_0_xxxxx_y, g_y_0_xxxxx_z, g_y_0_xxxxxx_x, g_y_0_xxxxxx_y, g_y_0_xxxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxx_x[k] = -g_y_0_xxxxx_x[k] * cd_x[k] + g_y_0_xxxxx_xx[k];

                g_y_0_xxxxxx_y[k] = -g_y_0_xxxxx_y[k] * cd_x[k] + g_y_0_xxxxx_xy[k];

                g_y_0_xxxxxx_z[k] = -g_y_0_xxxxx_z[k] * cd_x[k] + g_y_0_xxxxx_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxy_x = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 3);

            auto g_y_0_xxxxxy_y = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 4);

            auto g_y_0_xxxxxy_z = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxxy_x, g_y_0_xxxxxy_y, g_y_0_xxxxxy_z, g_y_0_xxxxy_x, g_y_0_xxxxy_xx, g_y_0_xxxxy_xy, g_y_0_xxxxy_xz, g_y_0_xxxxy_y, g_y_0_xxxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxy_x[k] = -g_y_0_xxxxy_x[k] * cd_x[k] + g_y_0_xxxxy_xx[k];

                g_y_0_xxxxxy_y[k] = -g_y_0_xxxxy_y[k] * cd_x[k] + g_y_0_xxxxy_xy[k];

                g_y_0_xxxxxy_z[k] = -g_y_0_xxxxy_z[k] * cd_x[k] + g_y_0_xxxxy_xz[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxz_x = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 6);

            auto g_y_0_xxxxxz_y = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 7);

            auto g_y_0_xxxxxz_z = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 8);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxxz_x, g_y_0_xxxxxz_y, g_y_0_xxxxxz_z, g_y_0_xxxxz_x, g_y_0_xxxxz_xx, g_y_0_xxxxz_xy, g_y_0_xxxxz_xz, g_y_0_xxxxz_y, g_y_0_xxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxz_x[k] = -g_y_0_xxxxz_x[k] * cd_x[k] + g_y_0_xxxxz_xx[k];

                g_y_0_xxxxxz_y[k] = -g_y_0_xxxxz_y[k] * cd_x[k] + g_y_0_xxxxz_xy[k];

                g_y_0_xxxxxz_z[k] = -g_y_0_xxxxz_z[k] * cd_x[k] + g_y_0_xxxxz_xz[k];
            }

            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxyy_x = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 9);

            auto g_y_0_xxxxyy_y = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 10);

            auto g_y_0_xxxxyy_z = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 11);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxyy_x, g_y_0_xxxxyy_y, g_y_0_xxxxyy_z, g_y_0_xxxyy_x, g_y_0_xxxyy_xx, g_y_0_xxxyy_xy, g_y_0_xxxyy_xz, g_y_0_xxxyy_y, g_y_0_xxxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxyy_x[k] = -g_y_0_xxxyy_x[k] * cd_x[k] + g_y_0_xxxyy_xx[k];

                g_y_0_xxxxyy_y[k] = -g_y_0_xxxyy_y[k] * cd_x[k] + g_y_0_xxxyy_xy[k];

                g_y_0_xxxxyy_z[k] = -g_y_0_xxxyy_z[k] * cd_x[k] + g_y_0_xxxyy_xz[k];
            }

            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxyz_x = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 12);

            auto g_y_0_xxxxyz_y = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 13);

            auto g_y_0_xxxxyz_z = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 14);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxyz_x, g_y_0_xxxxyz_y, g_y_0_xxxxyz_z, g_y_0_xxxyz_x, g_y_0_xxxyz_xx, g_y_0_xxxyz_xy, g_y_0_xxxyz_xz, g_y_0_xxxyz_y, g_y_0_xxxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxyz_x[k] = -g_y_0_xxxyz_x[k] * cd_x[k] + g_y_0_xxxyz_xx[k];

                g_y_0_xxxxyz_y[k] = -g_y_0_xxxyz_y[k] * cd_x[k] + g_y_0_xxxyz_xy[k];

                g_y_0_xxxxyz_z[k] = -g_y_0_xxxyz_z[k] * cd_x[k] + g_y_0_xxxyz_xz[k];
            }

            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxzz_x = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 15);

            auto g_y_0_xxxxzz_y = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 16);

            auto g_y_0_xxxxzz_z = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 17);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxzz_x, g_y_0_xxxxzz_y, g_y_0_xxxxzz_z, g_y_0_xxxzz_x, g_y_0_xxxzz_xx, g_y_0_xxxzz_xy, g_y_0_xxxzz_xz, g_y_0_xxxzz_y, g_y_0_xxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxzz_x[k] = -g_y_0_xxxzz_x[k] * cd_x[k] + g_y_0_xxxzz_xx[k];

                g_y_0_xxxxzz_y[k] = -g_y_0_xxxzz_y[k] * cd_x[k] + g_y_0_xxxzz_xy[k];

                g_y_0_xxxxzz_z[k] = -g_y_0_xxxzz_z[k] * cd_x[k] + g_y_0_xxxzz_xz[k];
            }

            /// Set up 18-21 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyyy_x = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 18);

            auto g_y_0_xxxyyy_y = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 19);

            auto g_y_0_xxxyyy_z = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 20);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyyy_x, g_y_0_xxxyyy_y, g_y_0_xxxyyy_z, g_y_0_xxyyy_x, g_y_0_xxyyy_xx, g_y_0_xxyyy_xy, g_y_0_xxyyy_xz, g_y_0_xxyyy_y, g_y_0_xxyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyyy_x[k] = -g_y_0_xxyyy_x[k] * cd_x[k] + g_y_0_xxyyy_xx[k];

                g_y_0_xxxyyy_y[k] = -g_y_0_xxyyy_y[k] * cd_x[k] + g_y_0_xxyyy_xy[k];

                g_y_0_xxxyyy_z[k] = -g_y_0_xxyyy_z[k] * cd_x[k] + g_y_0_xxyyy_xz[k];
            }

            /// Set up 21-24 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyyz_x = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 21);

            auto g_y_0_xxxyyz_y = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 22);

            auto g_y_0_xxxyyz_z = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 23);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyyz_x, g_y_0_xxxyyz_y, g_y_0_xxxyyz_z, g_y_0_xxyyz_x, g_y_0_xxyyz_xx, g_y_0_xxyyz_xy, g_y_0_xxyyz_xz, g_y_0_xxyyz_y, g_y_0_xxyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyyz_x[k] = -g_y_0_xxyyz_x[k] * cd_x[k] + g_y_0_xxyyz_xx[k];

                g_y_0_xxxyyz_y[k] = -g_y_0_xxyyz_y[k] * cd_x[k] + g_y_0_xxyyz_xy[k];

                g_y_0_xxxyyz_z[k] = -g_y_0_xxyyz_z[k] * cd_x[k] + g_y_0_xxyyz_xz[k];
            }

            /// Set up 24-27 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyzz_x = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 24);

            auto g_y_0_xxxyzz_y = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 25);

            auto g_y_0_xxxyzz_z = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 26);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyzz_x, g_y_0_xxxyzz_y, g_y_0_xxxyzz_z, g_y_0_xxyzz_x, g_y_0_xxyzz_xx, g_y_0_xxyzz_xy, g_y_0_xxyzz_xz, g_y_0_xxyzz_y, g_y_0_xxyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyzz_x[k] = -g_y_0_xxyzz_x[k] * cd_x[k] + g_y_0_xxyzz_xx[k];

                g_y_0_xxxyzz_y[k] = -g_y_0_xxyzz_y[k] * cd_x[k] + g_y_0_xxyzz_xy[k];

                g_y_0_xxxyzz_z[k] = -g_y_0_xxyzz_z[k] * cd_x[k] + g_y_0_xxyzz_xz[k];
            }

            /// Set up 27-30 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxzzz_x = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 27);

            auto g_y_0_xxxzzz_y = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 28);

            auto g_y_0_xxxzzz_z = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_x, g_y_0_xxxzzz_x, g_y_0_xxxzzz_y, g_y_0_xxxzzz_z, g_y_0_xxzzz_x, g_y_0_xxzzz_xx, g_y_0_xxzzz_xy, g_y_0_xxzzz_xz, g_y_0_xxzzz_y, g_y_0_xxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxzzz_x[k] = -g_y_0_xxzzz_x[k] * cd_x[k] + g_y_0_xxzzz_xx[k];

                g_y_0_xxxzzz_y[k] = -g_y_0_xxzzz_y[k] * cd_x[k] + g_y_0_xxzzz_xy[k];

                g_y_0_xxxzzz_z[k] = -g_y_0_xxzzz_z[k] * cd_x[k] + g_y_0_xxzzz_xz[k];
            }

            /// Set up 30-33 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyyy_x = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 30);

            auto g_y_0_xxyyyy_y = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 31);

            auto g_y_0_xxyyyy_z = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 32);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyyy_x, g_y_0_xxyyyy_y, g_y_0_xxyyyy_z, g_y_0_xyyyy_x, g_y_0_xyyyy_xx, g_y_0_xyyyy_xy, g_y_0_xyyyy_xz, g_y_0_xyyyy_y, g_y_0_xyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyyy_x[k] = -g_y_0_xyyyy_x[k] * cd_x[k] + g_y_0_xyyyy_xx[k];

                g_y_0_xxyyyy_y[k] = -g_y_0_xyyyy_y[k] * cd_x[k] + g_y_0_xyyyy_xy[k];

                g_y_0_xxyyyy_z[k] = -g_y_0_xyyyy_z[k] * cd_x[k] + g_y_0_xyyyy_xz[k];
            }

            /// Set up 33-36 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyyz_x = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 33);

            auto g_y_0_xxyyyz_y = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 34);

            auto g_y_0_xxyyyz_z = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 35);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyyz_x, g_y_0_xxyyyz_y, g_y_0_xxyyyz_z, g_y_0_xyyyz_x, g_y_0_xyyyz_xx, g_y_0_xyyyz_xy, g_y_0_xyyyz_xz, g_y_0_xyyyz_y, g_y_0_xyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyyz_x[k] = -g_y_0_xyyyz_x[k] * cd_x[k] + g_y_0_xyyyz_xx[k];

                g_y_0_xxyyyz_y[k] = -g_y_0_xyyyz_y[k] * cd_x[k] + g_y_0_xyyyz_xy[k];

                g_y_0_xxyyyz_z[k] = -g_y_0_xyyyz_z[k] * cd_x[k] + g_y_0_xyyyz_xz[k];
            }

            /// Set up 36-39 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyzz_x = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 36);

            auto g_y_0_xxyyzz_y = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 37);

            auto g_y_0_xxyyzz_z = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 38);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyzz_x, g_y_0_xxyyzz_y, g_y_0_xxyyzz_z, g_y_0_xyyzz_x, g_y_0_xyyzz_xx, g_y_0_xyyzz_xy, g_y_0_xyyzz_xz, g_y_0_xyyzz_y, g_y_0_xyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyzz_x[k] = -g_y_0_xyyzz_x[k] * cd_x[k] + g_y_0_xyyzz_xx[k];

                g_y_0_xxyyzz_y[k] = -g_y_0_xyyzz_y[k] * cd_x[k] + g_y_0_xyyzz_xy[k];

                g_y_0_xxyyzz_z[k] = -g_y_0_xyyzz_z[k] * cd_x[k] + g_y_0_xyyzz_xz[k];
            }

            /// Set up 39-42 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyzzz_x = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 39);

            auto g_y_0_xxyzzz_y = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 40);

            auto g_y_0_xxyzzz_z = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 41);

            #pragma omp simd aligned(cd_x, g_y_0_xxyzzz_x, g_y_0_xxyzzz_y, g_y_0_xxyzzz_z, g_y_0_xyzzz_x, g_y_0_xyzzz_xx, g_y_0_xyzzz_xy, g_y_0_xyzzz_xz, g_y_0_xyzzz_y, g_y_0_xyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyzzz_x[k] = -g_y_0_xyzzz_x[k] * cd_x[k] + g_y_0_xyzzz_xx[k];

                g_y_0_xxyzzz_y[k] = -g_y_0_xyzzz_y[k] * cd_x[k] + g_y_0_xyzzz_xy[k];

                g_y_0_xxyzzz_z[k] = -g_y_0_xyzzz_z[k] * cd_x[k] + g_y_0_xyzzz_xz[k];
            }

            /// Set up 42-45 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxzzzz_x = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 42);

            auto g_y_0_xxzzzz_y = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 43);

            auto g_y_0_xxzzzz_z = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 44);

            #pragma omp simd aligned(cd_x, g_y_0_xxzzzz_x, g_y_0_xxzzzz_y, g_y_0_xxzzzz_z, g_y_0_xzzzz_x, g_y_0_xzzzz_xx, g_y_0_xzzzz_xy, g_y_0_xzzzz_xz, g_y_0_xzzzz_y, g_y_0_xzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzzzz_x[k] = -g_y_0_xzzzz_x[k] * cd_x[k] + g_y_0_xzzzz_xx[k];

                g_y_0_xxzzzz_y[k] = -g_y_0_xzzzz_y[k] * cd_x[k] + g_y_0_xzzzz_xy[k];

                g_y_0_xxzzzz_z[k] = -g_y_0_xzzzz_z[k] * cd_x[k] + g_y_0_xzzzz_xz[k];
            }

            /// Set up 45-48 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyyy_x = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 45);

            auto g_y_0_xyyyyy_y = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 46);

            auto g_y_0_xyyyyy_z = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 47);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyyy_x, g_y_0_xyyyyy_y, g_y_0_xyyyyy_z, g_y_0_yyyyy_x, g_y_0_yyyyy_xx, g_y_0_yyyyy_xy, g_y_0_yyyyy_xz, g_y_0_yyyyy_y, g_y_0_yyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyyy_x[k] = -g_y_0_yyyyy_x[k] * cd_x[k] + g_y_0_yyyyy_xx[k];

                g_y_0_xyyyyy_y[k] = -g_y_0_yyyyy_y[k] * cd_x[k] + g_y_0_yyyyy_xy[k];

                g_y_0_xyyyyy_z[k] = -g_y_0_yyyyy_z[k] * cd_x[k] + g_y_0_yyyyy_xz[k];
            }

            /// Set up 48-51 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyyz_x = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 48);

            auto g_y_0_xyyyyz_y = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 49);

            auto g_y_0_xyyyyz_z = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 50);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyyz_x, g_y_0_xyyyyz_y, g_y_0_xyyyyz_z, g_y_0_yyyyz_x, g_y_0_yyyyz_xx, g_y_0_yyyyz_xy, g_y_0_yyyyz_xz, g_y_0_yyyyz_y, g_y_0_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyyz_x[k] = -g_y_0_yyyyz_x[k] * cd_x[k] + g_y_0_yyyyz_xx[k];

                g_y_0_xyyyyz_y[k] = -g_y_0_yyyyz_y[k] * cd_x[k] + g_y_0_yyyyz_xy[k];

                g_y_0_xyyyyz_z[k] = -g_y_0_yyyyz_z[k] * cd_x[k] + g_y_0_yyyyz_xz[k];
            }

            /// Set up 51-54 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyzz_x = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 51);

            auto g_y_0_xyyyzz_y = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 52);

            auto g_y_0_xyyyzz_z = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 53);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyzz_x, g_y_0_xyyyzz_y, g_y_0_xyyyzz_z, g_y_0_yyyzz_x, g_y_0_yyyzz_xx, g_y_0_yyyzz_xy, g_y_0_yyyzz_xz, g_y_0_yyyzz_y, g_y_0_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyzz_x[k] = -g_y_0_yyyzz_x[k] * cd_x[k] + g_y_0_yyyzz_xx[k];

                g_y_0_xyyyzz_y[k] = -g_y_0_yyyzz_y[k] * cd_x[k] + g_y_0_yyyzz_xy[k];

                g_y_0_xyyyzz_z[k] = -g_y_0_yyyzz_z[k] * cd_x[k] + g_y_0_yyyzz_xz[k];
            }

            /// Set up 54-57 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyzzz_x = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 54);

            auto g_y_0_xyyzzz_y = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 55);

            auto g_y_0_xyyzzz_z = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 56);

            #pragma omp simd aligned(cd_x, g_y_0_xyyzzz_x, g_y_0_xyyzzz_y, g_y_0_xyyzzz_z, g_y_0_yyzzz_x, g_y_0_yyzzz_xx, g_y_0_yyzzz_xy, g_y_0_yyzzz_xz, g_y_0_yyzzz_y, g_y_0_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyzzz_x[k] = -g_y_0_yyzzz_x[k] * cd_x[k] + g_y_0_yyzzz_xx[k];

                g_y_0_xyyzzz_y[k] = -g_y_0_yyzzz_y[k] * cd_x[k] + g_y_0_yyzzz_xy[k];

                g_y_0_xyyzzz_z[k] = -g_y_0_yyzzz_z[k] * cd_x[k] + g_y_0_yyzzz_xz[k];
            }

            /// Set up 57-60 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyzzzz_x = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 57);

            auto g_y_0_xyzzzz_y = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 58);

            auto g_y_0_xyzzzz_z = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 59);

            #pragma omp simd aligned(cd_x, g_y_0_xyzzzz_x, g_y_0_xyzzzz_y, g_y_0_xyzzzz_z, g_y_0_yzzzz_x, g_y_0_yzzzz_xx, g_y_0_yzzzz_xy, g_y_0_yzzzz_xz, g_y_0_yzzzz_y, g_y_0_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzzzz_x[k] = -g_y_0_yzzzz_x[k] * cd_x[k] + g_y_0_yzzzz_xx[k];

                g_y_0_xyzzzz_y[k] = -g_y_0_yzzzz_y[k] * cd_x[k] + g_y_0_yzzzz_xy[k];

                g_y_0_xyzzzz_z[k] = -g_y_0_yzzzz_z[k] * cd_x[k] + g_y_0_yzzzz_xz[k];
            }

            /// Set up 60-63 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzzzzz_x = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 60);

            auto g_y_0_xzzzzz_y = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 61);

            auto g_y_0_xzzzzz_z = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 62);

            #pragma omp simd aligned(cd_x, g_y_0_xzzzzz_x, g_y_0_xzzzzz_y, g_y_0_xzzzzz_z, g_y_0_zzzzz_x, g_y_0_zzzzz_xx, g_y_0_zzzzz_xy, g_y_0_zzzzz_xz, g_y_0_zzzzz_y, g_y_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzzzz_x[k] = -g_y_0_zzzzz_x[k] * cd_x[k] + g_y_0_zzzzz_xx[k];

                g_y_0_xzzzzz_y[k] = -g_y_0_zzzzz_y[k] * cd_x[k] + g_y_0_zzzzz_xy[k];

                g_y_0_xzzzzz_z[k] = -g_y_0_zzzzz_z[k] * cd_x[k] + g_y_0_zzzzz_xz[k];
            }

            /// Set up 63-66 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyyy_x = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 63);

            auto g_y_0_yyyyyy_y = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 64);

            auto g_y_0_yyyyyy_z = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 65);

            #pragma omp simd aligned(cd_y, g_y_0_yyyyy_x, g_y_0_yyyyy_xy, g_y_0_yyyyy_y, g_y_0_yyyyy_yy, g_y_0_yyyyy_yz, g_y_0_yyyyy_z, g_y_0_yyyyyy_x, g_y_0_yyyyyy_y, g_y_0_yyyyyy_z, g_yyyyy_x, g_yyyyy_y, g_yyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyyy_x[k] = -g_yyyyy_x[k] - g_y_0_yyyyy_x[k] * cd_y[k] + g_y_0_yyyyy_xy[k];

                g_y_0_yyyyyy_y[k] = -g_yyyyy_y[k] - g_y_0_yyyyy_y[k] * cd_y[k] + g_y_0_yyyyy_yy[k];

                g_y_0_yyyyyy_z[k] = -g_yyyyy_z[k] - g_y_0_yyyyy_z[k] * cd_y[k] + g_y_0_yyyyy_yz[k];
            }

            /// Set up 66-69 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyyz_x = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 66);

            auto g_y_0_yyyyyz_y = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 67);

            auto g_y_0_yyyyyz_z = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 68);

            #pragma omp simd aligned(cd_z, g_y_0_yyyyy_x, g_y_0_yyyyy_xz, g_y_0_yyyyy_y, g_y_0_yyyyy_yz, g_y_0_yyyyy_z, g_y_0_yyyyy_zz, g_y_0_yyyyyz_x, g_y_0_yyyyyz_y, g_y_0_yyyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyyz_x[k] = -g_y_0_yyyyy_x[k] * cd_z[k] + g_y_0_yyyyy_xz[k];

                g_y_0_yyyyyz_y[k] = -g_y_0_yyyyy_y[k] * cd_z[k] + g_y_0_yyyyy_yz[k];

                g_y_0_yyyyyz_z[k] = -g_y_0_yyyyy_z[k] * cd_z[k] + g_y_0_yyyyy_zz[k];
            }

            /// Set up 69-72 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyzz_x = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 69);

            auto g_y_0_yyyyzz_y = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 70);

            auto g_y_0_yyyyzz_z = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 71);

            #pragma omp simd aligned(cd_z, g_y_0_yyyyz_x, g_y_0_yyyyz_xz, g_y_0_yyyyz_y, g_y_0_yyyyz_yz, g_y_0_yyyyz_z, g_y_0_yyyyz_zz, g_y_0_yyyyzz_x, g_y_0_yyyyzz_y, g_y_0_yyyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyzz_x[k] = -g_y_0_yyyyz_x[k] * cd_z[k] + g_y_0_yyyyz_xz[k];

                g_y_0_yyyyzz_y[k] = -g_y_0_yyyyz_y[k] * cd_z[k] + g_y_0_yyyyz_yz[k];

                g_y_0_yyyyzz_z[k] = -g_y_0_yyyyz_z[k] * cd_z[k] + g_y_0_yyyyz_zz[k];
            }

            /// Set up 72-75 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyzzz_x = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 72);

            auto g_y_0_yyyzzz_y = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 73);

            auto g_y_0_yyyzzz_z = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 74);

            #pragma omp simd aligned(cd_z, g_y_0_yyyzz_x, g_y_0_yyyzz_xz, g_y_0_yyyzz_y, g_y_0_yyyzz_yz, g_y_0_yyyzz_z, g_y_0_yyyzz_zz, g_y_0_yyyzzz_x, g_y_0_yyyzzz_y, g_y_0_yyyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyzzz_x[k] = -g_y_0_yyyzz_x[k] * cd_z[k] + g_y_0_yyyzz_xz[k];

                g_y_0_yyyzzz_y[k] = -g_y_0_yyyzz_y[k] * cd_z[k] + g_y_0_yyyzz_yz[k];

                g_y_0_yyyzzz_z[k] = -g_y_0_yyyzz_z[k] * cd_z[k] + g_y_0_yyyzz_zz[k];
            }

            /// Set up 75-78 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyzzzz_x = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 75);

            auto g_y_0_yyzzzz_y = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 76);

            auto g_y_0_yyzzzz_z = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 77);

            #pragma omp simd aligned(cd_z, g_y_0_yyzzz_x, g_y_0_yyzzz_xz, g_y_0_yyzzz_y, g_y_0_yyzzz_yz, g_y_0_yyzzz_z, g_y_0_yyzzz_zz, g_y_0_yyzzzz_x, g_y_0_yyzzzz_y, g_y_0_yyzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzzzz_x[k] = -g_y_0_yyzzz_x[k] * cd_z[k] + g_y_0_yyzzz_xz[k];

                g_y_0_yyzzzz_y[k] = -g_y_0_yyzzz_y[k] * cd_z[k] + g_y_0_yyzzz_yz[k];

                g_y_0_yyzzzz_z[k] = -g_y_0_yyzzz_z[k] * cd_z[k] + g_y_0_yyzzz_zz[k];
            }

            /// Set up 78-81 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzzzzz_x = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 78);

            auto g_y_0_yzzzzz_y = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 79);

            auto g_y_0_yzzzzz_z = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 80);

            #pragma omp simd aligned(cd_z, g_y_0_yzzzz_x, g_y_0_yzzzz_xz, g_y_0_yzzzz_y, g_y_0_yzzzz_yz, g_y_0_yzzzz_z, g_y_0_yzzzz_zz, g_y_0_yzzzzz_x, g_y_0_yzzzzz_y, g_y_0_yzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzzzz_x[k] = -g_y_0_yzzzz_x[k] * cd_z[k] + g_y_0_yzzzz_xz[k];

                g_y_0_yzzzzz_y[k] = -g_y_0_yzzzz_y[k] * cd_z[k] + g_y_0_yzzzz_yz[k];

                g_y_0_yzzzzz_z[k] = -g_y_0_yzzzz_z[k] * cd_z[k] + g_y_0_yzzzz_zz[k];
            }

            /// Set up 81-84 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzzzzz_x = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 81);

            auto g_y_0_zzzzzz_y = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 82);

            auto g_y_0_zzzzzz_z = cbuffer.data(ip_geom_10_off + 84 * acomps * bcomps + 83);

            #pragma omp simd aligned(cd_z, g_y_0_zzzzz_x, g_y_0_zzzzz_xz, g_y_0_zzzzz_y, g_y_0_zzzzz_yz, g_y_0_zzzzz_z, g_y_0_zzzzz_zz, g_y_0_zzzzzz_x, g_y_0_zzzzzz_y, g_y_0_zzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzzzz_x[k] = -g_y_0_zzzzz_x[k] * cd_z[k] + g_y_0_zzzzz_xz[k];

                g_y_0_zzzzzz_y[k] = -g_y_0_zzzzz_y[k] * cd_z[k] + g_y_0_zzzzz_yz[k];

                g_y_0_zzzzzz_z[k] = -g_y_0_zzzzz_z[k] * cd_z[k] + g_y_0_zzzzz_zz[k];
            }
            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxx_x = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 0);

            auto g_z_0_xxxxxx_y = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 1);

            auto g_z_0_xxxxxx_z = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 2);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxx_x, g_z_0_xxxxx_xx, g_z_0_xxxxx_xy, g_z_0_xxxxx_xz, g_z_0_xxxxx_y, g_z_0_xxxxx_z, g_z_0_xxxxxx_x, g_z_0_xxxxxx_y, g_z_0_xxxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxx_x[k] = -g_z_0_xxxxx_x[k] * cd_x[k] + g_z_0_xxxxx_xx[k];

                g_z_0_xxxxxx_y[k] = -g_z_0_xxxxx_y[k] * cd_x[k] + g_z_0_xxxxx_xy[k];

                g_z_0_xxxxxx_z[k] = -g_z_0_xxxxx_z[k] * cd_x[k] + g_z_0_xxxxx_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxy_x = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 3);

            auto g_z_0_xxxxxy_y = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 4);

            auto g_z_0_xxxxxy_z = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxxy_x, g_z_0_xxxxxy_y, g_z_0_xxxxxy_z, g_z_0_xxxxy_x, g_z_0_xxxxy_xx, g_z_0_xxxxy_xy, g_z_0_xxxxy_xz, g_z_0_xxxxy_y, g_z_0_xxxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxy_x[k] = -g_z_0_xxxxy_x[k] * cd_x[k] + g_z_0_xxxxy_xx[k];

                g_z_0_xxxxxy_y[k] = -g_z_0_xxxxy_y[k] * cd_x[k] + g_z_0_xxxxy_xy[k];

                g_z_0_xxxxxy_z[k] = -g_z_0_xxxxy_z[k] * cd_x[k] + g_z_0_xxxxy_xz[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxz_x = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 6);

            auto g_z_0_xxxxxz_y = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 7);

            auto g_z_0_xxxxxz_z = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 8);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxxz_x, g_z_0_xxxxxz_y, g_z_0_xxxxxz_z, g_z_0_xxxxz_x, g_z_0_xxxxz_xx, g_z_0_xxxxz_xy, g_z_0_xxxxz_xz, g_z_0_xxxxz_y, g_z_0_xxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxz_x[k] = -g_z_0_xxxxz_x[k] * cd_x[k] + g_z_0_xxxxz_xx[k];

                g_z_0_xxxxxz_y[k] = -g_z_0_xxxxz_y[k] * cd_x[k] + g_z_0_xxxxz_xy[k];

                g_z_0_xxxxxz_z[k] = -g_z_0_xxxxz_z[k] * cd_x[k] + g_z_0_xxxxz_xz[k];
            }

            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxyy_x = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 9);

            auto g_z_0_xxxxyy_y = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 10);

            auto g_z_0_xxxxyy_z = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 11);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxyy_x, g_z_0_xxxxyy_y, g_z_0_xxxxyy_z, g_z_0_xxxyy_x, g_z_0_xxxyy_xx, g_z_0_xxxyy_xy, g_z_0_xxxyy_xz, g_z_0_xxxyy_y, g_z_0_xxxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxyy_x[k] = -g_z_0_xxxyy_x[k] * cd_x[k] + g_z_0_xxxyy_xx[k];

                g_z_0_xxxxyy_y[k] = -g_z_0_xxxyy_y[k] * cd_x[k] + g_z_0_xxxyy_xy[k];

                g_z_0_xxxxyy_z[k] = -g_z_0_xxxyy_z[k] * cd_x[k] + g_z_0_xxxyy_xz[k];
            }

            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxyz_x = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 12);

            auto g_z_0_xxxxyz_y = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 13);

            auto g_z_0_xxxxyz_z = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 14);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxyz_x, g_z_0_xxxxyz_y, g_z_0_xxxxyz_z, g_z_0_xxxyz_x, g_z_0_xxxyz_xx, g_z_0_xxxyz_xy, g_z_0_xxxyz_xz, g_z_0_xxxyz_y, g_z_0_xxxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxyz_x[k] = -g_z_0_xxxyz_x[k] * cd_x[k] + g_z_0_xxxyz_xx[k];

                g_z_0_xxxxyz_y[k] = -g_z_0_xxxyz_y[k] * cd_x[k] + g_z_0_xxxyz_xy[k];

                g_z_0_xxxxyz_z[k] = -g_z_0_xxxyz_z[k] * cd_x[k] + g_z_0_xxxyz_xz[k];
            }

            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxzz_x = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 15);

            auto g_z_0_xxxxzz_y = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 16);

            auto g_z_0_xxxxzz_z = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 17);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxzz_x, g_z_0_xxxxzz_y, g_z_0_xxxxzz_z, g_z_0_xxxzz_x, g_z_0_xxxzz_xx, g_z_0_xxxzz_xy, g_z_0_xxxzz_xz, g_z_0_xxxzz_y, g_z_0_xxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxzz_x[k] = -g_z_0_xxxzz_x[k] * cd_x[k] + g_z_0_xxxzz_xx[k];

                g_z_0_xxxxzz_y[k] = -g_z_0_xxxzz_y[k] * cd_x[k] + g_z_0_xxxzz_xy[k];

                g_z_0_xxxxzz_z[k] = -g_z_0_xxxzz_z[k] * cd_x[k] + g_z_0_xxxzz_xz[k];
            }

            /// Set up 18-21 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyyy_x = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 18);

            auto g_z_0_xxxyyy_y = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 19);

            auto g_z_0_xxxyyy_z = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 20);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyyy_x, g_z_0_xxxyyy_y, g_z_0_xxxyyy_z, g_z_0_xxyyy_x, g_z_0_xxyyy_xx, g_z_0_xxyyy_xy, g_z_0_xxyyy_xz, g_z_0_xxyyy_y, g_z_0_xxyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyyy_x[k] = -g_z_0_xxyyy_x[k] * cd_x[k] + g_z_0_xxyyy_xx[k];

                g_z_0_xxxyyy_y[k] = -g_z_0_xxyyy_y[k] * cd_x[k] + g_z_0_xxyyy_xy[k];

                g_z_0_xxxyyy_z[k] = -g_z_0_xxyyy_z[k] * cd_x[k] + g_z_0_xxyyy_xz[k];
            }

            /// Set up 21-24 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyyz_x = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 21);

            auto g_z_0_xxxyyz_y = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 22);

            auto g_z_0_xxxyyz_z = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 23);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyyz_x, g_z_0_xxxyyz_y, g_z_0_xxxyyz_z, g_z_0_xxyyz_x, g_z_0_xxyyz_xx, g_z_0_xxyyz_xy, g_z_0_xxyyz_xz, g_z_0_xxyyz_y, g_z_0_xxyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyyz_x[k] = -g_z_0_xxyyz_x[k] * cd_x[k] + g_z_0_xxyyz_xx[k];

                g_z_0_xxxyyz_y[k] = -g_z_0_xxyyz_y[k] * cd_x[k] + g_z_0_xxyyz_xy[k];

                g_z_0_xxxyyz_z[k] = -g_z_0_xxyyz_z[k] * cd_x[k] + g_z_0_xxyyz_xz[k];
            }

            /// Set up 24-27 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyzz_x = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 24);

            auto g_z_0_xxxyzz_y = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 25);

            auto g_z_0_xxxyzz_z = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 26);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyzz_x, g_z_0_xxxyzz_y, g_z_0_xxxyzz_z, g_z_0_xxyzz_x, g_z_0_xxyzz_xx, g_z_0_xxyzz_xy, g_z_0_xxyzz_xz, g_z_0_xxyzz_y, g_z_0_xxyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyzz_x[k] = -g_z_0_xxyzz_x[k] * cd_x[k] + g_z_0_xxyzz_xx[k];

                g_z_0_xxxyzz_y[k] = -g_z_0_xxyzz_y[k] * cd_x[k] + g_z_0_xxyzz_xy[k];

                g_z_0_xxxyzz_z[k] = -g_z_0_xxyzz_z[k] * cd_x[k] + g_z_0_xxyzz_xz[k];
            }

            /// Set up 27-30 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxzzz_x = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 27);

            auto g_z_0_xxxzzz_y = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 28);

            auto g_z_0_xxxzzz_z = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_x, g_z_0_xxxzzz_x, g_z_0_xxxzzz_y, g_z_0_xxxzzz_z, g_z_0_xxzzz_x, g_z_0_xxzzz_xx, g_z_0_xxzzz_xy, g_z_0_xxzzz_xz, g_z_0_xxzzz_y, g_z_0_xxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxzzz_x[k] = -g_z_0_xxzzz_x[k] * cd_x[k] + g_z_0_xxzzz_xx[k];

                g_z_0_xxxzzz_y[k] = -g_z_0_xxzzz_y[k] * cd_x[k] + g_z_0_xxzzz_xy[k];

                g_z_0_xxxzzz_z[k] = -g_z_0_xxzzz_z[k] * cd_x[k] + g_z_0_xxzzz_xz[k];
            }

            /// Set up 30-33 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyyy_x = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 30);

            auto g_z_0_xxyyyy_y = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 31);

            auto g_z_0_xxyyyy_z = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 32);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyyy_x, g_z_0_xxyyyy_y, g_z_0_xxyyyy_z, g_z_0_xyyyy_x, g_z_0_xyyyy_xx, g_z_0_xyyyy_xy, g_z_0_xyyyy_xz, g_z_0_xyyyy_y, g_z_0_xyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyyy_x[k] = -g_z_0_xyyyy_x[k] * cd_x[k] + g_z_0_xyyyy_xx[k];

                g_z_0_xxyyyy_y[k] = -g_z_0_xyyyy_y[k] * cd_x[k] + g_z_0_xyyyy_xy[k];

                g_z_0_xxyyyy_z[k] = -g_z_0_xyyyy_z[k] * cd_x[k] + g_z_0_xyyyy_xz[k];
            }

            /// Set up 33-36 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyyz_x = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 33);

            auto g_z_0_xxyyyz_y = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 34);

            auto g_z_0_xxyyyz_z = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 35);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyyz_x, g_z_0_xxyyyz_y, g_z_0_xxyyyz_z, g_z_0_xyyyz_x, g_z_0_xyyyz_xx, g_z_0_xyyyz_xy, g_z_0_xyyyz_xz, g_z_0_xyyyz_y, g_z_0_xyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyyz_x[k] = -g_z_0_xyyyz_x[k] * cd_x[k] + g_z_0_xyyyz_xx[k];

                g_z_0_xxyyyz_y[k] = -g_z_0_xyyyz_y[k] * cd_x[k] + g_z_0_xyyyz_xy[k];

                g_z_0_xxyyyz_z[k] = -g_z_0_xyyyz_z[k] * cd_x[k] + g_z_0_xyyyz_xz[k];
            }

            /// Set up 36-39 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyzz_x = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 36);

            auto g_z_0_xxyyzz_y = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 37);

            auto g_z_0_xxyyzz_z = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 38);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyzz_x, g_z_0_xxyyzz_y, g_z_0_xxyyzz_z, g_z_0_xyyzz_x, g_z_0_xyyzz_xx, g_z_0_xyyzz_xy, g_z_0_xyyzz_xz, g_z_0_xyyzz_y, g_z_0_xyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyzz_x[k] = -g_z_0_xyyzz_x[k] * cd_x[k] + g_z_0_xyyzz_xx[k];

                g_z_0_xxyyzz_y[k] = -g_z_0_xyyzz_y[k] * cd_x[k] + g_z_0_xyyzz_xy[k];

                g_z_0_xxyyzz_z[k] = -g_z_0_xyyzz_z[k] * cd_x[k] + g_z_0_xyyzz_xz[k];
            }

            /// Set up 39-42 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyzzz_x = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 39);

            auto g_z_0_xxyzzz_y = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 40);

            auto g_z_0_xxyzzz_z = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 41);

            #pragma omp simd aligned(cd_x, g_z_0_xxyzzz_x, g_z_0_xxyzzz_y, g_z_0_xxyzzz_z, g_z_0_xyzzz_x, g_z_0_xyzzz_xx, g_z_0_xyzzz_xy, g_z_0_xyzzz_xz, g_z_0_xyzzz_y, g_z_0_xyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyzzz_x[k] = -g_z_0_xyzzz_x[k] * cd_x[k] + g_z_0_xyzzz_xx[k];

                g_z_0_xxyzzz_y[k] = -g_z_0_xyzzz_y[k] * cd_x[k] + g_z_0_xyzzz_xy[k];

                g_z_0_xxyzzz_z[k] = -g_z_0_xyzzz_z[k] * cd_x[k] + g_z_0_xyzzz_xz[k];
            }

            /// Set up 42-45 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxzzzz_x = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 42);

            auto g_z_0_xxzzzz_y = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 43);

            auto g_z_0_xxzzzz_z = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 44);

            #pragma omp simd aligned(cd_x, g_z_0_xxzzzz_x, g_z_0_xxzzzz_y, g_z_0_xxzzzz_z, g_z_0_xzzzz_x, g_z_0_xzzzz_xx, g_z_0_xzzzz_xy, g_z_0_xzzzz_xz, g_z_0_xzzzz_y, g_z_0_xzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzzzz_x[k] = -g_z_0_xzzzz_x[k] * cd_x[k] + g_z_0_xzzzz_xx[k];

                g_z_0_xxzzzz_y[k] = -g_z_0_xzzzz_y[k] * cd_x[k] + g_z_0_xzzzz_xy[k];

                g_z_0_xxzzzz_z[k] = -g_z_0_xzzzz_z[k] * cd_x[k] + g_z_0_xzzzz_xz[k];
            }

            /// Set up 45-48 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyyy_x = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 45);

            auto g_z_0_xyyyyy_y = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 46);

            auto g_z_0_xyyyyy_z = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 47);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyyy_x, g_z_0_xyyyyy_y, g_z_0_xyyyyy_z, g_z_0_yyyyy_x, g_z_0_yyyyy_xx, g_z_0_yyyyy_xy, g_z_0_yyyyy_xz, g_z_0_yyyyy_y, g_z_0_yyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyyy_x[k] = -g_z_0_yyyyy_x[k] * cd_x[k] + g_z_0_yyyyy_xx[k];

                g_z_0_xyyyyy_y[k] = -g_z_0_yyyyy_y[k] * cd_x[k] + g_z_0_yyyyy_xy[k];

                g_z_0_xyyyyy_z[k] = -g_z_0_yyyyy_z[k] * cd_x[k] + g_z_0_yyyyy_xz[k];
            }

            /// Set up 48-51 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyyz_x = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 48);

            auto g_z_0_xyyyyz_y = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 49);

            auto g_z_0_xyyyyz_z = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 50);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyyz_x, g_z_0_xyyyyz_y, g_z_0_xyyyyz_z, g_z_0_yyyyz_x, g_z_0_yyyyz_xx, g_z_0_yyyyz_xy, g_z_0_yyyyz_xz, g_z_0_yyyyz_y, g_z_0_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyyz_x[k] = -g_z_0_yyyyz_x[k] * cd_x[k] + g_z_0_yyyyz_xx[k];

                g_z_0_xyyyyz_y[k] = -g_z_0_yyyyz_y[k] * cd_x[k] + g_z_0_yyyyz_xy[k];

                g_z_0_xyyyyz_z[k] = -g_z_0_yyyyz_z[k] * cd_x[k] + g_z_0_yyyyz_xz[k];
            }

            /// Set up 51-54 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyzz_x = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 51);

            auto g_z_0_xyyyzz_y = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 52);

            auto g_z_0_xyyyzz_z = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 53);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyzz_x, g_z_0_xyyyzz_y, g_z_0_xyyyzz_z, g_z_0_yyyzz_x, g_z_0_yyyzz_xx, g_z_0_yyyzz_xy, g_z_0_yyyzz_xz, g_z_0_yyyzz_y, g_z_0_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyzz_x[k] = -g_z_0_yyyzz_x[k] * cd_x[k] + g_z_0_yyyzz_xx[k];

                g_z_0_xyyyzz_y[k] = -g_z_0_yyyzz_y[k] * cd_x[k] + g_z_0_yyyzz_xy[k];

                g_z_0_xyyyzz_z[k] = -g_z_0_yyyzz_z[k] * cd_x[k] + g_z_0_yyyzz_xz[k];
            }

            /// Set up 54-57 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyzzz_x = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 54);

            auto g_z_0_xyyzzz_y = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 55);

            auto g_z_0_xyyzzz_z = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 56);

            #pragma omp simd aligned(cd_x, g_z_0_xyyzzz_x, g_z_0_xyyzzz_y, g_z_0_xyyzzz_z, g_z_0_yyzzz_x, g_z_0_yyzzz_xx, g_z_0_yyzzz_xy, g_z_0_yyzzz_xz, g_z_0_yyzzz_y, g_z_0_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyzzz_x[k] = -g_z_0_yyzzz_x[k] * cd_x[k] + g_z_0_yyzzz_xx[k];

                g_z_0_xyyzzz_y[k] = -g_z_0_yyzzz_y[k] * cd_x[k] + g_z_0_yyzzz_xy[k];

                g_z_0_xyyzzz_z[k] = -g_z_0_yyzzz_z[k] * cd_x[k] + g_z_0_yyzzz_xz[k];
            }

            /// Set up 57-60 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyzzzz_x = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 57);

            auto g_z_0_xyzzzz_y = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 58);

            auto g_z_0_xyzzzz_z = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 59);

            #pragma omp simd aligned(cd_x, g_z_0_xyzzzz_x, g_z_0_xyzzzz_y, g_z_0_xyzzzz_z, g_z_0_yzzzz_x, g_z_0_yzzzz_xx, g_z_0_yzzzz_xy, g_z_0_yzzzz_xz, g_z_0_yzzzz_y, g_z_0_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzzzz_x[k] = -g_z_0_yzzzz_x[k] * cd_x[k] + g_z_0_yzzzz_xx[k];

                g_z_0_xyzzzz_y[k] = -g_z_0_yzzzz_y[k] * cd_x[k] + g_z_0_yzzzz_xy[k];

                g_z_0_xyzzzz_z[k] = -g_z_0_yzzzz_z[k] * cd_x[k] + g_z_0_yzzzz_xz[k];
            }

            /// Set up 60-63 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzzzzz_x = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 60);

            auto g_z_0_xzzzzz_y = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 61);

            auto g_z_0_xzzzzz_z = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 62);

            #pragma omp simd aligned(cd_x, g_z_0_xzzzzz_x, g_z_0_xzzzzz_y, g_z_0_xzzzzz_z, g_z_0_zzzzz_x, g_z_0_zzzzz_xx, g_z_0_zzzzz_xy, g_z_0_zzzzz_xz, g_z_0_zzzzz_y, g_z_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzzzz_x[k] = -g_z_0_zzzzz_x[k] * cd_x[k] + g_z_0_zzzzz_xx[k];

                g_z_0_xzzzzz_y[k] = -g_z_0_zzzzz_y[k] * cd_x[k] + g_z_0_zzzzz_xy[k];

                g_z_0_xzzzzz_z[k] = -g_z_0_zzzzz_z[k] * cd_x[k] + g_z_0_zzzzz_xz[k];
            }

            /// Set up 63-66 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyyy_x = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 63);

            auto g_z_0_yyyyyy_y = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 64);

            auto g_z_0_yyyyyy_z = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 65);

            #pragma omp simd aligned(cd_y, g_z_0_yyyyy_x, g_z_0_yyyyy_xy, g_z_0_yyyyy_y, g_z_0_yyyyy_yy, g_z_0_yyyyy_yz, g_z_0_yyyyy_z, g_z_0_yyyyyy_x, g_z_0_yyyyyy_y, g_z_0_yyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyyy_x[k] = -g_z_0_yyyyy_x[k] * cd_y[k] + g_z_0_yyyyy_xy[k];

                g_z_0_yyyyyy_y[k] = -g_z_0_yyyyy_y[k] * cd_y[k] + g_z_0_yyyyy_yy[k];

                g_z_0_yyyyyy_z[k] = -g_z_0_yyyyy_z[k] * cd_y[k] + g_z_0_yyyyy_yz[k];
            }

            /// Set up 66-69 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyyz_x = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 66);

            auto g_z_0_yyyyyz_y = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 67);

            auto g_z_0_yyyyyz_z = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 68);

            #pragma omp simd aligned(cd_y, g_z_0_yyyyyz_x, g_z_0_yyyyyz_y, g_z_0_yyyyyz_z, g_z_0_yyyyz_x, g_z_0_yyyyz_xy, g_z_0_yyyyz_y, g_z_0_yyyyz_yy, g_z_0_yyyyz_yz, g_z_0_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyyz_x[k] = -g_z_0_yyyyz_x[k] * cd_y[k] + g_z_0_yyyyz_xy[k];

                g_z_0_yyyyyz_y[k] = -g_z_0_yyyyz_y[k] * cd_y[k] + g_z_0_yyyyz_yy[k];

                g_z_0_yyyyyz_z[k] = -g_z_0_yyyyz_z[k] * cd_y[k] + g_z_0_yyyyz_yz[k];
            }

            /// Set up 69-72 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyzz_x = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 69);

            auto g_z_0_yyyyzz_y = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 70);

            auto g_z_0_yyyyzz_z = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 71);

            #pragma omp simd aligned(cd_y, g_z_0_yyyyzz_x, g_z_0_yyyyzz_y, g_z_0_yyyyzz_z, g_z_0_yyyzz_x, g_z_0_yyyzz_xy, g_z_0_yyyzz_y, g_z_0_yyyzz_yy, g_z_0_yyyzz_yz, g_z_0_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyzz_x[k] = -g_z_0_yyyzz_x[k] * cd_y[k] + g_z_0_yyyzz_xy[k];

                g_z_0_yyyyzz_y[k] = -g_z_0_yyyzz_y[k] * cd_y[k] + g_z_0_yyyzz_yy[k];

                g_z_0_yyyyzz_z[k] = -g_z_0_yyyzz_z[k] * cd_y[k] + g_z_0_yyyzz_yz[k];
            }

            /// Set up 72-75 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyzzz_x = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 72);

            auto g_z_0_yyyzzz_y = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 73);

            auto g_z_0_yyyzzz_z = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 74);

            #pragma omp simd aligned(cd_y, g_z_0_yyyzzz_x, g_z_0_yyyzzz_y, g_z_0_yyyzzz_z, g_z_0_yyzzz_x, g_z_0_yyzzz_xy, g_z_0_yyzzz_y, g_z_0_yyzzz_yy, g_z_0_yyzzz_yz, g_z_0_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyzzz_x[k] = -g_z_0_yyzzz_x[k] * cd_y[k] + g_z_0_yyzzz_xy[k];

                g_z_0_yyyzzz_y[k] = -g_z_0_yyzzz_y[k] * cd_y[k] + g_z_0_yyzzz_yy[k];

                g_z_0_yyyzzz_z[k] = -g_z_0_yyzzz_z[k] * cd_y[k] + g_z_0_yyzzz_yz[k];
            }

            /// Set up 75-78 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyzzzz_x = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 75);

            auto g_z_0_yyzzzz_y = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 76);

            auto g_z_0_yyzzzz_z = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 77);

            #pragma omp simd aligned(cd_y, g_z_0_yyzzzz_x, g_z_0_yyzzzz_y, g_z_0_yyzzzz_z, g_z_0_yzzzz_x, g_z_0_yzzzz_xy, g_z_0_yzzzz_y, g_z_0_yzzzz_yy, g_z_0_yzzzz_yz, g_z_0_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzzzz_x[k] = -g_z_0_yzzzz_x[k] * cd_y[k] + g_z_0_yzzzz_xy[k];

                g_z_0_yyzzzz_y[k] = -g_z_0_yzzzz_y[k] * cd_y[k] + g_z_0_yzzzz_yy[k];

                g_z_0_yyzzzz_z[k] = -g_z_0_yzzzz_z[k] * cd_y[k] + g_z_0_yzzzz_yz[k];
            }

            /// Set up 78-81 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzzzzz_x = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 78);

            auto g_z_0_yzzzzz_y = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 79);

            auto g_z_0_yzzzzz_z = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 80);

            #pragma omp simd aligned(cd_y, g_z_0_yzzzzz_x, g_z_0_yzzzzz_y, g_z_0_yzzzzz_z, g_z_0_zzzzz_x, g_z_0_zzzzz_xy, g_z_0_zzzzz_y, g_z_0_zzzzz_yy, g_z_0_zzzzz_yz, g_z_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzzzz_x[k] = -g_z_0_zzzzz_x[k] * cd_y[k] + g_z_0_zzzzz_xy[k];

                g_z_0_yzzzzz_y[k] = -g_z_0_zzzzz_y[k] * cd_y[k] + g_z_0_zzzzz_yy[k];

                g_z_0_yzzzzz_z[k] = -g_z_0_zzzzz_z[k] * cd_y[k] + g_z_0_zzzzz_yz[k];
            }

            /// Set up 81-84 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzzzzz_x = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 81);

            auto g_z_0_zzzzzz_y = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 82);

            auto g_z_0_zzzzzz_z = cbuffer.data(ip_geom_10_off + 168 * acomps * bcomps + 83);

            #pragma omp simd aligned(cd_z, g_z_0_zzzzz_x, g_z_0_zzzzz_xz, g_z_0_zzzzz_y, g_z_0_zzzzz_yz, g_z_0_zzzzz_z, g_z_0_zzzzz_zz, g_z_0_zzzzzz_x, g_z_0_zzzzzz_y, g_z_0_zzzzzz_z, g_zzzzz_x, g_zzzzz_y, g_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzzzz_x[k] = -g_zzzzz_x[k] - g_z_0_zzzzz_x[k] * cd_z[k] + g_z_0_zzzzz_xz[k];

                g_z_0_zzzzzz_y[k] = -g_zzzzz_y[k] - g_z_0_zzzzz_y[k] * cd_z[k] + g_z_0_zzzzz_yz[k];

                g_z_0_zzzzzz_z[k] = -g_zzzzz_z[k] - g_z_0_zzzzz_z[k] * cd_z[k] + g_z_0_zzzzz_zz[k];
            }
        }
    }
}

} // erirec namespace

