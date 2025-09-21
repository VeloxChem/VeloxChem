#include "ElectronRepulsionGeom0010ContrRecXXPI.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxpi(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxpi,
                                            CSimdArray<double>& pbuffer,
                                            const size_t idx_xxsi,
                                            const size_t idx_geom_10_xxsi,
                                            const size_t idx_geom_10_xxsk,
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

            /// Set up components of auxilary buffer : SSSI

            const auto si_geom_10_off = idx_geom_10_xxsi + (i * bcomps + j) * 28;

            auto g_x_0_0_xxxxxx = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_0_xxxxxy = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_0_xxxxxz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_0_xxxxyy = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_0_xxxxyz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_0_xxxxzz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_0_xxxyyy = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_0_xxxyyz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_0_xxxyzz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_0_xxxzzz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_0_xxyyyy = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_0_xxyyyz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_0_xxyyzz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_0_xxyzzz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_0_xxzzzz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_0_xyyyyy = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_0_xyyyyz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_0_xyyyzz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_0_xyyzzz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_0_xyzzzz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_0_xzzzzz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_0_yyyyyy = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_0_yyyyyz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_0_yyyyzz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_0_yyyzzz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_0_yyzzzz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_0_yzzzzz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_0_zzzzzz = cbuffer.data(si_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_y_0_0_xxxxxx = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 0);

            auto g_y_0_0_xxxxxy = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 1);

            auto g_y_0_0_xxxxxz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 2);

            auto g_y_0_0_xxxxyy = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 3);

            auto g_y_0_0_xxxxyz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 4);

            auto g_y_0_0_xxxxzz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 5);

            auto g_y_0_0_xxxyyy = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 6);

            auto g_y_0_0_xxxyyz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 7);

            auto g_y_0_0_xxxyzz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 8);

            auto g_y_0_0_xxxzzz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 9);

            auto g_y_0_0_xxyyyy = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 10);

            auto g_y_0_0_xxyyyz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 11);

            auto g_y_0_0_xxyyzz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 12);

            auto g_y_0_0_xxyzzz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 13);

            auto g_y_0_0_xxzzzz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 14);

            auto g_y_0_0_xyyyyy = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 15);

            auto g_y_0_0_xyyyyz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 16);

            auto g_y_0_0_xyyyzz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 17);

            auto g_y_0_0_xyyzzz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 18);

            auto g_y_0_0_xyzzzz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 19);

            auto g_y_0_0_xzzzzz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 20);

            auto g_y_0_0_yyyyyy = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 21);

            auto g_y_0_0_yyyyyz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 22);

            auto g_y_0_0_yyyyzz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 23);

            auto g_y_0_0_yyyzzz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 24);

            auto g_y_0_0_yyzzzz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 25);

            auto g_y_0_0_yzzzzz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 26);

            auto g_y_0_0_zzzzzz = cbuffer.data(si_geom_10_off + 28 * acomps * bcomps + 27);

            auto g_z_0_0_xxxxxx = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 0);

            auto g_z_0_0_xxxxxy = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 1);

            auto g_z_0_0_xxxxxz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 2);

            auto g_z_0_0_xxxxyy = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 3);

            auto g_z_0_0_xxxxyz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 4);

            auto g_z_0_0_xxxxzz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 5);

            auto g_z_0_0_xxxyyy = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 6);

            auto g_z_0_0_xxxyyz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 7);

            auto g_z_0_0_xxxyzz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 8);

            auto g_z_0_0_xxxzzz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 9);

            auto g_z_0_0_xxyyyy = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 10);

            auto g_z_0_0_xxyyyz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 11);

            auto g_z_0_0_xxyyzz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 12);

            auto g_z_0_0_xxyzzz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 13);

            auto g_z_0_0_xxzzzz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 14);

            auto g_z_0_0_xyyyyy = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 15);

            auto g_z_0_0_xyyyyz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 16);

            auto g_z_0_0_xyyyzz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 17);

            auto g_z_0_0_xyyzzz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 18);

            auto g_z_0_0_xyzzzz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 19);

            auto g_z_0_0_xzzzzz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 20);

            auto g_z_0_0_yyyyyy = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 21);

            auto g_z_0_0_yyyyyz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 22);

            auto g_z_0_0_yyyyzz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 23);

            auto g_z_0_0_yyyzzz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 24);

            auto g_z_0_0_yyzzzz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 25);

            auto g_z_0_0_yzzzzz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 26);

            auto g_z_0_0_zzzzzz = cbuffer.data(si_geom_10_off + 56 * acomps * bcomps + 27);

            /// Set up components of auxilary buffer : SSSK

            const auto sk_geom_10_off = idx_geom_10_xxsk + (i * bcomps + j) * 36;

            auto g_x_0_0_xxxxxxx = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_0_xxxxxxy = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_0_xxxxxxz = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_0_xxxxxyy = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_0_xxxxxyz = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_0_xxxxxzz = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_0_xxxxyyy = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_0_xxxxyyz = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_0_xxxxyzz = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_0_xxxxzzz = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_0_xxxyyyy = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_0_xxxyyyz = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_0_xxxyyzz = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_0_xxxyzzz = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_0_xxxzzzz = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_0_xxyyyyy = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_0_xxyyyyz = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_0_xxyyyzz = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_0_xxyyzzz = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_0_xxyzzzz = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_0_xxzzzzz = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_0_xyyyyyy = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_0_xyyyyyz = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_0_xyyyyzz = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_0_xyyyzzz = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_0_xyyzzzz = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_0_xyzzzzz = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_0_xzzzzzz = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_0_yyyyyyy = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_0_yyyyyyz = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_0_yyyyyzz = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_0_yyyyzzz = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_0_yyyzzzz = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_0_yyzzzzz = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_0_yzzzzzz = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_0_zzzzzzz = cbuffer.data(sk_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_y_0_0_xxxxxxx = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 0);

            auto g_y_0_0_xxxxxxy = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 1);

            auto g_y_0_0_xxxxxxz = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 2);

            auto g_y_0_0_xxxxxyy = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 3);

            auto g_y_0_0_xxxxxyz = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 4);

            auto g_y_0_0_xxxxxzz = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 5);

            auto g_y_0_0_xxxxyyy = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 6);

            auto g_y_0_0_xxxxyyz = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 7);

            auto g_y_0_0_xxxxyzz = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 8);

            auto g_y_0_0_xxxxzzz = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 9);

            auto g_y_0_0_xxxyyyy = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 10);

            auto g_y_0_0_xxxyyyz = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 11);

            auto g_y_0_0_xxxyyzz = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 12);

            auto g_y_0_0_xxxyzzz = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 13);

            auto g_y_0_0_xxxzzzz = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 14);

            auto g_y_0_0_xxyyyyy = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 15);

            auto g_y_0_0_xxyyyyz = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 16);

            auto g_y_0_0_xxyyyzz = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 17);

            auto g_y_0_0_xxyyzzz = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 18);

            auto g_y_0_0_xxyzzzz = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 19);

            auto g_y_0_0_xxzzzzz = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 20);

            auto g_y_0_0_xyyyyyy = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 21);

            auto g_y_0_0_xyyyyyz = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 22);

            auto g_y_0_0_xyyyyzz = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 23);

            auto g_y_0_0_xyyyzzz = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 24);

            auto g_y_0_0_xyyzzzz = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 25);

            auto g_y_0_0_xyzzzzz = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 26);

            auto g_y_0_0_xzzzzzz = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 27);

            auto g_y_0_0_yyyyyyy = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 28);

            auto g_y_0_0_yyyyyyz = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 29);

            auto g_y_0_0_yyyyyzz = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 30);

            auto g_y_0_0_yyyyzzz = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 31);

            auto g_y_0_0_yyyzzzz = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 32);

            auto g_y_0_0_yyzzzzz = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 33);

            auto g_y_0_0_yzzzzzz = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 34);

            auto g_y_0_0_zzzzzzz = cbuffer.data(sk_geom_10_off + 36 * acomps * bcomps + 35);

            auto g_z_0_0_xxxxxxx = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 0);

            auto g_z_0_0_xxxxxxy = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 1);

            auto g_z_0_0_xxxxxxz = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 2);

            auto g_z_0_0_xxxxxyy = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 3);

            auto g_z_0_0_xxxxxyz = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 4);

            auto g_z_0_0_xxxxxzz = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 5);

            auto g_z_0_0_xxxxyyy = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 6);

            auto g_z_0_0_xxxxyyz = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 7);

            auto g_z_0_0_xxxxyzz = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 8);

            auto g_z_0_0_xxxxzzz = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 9);

            auto g_z_0_0_xxxyyyy = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 10);

            auto g_z_0_0_xxxyyyz = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 11);

            auto g_z_0_0_xxxyyzz = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 12);

            auto g_z_0_0_xxxyzzz = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 13);

            auto g_z_0_0_xxxzzzz = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 14);

            auto g_z_0_0_xxyyyyy = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 15);

            auto g_z_0_0_xxyyyyz = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 16);

            auto g_z_0_0_xxyyyzz = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 17);

            auto g_z_0_0_xxyyzzz = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 18);

            auto g_z_0_0_xxyzzzz = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 19);

            auto g_z_0_0_xxzzzzz = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 20);

            auto g_z_0_0_xyyyyyy = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 21);

            auto g_z_0_0_xyyyyyz = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 22);

            auto g_z_0_0_xyyyyzz = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 23);

            auto g_z_0_0_xyyyzzz = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 24);

            auto g_z_0_0_xyyzzzz = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 25);

            auto g_z_0_0_xyzzzzz = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 26);

            auto g_z_0_0_xzzzzzz = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 27);

            auto g_z_0_0_yyyyyyy = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 28);

            auto g_z_0_0_yyyyyyz = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 29);

            auto g_z_0_0_yyyyyzz = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 30);

            auto g_z_0_0_yyyyzzz = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 31);

            auto g_z_0_0_yyyzzzz = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 32);

            auto g_z_0_0_yyzzzzz = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 33);

            auto g_z_0_0_yzzzzzz = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 34);

            auto g_z_0_0_zzzzzzz = cbuffer.data(sk_geom_10_off + 72 * acomps * bcomps + 35);

            /// set up bra offset for contr_buffer_xxpi

            const auto pi_geom_10_off = idx_geom_10_xxpi + (i * bcomps + j) * 84;

            /// Set up 0-28 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_xxxxxx = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_x_xxxxxy = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_x_xxxxxz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_x_xxxxyy = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_x_xxxxyz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_x_xxxxzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_x_xxxyyy = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_x_xxxyyz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_x_xxxyzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_x_xxxzzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_x_xxyyyy = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_x_xxyyyz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_x_xxyyzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_x_xxyzzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_x_xxzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_x_xyyyyy = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_x_xyyyyz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_x_xyyyzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_x_xyyzzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_x_xyzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_x_xzzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_x_yyyyyy = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_x_yyyyyz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_x_yyyyzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_x_yyyzzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_x_yyzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_x_yzzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_x_zzzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 27);

            #pragma omp simd aligned(cd_x, g_0_xxxxxx, g_0_xxxxxy, g_0_xxxxxz, g_0_xxxxyy, g_0_xxxxyz, g_0_xxxxzz, g_0_xxxyyy, g_0_xxxyyz, g_0_xxxyzz, g_0_xxxzzz, g_0_xxyyyy, g_0_xxyyyz, g_0_xxyyzz, g_0_xxyzzz, g_0_xxzzzz, g_0_xyyyyy, g_0_xyyyyz, g_0_xyyyzz, g_0_xyyzzz, g_0_xyzzzz, g_0_xzzzzz, g_0_yyyyyy, g_0_yyyyyz, g_0_yyyyzz, g_0_yyyzzz, g_0_yyzzzz, g_0_yzzzzz, g_0_zzzzzz, g_x_0_0_xxxxxx, g_x_0_0_xxxxxxx, g_x_0_0_xxxxxxy, g_x_0_0_xxxxxxz, g_x_0_0_xxxxxy, g_x_0_0_xxxxxyy, g_x_0_0_xxxxxyz, g_x_0_0_xxxxxz, g_x_0_0_xxxxxzz, g_x_0_0_xxxxyy, g_x_0_0_xxxxyyy, g_x_0_0_xxxxyyz, g_x_0_0_xxxxyz, g_x_0_0_xxxxyzz, g_x_0_0_xxxxzz, g_x_0_0_xxxxzzz, g_x_0_0_xxxyyy, g_x_0_0_xxxyyyy, g_x_0_0_xxxyyyz, g_x_0_0_xxxyyz, g_x_0_0_xxxyyzz, g_x_0_0_xxxyzz, g_x_0_0_xxxyzzz, g_x_0_0_xxxzzz, g_x_0_0_xxxzzzz, g_x_0_0_xxyyyy, g_x_0_0_xxyyyyy, g_x_0_0_xxyyyyz, g_x_0_0_xxyyyz, g_x_0_0_xxyyyzz, g_x_0_0_xxyyzz, g_x_0_0_xxyyzzz, g_x_0_0_xxyzzz, g_x_0_0_xxyzzzz, g_x_0_0_xxzzzz, g_x_0_0_xxzzzzz, g_x_0_0_xyyyyy, g_x_0_0_xyyyyyy, g_x_0_0_xyyyyyz, g_x_0_0_xyyyyz, g_x_0_0_xyyyyzz, g_x_0_0_xyyyzz, g_x_0_0_xyyyzzz, g_x_0_0_xyyzzz, g_x_0_0_xyyzzzz, g_x_0_0_xyzzzz, g_x_0_0_xyzzzzz, g_x_0_0_xzzzzz, g_x_0_0_xzzzzzz, g_x_0_0_yyyyyy, g_x_0_0_yyyyyz, g_x_0_0_yyyyzz, g_x_0_0_yyyzzz, g_x_0_0_yyzzzz, g_x_0_0_yzzzzz, g_x_0_0_zzzzzz, g_x_0_x_xxxxxx, g_x_0_x_xxxxxy, g_x_0_x_xxxxxz, g_x_0_x_xxxxyy, g_x_0_x_xxxxyz, g_x_0_x_xxxxzz, g_x_0_x_xxxyyy, g_x_0_x_xxxyyz, g_x_0_x_xxxyzz, g_x_0_x_xxxzzz, g_x_0_x_xxyyyy, g_x_0_x_xxyyyz, g_x_0_x_xxyyzz, g_x_0_x_xxyzzz, g_x_0_x_xxzzzz, g_x_0_x_xyyyyy, g_x_0_x_xyyyyz, g_x_0_x_xyyyzz, g_x_0_x_xyyzzz, g_x_0_x_xyzzzz, g_x_0_x_xzzzzz, g_x_0_x_yyyyyy, g_x_0_x_yyyyyz, g_x_0_x_yyyyzz, g_x_0_x_yyyzzz, g_x_0_x_yyzzzz, g_x_0_x_yzzzzz, g_x_0_x_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_xxxxxx[k] = -g_0_xxxxxx[k] - g_x_0_0_xxxxxx[k] * cd_x[k] + g_x_0_0_xxxxxxx[k];

                g_x_0_x_xxxxxy[k] = -g_0_xxxxxy[k] - g_x_0_0_xxxxxy[k] * cd_x[k] + g_x_0_0_xxxxxxy[k];

                g_x_0_x_xxxxxz[k] = -g_0_xxxxxz[k] - g_x_0_0_xxxxxz[k] * cd_x[k] + g_x_0_0_xxxxxxz[k];

                g_x_0_x_xxxxyy[k] = -g_0_xxxxyy[k] - g_x_0_0_xxxxyy[k] * cd_x[k] + g_x_0_0_xxxxxyy[k];

                g_x_0_x_xxxxyz[k] = -g_0_xxxxyz[k] - g_x_0_0_xxxxyz[k] * cd_x[k] + g_x_0_0_xxxxxyz[k];

                g_x_0_x_xxxxzz[k] = -g_0_xxxxzz[k] - g_x_0_0_xxxxzz[k] * cd_x[k] + g_x_0_0_xxxxxzz[k];

                g_x_0_x_xxxyyy[k] = -g_0_xxxyyy[k] - g_x_0_0_xxxyyy[k] * cd_x[k] + g_x_0_0_xxxxyyy[k];

                g_x_0_x_xxxyyz[k] = -g_0_xxxyyz[k] - g_x_0_0_xxxyyz[k] * cd_x[k] + g_x_0_0_xxxxyyz[k];

                g_x_0_x_xxxyzz[k] = -g_0_xxxyzz[k] - g_x_0_0_xxxyzz[k] * cd_x[k] + g_x_0_0_xxxxyzz[k];

                g_x_0_x_xxxzzz[k] = -g_0_xxxzzz[k] - g_x_0_0_xxxzzz[k] * cd_x[k] + g_x_0_0_xxxxzzz[k];

                g_x_0_x_xxyyyy[k] = -g_0_xxyyyy[k] - g_x_0_0_xxyyyy[k] * cd_x[k] + g_x_0_0_xxxyyyy[k];

                g_x_0_x_xxyyyz[k] = -g_0_xxyyyz[k] - g_x_0_0_xxyyyz[k] * cd_x[k] + g_x_0_0_xxxyyyz[k];

                g_x_0_x_xxyyzz[k] = -g_0_xxyyzz[k] - g_x_0_0_xxyyzz[k] * cd_x[k] + g_x_0_0_xxxyyzz[k];

                g_x_0_x_xxyzzz[k] = -g_0_xxyzzz[k] - g_x_0_0_xxyzzz[k] * cd_x[k] + g_x_0_0_xxxyzzz[k];

                g_x_0_x_xxzzzz[k] = -g_0_xxzzzz[k] - g_x_0_0_xxzzzz[k] * cd_x[k] + g_x_0_0_xxxzzzz[k];

                g_x_0_x_xyyyyy[k] = -g_0_xyyyyy[k] - g_x_0_0_xyyyyy[k] * cd_x[k] + g_x_0_0_xxyyyyy[k];

                g_x_0_x_xyyyyz[k] = -g_0_xyyyyz[k] - g_x_0_0_xyyyyz[k] * cd_x[k] + g_x_0_0_xxyyyyz[k];

                g_x_0_x_xyyyzz[k] = -g_0_xyyyzz[k] - g_x_0_0_xyyyzz[k] * cd_x[k] + g_x_0_0_xxyyyzz[k];

                g_x_0_x_xyyzzz[k] = -g_0_xyyzzz[k] - g_x_0_0_xyyzzz[k] * cd_x[k] + g_x_0_0_xxyyzzz[k];

                g_x_0_x_xyzzzz[k] = -g_0_xyzzzz[k] - g_x_0_0_xyzzzz[k] * cd_x[k] + g_x_0_0_xxyzzzz[k];

                g_x_0_x_xzzzzz[k] = -g_0_xzzzzz[k] - g_x_0_0_xzzzzz[k] * cd_x[k] + g_x_0_0_xxzzzzz[k];

                g_x_0_x_yyyyyy[k] = -g_0_yyyyyy[k] - g_x_0_0_yyyyyy[k] * cd_x[k] + g_x_0_0_xyyyyyy[k];

                g_x_0_x_yyyyyz[k] = -g_0_yyyyyz[k] - g_x_0_0_yyyyyz[k] * cd_x[k] + g_x_0_0_xyyyyyz[k];

                g_x_0_x_yyyyzz[k] = -g_0_yyyyzz[k] - g_x_0_0_yyyyzz[k] * cd_x[k] + g_x_0_0_xyyyyzz[k];

                g_x_0_x_yyyzzz[k] = -g_0_yyyzzz[k] - g_x_0_0_yyyzzz[k] * cd_x[k] + g_x_0_0_xyyyzzz[k];

                g_x_0_x_yyzzzz[k] = -g_0_yyzzzz[k] - g_x_0_0_yyzzzz[k] * cd_x[k] + g_x_0_0_xyyzzzz[k];

                g_x_0_x_yzzzzz[k] = -g_0_yzzzzz[k] - g_x_0_0_yzzzzz[k] * cd_x[k] + g_x_0_0_xyzzzzz[k];

                g_x_0_x_zzzzzz[k] = -g_0_zzzzzz[k] - g_x_0_0_zzzzzz[k] * cd_x[k] + g_x_0_0_xzzzzzz[k];
            }

            /// Set up 28-56 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_xxxxxx = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_y_xxxxxy = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_y_xxxxxz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_y_xxxxyy = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_y_xxxxyz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_y_xxxxzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_y_xxxyyy = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_y_xxxyyz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_y_xxxyzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_y_xxxzzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_y_xxyyyy = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_y_xxyyyz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_y_xxyyzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_y_xxyzzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_y_xxzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_y_xyyyyy = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_y_xyyyyz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_y_xyyyzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_y_xyyzzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_y_xyzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_y_xzzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_y_yyyyyy = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_y_yyyyyz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_y_yyyyzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_y_yyyzzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_y_yyzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_y_yzzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_y_zzzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 55);

            #pragma omp simd aligned(cd_y, g_x_0_0_xxxxxx, g_x_0_0_xxxxxxy, g_x_0_0_xxxxxy, g_x_0_0_xxxxxyy, g_x_0_0_xxxxxyz, g_x_0_0_xxxxxz, g_x_0_0_xxxxyy, g_x_0_0_xxxxyyy, g_x_0_0_xxxxyyz, g_x_0_0_xxxxyz, g_x_0_0_xxxxyzz, g_x_0_0_xxxxzz, g_x_0_0_xxxyyy, g_x_0_0_xxxyyyy, g_x_0_0_xxxyyyz, g_x_0_0_xxxyyz, g_x_0_0_xxxyyzz, g_x_0_0_xxxyzz, g_x_0_0_xxxyzzz, g_x_0_0_xxxzzz, g_x_0_0_xxyyyy, g_x_0_0_xxyyyyy, g_x_0_0_xxyyyyz, g_x_0_0_xxyyyz, g_x_0_0_xxyyyzz, g_x_0_0_xxyyzz, g_x_0_0_xxyyzzz, g_x_0_0_xxyzzz, g_x_0_0_xxyzzzz, g_x_0_0_xxzzzz, g_x_0_0_xyyyyy, g_x_0_0_xyyyyyy, g_x_0_0_xyyyyyz, g_x_0_0_xyyyyz, g_x_0_0_xyyyyzz, g_x_0_0_xyyyzz, g_x_0_0_xyyyzzz, g_x_0_0_xyyzzz, g_x_0_0_xyyzzzz, g_x_0_0_xyzzzz, g_x_0_0_xyzzzzz, g_x_0_0_xzzzzz, g_x_0_0_yyyyyy, g_x_0_0_yyyyyyy, g_x_0_0_yyyyyyz, g_x_0_0_yyyyyz, g_x_0_0_yyyyyzz, g_x_0_0_yyyyzz, g_x_0_0_yyyyzzz, g_x_0_0_yyyzzz, g_x_0_0_yyyzzzz, g_x_0_0_yyzzzz, g_x_0_0_yyzzzzz, g_x_0_0_yzzzzz, g_x_0_0_yzzzzzz, g_x_0_0_zzzzzz, g_x_0_y_xxxxxx, g_x_0_y_xxxxxy, g_x_0_y_xxxxxz, g_x_0_y_xxxxyy, g_x_0_y_xxxxyz, g_x_0_y_xxxxzz, g_x_0_y_xxxyyy, g_x_0_y_xxxyyz, g_x_0_y_xxxyzz, g_x_0_y_xxxzzz, g_x_0_y_xxyyyy, g_x_0_y_xxyyyz, g_x_0_y_xxyyzz, g_x_0_y_xxyzzz, g_x_0_y_xxzzzz, g_x_0_y_xyyyyy, g_x_0_y_xyyyyz, g_x_0_y_xyyyzz, g_x_0_y_xyyzzz, g_x_0_y_xyzzzz, g_x_0_y_xzzzzz, g_x_0_y_yyyyyy, g_x_0_y_yyyyyz, g_x_0_y_yyyyzz, g_x_0_y_yyyzzz, g_x_0_y_yyzzzz, g_x_0_y_yzzzzz, g_x_0_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_xxxxxx[k] = -g_x_0_0_xxxxxx[k] * cd_y[k] + g_x_0_0_xxxxxxy[k];

                g_x_0_y_xxxxxy[k] = -g_x_0_0_xxxxxy[k] * cd_y[k] + g_x_0_0_xxxxxyy[k];

                g_x_0_y_xxxxxz[k] = -g_x_0_0_xxxxxz[k] * cd_y[k] + g_x_0_0_xxxxxyz[k];

                g_x_0_y_xxxxyy[k] = -g_x_0_0_xxxxyy[k] * cd_y[k] + g_x_0_0_xxxxyyy[k];

                g_x_0_y_xxxxyz[k] = -g_x_0_0_xxxxyz[k] * cd_y[k] + g_x_0_0_xxxxyyz[k];

                g_x_0_y_xxxxzz[k] = -g_x_0_0_xxxxzz[k] * cd_y[k] + g_x_0_0_xxxxyzz[k];

                g_x_0_y_xxxyyy[k] = -g_x_0_0_xxxyyy[k] * cd_y[k] + g_x_0_0_xxxyyyy[k];

                g_x_0_y_xxxyyz[k] = -g_x_0_0_xxxyyz[k] * cd_y[k] + g_x_0_0_xxxyyyz[k];

                g_x_0_y_xxxyzz[k] = -g_x_0_0_xxxyzz[k] * cd_y[k] + g_x_0_0_xxxyyzz[k];

                g_x_0_y_xxxzzz[k] = -g_x_0_0_xxxzzz[k] * cd_y[k] + g_x_0_0_xxxyzzz[k];

                g_x_0_y_xxyyyy[k] = -g_x_0_0_xxyyyy[k] * cd_y[k] + g_x_0_0_xxyyyyy[k];

                g_x_0_y_xxyyyz[k] = -g_x_0_0_xxyyyz[k] * cd_y[k] + g_x_0_0_xxyyyyz[k];

                g_x_0_y_xxyyzz[k] = -g_x_0_0_xxyyzz[k] * cd_y[k] + g_x_0_0_xxyyyzz[k];

                g_x_0_y_xxyzzz[k] = -g_x_0_0_xxyzzz[k] * cd_y[k] + g_x_0_0_xxyyzzz[k];

                g_x_0_y_xxzzzz[k] = -g_x_0_0_xxzzzz[k] * cd_y[k] + g_x_0_0_xxyzzzz[k];

                g_x_0_y_xyyyyy[k] = -g_x_0_0_xyyyyy[k] * cd_y[k] + g_x_0_0_xyyyyyy[k];

                g_x_0_y_xyyyyz[k] = -g_x_0_0_xyyyyz[k] * cd_y[k] + g_x_0_0_xyyyyyz[k];

                g_x_0_y_xyyyzz[k] = -g_x_0_0_xyyyzz[k] * cd_y[k] + g_x_0_0_xyyyyzz[k];

                g_x_0_y_xyyzzz[k] = -g_x_0_0_xyyzzz[k] * cd_y[k] + g_x_0_0_xyyyzzz[k];

                g_x_0_y_xyzzzz[k] = -g_x_0_0_xyzzzz[k] * cd_y[k] + g_x_0_0_xyyzzzz[k];

                g_x_0_y_xzzzzz[k] = -g_x_0_0_xzzzzz[k] * cd_y[k] + g_x_0_0_xyzzzzz[k];

                g_x_0_y_yyyyyy[k] = -g_x_0_0_yyyyyy[k] * cd_y[k] + g_x_0_0_yyyyyyy[k];

                g_x_0_y_yyyyyz[k] = -g_x_0_0_yyyyyz[k] * cd_y[k] + g_x_0_0_yyyyyyz[k];

                g_x_0_y_yyyyzz[k] = -g_x_0_0_yyyyzz[k] * cd_y[k] + g_x_0_0_yyyyyzz[k];

                g_x_0_y_yyyzzz[k] = -g_x_0_0_yyyzzz[k] * cd_y[k] + g_x_0_0_yyyyzzz[k];

                g_x_0_y_yyzzzz[k] = -g_x_0_0_yyzzzz[k] * cd_y[k] + g_x_0_0_yyyzzzz[k];

                g_x_0_y_yzzzzz[k] = -g_x_0_0_yzzzzz[k] * cd_y[k] + g_x_0_0_yyzzzzz[k];

                g_x_0_y_zzzzzz[k] = -g_x_0_0_zzzzzz[k] * cd_y[k] + g_x_0_0_yzzzzzz[k];
            }

            /// Set up 56-84 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_xxxxxx = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_z_xxxxxy = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_z_xxxxxz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_z_xxxxyy = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 59);

            auto g_x_0_z_xxxxyz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_z_xxxxzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_z_xxxyyy = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_x_0_z_xxxyyz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_z_xxxyzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_z_xxxzzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 65);

            auto g_x_0_z_xxyyyy = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_z_xxyyyz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_z_xxyyzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_z_xxyzzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_z_xxzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_z_xyyyyy = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_z_xyyyyz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_z_xyyyzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_z_xyyzzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_z_xyzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_z_xzzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_z_yyyyyy = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_z_yyyyyz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_z_yyyyzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_z_yyyzzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_z_yyzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_z_yzzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_z_zzzzzz = cbuffer.data(pi_geom_10_off + 0 * acomps * bcomps + 83);

            #pragma omp simd aligned(cd_z, g_x_0_0_xxxxxx, g_x_0_0_xxxxxxz, g_x_0_0_xxxxxy, g_x_0_0_xxxxxyz, g_x_0_0_xxxxxz, g_x_0_0_xxxxxzz, g_x_0_0_xxxxyy, g_x_0_0_xxxxyyz, g_x_0_0_xxxxyz, g_x_0_0_xxxxyzz, g_x_0_0_xxxxzz, g_x_0_0_xxxxzzz, g_x_0_0_xxxyyy, g_x_0_0_xxxyyyz, g_x_0_0_xxxyyz, g_x_0_0_xxxyyzz, g_x_0_0_xxxyzz, g_x_0_0_xxxyzzz, g_x_0_0_xxxzzz, g_x_0_0_xxxzzzz, g_x_0_0_xxyyyy, g_x_0_0_xxyyyyz, g_x_0_0_xxyyyz, g_x_0_0_xxyyyzz, g_x_0_0_xxyyzz, g_x_0_0_xxyyzzz, g_x_0_0_xxyzzz, g_x_0_0_xxyzzzz, g_x_0_0_xxzzzz, g_x_0_0_xxzzzzz, g_x_0_0_xyyyyy, g_x_0_0_xyyyyyz, g_x_0_0_xyyyyz, g_x_0_0_xyyyyzz, g_x_0_0_xyyyzz, g_x_0_0_xyyyzzz, g_x_0_0_xyyzzz, g_x_0_0_xyyzzzz, g_x_0_0_xyzzzz, g_x_0_0_xyzzzzz, g_x_0_0_xzzzzz, g_x_0_0_xzzzzzz, g_x_0_0_yyyyyy, g_x_0_0_yyyyyyz, g_x_0_0_yyyyyz, g_x_0_0_yyyyyzz, g_x_0_0_yyyyzz, g_x_0_0_yyyyzzz, g_x_0_0_yyyzzz, g_x_0_0_yyyzzzz, g_x_0_0_yyzzzz, g_x_0_0_yyzzzzz, g_x_0_0_yzzzzz, g_x_0_0_yzzzzzz, g_x_0_0_zzzzzz, g_x_0_0_zzzzzzz, g_x_0_z_xxxxxx, g_x_0_z_xxxxxy, g_x_0_z_xxxxxz, g_x_0_z_xxxxyy, g_x_0_z_xxxxyz, g_x_0_z_xxxxzz, g_x_0_z_xxxyyy, g_x_0_z_xxxyyz, g_x_0_z_xxxyzz, g_x_0_z_xxxzzz, g_x_0_z_xxyyyy, g_x_0_z_xxyyyz, g_x_0_z_xxyyzz, g_x_0_z_xxyzzz, g_x_0_z_xxzzzz, g_x_0_z_xyyyyy, g_x_0_z_xyyyyz, g_x_0_z_xyyyzz, g_x_0_z_xyyzzz, g_x_0_z_xyzzzz, g_x_0_z_xzzzzz, g_x_0_z_yyyyyy, g_x_0_z_yyyyyz, g_x_0_z_yyyyzz, g_x_0_z_yyyzzz, g_x_0_z_yyzzzz, g_x_0_z_yzzzzz, g_x_0_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_xxxxxx[k] = -g_x_0_0_xxxxxx[k] * cd_z[k] + g_x_0_0_xxxxxxz[k];

                g_x_0_z_xxxxxy[k] = -g_x_0_0_xxxxxy[k] * cd_z[k] + g_x_0_0_xxxxxyz[k];

                g_x_0_z_xxxxxz[k] = -g_x_0_0_xxxxxz[k] * cd_z[k] + g_x_0_0_xxxxxzz[k];

                g_x_0_z_xxxxyy[k] = -g_x_0_0_xxxxyy[k] * cd_z[k] + g_x_0_0_xxxxyyz[k];

                g_x_0_z_xxxxyz[k] = -g_x_0_0_xxxxyz[k] * cd_z[k] + g_x_0_0_xxxxyzz[k];

                g_x_0_z_xxxxzz[k] = -g_x_0_0_xxxxzz[k] * cd_z[k] + g_x_0_0_xxxxzzz[k];

                g_x_0_z_xxxyyy[k] = -g_x_0_0_xxxyyy[k] * cd_z[k] + g_x_0_0_xxxyyyz[k];

                g_x_0_z_xxxyyz[k] = -g_x_0_0_xxxyyz[k] * cd_z[k] + g_x_0_0_xxxyyzz[k];

                g_x_0_z_xxxyzz[k] = -g_x_0_0_xxxyzz[k] * cd_z[k] + g_x_0_0_xxxyzzz[k];

                g_x_0_z_xxxzzz[k] = -g_x_0_0_xxxzzz[k] * cd_z[k] + g_x_0_0_xxxzzzz[k];

                g_x_0_z_xxyyyy[k] = -g_x_0_0_xxyyyy[k] * cd_z[k] + g_x_0_0_xxyyyyz[k];

                g_x_0_z_xxyyyz[k] = -g_x_0_0_xxyyyz[k] * cd_z[k] + g_x_0_0_xxyyyzz[k];

                g_x_0_z_xxyyzz[k] = -g_x_0_0_xxyyzz[k] * cd_z[k] + g_x_0_0_xxyyzzz[k];

                g_x_0_z_xxyzzz[k] = -g_x_0_0_xxyzzz[k] * cd_z[k] + g_x_0_0_xxyzzzz[k];

                g_x_0_z_xxzzzz[k] = -g_x_0_0_xxzzzz[k] * cd_z[k] + g_x_0_0_xxzzzzz[k];

                g_x_0_z_xyyyyy[k] = -g_x_0_0_xyyyyy[k] * cd_z[k] + g_x_0_0_xyyyyyz[k];

                g_x_0_z_xyyyyz[k] = -g_x_0_0_xyyyyz[k] * cd_z[k] + g_x_0_0_xyyyyzz[k];

                g_x_0_z_xyyyzz[k] = -g_x_0_0_xyyyzz[k] * cd_z[k] + g_x_0_0_xyyyzzz[k];

                g_x_0_z_xyyzzz[k] = -g_x_0_0_xyyzzz[k] * cd_z[k] + g_x_0_0_xyyzzzz[k];

                g_x_0_z_xyzzzz[k] = -g_x_0_0_xyzzzz[k] * cd_z[k] + g_x_0_0_xyzzzzz[k];

                g_x_0_z_xzzzzz[k] = -g_x_0_0_xzzzzz[k] * cd_z[k] + g_x_0_0_xzzzzzz[k];

                g_x_0_z_yyyyyy[k] = -g_x_0_0_yyyyyy[k] * cd_z[k] + g_x_0_0_yyyyyyz[k];

                g_x_0_z_yyyyyz[k] = -g_x_0_0_yyyyyz[k] * cd_z[k] + g_x_0_0_yyyyyzz[k];

                g_x_0_z_yyyyzz[k] = -g_x_0_0_yyyyzz[k] * cd_z[k] + g_x_0_0_yyyyzzz[k];

                g_x_0_z_yyyzzz[k] = -g_x_0_0_yyyzzz[k] * cd_z[k] + g_x_0_0_yyyzzzz[k];

                g_x_0_z_yyzzzz[k] = -g_x_0_0_yyzzzz[k] * cd_z[k] + g_x_0_0_yyzzzzz[k];

                g_x_0_z_yzzzzz[k] = -g_x_0_0_yzzzzz[k] * cd_z[k] + g_x_0_0_yzzzzzz[k];

                g_x_0_z_zzzzzz[k] = -g_x_0_0_zzzzzz[k] * cd_z[k] + g_x_0_0_zzzzzzz[k];
            }
            /// Set up 0-28 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_xxxxxx = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 0);

            auto g_y_0_x_xxxxxy = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 1);

            auto g_y_0_x_xxxxxz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 2);

            auto g_y_0_x_xxxxyy = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 3);

            auto g_y_0_x_xxxxyz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 4);

            auto g_y_0_x_xxxxzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 5);

            auto g_y_0_x_xxxyyy = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 6);

            auto g_y_0_x_xxxyyz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 7);

            auto g_y_0_x_xxxyzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 8);

            auto g_y_0_x_xxxzzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 9);

            auto g_y_0_x_xxyyyy = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 10);

            auto g_y_0_x_xxyyyz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 11);

            auto g_y_0_x_xxyyzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 12);

            auto g_y_0_x_xxyzzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 13);

            auto g_y_0_x_xxzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 14);

            auto g_y_0_x_xyyyyy = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 15);

            auto g_y_0_x_xyyyyz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 16);

            auto g_y_0_x_xyyyzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 17);

            auto g_y_0_x_xyyzzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 18);

            auto g_y_0_x_xyzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 19);

            auto g_y_0_x_xzzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 20);

            auto g_y_0_x_yyyyyy = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 21);

            auto g_y_0_x_yyyyyz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 22);

            auto g_y_0_x_yyyyzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 23);

            auto g_y_0_x_yyyzzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 24);

            auto g_y_0_x_yyzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 25);

            auto g_y_0_x_yzzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 26);

            auto g_y_0_x_zzzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 27);

            #pragma omp simd aligned(cd_x, g_y_0_0_xxxxxx, g_y_0_0_xxxxxxx, g_y_0_0_xxxxxxy, g_y_0_0_xxxxxxz, g_y_0_0_xxxxxy, g_y_0_0_xxxxxyy, g_y_0_0_xxxxxyz, g_y_0_0_xxxxxz, g_y_0_0_xxxxxzz, g_y_0_0_xxxxyy, g_y_0_0_xxxxyyy, g_y_0_0_xxxxyyz, g_y_0_0_xxxxyz, g_y_0_0_xxxxyzz, g_y_0_0_xxxxzz, g_y_0_0_xxxxzzz, g_y_0_0_xxxyyy, g_y_0_0_xxxyyyy, g_y_0_0_xxxyyyz, g_y_0_0_xxxyyz, g_y_0_0_xxxyyzz, g_y_0_0_xxxyzz, g_y_0_0_xxxyzzz, g_y_0_0_xxxzzz, g_y_0_0_xxxzzzz, g_y_0_0_xxyyyy, g_y_0_0_xxyyyyy, g_y_0_0_xxyyyyz, g_y_0_0_xxyyyz, g_y_0_0_xxyyyzz, g_y_0_0_xxyyzz, g_y_0_0_xxyyzzz, g_y_0_0_xxyzzz, g_y_0_0_xxyzzzz, g_y_0_0_xxzzzz, g_y_0_0_xxzzzzz, g_y_0_0_xyyyyy, g_y_0_0_xyyyyyy, g_y_0_0_xyyyyyz, g_y_0_0_xyyyyz, g_y_0_0_xyyyyzz, g_y_0_0_xyyyzz, g_y_0_0_xyyyzzz, g_y_0_0_xyyzzz, g_y_0_0_xyyzzzz, g_y_0_0_xyzzzz, g_y_0_0_xyzzzzz, g_y_0_0_xzzzzz, g_y_0_0_xzzzzzz, g_y_0_0_yyyyyy, g_y_0_0_yyyyyz, g_y_0_0_yyyyzz, g_y_0_0_yyyzzz, g_y_0_0_yyzzzz, g_y_0_0_yzzzzz, g_y_0_0_zzzzzz, g_y_0_x_xxxxxx, g_y_0_x_xxxxxy, g_y_0_x_xxxxxz, g_y_0_x_xxxxyy, g_y_0_x_xxxxyz, g_y_0_x_xxxxzz, g_y_0_x_xxxyyy, g_y_0_x_xxxyyz, g_y_0_x_xxxyzz, g_y_0_x_xxxzzz, g_y_0_x_xxyyyy, g_y_0_x_xxyyyz, g_y_0_x_xxyyzz, g_y_0_x_xxyzzz, g_y_0_x_xxzzzz, g_y_0_x_xyyyyy, g_y_0_x_xyyyyz, g_y_0_x_xyyyzz, g_y_0_x_xyyzzz, g_y_0_x_xyzzzz, g_y_0_x_xzzzzz, g_y_0_x_yyyyyy, g_y_0_x_yyyyyz, g_y_0_x_yyyyzz, g_y_0_x_yyyzzz, g_y_0_x_yyzzzz, g_y_0_x_yzzzzz, g_y_0_x_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_xxxxxx[k] = -g_y_0_0_xxxxxx[k] * cd_x[k] + g_y_0_0_xxxxxxx[k];

                g_y_0_x_xxxxxy[k] = -g_y_0_0_xxxxxy[k] * cd_x[k] + g_y_0_0_xxxxxxy[k];

                g_y_0_x_xxxxxz[k] = -g_y_0_0_xxxxxz[k] * cd_x[k] + g_y_0_0_xxxxxxz[k];

                g_y_0_x_xxxxyy[k] = -g_y_0_0_xxxxyy[k] * cd_x[k] + g_y_0_0_xxxxxyy[k];

                g_y_0_x_xxxxyz[k] = -g_y_0_0_xxxxyz[k] * cd_x[k] + g_y_0_0_xxxxxyz[k];

                g_y_0_x_xxxxzz[k] = -g_y_0_0_xxxxzz[k] * cd_x[k] + g_y_0_0_xxxxxzz[k];

                g_y_0_x_xxxyyy[k] = -g_y_0_0_xxxyyy[k] * cd_x[k] + g_y_0_0_xxxxyyy[k];

                g_y_0_x_xxxyyz[k] = -g_y_0_0_xxxyyz[k] * cd_x[k] + g_y_0_0_xxxxyyz[k];

                g_y_0_x_xxxyzz[k] = -g_y_0_0_xxxyzz[k] * cd_x[k] + g_y_0_0_xxxxyzz[k];

                g_y_0_x_xxxzzz[k] = -g_y_0_0_xxxzzz[k] * cd_x[k] + g_y_0_0_xxxxzzz[k];

                g_y_0_x_xxyyyy[k] = -g_y_0_0_xxyyyy[k] * cd_x[k] + g_y_0_0_xxxyyyy[k];

                g_y_0_x_xxyyyz[k] = -g_y_0_0_xxyyyz[k] * cd_x[k] + g_y_0_0_xxxyyyz[k];

                g_y_0_x_xxyyzz[k] = -g_y_0_0_xxyyzz[k] * cd_x[k] + g_y_0_0_xxxyyzz[k];

                g_y_0_x_xxyzzz[k] = -g_y_0_0_xxyzzz[k] * cd_x[k] + g_y_0_0_xxxyzzz[k];

                g_y_0_x_xxzzzz[k] = -g_y_0_0_xxzzzz[k] * cd_x[k] + g_y_0_0_xxxzzzz[k];

                g_y_0_x_xyyyyy[k] = -g_y_0_0_xyyyyy[k] * cd_x[k] + g_y_0_0_xxyyyyy[k];

                g_y_0_x_xyyyyz[k] = -g_y_0_0_xyyyyz[k] * cd_x[k] + g_y_0_0_xxyyyyz[k];

                g_y_0_x_xyyyzz[k] = -g_y_0_0_xyyyzz[k] * cd_x[k] + g_y_0_0_xxyyyzz[k];

                g_y_0_x_xyyzzz[k] = -g_y_0_0_xyyzzz[k] * cd_x[k] + g_y_0_0_xxyyzzz[k];

                g_y_0_x_xyzzzz[k] = -g_y_0_0_xyzzzz[k] * cd_x[k] + g_y_0_0_xxyzzzz[k];

                g_y_0_x_xzzzzz[k] = -g_y_0_0_xzzzzz[k] * cd_x[k] + g_y_0_0_xxzzzzz[k];

                g_y_0_x_yyyyyy[k] = -g_y_0_0_yyyyyy[k] * cd_x[k] + g_y_0_0_xyyyyyy[k];

                g_y_0_x_yyyyyz[k] = -g_y_0_0_yyyyyz[k] * cd_x[k] + g_y_0_0_xyyyyyz[k];

                g_y_0_x_yyyyzz[k] = -g_y_0_0_yyyyzz[k] * cd_x[k] + g_y_0_0_xyyyyzz[k];

                g_y_0_x_yyyzzz[k] = -g_y_0_0_yyyzzz[k] * cd_x[k] + g_y_0_0_xyyyzzz[k];

                g_y_0_x_yyzzzz[k] = -g_y_0_0_yyzzzz[k] * cd_x[k] + g_y_0_0_xyyzzzz[k];

                g_y_0_x_yzzzzz[k] = -g_y_0_0_yzzzzz[k] * cd_x[k] + g_y_0_0_xyzzzzz[k];

                g_y_0_x_zzzzzz[k] = -g_y_0_0_zzzzzz[k] * cd_x[k] + g_y_0_0_xzzzzzz[k];
            }

            /// Set up 28-56 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_xxxxxx = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 28);

            auto g_y_0_y_xxxxxy = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 29);

            auto g_y_0_y_xxxxxz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 30);

            auto g_y_0_y_xxxxyy = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 31);

            auto g_y_0_y_xxxxyz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 32);

            auto g_y_0_y_xxxxzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 33);

            auto g_y_0_y_xxxyyy = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 34);

            auto g_y_0_y_xxxyyz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 35);

            auto g_y_0_y_xxxyzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 36);

            auto g_y_0_y_xxxzzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 37);

            auto g_y_0_y_xxyyyy = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 38);

            auto g_y_0_y_xxyyyz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 39);

            auto g_y_0_y_xxyyzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 40);

            auto g_y_0_y_xxyzzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 41);

            auto g_y_0_y_xxzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 42);

            auto g_y_0_y_xyyyyy = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 43);

            auto g_y_0_y_xyyyyz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 44);

            auto g_y_0_y_xyyyzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 45);

            auto g_y_0_y_xyyzzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 46);

            auto g_y_0_y_xyzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 47);

            auto g_y_0_y_xzzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 48);

            auto g_y_0_y_yyyyyy = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 49);

            auto g_y_0_y_yyyyyz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 50);

            auto g_y_0_y_yyyyzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 51);

            auto g_y_0_y_yyyzzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 52);

            auto g_y_0_y_yyzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 53);

            auto g_y_0_y_yzzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 54);

            auto g_y_0_y_zzzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 55);

            #pragma omp simd aligned(cd_y, g_0_xxxxxx, g_0_xxxxxy, g_0_xxxxxz, g_0_xxxxyy, g_0_xxxxyz, g_0_xxxxzz, g_0_xxxyyy, g_0_xxxyyz, g_0_xxxyzz, g_0_xxxzzz, g_0_xxyyyy, g_0_xxyyyz, g_0_xxyyzz, g_0_xxyzzz, g_0_xxzzzz, g_0_xyyyyy, g_0_xyyyyz, g_0_xyyyzz, g_0_xyyzzz, g_0_xyzzzz, g_0_xzzzzz, g_0_yyyyyy, g_0_yyyyyz, g_0_yyyyzz, g_0_yyyzzz, g_0_yyzzzz, g_0_yzzzzz, g_0_zzzzzz, g_y_0_0_xxxxxx, g_y_0_0_xxxxxxy, g_y_0_0_xxxxxy, g_y_0_0_xxxxxyy, g_y_0_0_xxxxxyz, g_y_0_0_xxxxxz, g_y_0_0_xxxxyy, g_y_0_0_xxxxyyy, g_y_0_0_xxxxyyz, g_y_0_0_xxxxyz, g_y_0_0_xxxxyzz, g_y_0_0_xxxxzz, g_y_0_0_xxxyyy, g_y_0_0_xxxyyyy, g_y_0_0_xxxyyyz, g_y_0_0_xxxyyz, g_y_0_0_xxxyyzz, g_y_0_0_xxxyzz, g_y_0_0_xxxyzzz, g_y_0_0_xxxzzz, g_y_0_0_xxyyyy, g_y_0_0_xxyyyyy, g_y_0_0_xxyyyyz, g_y_0_0_xxyyyz, g_y_0_0_xxyyyzz, g_y_0_0_xxyyzz, g_y_0_0_xxyyzzz, g_y_0_0_xxyzzz, g_y_0_0_xxyzzzz, g_y_0_0_xxzzzz, g_y_0_0_xyyyyy, g_y_0_0_xyyyyyy, g_y_0_0_xyyyyyz, g_y_0_0_xyyyyz, g_y_0_0_xyyyyzz, g_y_0_0_xyyyzz, g_y_0_0_xyyyzzz, g_y_0_0_xyyzzz, g_y_0_0_xyyzzzz, g_y_0_0_xyzzzz, g_y_0_0_xyzzzzz, g_y_0_0_xzzzzz, g_y_0_0_yyyyyy, g_y_0_0_yyyyyyy, g_y_0_0_yyyyyyz, g_y_0_0_yyyyyz, g_y_0_0_yyyyyzz, g_y_0_0_yyyyzz, g_y_0_0_yyyyzzz, g_y_0_0_yyyzzz, g_y_0_0_yyyzzzz, g_y_0_0_yyzzzz, g_y_0_0_yyzzzzz, g_y_0_0_yzzzzz, g_y_0_0_yzzzzzz, g_y_0_0_zzzzzz, g_y_0_y_xxxxxx, g_y_0_y_xxxxxy, g_y_0_y_xxxxxz, g_y_0_y_xxxxyy, g_y_0_y_xxxxyz, g_y_0_y_xxxxzz, g_y_0_y_xxxyyy, g_y_0_y_xxxyyz, g_y_0_y_xxxyzz, g_y_0_y_xxxzzz, g_y_0_y_xxyyyy, g_y_0_y_xxyyyz, g_y_0_y_xxyyzz, g_y_0_y_xxyzzz, g_y_0_y_xxzzzz, g_y_0_y_xyyyyy, g_y_0_y_xyyyyz, g_y_0_y_xyyyzz, g_y_0_y_xyyzzz, g_y_0_y_xyzzzz, g_y_0_y_xzzzzz, g_y_0_y_yyyyyy, g_y_0_y_yyyyyz, g_y_0_y_yyyyzz, g_y_0_y_yyyzzz, g_y_0_y_yyzzzz, g_y_0_y_yzzzzz, g_y_0_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_xxxxxx[k] = -g_0_xxxxxx[k] - g_y_0_0_xxxxxx[k] * cd_y[k] + g_y_0_0_xxxxxxy[k];

                g_y_0_y_xxxxxy[k] = -g_0_xxxxxy[k] - g_y_0_0_xxxxxy[k] * cd_y[k] + g_y_0_0_xxxxxyy[k];

                g_y_0_y_xxxxxz[k] = -g_0_xxxxxz[k] - g_y_0_0_xxxxxz[k] * cd_y[k] + g_y_0_0_xxxxxyz[k];

                g_y_0_y_xxxxyy[k] = -g_0_xxxxyy[k] - g_y_0_0_xxxxyy[k] * cd_y[k] + g_y_0_0_xxxxyyy[k];

                g_y_0_y_xxxxyz[k] = -g_0_xxxxyz[k] - g_y_0_0_xxxxyz[k] * cd_y[k] + g_y_0_0_xxxxyyz[k];

                g_y_0_y_xxxxzz[k] = -g_0_xxxxzz[k] - g_y_0_0_xxxxzz[k] * cd_y[k] + g_y_0_0_xxxxyzz[k];

                g_y_0_y_xxxyyy[k] = -g_0_xxxyyy[k] - g_y_0_0_xxxyyy[k] * cd_y[k] + g_y_0_0_xxxyyyy[k];

                g_y_0_y_xxxyyz[k] = -g_0_xxxyyz[k] - g_y_0_0_xxxyyz[k] * cd_y[k] + g_y_0_0_xxxyyyz[k];

                g_y_0_y_xxxyzz[k] = -g_0_xxxyzz[k] - g_y_0_0_xxxyzz[k] * cd_y[k] + g_y_0_0_xxxyyzz[k];

                g_y_0_y_xxxzzz[k] = -g_0_xxxzzz[k] - g_y_0_0_xxxzzz[k] * cd_y[k] + g_y_0_0_xxxyzzz[k];

                g_y_0_y_xxyyyy[k] = -g_0_xxyyyy[k] - g_y_0_0_xxyyyy[k] * cd_y[k] + g_y_0_0_xxyyyyy[k];

                g_y_0_y_xxyyyz[k] = -g_0_xxyyyz[k] - g_y_0_0_xxyyyz[k] * cd_y[k] + g_y_0_0_xxyyyyz[k];

                g_y_0_y_xxyyzz[k] = -g_0_xxyyzz[k] - g_y_0_0_xxyyzz[k] * cd_y[k] + g_y_0_0_xxyyyzz[k];

                g_y_0_y_xxyzzz[k] = -g_0_xxyzzz[k] - g_y_0_0_xxyzzz[k] * cd_y[k] + g_y_0_0_xxyyzzz[k];

                g_y_0_y_xxzzzz[k] = -g_0_xxzzzz[k] - g_y_0_0_xxzzzz[k] * cd_y[k] + g_y_0_0_xxyzzzz[k];

                g_y_0_y_xyyyyy[k] = -g_0_xyyyyy[k] - g_y_0_0_xyyyyy[k] * cd_y[k] + g_y_0_0_xyyyyyy[k];

                g_y_0_y_xyyyyz[k] = -g_0_xyyyyz[k] - g_y_0_0_xyyyyz[k] * cd_y[k] + g_y_0_0_xyyyyyz[k];

                g_y_0_y_xyyyzz[k] = -g_0_xyyyzz[k] - g_y_0_0_xyyyzz[k] * cd_y[k] + g_y_0_0_xyyyyzz[k];

                g_y_0_y_xyyzzz[k] = -g_0_xyyzzz[k] - g_y_0_0_xyyzzz[k] * cd_y[k] + g_y_0_0_xyyyzzz[k];

                g_y_0_y_xyzzzz[k] = -g_0_xyzzzz[k] - g_y_0_0_xyzzzz[k] * cd_y[k] + g_y_0_0_xyyzzzz[k];

                g_y_0_y_xzzzzz[k] = -g_0_xzzzzz[k] - g_y_0_0_xzzzzz[k] * cd_y[k] + g_y_0_0_xyzzzzz[k];

                g_y_0_y_yyyyyy[k] = -g_0_yyyyyy[k] - g_y_0_0_yyyyyy[k] * cd_y[k] + g_y_0_0_yyyyyyy[k];

                g_y_0_y_yyyyyz[k] = -g_0_yyyyyz[k] - g_y_0_0_yyyyyz[k] * cd_y[k] + g_y_0_0_yyyyyyz[k];

                g_y_0_y_yyyyzz[k] = -g_0_yyyyzz[k] - g_y_0_0_yyyyzz[k] * cd_y[k] + g_y_0_0_yyyyyzz[k];

                g_y_0_y_yyyzzz[k] = -g_0_yyyzzz[k] - g_y_0_0_yyyzzz[k] * cd_y[k] + g_y_0_0_yyyyzzz[k];

                g_y_0_y_yyzzzz[k] = -g_0_yyzzzz[k] - g_y_0_0_yyzzzz[k] * cd_y[k] + g_y_0_0_yyyzzzz[k];

                g_y_0_y_yzzzzz[k] = -g_0_yzzzzz[k] - g_y_0_0_yzzzzz[k] * cd_y[k] + g_y_0_0_yyzzzzz[k];

                g_y_0_y_zzzzzz[k] = -g_0_zzzzzz[k] - g_y_0_0_zzzzzz[k] * cd_y[k] + g_y_0_0_yzzzzzz[k];
            }

            /// Set up 56-84 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_xxxxxx = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 56);

            auto g_y_0_z_xxxxxy = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 57);

            auto g_y_0_z_xxxxxz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 58);

            auto g_y_0_z_xxxxyy = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 59);

            auto g_y_0_z_xxxxyz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 60);

            auto g_y_0_z_xxxxzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 61);

            auto g_y_0_z_xxxyyy = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 62);

            auto g_y_0_z_xxxyyz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 63);

            auto g_y_0_z_xxxyzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 64);

            auto g_y_0_z_xxxzzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 65);

            auto g_y_0_z_xxyyyy = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 66);

            auto g_y_0_z_xxyyyz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 67);

            auto g_y_0_z_xxyyzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 68);

            auto g_y_0_z_xxyzzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 69);

            auto g_y_0_z_xxzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 70);

            auto g_y_0_z_xyyyyy = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 71);

            auto g_y_0_z_xyyyyz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 72);

            auto g_y_0_z_xyyyzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 73);

            auto g_y_0_z_xyyzzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 74);

            auto g_y_0_z_xyzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 75);

            auto g_y_0_z_xzzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 76);

            auto g_y_0_z_yyyyyy = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 77);

            auto g_y_0_z_yyyyyz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 78);

            auto g_y_0_z_yyyyzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 79);

            auto g_y_0_z_yyyzzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 80);

            auto g_y_0_z_yyzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 81);

            auto g_y_0_z_yzzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 82);

            auto g_y_0_z_zzzzzz = cbuffer.data(pi_geom_10_off + 84 * acomps * bcomps + 83);

            #pragma omp simd aligned(cd_z, g_y_0_0_xxxxxx, g_y_0_0_xxxxxxz, g_y_0_0_xxxxxy, g_y_0_0_xxxxxyz, g_y_0_0_xxxxxz, g_y_0_0_xxxxxzz, g_y_0_0_xxxxyy, g_y_0_0_xxxxyyz, g_y_0_0_xxxxyz, g_y_0_0_xxxxyzz, g_y_0_0_xxxxzz, g_y_0_0_xxxxzzz, g_y_0_0_xxxyyy, g_y_0_0_xxxyyyz, g_y_0_0_xxxyyz, g_y_0_0_xxxyyzz, g_y_0_0_xxxyzz, g_y_0_0_xxxyzzz, g_y_0_0_xxxzzz, g_y_0_0_xxxzzzz, g_y_0_0_xxyyyy, g_y_0_0_xxyyyyz, g_y_0_0_xxyyyz, g_y_0_0_xxyyyzz, g_y_0_0_xxyyzz, g_y_0_0_xxyyzzz, g_y_0_0_xxyzzz, g_y_0_0_xxyzzzz, g_y_0_0_xxzzzz, g_y_0_0_xxzzzzz, g_y_0_0_xyyyyy, g_y_0_0_xyyyyyz, g_y_0_0_xyyyyz, g_y_0_0_xyyyyzz, g_y_0_0_xyyyzz, g_y_0_0_xyyyzzz, g_y_0_0_xyyzzz, g_y_0_0_xyyzzzz, g_y_0_0_xyzzzz, g_y_0_0_xyzzzzz, g_y_0_0_xzzzzz, g_y_0_0_xzzzzzz, g_y_0_0_yyyyyy, g_y_0_0_yyyyyyz, g_y_0_0_yyyyyz, g_y_0_0_yyyyyzz, g_y_0_0_yyyyzz, g_y_0_0_yyyyzzz, g_y_0_0_yyyzzz, g_y_0_0_yyyzzzz, g_y_0_0_yyzzzz, g_y_0_0_yyzzzzz, g_y_0_0_yzzzzz, g_y_0_0_yzzzzzz, g_y_0_0_zzzzzz, g_y_0_0_zzzzzzz, g_y_0_z_xxxxxx, g_y_0_z_xxxxxy, g_y_0_z_xxxxxz, g_y_0_z_xxxxyy, g_y_0_z_xxxxyz, g_y_0_z_xxxxzz, g_y_0_z_xxxyyy, g_y_0_z_xxxyyz, g_y_0_z_xxxyzz, g_y_0_z_xxxzzz, g_y_0_z_xxyyyy, g_y_0_z_xxyyyz, g_y_0_z_xxyyzz, g_y_0_z_xxyzzz, g_y_0_z_xxzzzz, g_y_0_z_xyyyyy, g_y_0_z_xyyyyz, g_y_0_z_xyyyzz, g_y_0_z_xyyzzz, g_y_0_z_xyzzzz, g_y_0_z_xzzzzz, g_y_0_z_yyyyyy, g_y_0_z_yyyyyz, g_y_0_z_yyyyzz, g_y_0_z_yyyzzz, g_y_0_z_yyzzzz, g_y_0_z_yzzzzz, g_y_0_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_xxxxxx[k] = -g_y_0_0_xxxxxx[k] * cd_z[k] + g_y_0_0_xxxxxxz[k];

                g_y_0_z_xxxxxy[k] = -g_y_0_0_xxxxxy[k] * cd_z[k] + g_y_0_0_xxxxxyz[k];

                g_y_0_z_xxxxxz[k] = -g_y_0_0_xxxxxz[k] * cd_z[k] + g_y_0_0_xxxxxzz[k];

                g_y_0_z_xxxxyy[k] = -g_y_0_0_xxxxyy[k] * cd_z[k] + g_y_0_0_xxxxyyz[k];

                g_y_0_z_xxxxyz[k] = -g_y_0_0_xxxxyz[k] * cd_z[k] + g_y_0_0_xxxxyzz[k];

                g_y_0_z_xxxxzz[k] = -g_y_0_0_xxxxzz[k] * cd_z[k] + g_y_0_0_xxxxzzz[k];

                g_y_0_z_xxxyyy[k] = -g_y_0_0_xxxyyy[k] * cd_z[k] + g_y_0_0_xxxyyyz[k];

                g_y_0_z_xxxyyz[k] = -g_y_0_0_xxxyyz[k] * cd_z[k] + g_y_0_0_xxxyyzz[k];

                g_y_0_z_xxxyzz[k] = -g_y_0_0_xxxyzz[k] * cd_z[k] + g_y_0_0_xxxyzzz[k];

                g_y_0_z_xxxzzz[k] = -g_y_0_0_xxxzzz[k] * cd_z[k] + g_y_0_0_xxxzzzz[k];

                g_y_0_z_xxyyyy[k] = -g_y_0_0_xxyyyy[k] * cd_z[k] + g_y_0_0_xxyyyyz[k];

                g_y_0_z_xxyyyz[k] = -g_y_0_0_xxyyyz[k] * cd_z[k] + g_y_0_0_xxyyyzz[k];

                g_y_0_z_xxyyzz[k] = -g_y_0_0_xxyyzz[k] * cd_z[k] + g_y_0_0_xxyyzzz[k];

                g_y_0_z_xxyzzz[k] = -g_y_0_0_xxyzzz[k] * cd_z[k] + g_y_0_0_xxyzzzz[k];

                g_y_0_z_xxzzzz[k] = -g_y_0_0_xxzzzz[k] * cd_z[k] + g_y_0_0_xxzzzzz[k];

                g_y_0_z_xyyyyy[k] = -g_y_0_0_xyyyyy[k] * cd_z[k] + g_y_0_0_xyyyyyz[k];

                g_y_0_z_xyyyyz[k] = -g_y_0_0_xyyyyz[k] * cd_z[k] + g_y_0_0_xyyyyzz[k];

                g_y_0_z_xyyyzz[k] = -g_y_0_0_xyyyzz[k] * cd_z[k] + g_y_0_0_xyyyzzz[k];

                g_y_0_z_xyyzzz[k] = -g_y_0_0_xyyzzz[k] * cd_z[k] + g_y_0_0_xyyzzzz[k];

                g_y_0_z_xyzzzz[k] = -g_y_0_0_xyzzzz[k] * cd_z[k] + g_y_0_0_xyzzzzz[k];

                g_y_0_z_xzzzzz[k] = -g_y_0_0_xzzzzz[k] * cd_z[k] + g_y_0_0_xzzzzzz[k];

                g_y_0_z_yyyyyy[k] = -g_y_0_0_yyyyyy[k] * cd_z[k] + g_y_0_0_yyyyyyz[k];

                g_y_0_z_yyyyyz[k] = -g_y_0_0_yyyyyz[k] * cd_z[k] + g_y_0_0_yyyyyzz[k];

                g_y_0_z_yyyyzz[k] = -g_y_0_0_yyyyzz[k] * cd_z[k] + g_y_0_0_yyyyzzz[k];

                g_y_0_z_yyyzzz[k] = -g_y_0_0_yyyzzz[k] * cd_z[k] + g_y_0_0_yyyzzzz[k];

                g_y_0_z_yyzzzz[k] = -g_y_0_0_yyzzzz[k] * cd_z[k] + g_y_0_0_yyzzzzz[k];

                g_y_0_z_yzzzzz[k] = -g_y_0_0_yzzzzz[k] * cd_z[k] + g_y_0_0_yzzzzzz[k];

                g_y_0_z_zzzzzz[k] = -g_y_0_0_zzzzzz[k] * cd_z[k] + g_y_0_0_zzzzzzz[k];
            }
            /// Set up 0-28 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_xxxxxx = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 0);

            auto g_z_0_x_xxxxxy = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 1);

            auto g_z_0_x_xxxxxz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 2);

            auto g_z_0_x_xxxxyy = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 3);

            auto g_z_0_x_xxxxyz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 4);

            auto g_z_0_x_xxxxzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 5);

            auto g_z_0_x_xxxyyy = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 6);

            auto g_z_0_x_xxxyyz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 7);

            auto g_z_0_x_xxxyzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 8);

            auto g_z_0_x_xxxzzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 9);

            auto g_z_0_x_xxyyyy = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 10);

            auto g_z_0_x_xxyyyz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 11);

            auto g_z_0_x_xxyyzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 12);

            auto g_z_0_x_xxyzzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 13);

            auto g_z_0_x_xxzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 14);

            auto g_z_0_x_xyyyyy = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 15);

            auto g_z_0_x_xyyyyz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 16);

            auto g_z_0_x_xyyyzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 17);

            auto g_z_0_x_xyyzzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 18);

            auto g_z_0_x_xyzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 19);

            auto g_z_0_x_xzzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 20);

            auto g_z_0_x_yyyyyy = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 21);

            auto g_z_0_x_yyyyyz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 22);

            auto g_z_0_x_yyyyzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 23);

            auto g_z_0_x_yyyzzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 24);

            auto g_z_0_x_yyzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 25);

            auto g_z_0_x_yzzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 26);

            auto g_z_0_x_zzzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 27);

            #pragma omp simd aligned(cd_x, g_z_0_0_xxxxxx, g_z_0_0_xxxxxxx, g_z_0_0_xxxxxxy, g_z_0_0_xxxxxxz, g_z_0_0_xxxxxy, g_z_0_0_xxxxxyy, g_z_0_0_xxxxxyz, g_z_0_0_xxxxxz, g_z_0_0_xxxxxzz, g_z_0_0_xxxxyy, g_z_0_0_xxxxyyy, g_z_0_0_xxxxyyz, g_z_0_0_xxxxyz, g_z_0_0_xxxxyzz, g_z_0_0_xxxxzz, g_z_0_0_xxxxzzz, g_z_0_0_xxxyyy, g_z_0_0_xxxyyyy, g_z_0_0_xxxyyyz, g_z_0_0_xxxyyz, g_z_0_0_xxxyyzz, g_z_0_0_xxxyzz, g_z_0_0_xxxyzzz, g_z_0_0_xxxzzz, g_z_0_0_xxxzzzz, g_z_0_0_xxyyyy, g_z_0_0_xxyyyyy, g_z_0_0_xxyyyyz, g_z_0_0_xxyyyz, g_z_0_0_xxyyyzz, g_z_0_0_xxyyzz, g_z_0_0_xxyyzzz, g_z_0_0_xxyzzz, g_z_0_0_xxyzzzz, g_z_0_0_xxzzzz, g_z_0_0_xxzzzzz, g_z_0_0_xyyyyy, g_z_0_0_xyyyyyy, g_z_0_0_xyyyyyz, g_z_0_0_xyyyyz, g_z_0_0_xyyyyzz, g_z_0_0_xyyyzz, g_z_0_0_xyyyzzz, g_z_0_0_xyyzzz, g_z_0_0_xyyzzzz, g_z_0_0_xyzzzz, g_z_0_0_xyzzzzz, g_z_0_0_xzzzzz, g_z_0_0_xzzzzzz, g_z_0_0_yyyyyy, g_z_0_0_yyyyyz, g_z_0_0_yyyyzz, g_z_0_0_yyyzzz, g_z_0_0_yyzzzz, g_z_0_0_yzzzzz, g_z_0_0_zzzzzz, g_z_0_x_xxxxxx, g_z_0_x_xxxxxy, g_z_0_x_xxxxxz, g_z_0_x_xxxxyy, g_z_0_x_xxxxyz, g_z_0_x_xxxxzz, g_z_0_x_xxxyyy, g_z_0_x_xxxyyz, g_z_0_x_xxxyzz, g_z_0_x_xxxzzz, g_z_0_x_xxyyyy, g_z_0_x_xxyyyz, g_z_0_x_xxyyzz, g_z_0_x_xxyzzz, g_z_0_x_xxzzzz, g_z_0_x_xyyyyy, g_z_0_x_xyyyyz, g_z_0_x_xyyyzz, g_z_0_x_xyyzzz, g_z_0_x_xyzzzz, g_z_0_x_xzzzzz, g_z_0_x_yyyyyy, g_z_0_x_yyyyyz, g_z_0_x_yyyyzz, g_z_0_x_yyyzzz, g_z_0_x_yyzzzz, g_z_0_x_yzzzzz, g_z_0_x_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_xxxxxx[k] = -g_z_0_0_xxxxxx[k] * cd_x[k] + g_z_0_0_xxxxxxx[k];

                g_z_0_x_xxxxxy[k] = -g_z_0_0_xxxxxy[k] * cd_x[k] + g_z_0_0_xxxxxxy[k];

                g_z_0_x_xxxxxz[k] = -g_z_0_0_xxxxxz[k] * cd_x[k] + g_z_0_0_xxxxxxz[k];

                g_z_0_x_xxxxyy[k] = -g_z_0_0_xxxxyy[k] * cd_x[k] + g_z_0_0_xxxxxyy[k];

                g_z_0_x_xxxxyz[k] = -g_z_0_0_xxxxyz[k] * cd_x[k] + g_z_0_0_xxxxxyz[k];

                g_z_0_x_xxxxzz[k] = -g_z_0_0_xxxxzz[k] * cd_x[k] + g_z_0_0_xxxxxzz[k];

                g_z_0_x_xxxyyy[k] = -g_z_0_0_xxxyyy[k] * cd_x[k] + g_z_0_0_xxxxyyy[k];

                g_z_0_x_xxxyyz[k] = -g_z_0_0_xxxyyz[k] * cd_x[k] + g_z_0_0_xxxxyyz[k];

                g_z_0_x_xxxyzz[k] = -g_z_0_0_xxxyzz[k] * cd_x[k] + g_z_0_0_xxxxyzz[k];

                g_z_0_x_xxxzzz[k] = -g_z_0_0_xxxzzz[k] * cd_x[k] + g_z_0_0_xxxxzzz[k];

                g_z_0_x_xxyyyy[k] = -g_z_0_0_xxyyyy[k] * cd_x[k] + g_z_0_0_xxxyyyy[k];

                g_z_0_x_xxyyyz[k] = -g_z_0_0_xxyyyz[k] * cd_x[k] + g_z_0_0_xxxyyyz[k];

                g_z_0_x_xxyyzz[k] = -g_z_0_0_xxyyzz[k] * cd_x[k] + g_z_0_0_xxxyyzz[k];

                g_z_0_x_xxyzzz[k] = -g_z_0_0_xxyzzz[k] * cd_x[k] + g_z_0_0_xxxyzzz[k];

                g_z_0_x_xxzzzz[k] = -g_z_0_0_xxzzzz[k] * cd_x[k] + g_z_0_0_xxxzzzz[k];

                g_z_0_x_xyyyyy[k] = -g_z_0_0_xyyyyy[k] * cd_x[k] + g_z_0_0_xxyyyyy[k];

                g_z_0_x_xyyyyz[k] = -g_z_0_0_xyyyyz[k] * cd_x[k] + g_z_0_0_xxyyyyz[k];

                g_z_0_x_xyyyzz[k] = -g_z_0_0_xyyyzz[k] * cd_x[k] + g_z_0_0_xxyyyzz[k];

                g_z_0_x_xyyzzz[k] = -g_z_0_0_xyyzzz[k] * cd_x[k] + g_z_0_0_xxyyzzz[k];

                g_z_0_x_xyzzzz[k] = -g_z_0_0_xyzzzz[k] * cd_x[k] + g_z_0_0_xxyzzzz[k];

                g_z_0_x_xzzzzz[k] = -g_z_0_0_xzzzzz[k] * cd_x[k] + g_z_0_0_xxzzzzz[k];

                g_z_0_x_yyyyyy[k] = -g_z_0_0_yyyyyy[k] * cd_x[k] + g_z_0_0_xyyyyyy[k];

                g_z_0_x_yyyyyz[k] = -g_z_0_0_yyyyyz[k] * cd_x[k] + g_z_0_0_xyyyyyz[k];

                g_z_0_x_yyyyzz[k] = -g_z_0_0_yyyyzz[k] * cd_x[k] + g_z_0_0_xyyyyzz[k];

                g_z_0_x_yyyzzz[k] = -g_z_0_0_yyyzzz[k] * cd_x[k] + g_z_0_0_xyyyzzz[k];

                g_z_0_x_yyzzzz[k] = -g_z_0_0_yyzzzz[k] * cd_x[k] + g_z_0_0_xyyzzzz[k];

                g_z_0_x_yzzzzz[k] = -g_z_0_0_yzzzzz[k] * cd_x[k] + g_z_0_0_xyzzzzz[k];

                g_z_0_x_zzzzzz[k] = -g_z_0_0_zzzzzz[k] * cd_x[k] + g_z_0_0_xzzzzzz[k];
            }

            /// Set up 28-56 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_xxxxxx = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 28);

            auto g_z_0_y_xxxxxy = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 29);

            auto g_z_0_y_xxxxxz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 30);

            auto g_z_0_y_xxxxyy = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 31);

            auto g_z_0_y_xxxxyz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 32);

            auto g_z_0_y_xxxxzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 33);

            auto g_z_0_y_xxxyyy = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 34);

            auto g_z_0_y_xxxyyz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 35);

            auto g_z_0_y_xxxyzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 36);

            auto g_z_0_y_xxxzzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 37);

            auto g_z_0_y_xxyyyy = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 38);

            auto g_z_0_y_xxyyyz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 39);

            auto g_z_0_y_xxyyzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 40);

            auto g_z_0_y_xxyzzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 41);

            auto g_z_0_y_xxzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 42);

            auto g_z_0_y_xyyyyy = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 43);

            auto g_z_0_y_xyyyyz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 44);

            auto g_z_0_y_xyyyzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 45);

            auto g_z_0_y_xyyzzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 46);

            auto g_z_0_y_xyzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 47);

            auto g_z_0_y_xzzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 48);

            auto g_z_0_y_yyyyyy = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 49);

            auto g_z_0_y_yyyyyz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 50);

            auto g_z_0_y_yyyyzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 51);

            auto g_z_0_y_yyyzzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 52);

            auto g_z_0_y_yyzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 53);

            auto g_z_0_y_yzzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 54);

            auto g_z_0_y_zzzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 55);

            #pragma omp simd aligned(cd_y, g_z_0_0_xxxxxx, g_z_0_0_xxxxxxy, g_z_0_0_xxxxxy, g_z_0_0_xxxxxyy, g_z_0_0_xxxxxyz, g_z_0_0_xxxxxz, g_z_0_0_xxxxyy, g_z_0_0_xxxxyyy, g_z_0_0_xxxxyyz, g_z_0_0_xxxxyz, g_z_0_0_xxxxyzz, g_z_0_0_xxxxzz, g_z_0_0_xxxyyy, g_z_0_0_xxxyyyy, g_z_0_0_xxxyyyz, g_z_0_0_xxxyyz, g_z_0_0_xxxyyzz, g_z_0_0_xxxyzz, g_z_0_0_xxxyzzz, g_z_0_0_xxxzzz, g_z_0_0_xxyyyy, g_z_0_0_xxyyyyy, g_z_0_0_xxyyyyz, g_z_0_0_xxyyyz, g_z_0_0_xxyyyzz, g_z_0_0_xxyyzz, g_z_0_0_xxyyzzz, g_z_0_0_xxyzzz, g_z_0_0_xxyzzzz, g_z_0_0_xxzzzz, g_z_0_0_xyyyyy, g_z_0_0_xyyyyyy, g_z_0_0_xyyyyyz, g_z_0_0_xyyyyz, g_z_0_0_xyyyyzz, g_z_0_0_xyyyzz, g_z_0_0_xyyyzzz, g_z_0_0_xyyzzz, g_z_0_0_xyyzzzz, g_z_0_0_xyzzzz, g_z_0_0_xyzzzzz, g_z_0_0_xzzzzz, g_z_0_0_yyyyyy, g_z_0_0_yyyyyyy, g_z_0_0_yyyyyyz, g_z_0_0_yyyyyz, g_z_0_0_yyyyyzz, g_z_0_0_yyyyzz, g_z_0_0_yyyyzzz, g_z_0_0_yyyzzz, g_z_0_0_yyyzzzz, g_z_0_0_yyzzzz, g_z_0_0_yyzzzzz, g_z_0_0_yzzzzz, g_z_0_0_yzzzzzz, g_z_0_0_zzzzzz, g_z_0_y_xxxxxx, g_z_0_y_xxxxxy, g_z_0_y_xxxxxz, g_z_0_y_xxxxyy, g_z_0_y_xxxxyz, g_z_0_y_xxxxzz, g_z_0_y_xxxyyy, g_z_0_y_xxxyyz, g_z_0_y_xxxyzz, g_z_0_y_xxxzzz, g_z_0_y_xxyyyy, g_z_0_y_xxyyyz, g_z_0_y_xxyyzz, g_z_0_y_xxyzzz, g_z_0_y_xxzzzz, g_z_0_y_xyyyyy, g_z_0_y_xyyyyz, g_z_0_y_xyyyzz, g_z_0_y_xyyzzz, g_z_0_y_xyzzzz, g_z_0_y_xzzzzz, g_z_0_y_yyyyyy, g_z_0_y_yyyyyz, g_z_0_y_yyyyzz, g_z_0_y_yyyzzz, g_z_0_y_yyzzzz, g_z_0_y_yzzzzz, g_z_0_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_xxxxxx[k] = -g_z_0_0_xxxxxx[k] * cd_y[k] + g_z_0_0_xxxxxxy[k];

                g_z_0_y_xxxxxy[k] = -g_z_0_0_xxxxxy[k] * cd_y[k] + g_z_0_0_xxxxxyy[k];

                g_z_0_y_xxxxxz[k] = -g_z_0_0_xxxxxz[k] * cd_y[k] + g_z_0_0_xxxxxyz[k];

                g_z_0_y_xxxxyy[k] = -g_z_0_0_xxxxyy[k] * cd_y[k] + g_z_0_0_xxxxyyy[k];

                g_z_0_y_xxxxyz[k] = -g_z_0_0_xxxxyz[k] * cd_y[k] + g_z_0_0_xxxxyyz[k];

                g_z_0_y_xxxxzz[k] = -g_z_0_0_xxxxzz[k] * cd_y[k] + g_z_0_0_xxxxyzz[k];

                g_z_0_y_xxxyyy[k] = -g_z_0_0_xxxyyy[k] * cd_y[k] + g_z_0_0_xxxyyyy[k];

                g_z_0_y_xxxyyz[k] = -g_z_0_0_xxxyyz[k] * cd_y[k] + g_z_0_0_xxxyyyz[k];

                g_z_0_y_xxxyzz[k] = -g_z_0_0_xxxyzz[k] * cd_y[k] + g_z_0_0_xxxyyzz[k];

                g_z_0_y_xxxzzz[k] = -g_z_0_0_xxxzzz[k] * cd_y[k] + g_z_0_0_xxxyzzz[k];

                g_z_0_y_xxyyyy[k] = -g_z_0_0_xxyyyy[k] * cd_y[k] + g_z_0_0_xxyyyyy[k];

                g_z_0_y_xxyyyz[k] = -g_z_0_0_xxyyyz[k] * cd_y[k] + g_z_0_0_xxyyyyz[k];

                g_z_0_y_xxyyzz[k] = -g_z_0_0_xxyyzz[k] * cd_y[k] + g_z_0_0_xxyyyzz[k];

                g_z_0_y_xxyzzz[k] = -g_z_0_0_xxyzzz[k] * cd_y[k] + g_z_0_0_xxyyzzz[k];

                g_z_0_y_xxzzzz[k] = -g_z_0_0_xxzzzz[k] * cd_y[k] + g_z_0_0_xxyzzzz[k];

                g_z_0_y_xyyyyy[k] = -g_z_0_0_xyyyyy[k] * cd_y[k] + g_z_0_0_xyyyyyy[k];

                g_z_0_y_xyyyyz[k] = -g_z_0_0_xyyyyz[k] * cd_y[k] + g_z_0_0_xyyyyyz[k];

                g_z_0_y_xyyyzz[k] = -g_z_0_0_xyyyzz[k] * cd_y[k] + g_z_0_0_xyyyyzz[k];

                g_z_0_y_xyyzzz[k] = -g_z_0_0_xyyzzz[k] * cd_y[k] + g_z_0_0_xyyyzzz[k];

                g_z_0_y_xyzzzz[k] = -g_z_0_0_xyzzzz[k] * cd_y[k] + g_z_0_0_xyyzzzz[k];

                g_z_0_y_xzzzzz[k] = -g_z_0_0_xzzzzz[k] * cd_y[k] + g_z_0_0_xyzzzzz[k];

                g_z_0_y_yyyyyy[k] = -g_z_0_0_yyyyyy[k] * cd_y[k] + g_z_0_0_yyyyyyy[k];

                g_z_0_y_yyyyyz[k] = -g_z_0_0_yyyyyz[k] * cd_y[k] + g_z_0_0_yyyyyyz[k];

                g_z_0_y_yyyyzz[k] = -g_z_0_0_yyyyzz[k] * cd_y[k] + g_z_0_0_yyyyyzz[k];

                g_z_0_y_yyyzzz[k] = -g_z_0_0_yyyzzz[k] * cd_y[k] + g_z_0_0_yyyyzzz[k];

                g_z_0_y_yyzzzz[k] = -g_z_0_0_yyzzzz[k] * cd_y[k] + g_z_0_0_yyyzzzz[k];

                g_z_0_y_yzzzzz[k] = -g_z_0_0_yzzzzz[k] * cd_y[k] + g_z_0_0_yyzzzzz[k];

                g_z_0_y_zzzzzz[k] = -g_z_0_0_zzzzzz[k] * cd_y[k] + g_z_0_0_yzzzzzz[k];
            }

            /// Set up 56-84 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_xxxxxx = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 56);

            auto g_z_0_z_xxxxxy = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 57);

            auto g_z_0_z_xxxxxz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 58);

            auto g_z_0_z_xxxxyy = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 59);

            auto g_z_0_z_xxxxyz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 60);

            auto g_z_0_z_xxxxzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 61);

            auto g_z_0_z_xxxyyy = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 62);

            auto g_z_0_z_xxxyyz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 63);

            auto g_z_0_z_xxxyzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 64);

            auto g_z_0_z_xxxzzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 65);

            auto g_z_0_z_xxyyyy = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 66);

            auto g_z_0_z_xxyyyz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 67);

            auto g_z_0_z_xxyyzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 68);

            auto g_z_0_z_xxyzzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 69);

            auto g_z_0_z_xxzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 70);

            auto g_z_0_z_xyyyyy = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 71);

            auto g_z_0_z_xyyyyz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 72);

            auto g_z_0_z_xyyyzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 73);

            auto g_z_0_z_xyyzzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 74);

            auto g_z_0_z_xyzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 75);

            auto g_z_0_z_xzzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 76);

            auto g_z_0_z_yyyyyy = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 77);

            auto g_z_0_z_yyyyyz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 78);

            auto g_z_0_z_yyyyzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 79);

            auto g_z_0_z_yyyzzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 80);

            auto g_z_0_z_yyzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 81);

            auto g_z_0_z_yzzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 82);

            auto g_z_0_z_zzzzzz = cbuffer.data(pi_geom_10_off + 168 * acomps * bcomps + 83);

            #pragma omp simd aligned(cd_z, g_0_xxxxxx, g_0_xxxxxy, g_0_xxxxxz, g_0_xxxxyy, g_0_xxxxyz, g_0_xxxxzz, g_0_xxxyyy, g_0_xxxyyz, g_0_xxxyzz, g_0_xxxzzz, g_0_xxyyyy, g_0_xxyyyz, g_0_xxyyzz, g_0_xxyzzz, g_0_xxzzzz, g_0_xyyyyy, g_0_xyyyyz, g_0_xyyyzz, g_0_xyyzzz, g_0_xyzzzz, g_0_xzzzzz, g_0_yyyyyy, g_0_yyyyyz, g_0_yyyyzz, g_0_yyyzzz, g_0_yyzzzz, g_0_yzzzzz, g_0_zzzzzz, g_z_0_0_xxxxxx, g_z_0_0_xxxxxxz, g_z_0_0_xxxxxy, g_z_0_0_xxxxxyz, g_z_0_0_xxxxxz, g_z_0_0_xxxxxzz, g_z_0_0_xxxxyy, g_z_0_0_xxxxyyz, g_z_0_0_xxxxyz, g_z_0_0_xxxxyzz, g_z_0_0_xxxxzz, g_z_0_0_xxxxzzz, g_z_0_0_xxxyyy, g_z_0_0_xxxyyyz, g_z_0_0_xxxyyz, g_z_0_0_xxxyyzz, g_z_0_0_xxxyzz, g_z_0_0_xxxyzzz, g_z_0_0_xxxzzz, g_z_0_0_xxxzzzz, g_z_0_0_xxyyyy, g_z_0_0_xxyyyyz, g_z_0_0_xxyyyz, g_z_0_0_xxyyyzz, g_z_0_0_xxyyzz, g_z_0_0_xxyyzzz, g_z_0_0_xxyzzz, g_z_0_0_xxyzzzz, g_z_0_0_xxzzzz, g_z_0_0_xxzzzzz, g_z_0_0_xyyyyy, g_z_0_0_xyyyyyz, g_z_0_0_xyyyyz, g_z_0_0_xyyyyzz, g_z_0_0_xyyyzz, g_z_0_0_xyyyzzz, g_z_0_0_xyyzzz, g_z_0_0_xyyzzzz, g_z_0_0_xyzzzz, g_z_0_0_xyzzzzz, g_z_0_0_xzzzzz, g_z_0_0_xzzzzzz, g_z_0_0_yyyyyy, g_z_0_0_yyyyyyz, g_z_0_0_yyyyyz, g_z_0_0_yyyyyzz, g_z_0_0_yyyyzz, g_z_0_0_yyyyzzz, g_z_0_0_yyyzzz, g_z_0_0_yyyzzzz, g_z_0_0_yyzzzz, g_z_0_0_yyzzzzz, g_z_0_0_yzzzzz, g_z_0_0_yzzzzzz, g_z_0_0_zzzzzz, g_z_0_0_zzzzzzz, g_z_0_z_xxxxxx, g_z_0_z_xxxxxy, g_z_0_z_xxxxxz, g_z_0_z_xxxxyy, g_z_0_z_xxxxyz, g_z_0_z_xxxxzz, g_z_0_z_xxxyyy, g_z_0_z_xxxyyz, g_z_0_z_xxxyzz, g_z_0_z_xxxzzz, g_z_0_z_xxyyyy, g_z_0_z_xxyyyz, g_z_0_z_xxyyzz, g_z_0_z_xxyzzz, g_z_0_z_xxzzzz, g_z_0_z_xyyyyy, g_z_0_z_xyyyyz, g_z_0_z_xyyyzz, g_z_0_z_xyyzzz, g_z_0_z_xyzzzz, g_z_0_z_xzzzzz, g_z_0_z_yyyyyy, g_z_0_z_yyyyyz, g_z_0_z_yyyyzz, g_z_0_z_yyyzzz, g_z_0_z_yyzzzz, g_z_0_z_yzzzzz, g_z_0_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_xxxxxx[k] = -g_0_xxxxxx[k] - g_z_0_0_xxxxxx[k] * cd_z[k] + g_z_0_0_xxxxxxz[k];

                g_z_0_z_xxxxxy[k] = -g_0_xxxxxy[k] - g_z_0_0_xxxxxy[k] * cd_z[k] + g_z_0_0_xxxxxyz[k];

                g_z_0_z_xxxxxz[k] = -g_0_xxxxxz[k] - g_z_0_0_xxxxxz[k] * cd_z[k] + g_z_0_0_xxxxxzz[k];

                g_z_0_z_xxxxyy[k] = -g_0_xxxxyy[k] - g_z_0_0_xxxxyy[k] * cd_z[k] + g_z_0_0_xxxxyyz[k];

                g_z_0_z_xxxxyz[k] = -g_0_xxxxyz[k] - g_z_0_0_xxxxyz[k] * cd_z[k] + g_z_0_0_xxxxyzz[k];

                g_z_0_z_xxxxzz[k] = -g_0_xxxxzz[k] - g_z_0_0_xxxxzz[k] * cd_z[k] + g_z_0_0_xxxxzzz[k];

                g_z_0_z_xxxyyy[k] = -g_0_xxxyyy[k] - g_z_0_0_xxxyyy[k] * cd_z[k] + g_z_0_0_xxxyyyz[k];

                g_z_0_z_xxxyyz[k] = -g_0_xxxyyz[k] - g_z_0_0_xxxyyz[k] * cd_z[k] + g_z_0_0_xxxyyzz[k];

                g_z_0_z_xxxyzz[k] = -g_0_xxxyzz[k] - g_z_0_0_xxxyzz[k] * cd_z[k] + g_z_0_0_xxxyzzz[k];

                g_z_0_z_xxxzzz[k] = -g_0_xxxzzz[k] - g_z_0_0_xxxzzz[k] * cd_z[k] + g_z_0_0_xxxzzzz[k];

                g_z_0_z_xxyyyy[k] = -g_0_xxyyyy[k] - g_z_0_0_xxyyyy[k] * cd_z[k] + g_z_0_0_xxyyyyz[k];

                g_z_0_z_xxyyyz[k] = -g_0_xxyyyz[k] - g_z_0_0_xxyyyz[k] * cd_z[k] + g_z_0_0_xxyyyzz[k];

                g_z_0_z_xxyyzz[k] = -g_0_xxyyzz[k] - g_z_0_0_xxyyzz[k] * cd_z[k] + g_z_0_0_xxyyzzz[k];

                g_z_0_z_xxyzzz[k] = -g_0_xxyzzz[k] - g_z_0_0_xxyzzz[k] * cd_z[k] + g_z_0_0_xxyzzzz[k];

                g_z_0_z_xxzzzz[k] = -g_0_xxzzzz[k] - g_z_0_0_xxzzzz[k] * cd_z[k] + g_z_0_0_xxzzzzz[k];

                g_z_0_z_xyyyyy[k] = -g_0_xyyyyy[k] - g_z_0_0_xyyyyy[k] * cd_z[k] + g_z_0_0_xyyyyyz[k];

                g_z_0_z_xyyyyz[k] = -g_0_xyyyyz[k] - g_z_0_0_xyyyyz[k] * cd_z[k] + g_z_0_0_xyyyyzz[k];

                g_z_0_z_xyyyzz[k] = -g_0_xyyyzz[k] - g_z_0_0_xyyyzz[k] * cd_z[k] + g_z_0_0_xyyyzzz[k];

                g_z_0_z_xyyzzz[k] = -g_0_xyyzzz[k] - g_z_0_0_xyyzzz[k] * cd_z[k] + g_z_0_0_xyyzzzz[k];

                g_z_0_z_xyzzzz[k] = -g_0_xyzzzz[k] - g_z_0_0_xyzzzz[k] * cd_z[k] + g_z_0_0_xyzzzzz[k];

                g_z_0_z_xzzzzz[k] = -g_0_xzzzzz[k] - g_z_0_0_xzzzzz[k] * cd_z[k] + g_z_0_0_xzzzzzz[k];

                g_z_0_z_yyyyyy[k] = -g_0_yyyyyy[k] - g_z_0_0_yyyyyy[k] * cd_z[k] + g_z_0_0_yyyyyyz[k];

                g_z_0_z_yyyyyz[k] = -g_0_yyyyyz[k] - g_z_0_0_yyyyyz[k] * cd_z[k] + g_z_0_0_yyyyyzz[k];

                g_z_0_z_yyyyzz[k] = -g_0_yyyyzz[k] - g_z_0_0_yyyyzz[k] * cd_z[k] + g_z_0_0_yyyyzzz[k];

                g_z_0_z_yyyzzz[k] = -g_0_yyyzzz[k] - g_z_0_0_yyyzzz[k] * cd_z[k] + g_z_0_0_yyyzzzz[k];

                g_z_0_z_yyzzzz[k] = -g_0_yyzzzz[k] - g_z_0_0_yyzzzz[k] * cd_z[k] + g_z_0_0_yyzzzzz[k];

                g_z_0_z_yzzzzz[k] = -g_0_yzzzzz[k] - g_z_0_0_yzzzzz[k] * cd_z[k] + g_z_0_0_yzzzzzz[k];

                g_z_0_z_zzzzzz[k] = -g_0_zzzzzz[k] - g_z_0_0_zzzzzz[k] * cd_z[k] + g_z_0_0_zzzzzzz[k];
            }
        }
    }
}

} // erirec namespace

