#include "ElectronRepulsionGeom0100ContrRecKDXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_kdxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_kdxx,
                                            const size_t idx_idxx,
                                            const size_t idx_geom_01_idxx,
                                            const size_t idx_geom_01_ifxx,
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
            /// Set up components of auxilary buffer : IDSS

            const auto id_off = idx_idxx + i * dcomps + j;

            auto g_xxxxxx_xx = cbuffer.data(id_off + 0 * ccomps * dcomps);

            auto g_xxxxxx_xy = cbuffer.data(id_off + 1 * ccomps * dcomps);

            auto g_xxxxxx_xz = cbuffer.data(id_off + 2 * ccomps * dcomps);

            auto g_xxxxxx_yy = cbuffer.data(id_off + 3 * ccomps * dcomps);

            auto g_xxxxxx_yz = cbuffer.data(id_off + 4 * ccomps * dcomps);

            auto g_xxxxxx_zz = cbuffer.data(id_off + 5 * ccomps * dcomps);

            auto g_xxxxxy_xx = cbuffer.data(id_off + 6 * ccomps * dcomps);

            auto g_xxxxxy_xy = cbuffer.data(id_off + 7 * ccomps * dcomps);

            auto g_xxxxxy_xz = cbuffer.data(id_off + 8 * ccomps * dcomps);

            auto g_xxxxxy_yy = cbuffer.data(id_off + 9 * ccomps * dcomps);

            auto g_xxxxxy_yz = cbuffer.data(id_off + 10 * ccomps * dcomps);

            auto g_xxxxxy_zz = cbuffer.data(id_off + 11 * ccomps * dcomps);

            auto g_xxxxxz_xx = cbuffer.data(id_off + 12 * ccomps * dcomps);

            auto g_xxxxxz_xy = cbuffer.data(id_off + 13 * ccomps * dcomps);

            auto g_xxxxxz_xz = cbuffer.data(id_off + 14 * ccomps * dcomps);

            auto g_xxxxxz_yy = cbuffer.data(id_off + 15 * ccomps * dcomps);

            auto g_xxxxxz_yz = cbuffer.data(id_off + 16 * ccomps * dcomps);

            auto g_xxxxxz_zz = cbuffer.data(id_off + 17 * ccomps * dcomps);

            auto g_xxxxyy_xx = cbuffer.data(id_off + 18 * ccomps * dcomps);

            auto g_xxxxyy_xy = cbuffer.data(id_off + 19 * ccomps * dcomps);

            auto g_xxxxyy_xz = cbuffer.data(id_off + 20 * ccomps * dcomps);

            auto g_xxxxyy_yy = cbuffer.data(id_off + 21 * ccomps * dcomps);

            auto g_xxxxyy_yz = cbuffer.data(id_off + 22 * ccomps * dcomps);

            auto g_xxxxyy_zz = cbuffer.data(id_off + 23 * ccomps * dcomps);

            auto g_xxxxyz_xx = cbuffer.data(id_off + 24 * ccomps * dcomps);

            auto g_xxxxyz_xy = cbuffer.data(id_off + 25 * ccomps * dcomps);

            auto g_xxxxyz_xz = cbuffer.data(id_off + 26 * ccomps * dcomps);

            auto g_xxxxyz_yy = cbuffer.data(id_off + 27 * ccomps * dcomps);

            auto g_xxxxyz_yz = cbuffer.data(id_off + 28 * ccomps * dcomps);

            auto g_xxxxyz_zz = cbuffer.data(id_off + 29 * ccomps * dcomps);

            auto g_xxxxzz_xx = cbuffer.data(id_off + 30 * ccomps * dcomps);

            auto g_xxxxzz_xy = cbuffer.data(id_off + 31 * ccomps * dcomps);

            auto g_xxxxzz_xz = cbuffer.data(id_off + 32 * ccomps * dcomps);

            auto g_xxxxzz_yy = cbuffer.data(id_off + 33 * ccomps * dcomps);

            auto g_xxxxzz_yz = cbuffer.data(id_off + 34 * ccomps * dcomps);

            auto g_xxxxzz_zz = cbuffer.data(id_off + 35 * ccomps * dcomps);

            auto g_xxxyyy_xx = cbuffer.data(id_off + 36 * ccomps * dcomps);

            auto g_xxxyyy_xy = cbuffer.data(id_off + 37 * ccomps * dcomps);

            auto g_xxxyyy_xz = cbuffer.data(id_off + 38 * ccomps * dcomps);

            auto g_xxxyyy_yy = cbuffer.data(id_off + 39 * ccomps * dcomps);

            auto g_xxxyyy_yz = cbuffer.data(id_off + 40 * ccomps * dcomps);

            auto g_xxxyyy_zz = cbuffer.data(id_off + 41 * ccomps * dcomps);

            auto g_xxxyyz_xx = cbuffer.data(id_off + 42 * ccomps * dcomps);

            auto g_xxxyyz_xy = cbuffer.data(id_off + 43 * ccomps * dcomps);

            auto g_xxxyyz_xz = cbuffer.data(id_off + 44 * ccomps * dcomps);

            auto g_xxxyyz_yy = cbuffer.data(id_off + 45 * ccomps * dcomps);

            auto g_xxxyyz_yz = cbuffer.data(id_off + 46 * ccomps * dcomps);

            auto g_xxxyyz_zz = cbuffer.data(id_off + 47 * ccomps * dcomps);

            auto g_xxxyzz_xx = cbuffer.data(id_off + 48 * ccomps * dcomps);

            auto g_xxxyzz_xy = cbuffer.data(id_off + 49 * ccomps * dcomps);

            auto g_xxxyzz_xz = cbuffer.data(id_off + 50 * ccomps * dcomps);

            auto g_xxxyzz_yy = cbuffer.data(id_off + 51 * ccomps * dcomps);

            auto g_xxxyzz_yz = cbuffer.data(id_off + 52 * ccomps * dcomps);

            auto g_xxxyzz_zz = cbuffer.data(id_off + 53 * ccomps * dcomps);

            auto g_xxxzzz_xx = cbuffer.data(id_off + 54 * ccomps * dcomps);

            auto g_xxxzzz_xy = cbuffer.data(id_off + 55 * ccomps * dcomps);

            auto g_xxxzzz_xz = cbuffer.data(id_off + 56 * ccomps * dcomps);

            auto g_xxxzzz_yy = cbuffer.data(id_off + 57 * ccomps * dcomps);

            auto g_xxxzzz_yz = cbuffer.data(id_off + 58 * ccomps * dcomps);

            auto g_xxxzzz_zz = cbuffer.data(id_off + 59 * ccomps * dcomps);

            auto g_xxyyyy_xx = cbuffer.data(id_off + 60 * ccomps * dcomps);

            auto g_xxyyyy_xy = cbuffer.data(id_off + 61 * ccomps * dcomps);

            auto g_xxyyyy_xz = cbuffer.data(id_off + 62 * ccomps * dcomps);

            auto g_xxyyyy_yy = cbuffer.data(id_off + 63 * ccomps * dcomps);

            auto g_xxyyyy_yz = cbuffer.data(id_off + 64 * ccomps * dcomps);

            auto g_xxyyyy_zz = cbuffer.data(id_off + 65 * ccomps * dcomps);

            auto g_xxyyyz_xx = cbuffer.data(id_off + 66 * ccomps * dcomps);

            auto g_xxyyyz_xy = cbuffer.data(id_off + 67 * ccomps * dcomps);

            auto g_xxyyyz_xz = cbuffer.data(id_off + 68 * ccomps * dcomps);

            auto g_xxyyyz_yy = cbuffer.data(id_off + 69 * ccomps * dcomps);

            auto g_xxyyyz_yz = cbuffer.data(id_off + 70 * ccomps * dcomps);

            auto g_xxyyyz_zz = cbuffer.data(id_off + 71 * ccomps * dcomps);

            auto g_xxyyzz_xx = cbuffer.data(id_off + 72 * ccomps * dcomps);

            auto g_xxyyzz_xy = cbuffer.data(id_off + 73 * ccomps * dcomps);

            auto g_xxyyzz_xz = cbuffer.data(id_off + 74 * ccomps * dcomps);

            auto g_xxyyzz_yy = cbuffer.data(id_off + 75 * ccomps * dcomps);

            auto g_xxyyzz_yz = cbuffer.data(id_off + 76 * ccomps * dcomps);

            auto g_xxyyzz_zz = cbuffer.data(id_off + 77 * ccomps * dcomps);

            auto g_xxyzzz_xx = cbuffer.data(id_off + 78 * ccomps * dcomps);

            auto g_xxyzzz_xy = cbuffer.data(id_off + 79 * ccomps * dcomps);

            auto g_xxyzzz_xz = cbuffer.data(id_off + 80 * ccomps * dcomps);

            auto g_xxyzzz_yy = cbuffer.data(id_off + 81 * ccomps * dcomps);

            auto g_xxyzzz_yz = cbuffer.data(id_off + 82 * ccomps * dcomps);

            auto g_xxyzzz_zz = cbuffer.data(id_off + 83 * ccomps * dcomps);

            auto g_xxzzzz_xx = cbuffer.data(id_off + 84 * ccomps * dcomps);

            auto g_xxzzzz_xy = cbuffer.data(id_off + 85 * ccomps * dcomps);

            auto g_xxzzzz_xz = cbuffer.data(id_off + 86 * ccomps * dcomps);

            auto g_xxzzzz_yy = cbuffer.data(id_off + 87 * ccomps * dcomps);

            auto g_xxzzzz_yz = cbuffer.data(id_off + 88 * ccomps * dcomps);

            auto g_xxzzzz_zz = cbuffer.data(id_off + 89 * ccomps * dcomps);

            auto g_xyyyyy_xx = cbuffer.data(id_off + 90 * ccomps * dcomps);

            auto g_xyyyyy_xy = cbuffer.data(id_off + 91 * ccomps * dcomps);

            auto g_xyyyyy_xz = cbuffer.data(id_off + 92 * ccomps * dcomps);

            auto g_xyyyyy_yy = cbuffer.data(id_off + 93 * ccomps * dcomps);

            auto g_xyyyyy_yz = cbuffer.data(id_off + 94 * ccomps * dcomps);

            auto g_xyyyyy_zz = cbuffer.data(id_off + 95 * ccomps * dcomps);

            auto g_xyyyyz_xx = cbuffer.data(id_off + 96 * ccomps * dcomps);

            auto g_xyyyyz_xy = cbuffer.data(id_off + 97 * ccomps * dcomps);

            auto g_xyyyyz_xz = cbuffer.data(id_off + 98 * ccomps * dcomps);

            auto g_xyyyyz_yy = cbuffer.data(id_off + 99 * ccomps * dcomps);

            auto g_xyyyyz_yz = cbuffer.data(id_off + 100 * ccomps * dcomps);

            auto g_xyyyyz_zz = cbuffer.data(id_off + 101 * ccomps * dcomps);

            auto g_xyyyzz_xx = cbuffer.data(id_off + 102 * ccomps * dcomps);

            auto g_xyyyzz_xy = cbuffer.data(id_off + 103 * ccomps * dcomps);

            auto g_xyyyzz_xz = cbuffer.data(id_off + 104 * ccomps * dcomps);

            auto g_xyyyzz_yy = cbuffer.data(id_off + 105 * ccomps * dcomps);

            auto g_xyyyzz_yz = cbuffer.data(id_off + 106 * ccomps * dcomps);

            auto g_xyyyzz_zz = cbuffer.data(id_off + 107 * ccomps * dcomps);

            auto g_xyyzzz_xx = cbuffer.data(id_off + 108 * ccomps * dcomps);

            auto g_xyyzzz_xy = cbuffer.data(id_off + 109 * ccomps * dcomps);

            auto g_xyyzzz_xz = cbuffer.data(id_off + 110 * ccomps * dcomps);

            auto g_xyyzzz_yy = cbuffer.data(id_off + 111 * ccomps * dcomps);

            auto g_xyyzzz_yz = cbuffer.data(id_off + 112 * ccomps * dcomps);

            auto g_xyyzzz_zz = cbuffer.data(id_off + 113 * ccomps * dcomps);

            auto g_xyzzzz_xx = cbuffer.data(id_off + 114 * ccomps * dcomps);

            auto g_xyzzzz_xy = cbuffer.data(id_off + 115 * ccomps * dcomps);

            auto g_xyzzzz_xz = cbuffer.data(id_off + 116 * ccomps * dcomps);

            auto g_xyzzzz_yy = cbuffer.data(id_off + 117 * ccomps * dcomps);

            auto g_xyzzzz_yz = cbuffer.data(id_off + 118 * ccomps * dcomps);

            auto g_xyzzzz_zz = cbuffer.data(id_off + 119 * ccomps * dcomps);

            auto g_xzzzzz_xx = cbuffer.data(id_off + 120 * ccomps * dcomps);

            auto g_xzzzzz_xy = cbuffer.data(id_off + 121 * ccomps * dcomps);

            auto g_xzzzzz_xz = cbuffer.data(id_off + 122 * ccomps * dcomps);

            auto g_xzzzzz_yy = cbuffer.data(id_off + 123 * ccomps * dcomps);

            auto g_xzzzzz_yz = cbuffer.data(id_off + 124 * ccomps * dcomps);

            auto g_xzzzzz_zz = cbuffer.data(id_off + 125 * ccomps * dcomps);

            auto g_yyyyyy_xx = cbuffer.data(id_off + 126 * ccomps * dcomps);

            auto g_yyyyyy_xy = cbuffer.data(id_off + 127 * ccomps * dcomps);

            auto g_yyyyyy_xz = cbuffer.data(id_off + 128 * ccomps * dcomps);

            auto g_yyyyyy_yy = cbuffer.data(id_off + 129 * ccomps * dcomps);

            auto g_yyyyyy_yz = cbuffer.data(id_off + 130 * ccomps * dcomps);

            auto g_yyyyyy_zz = cbuffer.data(id_off + 131 * ccomps * dcomps);

            auto g_yyyyyz_xx = cbuffer.data(id_off + 132 * ccomps * dcomps);

            auto g_yyyyyz_xy = cbuffer.data(id_off + 133 * ccomps * dcomps);

            auto g_yyyyyz_xz = cbuffer.data(id_off + 134 * ccomps * dcomps);

            auto g_yyyyyz_yy = cbuffer.data(id_off + 135 * ccomps * dcomps);

            auto g_yyyyyz_yz = cbuffer.data(id_off + 136 * ccomps * dcomps);

            auto g_yyyyyz_zz = cbuffer.data(id_off + 137 * ccomps * dcomps);

            auto g_yyyyzz_xx = cbuffer.data(id_off + 138 * ccomps * dcomps);

            auto g_yyyyzz_xy = cbuffer.data(id_off + 139 * ccomps * dcomps);

            auto g_yyyyzz_xz = cbuffer.data(id_off + 140 * ccomps * dcomps);

            auto g_yyyyzz_yy = cbuffer.data(id_off + 141 * ccomps * dcomps);

            auto g_yyyyzz_yz = cbuffer.data(id_off + 142 * ccomps * dcomps);

            auto g_yyyyzz_zz = cbuffer.data(id_off + 143 * ccomps * dcomps);

            auto g_yyyzzz_xx = cbuffer.data(id_off + 144 * ccomps * dcomps);

            auto g_yyyzzz_xy = cbuffer.data(id_off + 145 * ccomps * dcomps);

            auto g_yyyzzz_xz = cbuffer.data(id_off + 146 * ccomps * dcomps);

            auto g_yyyzzz_yy = cbuffer.data(id_off + 147 * ccomps * dcomps);

            auto g_yyyzzz_yz = cbuffer.data(id_off + 148 * ccomps * dcomps);

            auto g_yyyzzz_zz = cbuffer.data(id_off + 149 * ccomps * dcomps);

            auto g_yyzzzz_xx = cbuffer.data(id_off + 150 * ccomps * dcomps);

            auto g_yyzzzz_xy = cbuffer.data(id_off + 151 * ccomps * dcomps);

            auto g_yyzzzz_xz = cbuffer.data(id_off + 152 * ccomps * dcomps);

            auto g_yyzzzz_yy = cbuffer.data(id_off + 153 * ccomps * dcomps);

            auto g_yyzzzz_yz = cbuffer.data(id_off + 154 * ccomps * dcomps);

            auto g_yyzzzz_zz = cbuffer.data(id_off + 155 * ccomps * dcomps);

            auto g_yzzzzz_xx = cbuffer.data(id_off + 156 * ccomps * dcomps);

            auto g_yzzzzz_xy = cbuffer.data(id_off + 157 * ccomps * dcomps);

            auto g_yzzzzz_xz = cbuffer.data(id_off + 158 * ccomps * dcomps);

            auto g_yzzzzz_yy = cbuffer.data(id_off + 159 * ccomps * dcomps);

            auto g_yzzzzz_yz = cbuffer.data(id_off + 160 * ccomps * dcomps);

            auto g_yzzzzz_zz = cbuffer.data(id_off + 161 * ccomps * dcomps);

            auto g_zzzzzz_xx = cbuffer.data(id_off + 162 * ccomps * dcomps);

            auto g_zzzzzz_xy = cbuffer.data(id_off + 163 * ccomps * dcomps);

            auto g_zzzzzz_xz = cbuffer.data(id_off + 164 * ccomps * dcomps);

            auto g_zzzzzz_yy = cbuffer.data(id_off + 165 * ccomps * dcomps);

            auto g_zzzzzz_yz = cbuffer.data(id_off + 166 * ccomps * dcomps);

            auto g_zzzzzz_zz = cbuffer.data(id_off + 167 * ccomps * dcomps);

            /// Set up components of auxilary buffer : IDSS

            const auto id_geom_01_off = idx_geom_01_idxx + i * dcomps + j;

            auto g_0_x_xxxxxx_xx = cbuffer.data(id_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xy = cbuffer.data(id_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xz = cbuffer.data(id_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxxxx_yy = cbuffer.data(id_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxxxx_yz = cbuffer.data(id_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxxxx_zz = cbuffer.data(id_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xx = cbuffer.data(id_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xy = cbuffer.data(id_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xz = cbuffer.data(id_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxxxxy_yy = cbuffer.data(id_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxxxxy_yz = cbuffer.data(id_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxxxy_zz = cbuffer.data(id_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xx = cbuffer.data(id_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xy = cbuffer.data(id_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xz = cbuffer.data(id_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxxxxz_yy = cbuffer.data(id_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxxxz_yz = cbuffer.data(id_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxxxz_zz = cbuffer.data(id_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xx = cbuffer.data(id_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xy = cbuffer.data(id_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xz = cbuffer.data(id_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxxxyy_yy = cbuffer.data(id_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxxxyy_yz = cbuffer.data(id_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxxxyy_zz = cbuffer.data(id_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xx = cbuffer.data(id_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xy = cbuffer.data(id_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xz = cbuffer.data(id_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxxxyz_yy = cbuffer.data(id_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxxxyz_yz = cbuffer.data(id_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxxxyz_zz = cbuffer.data(id_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xx = cbuffer.data(id_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xy = cbuffer.data(id_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xz = cbuffer.data(id_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxxxzz_yy = cbuffer.data(id_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxxxzz_yz = cbuffer.data(id_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxxxzz_zz = cbuffer.data(id_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xx = cbuffer.data(id_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xy = cbuffer.data(id_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xz = cbuffer.data(id_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxxyyy_yy = cbuffer.data(id_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxxyyy_yz = cbuffer.data(id_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxxyyy_zz = cbuffer.data(id_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xx = cbuffer.data(id_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xy = cbuffer.data(id_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xz = cbuffer.data(id_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xxxyyz_yy = cbuffer.data(id_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xxxyyz_yz = cbuffer.data(id_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xxxyyz_zz = cbuffer.data(id_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xx = cbuffer.data(id_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xy = cbuffer.data(id_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xz = cbuffer.data(id_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xxxyzz_yy = cbuffer.data(id_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xxxyzz_yz = cbuffer.data(id_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xxxyzz_zz = cbuffer.data(id_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xx = cbuffer.data(id_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xy = cbuffer.data(id_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xz = cbuffer.data(id_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xxxzzz_yy = cbuffer.data(id_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xxxzzz_yz = cbuffer.data(id_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xxxzzz_zz = cbuffer.data(id_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xx = cbuffer.data(id_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xy = cbuffer.data(id_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xz = cbuffer.data(id_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xxyyyy_yy = cbuffer.data(id_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xxyyyy_yz = cbuffer.data(id_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xxyyyy_zz = cbuffer.data(id_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xx = cbuffer.data(id_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xy = cbuffer.data(id_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xz = cbuffer.data(id_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xxyyyz_yy = cbuffer.data(id_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xxyyyz_yz = cbuffer.data(id_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xxyyyz_zz = cbuffer.data(id_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xx = cbuffer.data(id_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xy = cbuffer.data(id_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xz = cbuffer.data(id_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xxyyzz_yy = cbuffer.data(id_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xxyyzz_yz = cbuffer.data(id_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xxyyzz_zz = cbuffer.data(id_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xx = cbuffer.data(id_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xy = cbuffer.data(id_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xz = cbuffer.data(id_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xxyzzz_yy = cbuffer.data(id_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xxyzzz_yz = cbuffer.data(id_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xxyzzz_zz = cbuffer.data(id_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xx = cbuffer.data(id_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xy = cbuffer.data(id_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xz = cbuffer.data(id_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xxzzzz_yy = cbuffer.data(id_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xxzzzz_yz = cbuffer.data(id_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xxzzzz_zz = cbuffer.data(id_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xx = cbuffer.data(id_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xy = cbuffer.data(id_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xz = cbuffer.data(id_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_xyyyyy_yy = cbuffer.data(id_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_xyyyyy_yz = cbuffer.data(id_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_xyyyyy_zz = cbuffer.data(id_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xx = cbuffer.data(id_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xy = cbuffer.data(id_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xz = cbuffer.data(id_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_xyyyyz_yy = cbuffer.data(id_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_xyyyyz_yz = cbuffer.data(id_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_xyyyyz_zz = cbuffer.data(id_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xx = cbuffer.data(id_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xy = cbuffer.data(id_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xz = cbuffer.data(id_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_xyyyzz_yy = cbuffer.data(id_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_xyyyzz_yz = cbuffer.data(id_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_xyyyzz_zz = cbuffer.data(id_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xx = cbuffer.data(id_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xy = cbuffer.data(id_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xz = cbuffer.data(id_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_xyyzzz_yy = cbuffer.data(id_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_xyyzzz_yz = cbuffer.data(id_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_xyyzzz_zz = cbuffer.data(id_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xx = cbuffer.data(id_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xy = cbuffer.data(id_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xz = cbuffer.data(id_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_xyzzzz_yy = cbuffer.data(id_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_xyzzzz_yz = cbuffer.data(id_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_xyzzzz_zz = cbuffer.data(id_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xx = cbuffer.data(id_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xy = cbuffer.data(id_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xz = cbuffer.data(id_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_xzzzzz_yy = cbuffer.data(id_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_xzzzzz_yz = cbuffer.data(id_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_xzzzzz_zz = cbuffer.data(id_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xx = cbuffer.data(id_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xy = cbuffer.data(id_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xz = cbuffer.data(id_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_yyyyyy_yy = cbuffer.data(id_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_yyyyyy_yz = cbuffer.data(id_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_yyyyyy_zz = cbuffer.data(id_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xx = cbuffer.data(id_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xy = cbuffer.data(id_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xz = cbuffer.data(id_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_x_yyyyyz_yy = cbuffer.data(id_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_yyyyyz_yz = cbuffer.data(id_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_yyyyyz_zz = cbuffer.data(id_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xx = cbuffer.data(id_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xy = cbuffer.data(id_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xz = cbuffer.data(id_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_yyyyzz_yy = cbuffer.data(id_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_yyyyzz_yz = cbuffer.data(id_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_yyyyzz_zz = cbuffer.data(id_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xx = cbuffer.data(id_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xy = cbuffer.data(id_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xz = cbuffer.data(id_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_yyyzzz_yy = cbuffer.data(id_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_yyyzzz_yz = cbuffer.data(id_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_yyyzzz_zz = cbuffer.data(id_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xx = cbuffer.data(id_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xy = cbuffer.data(id_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xz = cbuffer.data(id_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_x_yyzzzz_yy = cbuffer.data(id_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_x_yyzzzz_yz = cbuffer.data(id_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_x_yyzzzz_zz = cbuffer.data(id_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xx = cbuffer.data(id_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xy = cbuffer.data(id_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xz = cbuffer.data(id_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_x_yzzzzz_yy = cbuffer.data(id_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_x_yzzzzz_yz = cbuffer.data(id_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_x_yzzzzz_zz = cbuffer.data(id_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xx = cbuffer.data(id_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xy = cbuffer.data(id_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xz = cbuffer.data(id_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_x_zzzzzz_yy = cbuffer.data(id_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_x_zzzzzz_yz = cbuffer.data(id_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_x_zzzzzz_zz = cbuffer.data(id_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xx = cbuffer.data(id_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xy = cbuffer.data(id_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xz = cbuffer.data(id_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_y_xxxxxx_yy = cbuffer.data(id_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_y_xxxxxx_yz = cbuffer.data(id_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_y_xxxxxx_zz = cbuffer.data(id_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xx = cbuffer.data(id_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xy = cbuffer.data(id_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xz = cbuffer.data(id_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_y_xxxxxy_yy = cbuffer.data(id_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_y_xxxxxy_yz = cbuffer.data(id_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_y_xxxxxy_zz = cbuffer.data(id_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xx = cbuffer.data(id_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xy = cbuffer.data(id_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xz = cbuffer.data(id_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_y_xxxxxz_yy = cbuffer.data(id_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_y_xxxxxz_yz = cbuffer.data(id_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_y_xxxxxz_zz = cbuffer.data(id_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xx = cbuffer.data(id_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xy = cbuffer.data(id_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xz = cbuffer.data(id_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_y_xxxxyy_yy = cbuffer.data(id_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_y_xxxxyy_yz = cbuffer.data(id_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_y_xxxxyy_zz = cbuffer.data(id_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xx = cbuffer.data(id_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xy = cbuffer.data(id_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xz = cbuffer.data(id_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_y_xxxxyz_yy = cbuffer.data(id_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_y_xxxxyz_yz = cbuffer.data(id_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_y_xxxxyz_zz = cbuffer.data(id_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xx = cbuffer.data(id_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xy = cbuffer.data(id_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xz = cbuffer.data(id_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_y_xxxxzz_yy = cbuffer.data(id_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_y_xxxxzz_yz = cbuffer.data(id_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_y_xxxxzz_zz = cbuffer.data(id_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xx = cbuffer.data(id_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xy = cbuffer.data(id_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xz = cbuffer.data(id_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_y_xxxyyy_yy = cbuffer.data(id_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_y_xxxyyy_yz = cbuffer.data(id_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_y_xxxyyy_zz = cbuffer.data(id_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xx = cbuffer.data(id_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xy = cbuffer.data(id_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xz = cbuffer.data(id_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_y_xxxyyz_yy = cbuffer.data(id_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_y_xxxyyz_yz = cbuffer.data(id_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_y_xxxyyz_zz = cbuffer.data(id_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xx = cbuffer.data(id_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xy = cbuffer.data(id_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xz = cbuffer.data(id_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_y_xxxyzz_yy = cbuffer.data(id_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_y_xxxyzz_yz = cbuffer.data(id_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_y_xxxyzz_zz = cbuffer.data(id_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xx = cbuffer.data(id_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xy = cbuffer.data(id_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xz = cbuffer.data(id_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_y_xxxzzz_yy = cbuffer.data(id_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_y_xxxzzz_yz = cbuffer.data(id_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_y_xxxzzz_zz = cbuffer.data(id_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xx = cbuffer.data(id_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xy = cbuffer.data(id_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xz = cbuffer.data(id_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_y_xxyyyy_yy = cbuffer.data(id_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_y_xxyyyy_yz = cbuffer.data(id_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_y_xxyyyy_zz = cbuffer.data(id_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xx = cbuffer.data(id_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xy = cbuffer.data(id_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xz = cbuffer.data(id_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_y_xxyyyz_yy = cbuffer.data(id_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_y_xxyyyz_yz = cbuffer.data(id_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_y_xxyyyz_zz = cbuffer.data(id_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xx = cbuffer.data(id_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xy = cbuffer.data(id_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xz = cbuffer.data(id_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_y_xxyyzz_yy = cbuffer.data(id_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_y_xxyyzz_yz = cbuffer.data(id_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_y_xxyyzz_zz = cbuffer.data(id_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xx = cbuffer.data(id_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xy = cbuffer.data(id_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xz = cbuffer.data(id_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_y_xxyzzz_yy = cbuffer.data(id_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_y_xxyzzz_yz = cbuffer.data(id_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_y_xxyzzz_zz = cbuffer.data(id_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xx = cbuffer.data(id_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xy = cbuffer.data(id_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xz = cbuffer.data(id_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_y_xxzzzz_yy = cbuffer.data(id_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_y_xxzzzz_yz = cbuffer.data(id_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_y_xxzzzz_zz = cbuffer.data(id_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xx = cbuffer.data(id_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xy = cbuffer.data(id_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xz = cbuffer.data(id_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_y_xyyyyy_yy = cbuffer.data(id_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_y_xyyyyy_yz = cbuffer.data(id_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_y_xyyyyy_zz = cbuffer.data(id_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xx = cbuffer.data(id_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xy = cbuffer.data(id_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xz = cbuffer.data(id_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_y_xyyyyz_yy = cbuffer.data(id_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_y_xyyyyz_yz = cbuffer.data(id_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_y_xyyyyz_zz = cbuffer.data(id_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xx = cbuffer.data(id_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xy = cbuffer.data(id_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xz = cbuffer.data(id_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_y_xyyyzz_yy = cbuffer.data(id_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_y_xyyyzz_yz = cbuffer.data(id_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_y_xyyyzz_zz = cbuffer.data(id_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xx = cbuffer.data(id_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xy = cbuffer.data(id_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xz = cbuffer.data(id_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_y_xyyzzz_yy = cbuffer.data(id_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_y_xyyzzz_yz = cbuffer.data(id_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_y_xyyzzz_zz = cbuffer.data(id_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xx = cbuffer.data(id_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xy = cbuffer.data(id_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xz = cbuffer.data(id_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_y_xyzzzz_yy = cbuffer.data(id_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_y_xyzzzz_yz = cbuffer.data(id_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_y_xyzzzz_zz = cbuffer.data(id_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xx = cbuffer.data(id_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xy = cbuffer.data(id_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xz = cbuffer.data(id_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_y_xzzzzz_yy = cbuffer.data(id_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_y_xzzzzz_yz = cbuffer.data(id_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_y_xzzzzz_zz = cbuffer.data(id_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xx = cbuffer.data(id_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xy = cbuffer.data(id_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xz = cbuffer.data(id_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_y_yyyyyy_yy = cbuffer.data(id_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_y_yyyyyy_yz = cbuffer.data(id_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_y_yyyyyy_zz = cbuffer.data(id_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xx = cbuffer.data(id_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xy = cbuffer.data(id_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xz = cbuffer.data(id_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_y_yyyyyz_yy = cbuffer.data(id_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_y_yyyyyz_yz = cbuffer.data(id_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_y_yyyyyz_zz = cbuffer.data(id_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xx = cbuffer.data(id_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xy = cbuffer.data(id_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xz = cbuffer.data(id_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_y_yyyyzz_yy = cbuffer.data(id_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_y_yyyyzz_yz = cbuffer.data(id_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_y_yyyyzz_zz = cbuffer.data(id_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xx = cbuffer.data(id_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xy = cbuffer.data(id_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xz = cbuffer.data(id_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_y_yyyzzz_yy = cbuffer.data(id_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_y_yyyzzz_yz = cbuffer.data(id_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_y_yyyzzz_zz = cbuffer.data(id_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xx = cbuffer.data(id_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xy = cbuffer.data(id_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xz = cbuffer.data(id_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_y_yyzzzz_yy = cbuffer.data(id_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_y_yyzzzz_yz = cbuffer.data(id_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_y_yyzzzz_zz = cbuffer.data(id_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xx = cbuffer.data(id_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xy = cbuffer.data(id_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xz = cbuffer.data(id_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_y_yzzzzz_yy = cbuffer.data(id_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_y_yzzzzz_yz = cbuffer.data(id_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_y_yzzzzz_zz = cbuffer.data(id_geom_01_off + 329 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xx = cbuffer.data(id_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xy = cbuffer.data(id_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xz = cbuffer.data(id_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_y_zzzzzz_yy = cbuffer.data(id_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_y_zzzzzz_yz = cbuffer.data(id_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_y_zzzzzz_zz = cbuffer.data(id_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xx = cbuffer.data(id_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xy = cbuffer.data(id_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xz = cbuffer.data(id_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_z_xxxxxx_yy = cbuffer.data(id_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_z_xxxxxx_yz = cbuffer.data(id_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_z_xxxxxx_zz = cbuffer.data(id_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xx = cbuffer.data(id_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xy = cbuffer.data(id_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xz = cbuffer.data(id_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_z_xxxxxy_yy = cbuffer.data(id_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_z_xxxxxy_yz = cbuffer.data(id_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_z_xxxxxy_zz = cbuffer.data(id_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xx = cbuffer.data(id_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xy = cbuffer.data(id_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xz = cbuffer.data(id_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_z_xxxxxz_yy = cbuffer.data(id_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_z_xxxxxz_yz = cbuffer.data(id_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_z_xxxxxz_zz = cbuffer.data(id_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xx = cbuffer.data(id_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xy = cbuffer.data(id_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xz = cbuffer.data(id_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_z_xxxxyy_yy = cbuffer.data(id_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_z_xxxxyy_yz = cbuffer.data(id_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_z_xxxxyy_zz = cbuffer.data(id_geom_01_off + 359 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xx = cbuffer.data(id_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xy = cbuffer.data(id_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xz = cbuffer.data(id_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_z_xxxxyz_yy = cbuffer.data(id_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_z_xxxxyz_yz = cbuffer.data(id_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_z_xxxxyz_zz = cbuffer.data(id_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xx = cbuffer.data(id_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xy = cbuffer.data(id_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xz = cbuffer.data(id_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_z_xxxxzz_yy = cbuffer.data(id_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_z_xxxxzz_yz = cbuffer.data(id_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_z_xxxxzz_zz = cbuffer.data(id_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xx = cbuffer.data(id_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xy = cbuffer.data(id_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xz = cbuffer.data(id_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_z_xxxyyy_yy = cbuffer.data(id_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_z_xxxyyy_yz = cbuffer.data(id_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_z_xxxyyy_zz = cbuffer.data(id_geom_01_off + 377 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xx = cbuffer.data(id_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xy = cbuffer.data(id_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xz = cbuffer.data(id_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_z_xxxyyz_yy = cbuffer.data(id_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_z_xxxyyz_yz = cbuffer.data(id_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_z_xxxyyz_zz = cbuffer.data(id_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xx = cbuffer.data(id_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xy = cbuffer.data(id_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xz = cbuffer.data(id_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_z_xxxyzz_yy = cbuffer.data(id_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_z_xxxyzz_yz = cbuffer.data(id_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_z_xxxyzz_zz = cbuffer.data(id_geom_01_off + 389 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xx = cbuffer.data(id_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xy = cbuffer.data(id_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xz = cbuffer.data(id_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_z_xxxzzz_yy = cbuffer.data(id_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_z_xxxzzz_yz = cbuffer.data(id_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_z_xxxzzz_zz = cbuffer.data(id_geom_01_off + 395 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xx = cbuffer.data(id_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xy = cbuffer.data(id_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xz = cbuffer.data(id_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_z_xxyyyy_yy = cbuffer.data(id_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_z_xxyyyy_yz = cbuffer.data(id_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_z_xxyyyy_zz = cbuffer.data(id_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xx = cbuffer.data(id_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xy = cbuffer.data(id_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xz = cbuffer.data(id_geom_01_off + 404 * ccomps * dcomps);

            auto g_0_z_xxyyyz_yy = cbuffer.data(id_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_z_xxyyyz_yz = cbuffer.data(id_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_z_xxyyyz_zz = cbuffer.data(id_geom_01_off + 407 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xx = cbuffer.data(id_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xy = cbuffer.data(id_geom_01_off + 409 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xz = cbuffer.data(id_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_z_xxyyzz_yy = cbuffer.data(id_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_z_xxyyzz_yz = cbuffer.data(id_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_z_xxyyzz_zz = cbuffer.data(id_geom_01_off + 413 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xx = cbuffer.data(id_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xy = cbuffer.data(id_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xz = cbuffer.data(id_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_z_xxyzzz_yy = cbuffer.data(id_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_z_xxyzzz_yz = cbuffer.data(id_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_z_xxyzzz_zz = cbuffer.data(id_geom_01_off + 419 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xx = cbuffer.data(id_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xy = cbuffer.data(id_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xz = cbuffer.data(id_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_z_xxzzzz_yy = cbuffer.data(id_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_z_xxzzzz_yz = cbuffer.data(id_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_z_xxzzzz_zz = cbuffer.data(id_geom_01_off + 425 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xx = cbuffer.data(id_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xy = cbuffer.data(id_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xz = cbuffer.data(id_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_z_xyyyyy_yy = cbuffer.data(id_geom_01_off + 429 * ccomps * dcomps);

            auto g_0_z_xyyyyy_yz = cbuffer.data(id_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_z_xyyyyy_zz = cbuffer.data(id_geom_01_off + 431 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xx = cbuffer.data(id_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xy = cbuffer.data(id_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xz = cbuffer.data(id_geom_01_off + 434 * ccomps * dcomps);

            auto g_0_z_xyyyyz_yy = cbuffer.data(id_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_z_xyyyyz_yz = cbuffer.data(id_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_z_xyyyyz_zz = cbuffer.data(id_geom_01_off + 437 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xx = cbuffer.data(id_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xy = cbuffer.data(id_geom_01_off + 439 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xz = cbuffer.data(id_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_z_xyyyzz_yy = cbuffer.data(id_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_z_xyyyzz_yz = cbuffer.data(id_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_z_xyyyzz_zz = cbuffer.data(id_geom_01_off + 443 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xx = cbuffer.data(id_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xy = cbuffer.data(id_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xz = cbuffer.data(id_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_z_xyyzzz_yy = cbuffer.data(id_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_z_xyyzzz_yz = cbuffer.data(id_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_z_xyyzzz_zz = cbuffer.data(id_geom_01_off + 449 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xx = cbuffer.data(id_geom_01_off + 450 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xy = cbuffer.data(id_geom_01_off + 451 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xz = cbuffer.data(id_geom_01_off + 452 * ccomps * dcomps);

            auto g_0_z_xyzzzz_yy = cbuffer.data(id_geom_01_off + 453 * ccomps * dcomps);

            auto g_0_z_xyzzzz_yz = cbuffer.data(id_geom_01_off + 454 * ccomps * dcomps);

            auto g_0_z_xyzzzz_zz = cbuffer.data(id_geom_01_off + 455 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xx = cbuffer.data(id_geom_01_off + 456 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xy = cbuffer.data(id_geom_01_off + 457 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xz = cbuffer.data(id_geom_01_off + 458 * ccomps * dcomps);

            auto g_0_z_xzzzzz_yy = cbuffer.data(id_geom_01_off + 459 * ccomps * dcomps);

            auto g_0_z_xzzzzz_yz = cbuffer.data(id_geom_01_off + 460 * ccomps * dcomps);

            auto g_0_z_xzzzzz_zz = cbuffer.data(id_geom_01_off + 461 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xx = cbuffer.data(id_geom_01_off + 462 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xy = cbuffer.data(id_geom_01_off + 463 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xz = cbuffer.data(id_geom_01_off + 464 * ccomps * dcomps);

            auto g_0_z_yyyyyy_yy = cbuffer.data(id_geom_01_off + 465 * ccomps * dcomps);

            auto g_0_z_yyyyyy_yz = cbuffer.data(id_geom_01_off + 466 * ccomps * dcomps);

            auto g_0_z_yyyyyy_zz = cbuffer.data(id_geom_01_off + 467 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xx = cbuffer.data(id_geom_01_off + 468 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xy = cbuffer.data(id_geom_01_off + 469 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xz = cbuffer.data(id_geom_01_off + 470 * ccomps * dcomps);

            auto g_0_z_yyyyyz_yy = cbuffer.data(id_geom_01_off + 471 * ccomps * dcomps);

            auto g_0_z_yyyyyz_yz = cbuffer.data(id_geom_01_off + 472 * ccomps * dcomps);

            auto g_0_z_yyyyyz_zz = cbuffer.data(id_geom_01_off + 473 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xx = cbuffer.data(id_geom_01_off + 474 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xy = cbuffer.data(id_geom_01_off + 475 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xz = cbuffer.data(id_geom_01_off + 476 * ccomps * dcomps);

            auto g_0_z_yyyyzz_yy = cbuffer.data(id_geom_01_off + 477 * ccomps * dcomps);

            auto g_0_z_yyyyzz_yz = cbuffer.data(id_geom_01_off + 478 * ccomps * dcomps);

            auto g_0_z_yyyyzz_zz = cbuffer.data(id_geom_01_off + 479 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xx = cbuffer.data(id_geom_01_off + 480 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xy = cbuffer.data(id_geom_01_off + 481 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xz = cbuffer.data(id_geom_01_off + 482 * ccomps * dcomps);

            auto g_0_z_yyyzzz_yy = cbuffer.data(id_geom_01_off + 483 * ccomps * dcomps);

            auto g_0_z_yyyzzz_yz = cbuffer.data(id_geom_01_off + 484 * ccomps * dcomps);

            auto g_0_z_yyyzzz_zz = cbuffer.data(id_geom_01_off + 485 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xx = cbuffer.data(id_geom_01_off + 486 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xy = cbuffer.data(id_geom_01_off + 487 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xz = cbuffer.data(id_geom_01_off + 488 * ccomps * dcomps);

            auto g_0_z_yyzzzz_yy = cbuffer.data(id_geom_01_off + 489 * ccomps * dcomps);

            auto g_0_z_yyzzzz_yz = cbuffer.data(id_geom_01_off + 490 * ccomps * dcomps);

            auto g_0_z_yyzzzz_zz = cbuffer.data(id_geom_01_off + 491 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xx = cbuffer.data(id_geom_01_off + 492 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xy = cbuffer.data(id_geom_01_off + 493 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xz = cbuffer.data(id_geom_01_off + 494 * ccomps * dcomps);

            auto g_0_z_yzzzzz_yy = cbuffer.data(id_geom_01_off + 495 * ccomps * dcomps);

            auto g_0_z_yzzzzz_yz = cbuffer.data(id_geom_01_off + 496 * ccomps * dcomps);

            auto g_0_z_yzzzzz_zz = cbuffer.data(id_geom_01_off + 497 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xx = cbuffer.data(id_geom_01_off + 498 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xy = cbuffer.data(id_geom_01_off + 499 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xz = cbuffer.data(id_geom_01_off + 500 * ccomps * dcomps);

            auto g_0_z_zzzzzz_yy = cbuffer.data(id_geom_01_off + 501 * ccomps * dcomps);

            auto g_0_z_zzzzzz_yz = cbuffer.data(id_geom_01_off + 502 * ccomps * dcomps);

            auto g_0_z_zzzzzz_zz = cbuffer.data(id_geom_01_off + 503 * ccomps * dcomps);

            /// Set up components of auxilary buffer : IFSS

            const auto if_geom_01_off = idx_geom_01_ifxx + i * dcomps + j;

            auto g_0_x_xxxxxx_xxx = cbuffer.data(if_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xxy = cbuffer.data(if_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xxz = cbuffer.data(if_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xyy = cbuffer.data(if_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xyz = cbuffer.data(if_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xzz = cbuffer.data(if_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxxxxx_yyy = cbuffer.data(if_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxxxx_yyz = cbuffer.data(if_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxxxx_yzz = cbuffer.data(if_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxxxxx_zzz = cbuffer.data(if_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xxx = cbuffer.data(if_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xxy = cbuffer.data(if_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xxz = cbuffer.data(if_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xyy = cbuffer.data(if_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xyz = cbuffer.data(if_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xzz = cbuffer.data(if_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxxxy_yyy = cbuffer.data(if_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxxxy_yyz = cbuffer.data(if_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxxxxy_yzz = cbuffer.data(if_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxxxxy_zzz = cbuffer.data(if_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xxx = cbuffer.data(if_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xxy = cbuffer.data(if_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xxz = cbuffer.data(if_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xyy = cbuffer.data(if_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xyz = cbuffer.data(if_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xzz = cbuffer.data(if_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxxxxz_yyy = cbuffer.data(if_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxxxxz_yyz = cbuffer.data(if_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxxxxz_yzz = cbuffer.data(if_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxxxxz_zzz = cbuffer.data(if_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xxx = cbuffer.data(if_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xxy = cbuffer.data(if_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xxz = cbuffer.data(if_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xyy = cbuffer.data(if_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xyz = cbuffer.data(if_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xzz = cbuffer.data(if_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xxxxyy_yyy = cbuffer.data(if_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxxxyy_yyz = cbuffer.data(if_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxxxyy_yzz = cbuffer.data(if_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxxxyy_zzz = cbuffer.data(if_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xxx = cbuffer.data(if_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xxy = cbuffer.data(if_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xxz = cbuffer.data(if_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xyy = cbuffer.data(if_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xyz = cbuffer.data(if_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xzz = cbuffer.data(if_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xxxxyz_yyy = cbuffer.data(if_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xxxxyz_yyz = cbuffer.data(if_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xxxxyz_yzz = cbuffer.data(if_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xxxxyz_zzz = cbuffer.data(if_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xxx = cbuffer.data(if_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xxy = cbuffer.data(if_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xxz = cbuffer.data(if_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xyy = cbuffer.data(if_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xyz = cbuffer.data(if_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xzz = cbuffer.data(if_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xxxxzz_yyy = cbuffer.data(if_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xxxxzz_yyz = cbuffer.data(if_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xxxxzz_yzz = cbuffer.data(if_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xxxxzz_zzz = cbuffer.data(if_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xxx = cbuffer.data(if_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xxy = cbuffer.data(if_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xxz = cbuffer.data(if_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xyy = cbuffer.data(if_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xyz = cbuffer.data(if_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xzz = cbuffer.data(if_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xxxyyy_yyy = cbuffer.data(if_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xxxyyy_yyz = cbuffer.data(if_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xxxyyy_yzz = cbuffer.data(if_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xxxyyy_zzz = cbuffer.data(if_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xxx = cbuffer.data(if_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xxy = cbuffer.data(if_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xxz = cbuffer.data(if_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xyy = cbuffer.data(if_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xyz = cbuffer.data(if_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xzz = cbuffer.data(if_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xxxyyz_yyy = cbuffer.data(if_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xxxyyz_yyz = cbuffer.data(if_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xxxyyz_yzz = cbuffer.data(if_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xxxyyz_zzz = cbuffer.data(if_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xxx = cbuffer.data(if_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xxy = cbuffer.data(if_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xxz = cbuffer.data(if_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xyy = cbuffer.data(if_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xyz = cbuffer.data(if_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xzz = cbuffer.data(if_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xxxyzz_yyy = cbuffer.data(if_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xxxyzz_yyz = cbuffer.data(if_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xxxyzz_yzz = cbuffer.data(if_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xxxyzz_zzz = cbuffer.data(if_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xxx = cbuffer.data(if_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xxy = cbuffer.data(if_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xxz = cbuffer.data(if_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xyy = cbuffer.data(if_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xyz = cbuffer.data(if_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xzz = cbuffer.data(if_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_xxxzzz_yyy = cbuffer.data(if_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_xxxzzz_yyz = cbuffer.data(if_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_xxxzzz_yzz = cbuffer.data(if_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_xxxzzz_zzz = cbuffer.data(if_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xxx = cbuffer.data(if_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xxy = cbuffer.data(if_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xxz = cbuffer.data(if_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xyy = cbuffer.data(if_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xyz = cbuffer.data(if_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xzz = cbuffer.data(if_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_xxyyyy_yyy = cbuffer.data(if_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_xxyyyy_yyz = cbuffer.data(if_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_xxyyyy_yzz = cbuffer.data(if_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_xxyyyy_zzz = cbuffer.data(if_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xxx = cbuffer.data(if_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xxy = cbuffer.data(if_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xxz = cbuffer.data(if_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xyy = cbuffer.data(if_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xyz = cbuffer.data(if_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xzz = cbuffer.data(if_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_xxyyyz_yyy = cbuffer.data(if_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_xxyyyz_yyz = cbuffer.data(if_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_xxyyyz_yzz = cbuffer.data(if_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_xxyyyz_zzz = cbuffer.data(if_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xxx = cbuffer.data(if_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xxy = cbuffer.data(if_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xxz = cbuffer.data(if_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xyy = cbuffer.data(if_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xyz = cbuffer.data(if_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xzz = cbuffer.data(if_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_x_xxyyzz_yyy = cbuffer.data(if_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_xxyyzz_yyz = cbuffer.data(if_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_xxyyzz_yzz = cbuffer.data(if_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_xxyyzz_zzz = cbuffer.data(if_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xxx = cbuffer.data(if_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xxy = cbuffer.data(if_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xxz = cbuffer.data(if_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xyy = cbuffer.data(if_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xyz = cbuffer.data(if_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xzz = cbuffer.data(if_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_xxyzzz_yyy = cbuffer.data(if_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_xxyzzz_yyz = cbuffer.data(if_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_x_xxyzzz_yzz = cbuffer.data(if_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_xxyzzz_zzz = cbuffer.data(if_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xxx = cbuffer.data(if_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xxy = cbuffer.data(if_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xxz = cbuffer.data(if_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xyy = cbuffer.data(if_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xyz = cbuffer.data(if_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xzz = cbuffer.data(if_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_xxzzzz_yyy = cbuffer.data(if_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_xxzzzz_yyz = cbuffer.data(if_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_xxzzzz_yzz = cbuffer.data(if_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_xxzzzz_zzz = cbuffer.data(if_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xxx = cbuffer.data(if_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xxy = cbuffer.data(if_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xxz = cbuffer.data(if_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xyy = cbuffer.data(if_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xyz = cbuffer.data(if_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xzz = cbuffer.data(if_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_x_xyyyyy_yyy = cbuffer.data(if_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_x_xyyyyy_yyz = cbuffer.data(if_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_x_xyyyyy_yzz = cbuffer.data(if_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_x_xyyyyy_zzz = cbuffer.data(if_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xxx = cbuffer.data(if_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xxy = cbuffer.data(if_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xxz = cbuffer.data(if_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xyy = cbuffer.data(if_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xyz = cbuffer.data(if_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xzz = cbuffer.data(if_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_x_xyyyyz_yyy = cbuffer.data(if_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_x_xyyyyz_yyz = cbuffer.data(if_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_x_xyyyyz_yzz = cbuffer.data(if_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_x_xyyyyz_zzz = cbuffer.data(if_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xxx = cbuffer.data(if_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xxy = cbuffer.data(if_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xxz = cbuffer.data(if_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xyy = cbuffer.data(if_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xyz = cbuffer.data(if_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xzz = cbuffer.data(if_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_x_xyyyzz_yyy = cbuffer.data(if_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_x_xyyyzz_yyz = cbuffer.data(if_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_x_xyyyzz_yzz = cbuffer.data(if_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_x_xyyyzz_zzz = cbuffer.data(if_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xxx = cbuffer.data(if_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xxy = cbuffer.data(if_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xxz = cbuffer.data(if_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xyy = cbuffer.data(if_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xyz = cbuffer.data(if_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xzz = cbuffer.data(if_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_x_xyyzzz_yyy = cbuffer.data(if_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_x_xyyzzz_yyz = cbuffer.data(if_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_x_xyyzzz_yzz = cbuffer.data(if_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_x_xyyzzz_zzz = cbuffer.data(if_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xxx = cbuffer.data(if_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xxy = cbuffer.data(if_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xxz = cbuffer.data(if_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xyy = cbuffer.data(if_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xyz = cbuffer.data(if_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xzz = cbuffer.data(if_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_x_xyzzzz_yyy = cbuffer.data(if_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_x_xyzzzz_yyz = cbuffer.data(if_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_x_xyzzzz_yzz = cbuffer.data(if_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_x_xyzzzz_zzz = cbuffer.data(if_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xxx = cbuffer.data(if_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xxy = cbuffer.data(if_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xxz = cbuffer.data(if_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xyy = cbuffer.data(if_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xyz = cbuffer.data(if_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xzz = cbuffer.data(if_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_x_xzzzzz_yyy = cbuffer.data(if_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_x_xzzzzz_yyz = cbuffer.data(if_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_x_xzzzzz_yzz = cbuffer.data(if_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_x_xzzzzz_zzz = cbuffer.data(if_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xxx = cbuffer.data(if_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xxy = cbuffer.data(if_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xxz = cbuffer.data(if_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xyy = cbuffer.data(if_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xyz = cbuffer.data(if_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xzz = cbuffer.data(if_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_x_yyyyyy_yyy = cbuffer.data(if_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_x_yyyyyy_yyz = cbuffer.data(if_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_x_yyyyyy_yzz = cbuffer.data(if_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_x_yyyyyy_zzz = cbuffer.data(if_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xxx = cbuffer.data(if_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xxy = cbuffer.data(if_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xxz = cbuffer.data(if_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xyy = cbuffer.data(if_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xyz = cbuffer.data(if_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xzz = cbuffer.data(if_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_x_yyyyyz_yyy = cbuffer.data(if_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_x_yyyyyz_yyz = cbuffer.data(if_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_x_yyyyyz_yzz = cbuffer.data(if_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_x_yyyyyz_zzz = cbuffer.data(if_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xxx = cbuffer.data(if_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xxy = cbuffer.data(if_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xxz = cbuffer.data(if_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xyy = cbuffer.data(if_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xyz = cbuffer.data(if_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xzz = cbuffer.data(if_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_x_yyyyzz_yyy = cbuffer.data(if_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_x_yyyyzz_yyz = cbuffer.data(if_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_x_yyyyzz_yzz = cbuffer.data(if_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_x_yyyyzz_zzz = cbuffer.data(if_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xxx = cbuffer.data(if_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xxy = cbuffer.data(if_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xxz = cbuffer.data(if_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xyy = cbuffer.data(if_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xyz = cbuffer.data(if_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xzz = cbuffer.data(if_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_x_yyyzzz_yyy = cbuffer.data(if_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_x_yyyzzz_yyz = cbuffer.data(if_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_x_yyyzzz_yzz = cbuffer.data(if_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_x_yyyzzz_zzz = cbuffer.data(if_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xxx = cbuffer.data(if_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xxy = cbuffer.data(if_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xxz = cbuffer.data(if_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xyy = cbuffer.data(if_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xyz = cbuffer.data(if_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xzz = cbuffer.data(if_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_x_yyzzzz_yyy = cbuffer.data(if_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_x_yyzzzz_yyz = cbuffer.data(if_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_x_yyzzzz_yzz = cbuffer.data(if_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_x_yyzzzz_zzz = cbuffer.data(if_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xxx = cbuffer.data(if_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xxy = cbuffer.data(if_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xxz = cbuffer.data(if_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xyy = cbuffer.data(if_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xyz = cbuffer.data(if_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xzz = cbuffer.data(if_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_x_yzzzzz_yyy = cbuffer.data(if_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_x_yzzzzz_yyz = cbuffer.data(if_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_x_yzzzzz_yzz = cbuffer.data(if_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_x_yzzzzz_zzz = cbuffer.data(if_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xxx = cbuffer.data(if_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xxy = cbuffer.data(if_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xxz = cbuffer.data(if_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xyy = cbuffer.data(if_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xyz = cbuffer.data(if_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xzz = cbuffer.data(if_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_x_zzzzzz_yyy = cbuffer.data(if_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_x_zzzzzz_yyz = cbuffer.data(if_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_x_zzzzzz_yzz = cbuffer.data(if_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_x_zzzzzz_zzz = cbuffer.data(if_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xxx = cbuffer.data(if_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xxy = cbuffer.data(if_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xxz = cbuffer.data(if_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xyy = cbuffer.data(if_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xyz = cbuffer.data(if_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xzz = cbuffer.data(if_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_y_xxxxxx_yyy = cbuffer.data(if_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_y_xxxxxx_yyz = cbuffer.data(if_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_y_xxxxxx_yzz = cbuffer.data(if_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_y_xxxxxx_zzz = cbuffer.data(if_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xxx = cbuffer.data(if_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xxy = cbuffer.data(if_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xxz = cbuffer.data(if_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xyy = cbuffer.data(if_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xyz = cbuffer.data(if_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xzz = cbuffer.data(if_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_y_xxxxxy_yyy = cbuffer.data(if_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_y_xxxxxy_yyz = cbuffer.data(if_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_y_xxxxxy_yzz = cbuffer.data(if_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_y_xxxxxy_zzz = cbuffer.data(if_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xxx = cbuffer.data(if_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xxy = cbuffer.data(if_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xxz = cbuffer.data(if_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xyy = cbuffer.data(if_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xyz = cbuffer.data(if_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xzz = cbuffer.data(if_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_y_xxxxxz_yyy = cbuffer.data(if_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_y_xxxxxz_yyz = cbuffer.data(if_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_y_xxxxxz_yzz = cbuffer.data(if_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_y_xxxxxz_zzz = cbuffer.data(if_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xxx = cbuffer.data(if_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xxy = cbuffer.data(if_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xxz = cbuffer.data(if_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xyy = cbuffer.data(if_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xyz = cbuffer.data(if_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xzz = cbuffer.data(if_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_y_xxxxyy_yyy = cbuffer.data(if_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_y_xxxxyy_yyz = cbuffer.data(if_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_y_xxxxyy_yzz = cbuffer.data(if_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_y_xxxxyy_zzz = cbuffer.data(if_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xxx = cbuffer.data(if_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xxy = cbuffer.data(if_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xxz = cbuffer.data(if_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xyy = cbuffer.data(if_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xyz = cbuffer.data(if_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xzz = cbuffer.data(if_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_y_xxxxyz_yyy = cbuffer.data(if_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_y_xxxxyz_yyz = cbuffer.data(if_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_y_xxxxyz_yzz = cbuffer.data(if_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_y_xxxxyz_zzz = cbuffer.data(if_geom_01_off + 329 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xxx = cbuffer.data(if_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xxy = cbuffer.data(if_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xxz = cbuffer.data(if_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xyy = cbuffer.data(if_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xyz = cbuffer.data(if_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xzz = cbuffer.data(if_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_y_xxxxzz_yyy = cbuffer.data(if_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_y_xxxxzz_yyz = cbuffer.data(if_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_y_xxxxzz_yzz = cbuffer.data(if_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_y_xxxxzz_zzz = cbuffer.data(if_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xxx = cbuffer.data(if_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xxy = cbuffer.data(if_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xxz = cbuffer.data(if_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xyy = cbuffer.data(if_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xyz = cbuffer.data(if_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xzz = cbuffer.data(if_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_y_xxxyyy_yyy = cbuffer.data(if_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_y_xxxyyy_yyz = cbuffer.data(if_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_y_xxxyyy_yzz = cbuffer.data(if_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_y_xxxyyy_zzz = cbuffer.data(if_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xxx = cbuffer.data(if_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xxy = cbuffer.data(if_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xxz = cbuffer.data(if_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xyy = cbuffer.data(if_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xyz = cbuffer.data(if_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xzz = cbuffer.data(if_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_y_xxxyyz_yyy = cbuffer.data(if_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_y_xxxyyz_yyz = cbuffer.data(if_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_y_xxxyyz_yzz = cbuffer.data(if_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_y_xxxyyz_zzz = cbuffer.data(if_geom_01_off + 359 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xxx = cbuffer.data(if_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xxy = cbuffer.data(if_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xxz = cbuffer.data(if_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xyy = cbuffer.data(if_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xyz = cbuffer.data(if_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xzz = cbuffer.data(if_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_y_xxxyzz_yyy = cbuffer.data(if_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_y_xxxyzz_yyz = cbuffer.data(if_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_y_xxxyzz_yzz = cbuffer.data(if_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_y_xxxyzz_zzz = cbuffer.data(if_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xxx = cbuffer.data(if_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xxy = cbuffer.data(if_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xxz = cbuffer.data(if_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xyy = cbuffer.data(if_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xyz = cbuffer.data(if_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xzz = cbuffer.data(if_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_y_xxxzzz_yyy = cbuffer.data(if_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_y_xxxzzz_yyz = cbuffer.data(if_geom_01_off + 377 * ccomps * dcomps);

            auto g_0_y_xxxzzz_yzz = cbuffer.data(if_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_y_xxxzzz_zzz = cbuffer.data(if_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xxx = cbuffer.data(if_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xxy = cbuffer.data(if_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xxz = cbuffer.data(if_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xyy = cbuffer.data(if_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xyz = cbuffer.data(if_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xzz = cbuffer.data(if_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_y_xxyyyy_yyy = cbuffer.data(if_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_y_xxyyyy_yyz = cbuffer.data(if_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_y_xxyyyy_yzz = cbuffer.data(if_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_y_xxyyyy_zzz = cbuffer.data(if_geom_01_off + 389 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xxx = cbuffer.data(if_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xxy = cbuffer.data(if_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xxz = cbuffer.data(if_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xyy = cbuffer.data(if_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xyz = cbuffer.data(if_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xzz = cbuffer.data(if_geom_01_off + 395 * ccomps * dcomps);

            auto g_0_y_xxyyyz_yyy = cbuffer.data(if_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_y_xxyyyz_yyz = cbuffer.data(if_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_y_xxyyyz_yzz = cbuffer.data(if_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_y_xxyyyz_zzz = cbuffer.data(if_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xxx = cbuffer.data(if_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xxy = cbuffer.data(if_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xxz = cbuffer.data(if_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xyy = cbuffer.data(if_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xyz = cbuffer.data(if_geom_01_off + 404 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xzz = cbuffer.data(if_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_y_xxyyzz_yyy = cbuffer.data(if_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_y_xxyyzz_yyz = cbuffer.data(if_geom_01_off + 407 * ccomps * dcomps);

            auto g_0_y_xxyyzz_yzz = cbuffer.data(if_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_y_xxyyzz_zzz = cbuffer.data(if_geom_01_off + 409 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xxx = cbuffer.data(if_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xxy = cbuffer.data(if_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xxz = cbuffer.data(if_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xyy = cbuffer.data(if_geom_01_off + 413 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xyz = cbuffer.data(if_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xzz = cbuffer.data(if_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_y_xxyzzz_yyy = cbuffer.data(if_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_y_xxyzzz_yyz = cbuffer.data(if_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_y_xxyzzz_yzz = cbuffer.data(if_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_y_xxyzzz_zzz = cbuffer.data(if_geom_01_off + 419 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xxx = cbuffer.data(if_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xxy = cbuffer.data(if_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xxz = cbuffer.data(if_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xyy = cbuffer.data(if_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xyz = cbuffer.data(if_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xzz = cbuffer.data(if_geom_01_off + 425 * ccomps * dcomps);

            auto g_0_y_xxzzzz_yyy = cbuffer.data(if_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_y_xxzzzz_yyz = cbuffer.data(if_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_y_xxzzzz_yzz = cbuffer.data(if_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_y_xxzzzz_zzz = cbuffer.data(if_geom_01_off + 429 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xxx = cbuffer.data(if_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xxy = cbuffer.data(if_geom_01_off + 431 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xxz = cbuffer.data(if_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xyy = cbuffer.data(if_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xyz = cbuffer.data(if_geom_01_off + 434 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xzz = cbuffer.data(if_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_y_xyyyyy_yyy = cbuffer.data(if_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_y_xyyyyy_yyz = cbuffer.data(if_geom_01_off + 437 * ccomps * dcomps);

            auto g_0_y_xyyyyy_yzz = cbuffer.data(if_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_y_xyyyyy_zzz = cbuffer.data(if_geom_01_off + 439 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xxx = cbuffer.data(if_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xxy = cbuffer.data(if_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xxz = cbuffer.data(if_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xyy = cbuffer.data(if_geom_01_off + 443 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xyz = cbuffer.data(if_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xzz = cbuffer.data(if_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_y_xyyyyz_yyy = cbuffer.data(if_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_y_xyyyyz_yyz = cbuffer.data(if_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_y_xyyyyz_yzz = cbuffer.data(if_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_y_xyyyyz_zzz = cbuffer.data(if_geom_01_off + 449 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xxx = cbuffer.data(if_geom_01_off + 450 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xxy = cbuffer.data(if_geom_01_off + 451 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xxz = cbuffer.data(if_geom_01_off + 452 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xyy = cbuffer.data(if_geom_01_off + 453 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xyz = cbuffer.data(if_geom_01_off + 454 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xzz = cbuffer.data(if_geom_01_off + 455 * ccomps * dcomps);

            auto g_0_y_xyyyzz_yyy = cbuffer.data(if_geom_01_off + 456 * ccomps * dcomps);

            auto g_0_y_xyyyzz_yyz = cbuffer.data(if_geom_01_off + 457 * ccomps * dcomps);

            auto g_0_y_xyyyzz_yzz = cbuffer.data(if_geom_01_off + 458 * ccomps * dcomps);

            auto g_0_y_xyyyzz_zzz = cbuffer.data(if_geom_01_off + 459 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xxx = cbuffer.data(if_geom_01_off + 460 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xxy = cbuffer.data(if_geom_01_off + 461 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xxz = cbuffer.data(if_geom_01_off + 462 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xyy = cbuffer.data(if_geom_01_off + 463 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xyz = cbuffer.data(if_geom_01_off + 464 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xzz = cbuffer.data(if_geom_01_off + 465 * ccomps * dcomps);

            auto g_0_y_xyyzzz_yyy = cbuffer.data(if_geom_01_off + 466 * ccomps * dcomps);

            auto g_0_y_xyyzzz_yyz = cbuffer.data(if_geom_01_off + 467 * ccomps * dcomps);

            auto g_0_y_xyyzzz_yzz = cbuffer.data(if_geom_01_off + 468 * ccomps * dcomps);

            auto g_0_y_xyyzzz_zzz = cbuffer.data(if_geom_01_off + 469 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xxx = cbuffer.data(if_geom_01_off + 470 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xxy = cbuffer.data(if_geom_01_off + 471 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xxz = cbuffer.data(if_geom_01_off + 472 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xyy = cbuffer.data(if_geom_01_off + 473 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xyz = cbuffer.data(if_geom_01_off + 474 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xzz = cbuffer.data(if_geom_01_off + 475 * ccomps * dcomps);

            auto g_0_y_xyzzzz_yyy = cbuffer.data(if_geom_01_off + 476 * ccomps * dcomps);

            auto g_0_y_xyzzzz_yyz = cbuffer.data(if_geom_01_off + 477 * ccomps * dcomps);

            auto g_0_y_xyzzzz_yzz = cbuffer.data(if_geom_01_off + 478 * ccomps * dcomps);

            auto g_0_y_xyzzzz_zzz = cbuffer.data(if_geom_01_off + 479 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xxx = cbuffer.data(if_geom_01_off + 480 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xxy = cbuffer.data(if_geom_01_off + 481 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xxz = cbuffer.data(if_geom_01_off + 482 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xyy = cbuffer.data(if_geom_01_off + 483 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xyz = cbuffer.data(if_geom_01_off + 484 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xzz = cbuffer.data(if_geom_01_off + 485 * ccomps * dcomps);

            auto g_0_y_xzzzzz_yyy = cbuffer.data(if_geom_01_off + 486 * ccomps * dcomps);

            auto g_0_y_xzzzzz_yyz = cbuffer.data(if_geom_01_off + 487 * ccomps * dcomps);

            auto g_0_y_xzzzzz_yzz = cbuffer.data(if_geom_01_off + 488 * ccomps * dcomps);

            auto g_0_y_xzzzzz_zzz = cbuffer.data(if_geom_01_off + 489 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xxx = cbuffer.data(if_geom_01_off + 490 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xxy = cbuffer.data(if_geom_01_off + 491 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xxz = cbuffer.data(if_geom_01_off + 492 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xyy = cbuffer.data(if_geom_01_off + 493 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xyz = cbuffer.data(if_geom_01_off + 494 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xzz = cbuffer.data(if_geom_01_off + 495 * ccomps * dcomps);

            auto g_0_y_yyyyyy_yyy = cbuffer.data(if_geom_01_off + 496 * ccomps * dcomps);

            auto g_0_y_yyyyyy_yyz = cbuffer.data(if_geom_01_off + 497 * ccomps * dcomps);

            auto g_0_y_yyyyyy_yzz = cbuffer.data(if_geom_01_off + 498 * ccomps * dcomps);

            auto g_0_y_yyyyyy_zzz = cbuffer.data(if_geom_01_off + 499 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xxx = cbuffer.data(if_geom_01_off + 500 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xxy = cbuffer.data(if_geom_01_off + 501 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xxz = cbuffer.data(if_geom_01_off + 502 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xyy = cbuffer.data(if_geom_01_off + 503 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xyz = cbuffer.data(if_geom_01_off + 504 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xzz = cbuffer.data(if_geom_01_off + 505 * ccomps * dcomps);

            auto g_0_y_yyyyyz_yyy = cbuffer.data(if_geom_01_off + 506 * ccomps * dcomps);

            auto g_0_y_yyyyyz_yyz = cbuffer.data(if_geom_01_off + 507 * ccomps * dcomps);

            auto g_0_y_yyyyyz_yzz = cbuffer.data(if_geom_01_off + 508 * ccomps * dcomps);

            auto g_0_y_yyyyyz_zzz = cbuffer.data(if_geom_01_off + 509 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xxx = cbuffer.data(if_geom_01_off + 510 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xxy = cbuffer.data(if_geom_01_off + 511 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xxz = cbuffer.data(if_geom_01_off + 512 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xyy = cbuffer.data(if_geom_01_off + 513 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xyz = cbuffer.data(if_geom_01_off + 514 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xzz = cbuffer.data(if_geom_01_off + 515 * ccomps * dcomps);

            auto g_0_y_yyyyzz_yyy = cbuffer.data(if_geom_01_off + 516 * ccomps * dcomps);

            auto g_0_y_yyyyzz_yyz = cbuffer.data(if_geom_01_off + 517 * ccomps * dcomps);

            auto g_0_y_yyyyzz_yzz = cbuffer.data(if_geom_01_off + 518 * ccomps * dcomps);

            auto g_0_y_yyyyzz_zzz = cbuffer.data(if_geom_01_off + 519 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xxx = cbuffer.data(if_geom_01_off + 520 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xxy = cbuffer.data(if_geom_01_off + 521 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xxz = cbuffer.data(if_geom_01_off + 522 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xyy = cbuffer.data(if_geom_01_off + 523 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xyz = cbuffer.data(if_geom_01_off + 524 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xzz = cbuffer.data(if_geom_01_off + 525 * ccomps * dcomps);

            auto g_0_y_yyyzzz_yyy = cbuffer.data(if_geom_01_off + 526 * ccomps * dcomps);

            auto g_0_y_yyyzzz_yyz = cbuffer.data(if_geom_01_off + 527 * ccomps * dcomps);

            auto g_0_y_yyyzzz_yzz = cbuffer.data(if_geom_01_off + 528 * ccomps * dcomps);

            auto g_0_y_yyyzzz_zzz = cbuffer.data(if_geom_01_off + 529 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xxx = cbuffer.data(if_geom_01_off + 530 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xxy = cbuffer.data(if_geom_01_off + 531 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xxz = cbuffer.data(if_geom_01_off + 532 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xyy = cbuffer.data(if_geom_01_off + 533 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xyz = cbuffer.data(if_geom_01_off + 534 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xzz = cbuffer.data(if_geom_01_off + 535 * ccomps * dcomps);

            auto g_0_y_yyzzzz_yyy = cbuffer.data(if_geom_01_off + 536 * ccomps * dcomps);

            auto g_0_y_yyzzzz_yyz = cbuffer.data(if_geom_01_off + 537 * ccomps * dcomps);

            auto g_0_y_yyzzzz_yzz = cbuffer.data(if_geom_01_off + 538 * ccomps * dcomps);

            auto g_0_y_yyzzzz_zzz = cbuffer.data(if_geom_01_off + 539 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xxx = cbuffer.data(if_geom_01_off + 540 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xxy = cbuffer.data(if_geom_01_off + 541 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xxz = cbuffer.data(if_geom_01_off + 542 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xyy = cbuffer.data(if_geom_01_off + 543 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xyz = cbuffer.data(if_geom_01_off + 544 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xzz = cbuffer.data(if_geom_01_off + 545 * ccomps * dcomps);

            auto g_0_y_yzzzzz_yyy = cbuffer.data(if_geom_01_off + 546 * ccomps * dcomps);

            auto g_0_y_yzzzzz_yyz = cbuffer.data(if_geom_01_off + 547 * ccomps * dcomps);

            auto g_0_y_yzzzzz_yzz = cbuffer.data(if_geom_01_off + 548 * ccomps * dcomps);

            auto g_0_y_yzzzzz_zzz = cbuffer.data(if_geom_01_off + 549 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xxx = cbuffer.data(if_geom_01_off + 550 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xxy = cbuffer.data(if_geom_01_off + 551 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xxz = cbuffer.data(if_geom_01_off + 552 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xyy = cbuffer.data(if_geom_01_off + 553 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xyz = cbuffer.data(if_geom_01_off + 554 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xzz = cbuffer.data(if_geom_01_off + 555 * ccomps * dcomps);

            auto g_0_y_zzzzzz_yyy = cbuffer.data(if_geom_01_off + 556 * ccomps * dcomps);

            auto g_0_y_zzzzzz_yyz = cbuffer.data(if_geom_01_off + 557 * ccomps * dcomps);

            auto g_0_y_zzzzzz_yzz = cbuffer.data(if_geom_01_off + 558 * ccomps * dcomps);

            auto g_0_y_zzzzzz_zzz = cbuffer.data(if_geom_01_off + 559 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xxx = cbuffer.data(if_geom_01_off + 560 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xxy = cbuffer.data(if_geom_01_off + 561 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xxz = cbuffer.data(if_geom_01_off + 562 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xyy = cbuffer.data(if_geom_01_off + 563 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xyz = cbuffer.data(if_geom_01_off + 564 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xzz = cbuffer.data(if_geom_01_off + 565 * ccomps * dcomps);

            auto g_0_z_xxxxxx_yyy = cbuffer.data(if_geom_01_off + 566 * ccomps * dcomps);

            auto g_0_z_xxxxxx_yyz = cbuffer.data(if_geom_01_off + 567 * ccomps * dcomps);

            auto g_0_z_xxxxxx_yzz = cbuffer.data(if_geom_01_off + 568 * ccomps * dcomps);

            auto g_0_z_xxxxxx_zzz = cbuffer.data(if_geom_01_off + 569 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xxx = cbuffer.data(if_geom_01_off + 570 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xxy = cbuffer.data(if_geom_01_off + 571 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xxz = cbuffer.data(if_geom_01_off + 572 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xyy = cbuffer.data(if_geom_01_off + 573 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xyz = cbuffer.data(if_geom_01_off + 574 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xzz = cbuffer.data(if_geom_01_off + 575 * ccomps * dcomps);

            auto g_0_z_xxxxxy_yyy = cbuffer.data(if_geom_01_off + 576 * ccomps * dcomps);

            auto g_0_z_xxxxxy_yyz = cbuffer.data(if_geom_01_off + 577 * ccomps * dcomps);

            auto g_0_z_xxxxxy_yzz = cbuffer.data(if_geom_01_off + 578 * ccomps * dcomps);

            auto g_0_z_xxxxxy_zzz = cbuffer.data(if_geom_01_off + 579 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xxx = cbuffer.data(if_geom_01_off + 580 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xxy = cbuffer.data(if_geom_01_off + 581 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xxz = cbuffer.data(if_geom_01_off + 582 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xyy = cbuffer.data(if_geom_01_off + 583 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xyz = cbuffer.data(if_geom_01_off + 584 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xzz = cbuffer.data(if_geom_01_off + 585 * ccomps * dcomps);

            auto g_0_z_xxxxxz_yyy = cbuffer.data(if_geom_01_off + 586 * ccomps * dcomps);

            auto g_0_z_xxxxxz_yyz = cbuffer.data(if_geom_01_off + 587 * ccomps * dcomps);

            auto g_0_z_xxxxxz_yzz = cbuffer.data(if_geom_01_off + 588 * ccomps * dcomps);

            auto g_0_z_xxxxxz_zzz = cbuffer.data(if_geom_01_off + 589 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xxx = cbuffer.data(if_geom_01_off + 590 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xxy = cbuffer.data(if_geom_01_off + 591 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xxz = cbuffer.data(if_geom_01_off + 592 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xyy = cbuffer.data(if_geom_01_off + 593 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xyz = cbuffer.data(if_geom_01_off + 594 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xzz = cbuffer.data(if_geom_01_off + 595 * ccomps * dcomps);

            auto g_0_z_xxxxyy_yyy = cbuffer.data(if_geom_01_off + 596 * ccomps * dcomps);

            auto g_0_z_xxxxyy_yyz = cbuffer.data(if_geom_01_off + 597 * ccomps * dcomps);

            auto g_0_z_xxxxyy_yzz = cbuffer.data(if_geom_01_off + 598 * ccomps * dcomps);

            auto g_0_z_xxxxyy_zzz = cbuffer.data(if_geom_01_off + 599 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xxx = cbuffer.data(if_geom_01_off + 600 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xxy = cbuffer.data(if_geom_01_off + 601 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xxz = cbuffer.data(if_geom_01_off + 602 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xyy = cbuffer.data(if_geom_01_off + 603 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xyz = cbuffer.data(if_geom_01_off + 604 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xzz = cbuffer.data(if_geom_01_off + 605 * ccomps * dcomps);

            auto g_0_z_xxxxyz_yyy = cbuffer.data(if_geom_01_off + 606 * ccomps * dcomps);

            auto g_0_z_xxxxyz_yyz = cbuffer.data(if_geom_01_off + 607 * ccomps * dcomps);

            auto g_0_z_xxxxyz_yzz = cbuffer.data(if_geom_01_off + 608 * ccomps * dcomps);

            auto g_0_z_xxxxyz_zzz = cbuffer.data(if_geom_01_off + 609 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xxx = cbuffer.data(if_geom_01_off + 610 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xxy = cbuffer.data(if_geom_01_off + 611 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xxz = cbuffer.data(if_geom_01_off + 612 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xyy = cbuffer.data(if_geom_01_off + 613 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xyz = cbuffer.data(if_geom_01_off + 614 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xzz = cbuffer.data(if_geom_01_off + 615 * ccomps * dcomps);

            auto g_0_z_xxxxzz_yyy = cbuffer.data(if_geom_01_off + 616 * ccomps * dcomps);

            auto g_0_z_xxxxzz_yyz = cbuffer.data(if_geom_01_off + 617 * ccomps * dcomps);

            auto g_0_z_xxxxzz_yzz = cbuffer.data(if_geom_01_off + 618 * ccomps * dcomps);

            auto g_0_z_xxxxzz_zzz = cbuffer.data(if_geom_01_off + 619 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xxx = cbuffer.data(if_geom_01_off + 620 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xxy = cbuffer.data(if_geom_01_off + 621 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xxz = cbuffer.data(if_geom_01_off + 622 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xyy = cbuffer.data(if_geom_01_off + 623 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xyz = cbuffer.data(if_geom_01_off + 624 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xzz = cbuffer.data(if_geom_01_off + 625 * ccomps * dcomps);

            auto g_0_z_xxxyyy_yyy = cbuffer.data(if_geom_01_off + 626 * ccomps * dcomps);

            auto g_0_z_xxxyyy_yyz = cbuffer.data(if_geom_01_off + 627 * ccomps * dcomps);

            auto g_0_z_xxxyyy_yzz = cbuffer.data(if_geom_01_off + 628 * ccomps * dcomps);

            auto g_0_z_xxxyyy_zzz = cbuffer.data(if_geom_01_off + 629 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xxx = cbuffer.data(if_geom_01_off + 630 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xxy = cbuffer.data(if_geom_01_off + 631 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xxz = cbuffer.data(if_geom_01_off + 632 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xyy = cbuffer.data(if_geom_01_off + 633 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xyz = cbuffer.data(if_geom_01_off + 634 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xzz = cbuffer.data(if_geom_01_off + 635 * ccomps * dcomps);

            auto g_0_z_xxxyyz_yyy = cbuffer.data(if_geom_01_off + 636 * ccomps * dcomps);

            auto g_0_z_xxxyyz_yyz = cbuffer.data(if_geom_01_off + 637 * ccomps * dcomps);

            auto g_0_z_xxxyyz_yzz = cbuffer.data(if_geom_01_off + 638 * ccomps * dcomps);

            auto g_0_z_xxxyyz_zzz = cbuffer.data(if_geom_01_off + 639 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xxx = cbuffer.data(if_geom_01_off + 640 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xxy = cbuffer.data(if_geom_01_off + 641 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xxz = cbuffer.data(if_geom_01_off + 642 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xyy = cbuffer.data(if_geom_01_off + 643 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xyz = cbuffer.data(if_geom_01_off + 644 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xzz = cbuffer.data(if_geom_01_off + 645 * ccomps * dcomps);

            auto g_0_z_xxxyzz_yyy = cbuffer.data(if_geom_01_off + 646 * ccomps * dcomps);

            auto g_0_z_xxxyzz_yyz = cbuffer.data(if_geom_01_off + 647 * ccomps * dcomps);

            auto g_0_z_xxxyzz_yzz = cbuffer.data(if_geom_01_off + 648 * ccomps * dcomps);

            auto g_0_z_xxxyzz_zzz = cbuffer.data(if_geom_01_off + 649 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xxx = cbuffer.data(if_geom_01_off + 650 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xxy = cbuffer.data(if_geom_01_off + 651 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xxz = cbuffer.data(if_geom_01_off + 652 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xyy = cbuffer.data(if_geom_01_off + 653 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xyz = cbuffer.data(if_geom_01_off + 654 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xzz = cbuffer.data(if_geom_01_off + 655 * ccomps * dcomps);

            auto g_0_z_xxxzzz_yyy = cbuffer.data(if_geom_01_off + 656 * ccomps * dcomps);

            auto g_0_z_xxxzzz_yyz = cbuffer.data(if_geom_01_off + 657 * ccomps * dcomps);

            auto g_0_z_xxxzzz_yzz = cbuffer.data(if_geom_01_off + 658 * ccomps * dcomps);

            auto g_0_z_xxxzzz_zzz = cbuffer.data(if_geom_01_off + 659 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xxx = cbuffer.data(if_geom_01_off + 660 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xxy = cbuffer.data(if_geom_01_off + 661 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xxz = cbuffer.data(if_geom_01_off + 662 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xyy = cbuffer.data(if_geom_01_off + 663 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xyz = cbuffer.data(if_geom_01_off + 664 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xzz = cbuffer.data(if_geom_01_off + 665 * ccomps * dcomps);

            auto g_0_z_xxyyyy_yyy = cbuffer.data(if_geom_01_off + 666 * ccomps * dcomps);

            auto g_0_z_xxyyyy_yyz = cbuffer.data(if_geom_01_off + 667 * ccomps * dcomps);

            auto g_0_z_xxyyyy_yzz = cbuffer.data(if_geom_01_off + 668 * ccomps * dcomps);

            auto g_0_z_xxyyyy_zzz = cbuffer.data(if_geom_01_off + 669 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xxx = cbuffer.data(if_geom_01_off + 670 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xxy = cbuffer.data(if_geom_01_off + 671 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xxz = cbuffer.data(if_geom_01_off + 672 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xyy = cbuffer.data(if_geom_01_off + 673 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xyz = cbuffer.data(if_geom_01_off + 674 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xzz = cbuffer.data(if_geom_01_off + 675 * ccomps * dcomps);

            auto g_0_z_xxyyyz_yyy = cbuffer.data(if_geom_01_off + 676 * ccomps * dcomps);

            auto g_0_z_xxyyyz_yyz = cbuffer.data(if_geom_01_off + 677 * ccomps * dcomps);

            auto g_0_z_xxyyyz_yzz = cbuffer.data(if_geom_01_off + 678 * ccomps * dcomps);

            auto g_0_z_xxyyyz_zzz = cbuffer.data(if_geom_01_off + 679 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xxx = cbuffer.data(if_geom_01_off + 680 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xxy = cbuffer.data(if_geom_01_off + 681 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xxz = cbuffer.data(if_geom_01_off + 682 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xyy = cbuffer.data(if_geom_01_off + 683 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xyz = cbuffer.data(if_geom_01_off + 684 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xzz = cbuffer.data(if_geom_01_off + 685 * ccomps * dcomps);

            auto g_0_z_xxyyzz_yyy = cbuffer.data(if_geom_01_off + 686 * ccomps * dcomps);

            auto g_0_z_xxyyzz_yyz = cbuffer.data(if_geom_01_off + 687 * ccomps * dcomps);

            auto g_0_z_xxyyzz_yzz = cbuffer.data(if_geom_01_off + 688 * ccomps * dcomps);

            auto g_0_z_xxyyzz_zzz = cbuffer.data(if_geom_01_off + 689 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xxx = cbuffer.data(if_geom_01_off + 690 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xxy = cbuffer.data(if_geom_01_off + 691 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xxz = cbuffer.data(if_geom_01_off + 692 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xyy = cbuffer.data(if_geom_01_off + 693 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xyz = cbuffer.data(if_geom_01_off + 694 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xzz = cbuffer.data(if_geom_01_off + 695 * ccomps * dcomps);

            auto g_0_z_xxyzzz_yyy = cbuffer.data(if_geom_01_off + 696 * ccomps * dcomps);

            auto g_0_z_xxyzzz_yyz = cbuffer.data(if_geom_01_off + 697 * ccomps * dcomps);

            auto g_0_z_xxyzzz_yzz = cbuffer.data(if_geom_01_off + 698 * ccomps * dcomps);

            auto g_0_z_xxyzzz_zzz = cbuffer.data(if_geom_01_off + 699 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xxx = cbuffer.data(if_geom_01_off + 700 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xxy = cbuffer.data(if_geom_01_off + 701 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xxz = cbuffer.data(if_geom_01_off + 702 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xyy = cbuffer.data(if_geom_01_off + 703 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xyz = cbuffer.data(if_geom_01_off + 704 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xzz = cbuffer.data(if_geom_01_off + 705 * ccomps * dcomps);

            auto g_0_z_xxzzzz_yyy = cbuffer.data(if_geom_01_off + 706 * ccomps * dcomps);

            auto g_0_z_xxzzzz_yyz = cbuffer.data(if_geom_01_off + 707 * ccomps * dcomps);

            auto g_0_z_xxzzzz_yzz = cbuffer.data(if_geom_01_off + 708 * ccomps * dcomps);

            auto g_0_z_xxzzzz_zzz = cbuffer.data(if_geom_01_off + 709 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xxx = cbuffer.data(if_geom_01_off + 710 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xxy = cbuffer.data(if_geom_01_off + 711 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xxz = cbuffer.data(if_geom_01_off + 712 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xyy = cbuffer.data(if_geom_01_off + 713 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xyz = cbuffer.data(if_geom_01_off + 714 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xzz = cbuffer.data(if_geom_01_off + 715 * ccomps * dcomps);

            auto g_0_z_xyyyyy_yyy = cbuffer.data(if_geom_01_off + 716 * ccomps * dcomps);

            auto g_0_z_xyyyyy_yyz = cbuffer.data(if_geom_01_off + 717 * ccomps * dcomps);

            auto g_0_z_xyyyyy_yzz = cbuffer.data(if_geom_01_off + 718 * ccomps * dcomps);

            auto g_0_z_xyyyyy_zzz = cbuffer.data(if_geom_01_off + 719 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xxx = cbuffer.data(if_geom_01_off + 720 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xxy = cbuffer.data(if_geom_01_off + 721 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xxz = cbuffer.data(if_geom_01_off + 722 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xyy = cbuffer.data(if_geom_01_off + 723 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xyz = cbuffer.data(if_geom_01_off + 724 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xzz = cbuffer.data(if_geom_01_off + 725 * ccomps * dcomps);

            auto g_0_z_xyyyyz_yyy = cbuffer.data(if_geom_01_off + 726 * ccomps * dcomps);

            auto g_0_z_xyyyyz_yyz = cbuffer.data(if_geom_01_off + 727 * ccomps * dcomps);

            auto g_0_z_xyyyyz_yzz = cbuffer.data(if_geom_01_off + 728 * ccomps * dcomps);

            auto g_0_z_xyyyyz_zzz = cbuffer.data(if_geom_01_off + 729 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xxx = cbuffer.data(if_geom_01_off + 730 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xxy = cbuffer.data(if_geom_01_off + 731 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xxz = cbuffer.data(if_geom_01_off + 732 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xyy = cbuffer.data(if_geom_01_off + 733 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xyz = cbuffer.data(if_geom_01_off + 734 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xzz = cbuffer.data(if_geom_01_off + 735 * ccomps * dcomps);

            auto g_0_z_xyyyzz_yyy = cbuffer.data(if_geom_01_off + 736 * ccomps * dcomps);

            auto g_0_z_xyyyzz_yyz = cbuffer.data(if_geom_01_off + 737 * ccomps * dcomps);

            auto g_0_z_xyyyzz_yzz = cbuffer.data(if_geom_01_off + 738 * ccomps * dcomps);

            auto g_0_z_xyyyzz_zzz = cbuffer.data(if_geom_01_off + 739 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xxx = cbuffer.data(if_geom_01_off + 740 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xxy = cbuffer.data(if_geom_01_off + 741 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xxz = cbuffer.data(if_geom_01_off + 742 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xyy = cbuffer.data(if_geom_01_off + 743 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xyz = cbuffer.data(if_geom_01_off + 744 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xzz = cbuffer.data(if_geom_01_off + 745 * ccomps * dcomps);

            auto g_0_z_xyyzzz_yyy = cbuffer.data(if_geom_01_off + 746 * ccomps * dcomps);

            auto g_0_z_xyyzzz_yyz = cbuffer.data(if_geom_01_off + 747 * ccomps * dcomps);

            auto g_0_z_xyyzzz_yzz = cbuffer.data(if_geom_01_off + 748 * ccomps * dcomps);

            auto g_0_z_xyyzzz_zzz = cbuffer.data(if_geom_01_off + 749 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xxx = cbuffer.data(if_geom_01_off + 750 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xxy = cbuffer.data(if_geom_01_off + 751 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xxz = cbuffer.data(if_geom_01_off + 752 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xyy = cbuffer.data(if_geom_01_off + 753 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xyz = cbuffer.data(if_geom_01_off + 754 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xzz = cbuffer.data(if_geom_01_off + 755 * ccomps * dcomps);

            auto g_0_z_xyzzzz_yyy = cbuffer.data(if_geom_01_off + 756 * ccomps * dcomps);

            auto g_0_z_xyzzzz_yyz = cbuffer.data(if_geom_01_off + 757 * ccomps * dcomps);

            auto g_0_z_xyzzzz_yzz = cbuffer.data(if_geom_01_off + 758 * ccomps * dcomps);

            auto g_0_z_xyzzzz_zzz = cbuffer.data(if_geom_01_off + 759 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xxx = cbuffer.data(if_geom_01_off + 760 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xxy = cbuffer.data(if_geom_01_off + 761 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xxz = cbuffer.data(if_geom_01_off + 762 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xyy = cbuffer.data(if_geom_01_off + 763 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xyz = cbuffer.data(if_geom_01_off + 764 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xzz = cbuffer.data(if_geom_01_off + 765 * ccomps * dcomps);

            auto g_0_z_xzzzzz_yyy = cbuffer.data(if_geom_01_off + 766 * ccomps * dcomps);

            auto g_0_z_xzzzzz_yyz = cbuffer.data(if_geom_01_off + 767 * ccomps * dcomps);

            auto g_0_z_xzzzzz_yzz = cbuffer.data(if_geom_01_off + 768 * ccomps * dcomps);

            auto g_0_z_xzzzzz_zzz = cbuffer.data(if_geom_01_off + 769 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xxx = cbuffer.data(if_geom_01_off + 770 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xxy = cbuffer.data(if_geom_01_off + 771 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xxz = cbuffer.data(if_geom_01_off + 772 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xyy = cbuffer.data(if_geom_01_off + 773 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xyz = cbuffer.data(if_geom_01_off + 774 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xzz = cbuffer.data(if_geom_01_off + 775 * ccomps * dcomps);

            auto g_0_z_yyyyyy_yyy = cbuffer.data(if_geom_01_off + 776 * ccomps * dcomps);

            auto g_0_z_yyyyyy_yyz = cbuffer.data(if_geom_01_off + 777 * ccomps * dcomps);

            auto g_0_z_yyyyyy_yzz = cbuffer.data(if_geom_01_off + 778 * ccomps * dcomps);

            auto g_0_z_yyyyyy_zzz = cbuffer.data(if_geom_01_off + 779 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xxx = cbuffer.data(if_geom_01_off + 780 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xxy = cbuffer.data(if_geom_01_off + 781 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xxz = cbuffer.data(if_geom_01_off + 782 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xyy = cbuffer.data(if_geom_01_off + 783 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xyz = cbuffer.data(if_geom_01_off + 784 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xzz = cbuffer.data(if_geom_01_off + 785 * ccomps * dcomps);

            auto g_0_z_yyyyyz_yyy = cbuffer.data(if_geom_01_off + 786 * ccomps * dcomps);

            auto g_0_z_yyyyyz_yyz = cbuffer.data(if_geom_01_off + 787 * ccomps * dcomps);

            auto g_0_z_yyyyyz_yzz = cbuffer.data(if_geom_01_off + 788 * ccomps * dcomps);

            auto g_0_z_yyyyyz_zzz = cbuffer.data(if_geom_01_off + 789 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xxx = cbuffer.data(if_geom_01_off + 790 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xxy = cbuffer.data(if_geom_01_off + 791 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xxz = cbuffer.data(if_geom_01_off + 792 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xyy = cbuffer.data(if_geom_01_off + 793 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xyz = cbuffer.data(if_geom_01_off + 794 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xzz = cbuffer.data(if_geom_01_off + 795 * ccomps * dcomps);

            auto g_0_z_yyyyzz_yyy = cbuffer.data(if_geom_01_off + 796 * ccomps * dcomps);

            auto g_0_z_yyyyzz_yyz = cbuffer.data(if_geom_01_off + 797 * ccomps * dcomps);

            auto g_0_z_yyyyzz_yzz = cbuffer.data(if_geom_01_off + 798 * ccomps * dcomps);

            auto g_0_z_yyyyzz_zzz = cbuffer.data(if_geom_01_off + 799 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xxx = cbuffer.data(if_geom_01_off + 800 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xxy = cbuffer.data(if_geom_01_off + 801 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xxz = cbuffer.data(if_geom_01_off + 802 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xyy = cbuffer.data(if_geom_01_off + 803 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xyz = cbuffer.data(if_geom_01_off + 804 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xzz = cbuffer.data(if_geom_01_off + 805 * ccomps * dcomps);

            auto g_0_z_yyyzzz_yyy = cbuffer.data(if_geom_01_off + 806 * ccomps * dcomps);

            auto g_0_z_yyyzzz_yyz = cbuffer.data(if_geom_01_off + 807 * ccomps * dcomps);

            auto g_0_z_yyyzzz_yzz = cbuffer.data(if_geom_01_off + 808 * ccomps * dcomps);

            auto g_0_z_yyyzzz_zzz = cbuffer.data(if_geom_01_off + 809 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xxx = cbuffer.data(if_geom_01_off + 810 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xxy = cbuffer.data(if_geom_01_off + 811 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xxz = cbuffer.data(if_geom_01_off + 812 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xyy = cbuffer.data(if_geom_01_off + 813 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xyz = cbuffer.data(if_geom_01_off + 814 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xzz = cbuffer.data(if_geom_01_off + 815 * ccomps * dcomps);

            auto g_0_z_yyzzzz_yyy = cbuffer.data(if_geom_01_off + 816 * ccomps * dcomps);

            auto g_0_z_yyzzzz_yyz = cbuffer.data(if_geom_01_off + 817 * ccomps * dcomps);

            auto g_0_z_yyzzzz_yzz = cbuffer.data(if_geom_01_off + 818 * ccomps * dcomps);

            auto g_0_z_yyzzzz_zzz = cbuffer.data(if_geom_01_off + 819 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xxx = cbuffer.data(if_geom_01_off + 820 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xxy = cbuffer.data(if_geom_01_off + 821 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xxz = cbuffer.data(if_geom_01_off + 822 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xyy = cbuffer.data(if_geom_01_off + 823 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xyz = cbuffer.data(if_geom_01_off + 824 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xzz = cbuffer.data(if_geom_01_off + 825 * ccomps * dcomps);

            auto g_0_z_yzzzzz_yyy = cbuffer.data(if_geom_01_off + 826 * ccomps * dcomps);

            auto g_0_z_yzzzzz_yyz = cbuffer.data(if_geom_01_off + 827 * ccomps * dcomps);

            auto g_0_z_yzzzzz_yzz = cbuffer.data(if_geom_01_off + 828 * ccomps * dcomps);

            auto g_0_z_yzzzzz_zzz = cbuffer.data(if_geom_01_off + 829 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xxx = cbuffer.data(if_geom_01_off + 830 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xxy = cbuffer.data(if_geom_01_off + 831 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xxz = cbuffer.data(if_geom_01_off + 832 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xyy = cbuffer.data(if_geom_01_off + 833 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xyz = cbuffer.data(if_geom_01_off + 834 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xzz = cbuffer.data(if_geom_01_off + 835 * ccomps * dcomps);

            auto g_0_z_zzzzzz_yyy = cbuffer.data(if_geom_01_off + 836 * ccomps * dcomps);

            auto g_0_z_zzzzzz_yyz = cbuffer.data(if_geom_01_off + 837 * ccomps * dcomps);

            auto g_0_z_zzzzzz_yzz = cbuffer.data(if_geom_01_off + 838 * ccomps * dcomps);

            auto g_0_z_zzzzzz_zzz = cbuffer.data(if_geom_01_off + 839 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_kdxx

            const auto kd_geom_01_off = idx_geom_01_kdxx + i * dcomps + j;

            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxxx_xx = cbuffer.data(kd_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxxxxx_xy = cbuffer.data(kd_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxxxxx_xz = cbuffer.data(kd_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxxxxx_yy = cbuffer.data(kd_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxxxxx_yz = cbuffer.data(kd_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxxxxx_zz = cbuffer.data(kd_geom_01_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxxx_xx, g_0_x_xxxxxx_xxx, g_0_x_xxxxxx_xxy, g_0_x_xxxxxx_xxz, g_0_x_xxxxxx_xy, g_0_x_xxxxxx_xyy, g_0_x_xxxxxx_xyz, g_0_x_xxxxxx_xz, g_0_x_xxxxxx_xzz, g_0_x_xxxxxx_yy, g_0_x_xxxxxx_yz, g_0_x_xxxxxx_zz, g_0_x_xxxxxxx_xx, g_0_x_xxxxxxx_xy, g_0_x_xxxxxxx_xz, g_0_x_xxxxxxx_yy, g_0_x_xxxxxxx_yz, g_0_x_xxxxxxx_zz, g_xxxxxx_xx, g_xxxxxx_xy, g_xxxxxx_xz, g_xxxxxx_yy, g_xxxxxx_yz, g_xxxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxxx_xx[k] = g_xxxxxx_xx[k] - g_0_x_xxxxxx_xx[k] * ab_x + g_0_x_xxxxxx_xxx[k];

                g_0_x_xxxxxxx_xy[k] = g_xxxxxx_xy[k] - g_0_x_xxxxxx_xy[k] * ab_x + g_0_x_xxxxxx_xxy[k];

                g_0_x_xxxxxxx_xz[k] = g_xxxxxx_xz[k] - g_0_x_xxxxxx_xz[k] * ab_x + g_0_x_xxxxxx_xxz[k];

                g_0_x_xxxxxxx_yy[k] = g_xxxxxx_yy[k] - g_0_x_xxxxxx_yy[k] * ab_x + g_0_x_xxxxxx_xyy[k];

                g_0_x_xxxxxxx_yz[k] = g_xxxxxx_yz[k] - g_0_x_xxxxxx_yz[k] * ab_x + g_0_x_xxxxxx_xyz[k];

                g_0_x_xxxxxxx_zz[k] = g_xxxxxx_zz[k] - g_0_x_xxxxxx_zz[k] * ab_x + g_0_x_xxxxxx_xzz[k];
            }

            /// Set up 6-12 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxxy_xx = cbuffer.data(kd_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxxxxy_xy = cbuffer.data(kd_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxxxxy_xz = cbuffer.data(kd_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxxxxxy_yy = cbuffer.data(kd_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxxxxxy_yz = cbuffer.data(kd_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxxxxy_zz = cbuffer.data(kd_geom_01_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxxx_xx, g_0_x_xxxxxx_xxy, g_0_x_xxxxxx_xy, g_0_x_xxxxxx_xyy, g_0_x_xxxxxx_xyz, g_0_x_xxxxxx_xz, g_0_x_xxxxxx_yy, g_0_x_xxxxxx_yyy, g_0_x_xxxxxx_yyz, g_0_x_xxxxxx_yz, g_0_x_xxxxxx_yzz, g_0_x_xxxxxx_zz, g_0_x_xxxxxxy_xx, g_0_x_xxxxxxy_xy, g_0_x_xxxxxxy_xz, g_0_x_xxxxxxy_yy, g_0_x_xxxxxxy_yz, g_0_x_xxxxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxxy_xx[k] = -g_0_x_xxxxxx_xx[k] * ab_y + g_0_x_xxxxxx_xxy[k];

                g_0_x_xxxxxxy_xy[k] = -g_0_x_xxxxxx_xy[k] * ab_y + g_0_x_xxxxxx_xyy[k];

                g_0_x_xxxxxxy_xz[k] = -g_0_x_xxxxxx_xz[k] * ab_y + g_0_x_xxxxxx_xyz[k];

                g_0_x_xxxxxxy_yy[k] = -g_0_x_xxxxxx_yy[k] * ab_y + g_0_x_xxxxxx_yyy[k];

                g_0_x_xxxxxxy_yz[k] = -g_0_x_xxxxxx_yz[k] * ab_y + g_0_x_xxxxxx_yyz[k];

                g_0_x_xxxxxxy_zz[k] = -g_0_x_xxxxxx_zz[k] * ab_y + g_0_x_xxxxxx_yzz[k];
            }

            /// Set up 12-18 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxxz_xx = cbuffer.data(kd_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxxxxz_xy = cbuffer.data(kd_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxxxxz_xz = cbuffer.data(kd_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxxxxxz_yy = cbuffer.data(kd_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxxxxz_yz = cbuffer.data(kd_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxxxxz_zz = cbuffer.data(kd_geom_01_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxxx_xx, g_0_x_xxxxxx_xxz, g_0_x_xxxxxx_xy, g_0_x_xxxxxx_xyz, g_0_x_xxxxxx_xz, g_0_x_xxxxxx_xzz, g_0_x_xxxxxx_yy, g_0_x_xxxxxx_yyz, g_0_x_xxxxxx_yz, g_0_x_xxxxxx_yzz, g_0_x_xxxxxx_zz, g_0_x_xxxxxx_zzz, g_0_x_xxxxxxz_xx, g_0_x_xxxxxxz_xy, g_0_x_xxxxxxz_xz, g_0_x_xxxxxxz_yy, g_0_x_xxxxxxz_yz, g_0_x_xxxxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxxz_xx[k] = -g_0_x_xxxxxx_xx[k] * ab_z + g_0_x_xxxxxx_xxz[k];

                g_0_x_xxxxxxz_xy[k] = -g_0_x_xxxxxx_xy[k] * ab_z + g_0_x_xxxxxx_xyz[k];

                g_0_x_xxxxxxz_xz[k] = -g_0_x_xxxxxx_xz[k] * ab_z + g_0_x_xxxxxx_xzz[k];

                g_0_x_xxxxxxz_yy[k] = -g_0_x_xxxxxx_yy[k] * ab_z + g_0_x_xxxxxx_yyz[k];

                g_0_x_xxxxxxz_yz[k] = -g_0_x_xxxxxx_yz[k] * ab_z + g_0_x_xxxxxx_yzz[k];

                g_0_x_xxxxxxz_zz[k] = -g_0_x_xxxxxx_zz[k] * ab_z + g_0_x_xxxxxx_zzz[k];
            }

            /// Set up 18-24 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxyy_xx = cbuffer.data(kd_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxxxxyy_xy = cbuffer.data(kd_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxxxxyy_xz = cbuffer.data(kd_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxxxxyy_yy = cbuffer.data(kd_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxxxxyy_yz = cbuffer.data(kd_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxxxxyy_zz = cbuffer.data(kd_geom_01_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxxy_xx, g_0_x_xxxxxy_xxy, g_0_x_xxxxxy_xy, g_0_x_xxxxxy_xyy, g_0_x_xxxxxy_xyz, g_0_x_xxxxxy_xz, g_0_x_xxxxxy_yy, g_0_x_xxxxxy_yyy, g_0_x_xxxxxy_yyz, g_0_x_xxxxxy_yz, g_0_x_xxxxxy_yzz, g_0_x_xxxxxy_zz, g_0_x_xxxxxyy_xx, g_0_x_xxxxxyy_xy, g_0_x_xxxxxyy_xz, g_0_x_xxxxxyy_yy, g_0_x_xxxxxyy_yz, g_0_x_xxxxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxyy_xx[k] = -g_0_x_xxxxxy_xx[k] * ab_y + g_0_x_xxxxxy_xxy[k];

                g_0_x_xxxxxyy_xy[k] = -g_0_x_xxxxxy_xy[k] * ab_y + g_0_x_xxxxxy_xyy[k];

                g_0_x_xxxxxyy_xz[k] = -g_0_x_xxxxxy_xz[k] * ab_y + g_0_x_xxxxxy_xyz[k];

                g_0_x_xxxxxyy_yy[k] = -g_0_x_xxxxxy_yy[k] * ab_y + g_0_x_xxxxxy_yyy[k];

                g_0_x_xxxxxyy_yz[k] = -g_0_x_xxxxxy_yz[k] * ab_y + g_0_x_xxxxxy_yyz[k];

                g_0_x_xxxxxyy_zz[k] = -g_0_x_xxxxxy_zz[k] * ab_y + g_0_x_xxxxxy_yzz[k];
            }

            /// Set up 24-30 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxyz_xx = cbuffer.data(kd_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxxxxyz_xy = cbuffer.data(kd_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxxxxyz_xz = cbuffer.data(kd_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxxxxyz_yy = cbuffer.data(kd_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxxxxyz_yz = cbuffer.data(kd_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxxxxyz_zz = cbuffer.data(kd_geom_01_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxxyz_xx, g_0_x_xxxxxyz_xy, g_0_x_xxxxxyz_xz, g_0_x_xxxxxyz_yy, g_0_x_xxxxxyz_yz, g_0_x_xxxxxyz_zz, g_0_x_xxxxxz_xx, g_0_x_xxxxxz_xxy, g_0_x_xxxxxz_xy, g_0_x_xxxxxz_xyy, g_0_x_xxxxxz_xyz, g_0_x_xxxxxz_xz, g_0_x_xxxxxz_yy, g_0_x_xxxxxz_yyy, g_0_x_xxxxxz_yyz, g_0_x_xxxxxz_yz, g_0_x_xxxxxz_yzz, g_0_x_xxxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxyz_xx[k] = -g_0_x_xxxxxz_xx[k] * ab_y + g_0_x_xxxxxz_xxy[k];

                g_0_x_xxxxxyz_xy[k] = -g_0_x_xxxxxz_xy[k] * ab_y + g_0_x_xxxxxz_xyy[k];

                g_0_x_xxxxxyz_xz[k] = -g_0_x_xxxxxz_xz[k] * ab_y + g_0_x_xxxxxz_xyz[k];

                g_0_x_xxxxxyz_yy[k] = -g_0_x_xxxxxz_yy[k] * ab_y + g_0_x_xxxxxz_yyy[k];

                g_0_x_xxxxxyz_yz[k] = -g_0_x_xxxxxz_yz[k] * ab_y + g_0_x_xxxxxz_yyz[k];

                g_0_x_xxxxxyz_zz[k] = -g_0_x_xxxxxz_zz[k] * ab_y + g_0_x_xxxxxz_yzz[k];
            }

            /// Set up 30-36 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxzz_xx = cbuffer.data(kd_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxxxxzz_xy = cbuffer.data(kd_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxxxxzz_xz = cbuffer.data(kd_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxxxxzz_yy = cbuffer.data(kd_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxxxxzz_yz = cbuffer.data(kd_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxxxxzz_zz = cbuffer.data(kd_geom_01_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxxz_xx, g_0_x_xxxxxz_xxz, g_0_x_xxxxxz_xy, g_0_x_xxxxxz_xyz, g_0_x_xxxxxz_xz, g_0_x_xxxxxz_xzz, g_0_x_xxxxxz_yy, g_0_x_xxxxxz_yyz, g_0_x_xxxxxz_yz, g_0_x_xxxxxz_yzz, g_0_x_xxxxxz_zz, g_0_x_xxxxxz_zzz, g_0_x_xxxxxzz_xx, g_0_x_xxxxxzz_xy, g_0_x_xxxxxzz_xz, g_0_x_xxxxxzz_yy, g_0_x_xxxxxzz_yz, g_0_x_xxxxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxzz_xx[k] = -g_0_x_xxxxxz_xx[k] * ab_z + g_0_x_xxxxxz_xxz[k];

                g_0_x_xxxxxzz_xy[k] = -g_0_x_xxxxxz_xy[k] * ab_z + g_0_x_xxxxxz_xyz[k];

                g_0_x_xxxxxzz_xz[k] = -g_0_x_xxxxxz_xz[k] * ab_z + g_0_x_xxxxxz_xzz[k];

                g_0_x_xxxxxzz_yy[k] = -g_0_x_xxxxxz_yy[k] * ab_z + g_0_x_xxxxxz_yyz[k];

                g_0_x_xxxxxzz_yz[k] = -g_0_x_xxxxxz_yz[k] * ab_z + g_0_x_xxxxxz_yzz[k];

                g_0_x_xxxxxzz_zz[k] = -g_0_x_xxxxxz_zz[k] * ab_z + g_0_x_xxxxxz_zzz[k];
            }

            /// Set up 36-42 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxyyy_xx = cbuffer.data(kd_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxxxyyy_xy = cbuffer.data(kd_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxxxyyy_xz = cbuffer.data(kd_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxxxyyy_yy = cbuffer.data(kd_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxxxyyy_yz = cbuffer.data(kd_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxxxyyy_zz = cbuffer.data(kd_geom_01_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxyy_xx, g_0_x_xxxxyy_xxy, g_0_x_xxxxyy_xy, g_0_x_xxxxyy_xyy, g_0_x_xxxxyy_xyz, g_0_x_xxxxyy_xz, g_0_x_xxxxyy_yy, g_0_x_xxxxyy_yyy, g_0_x_xxxxyy_yyz, g_0_x_xxxxyy_yz, g_0_x_xxxxyy_yzz, g_0_x_xxxxyy_zz, g_0_x_xxxxyyy_xx, g_0_x_xxxxyyy_xy, g_0_x_xxxxyyy_xz, g_0_x_xxxxyyy_yy, g_0_x_xxxxyyy_yz, g_0_x_xxxxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxyyy_xx[k] = -g_0_x_xxxxyy_xx[k] * ab_y + g_0_x_xxxxyy_xxy[k];

                g_0_x_xxxxyyy_xy[k] = -g_0_x_xxxxyy_xy[k] * ab_y + g_0_x_xxxxyy_xyy[k];

                g_0_x_xxxxyyy_xz[k] = -g_0_x_xxxxyy_xz[k] * ab_y + g_0_x_xxxxyy_xyz[k];

                g_0_x_xxxxyyy_yy[k] = -g_0_x_xxxxyy_yy[k] * ab_y + g_0_x_xxxxyy_yyy[k];

                g_0_x_xxxxyyy_yz[k] = -g_0_x_xxxxyy_yz[k] * ab_y + g_0_x_xxxxyy_yyz[k];

                g_0_x_xxxxyyy_zz[k] = -g_0_x_xxxxyy_zz[k] * ab_y + g_0_x_xxxxyy_yzz[k];
            }

            /// Set up 42-48 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxyyz_xx = cbuffer.data(kd_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxxxyyz_xy = cbuffer.data(kd_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxxxyyz_xz = cbuffer.data(kd_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xxxxyyz_yy = cbuffer.data(kd_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xxxxyyz_yz = cbuffer.data(kd_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xxxxyyz_zz = cbuffer.data(kd_geom_01_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxyyz_xx, g_0_x_xxxxyyz_xy, g_0_x_xxxxyyz_xz, g_0_x_xxxxyyz_yy, g_0_x_xxxxyyz_yz, g_0_x_xxxxyyz_zz, g_0_x_xxxxyz_xx, g_0_x_xxxxyz_xxy, g_0_x_xxxxyz_xy, g_0_x_xxxxyz_xyy, g_0_x_xxxxyz_xyz, g_0_x_xxxxyz_xz, g_0_x_xxxxyz_yy, g_0_x_xxxxyz_yyy, g_0_x_xxxxyz_yyz, g_0_x_xxxxyz_yz, g_0_x_xxxxyz_yzz, g_0_x_xxxxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxyyz_xx[k] = -g_0_x_xxxxyz_xx[k] * ab_y + g_0_x_xxxxyz_xxy[k];

                g_0_x_xxxxyyz_xy[k] = -g_0_x_xxxxyz_xy[k] * ab_y + g_0_x_xxxxyz_xyy[k];

                g_0_x_xxxxyyz_xz[k] = -g_0_x_xxxxyz_xz[k] * ab_y + g_0_x_xxxxyz_xyz[k];

                g_0_x_xxxxyyz_yy[k] = -g_0_x_xxxxyz_yy[k] * ab_y + g_0_x_xxxxyz_yyy[k];

                g_0_x_xxxxyyz_yz[k] = -g_0_x_xxxxyz_yz[k] * ab_y + g_0_x_xxxxyz_yyz[k];

                g_0_x_xxxxyyz_zz[k] = -g_0_x_xxxxyz_zz[k] * ab_y + g_0_x_xxxxyz_yzz[k];
            }

            /// Set up 48-54 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxyzz_xx = cbuffer.data(kd_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xxxxyzz_xy = cbuffer.data(kd_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xxxxyzz_xz = cbuffer.data(kd_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xxxxyzz_yy = cbuffer.data(kd_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xxxxyzz_yz = cbuffer.data(kd_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xxxxyzz_zz = cbuffer.data(kd_geom_01_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxyzz_xx, g_0_x_xxxxyzz_xy, g_0_x_xxxxyzz_xz, g_0_x_xxxxyzz_yy, g_0_x_xxxxyzz_yz, g_0_x_xxxxyzz_zz, g_0_x_xxxxzz_xx, g_0_x_xxxxzz_xxy, g_0_x_xxxxzz_xy, g_0_x_xxxxzz_xyy, g_0_x_xxxxzz_xyz, g_0_x_xxxxzz_xz, g_0_x_xxxxzz_yy, g_0_x_xxxxzz_yyy, g_0_x_xxxxzz_yyz, g_0_x_xxxxzz_yz, g_0_x_xxxxzz_yzz, g_0_x_xxxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxyzz_xx[k] = -g_0_x_xxxxzz_xx[k] * ab_y + g_0_x_xxxxzz_xxy[k];

                g_0_x_xxxxyzz_xy[k] = -g_0_x_xxxxzz_xy[k] * ab_y + g_0_x_xxxxzz_xyy[k];

                g_0_x_xxxxyzz_xz[k] = -g_0_x_xxxxzz_xz[k] * ab_y + g_0_x_xxxxzz_xyz[k];

                g_0_x_xxxxyzz_yy[k] = -g_0_x_xxxxzz_yy[k] * ab_y + g_0_x_xxxxzz_yyy[k];

                g_0_x_xxxxyzz_yz[k] = -g_0_x_xxxxzz_yz[k] * ab_y + g_0_x_xxxxzz_yyz[k];

                g_0_x_xxxxyzz_zz[k] = -g_0_x_xxxxzz_zz[k] * ab_y + g_0_x_xxxxzz_yzz[k];
            }

            /// Set up 54-60 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxzzz_xx = cbuffer.data(kd_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xxxxzzz_xy = cbuffer.data(kd_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xxxxzzz_xz = cbuffer.data(kd_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xxxxzzz_yy = cbuffer.data(kd_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xxxxzzz_yz = cbuffer.data(kd_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xxxxzzz_zz = cbuffer.data(kd_geom_01_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxzz_xx, g_0_x_xxxxzz_xxz, g_0_x_xxxxzz_xy, g_0_x_xxxxzz_xyz, g_0_x_xxxxzz_xz, g_0_x_xxxxzz_xzz, g_0_x_xxxxzz_yy, g_0_x_xxxxzz_yyz, g_0_x_xxxxzz_yz, g_0_x_xxxxzz_yzz, g_0_x_xxxxzz_zz, g_0_x_xxxxzz_zzz, g_0_x_xxxxzzz_xx, g_0_x_xxxxzzz_xy, g_0_x_xxxxzzz_xz, g_0_x_xxxxzzz_yy, g_0_x_xxxxzzz_yz, g_0_x_xxxxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxzzz_xx[k] = -g_0_x_xxxxzz_xx[k] * ab_z + g_0_x_xxxxzz_xxz[k];

                g_0_x_xxxxzzz_xy[k] = -g_0_x_xxxxzz_xy[k] * ab_z + g_0_x_xxxxzz_xyz[k];

                g_0_x_xxxxzzz_xz[k] = -g_0_x_xxxxzz_xz[k] * ab_z + g_0_x_xxxxzz_xzz[k];

                g_0_x_xxxxzzz_yy[k] = -g_0_x_xxxxzz_yy[k] * ab_z + g_0_x_xxxxzz_yyz[k];

                g_0_x_xxxxzzz_yz[k] = -g_0_x_xxxxzz_yz[k] * ab_z + g_0_x_xxxxzz_yzz[k];

                g_0_x_xxxxzzz_zz[k] = -g_0_x_xxxxzz_zz[k] * ab_z + g_0_x_xxxxzz_zzz[k];
            }

            /// Set up 60-66 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyyyy_xx = cbuffer.data(kd_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xxxyyyy_xy = cbuffer.data(kd_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xxxyyyy_xz = cbuffer.data(kd_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xxxyyyy_yy = cbuffer.data(kd_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xxxyyyy_yz = cbuffer.data(kd_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xxxyyyy_zz = cbuffer.data(kd_geom_01_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyyy_xx, g_0_x_xxxyyy_xxy, g_0_x_xxxyyy_xy, g_0_x_xxxyyy_xyy, g_0_x_xxxyyy_xyz, g_0_x_xxxyyy_xz, g_0_x_xxxyyy_yy, g_0_x_xxxyyy_yyy, g_0_x_xxxyyy_yyz, g_0_x_xxxyyy_yz, g_0_x_xxxyyy_yzz, g_0_x_xxxyyy_zz, g_0_x_xxxyyyy_xx, g_0_x_xxxyyyy_xy, g_0_x_xxxyyyy_xz, g_0_x_xxxyyyy_yy, g_0_x_xxxyyyy_yz, g_0_x_xxxyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyyyy_xx[k] = -g_0_x_xxxyyy_xx[k] * ab_y + g_0_x_xxxyyy_xxy[k];

                g_0_x_xxxyyyy_xy[k] = -g_0_x_xxxyyy_xy[k] * ab_y + g_0_x_xxxyyy_xyy[k];

                g_0_x_xxxyyyy_xz[k] = -g_0_x_xxxyyy_xz[k] * ab_y + g_0_x_xxxyyy_xyz[k];

                g_0_x_xxxyyyy_yy[k] = -g_0_x_xxxyyy_yy[k] * ab_y + g_0_x_xxxyyy_yyy[k];

                g_0_x_xxxyyyy_yz[k] = -g_0_x_xxxyyy_yz[k] * ab_y + g_0_x_xxxyyy_yyz[k];

                g_0_x_xxxyyyy_zz[k] = -g_0_x_xxxyyy_zz[k] * ab_y + g_0_x_xxxyyy_yzz[k];
            }

            /// Set up 66-72 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyyyz_xx = cbuffer.data(kd_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xxxyyyz_xy = cbuffer.data(kd_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xxxyyyz_xz = cbuffer.data(kd_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xxxyyyz_yy = cbuffer.data(kd_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xxxyyyz_yz = cbuffer.data(kd_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xxxyyyz_zz = cbuffer.data(kd_geom_01_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyyyz_xx, g_0_x_xxxyyyz_xy, g_0_x_xxxyyyz_xz, g_0_x_xxxyyyz_yy, g_0_x_xxxyyyz_yz, g_0_x_xxxyyyz_zz, g_0_x_xxxyyz_xx, g_0_x_xxxyyz_xxy, g_0_x_xxxyyz_xy, g_0_x_xxxyyz_xyy, g_0_x_xxxyyz_xyz, g_0_x_xxxyyz_xz, g_0_x_xxxyyz_yy, g_0_x_xxxyyz_yyy, g_0_x_xxxyyz_yyz, g_0_x_xxxyyz_yz, g_0_x_xxxyyz_yzz, g_0_x_xxxyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyyyz_xx[k] = -g_0_x_xxxyyz_xx[k] * ab_y + g_0_x_xxxyyz_xxy[k];

                g_0_x_xxxyyyz_xy[k] = -g_0_x_xxxyyz_xy[k] * ab_y + g_0_x_xxxyyz_xyy[k];

                g_0_x_xxxyyyz_xz[k] = -g_0_x_xxxyyz_xz[k] * ab_y + g_0_x_xxxyyz_xyz[k];

                g_0_x_xxxyyyz_yy[k] = -g_0_x_xxxyyz_yy[k] * ab_y + g_0_x_xxxyyz_yyy[k];

                g_0_x_xxxyyyz_yz[k] = -g_0_x_xxxyyz_yz[k] * ab_y + g_0_x_xxxyyz_yyz[k];

                g_0_x_xxxyyyz_zz[k] = -g_0_x_xxxyyz_zz[k] * ab_y + g_0_x_xxxyyz_yzz[k];
            }

            /// Set up 72-78 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyyzz_xx = cbuffer.data(kd_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xxxyyzz_xy = cbuffer.data(kd_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xxxyyzz_xz = cbuffer.data(kd_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xxxyyzz_yy = cbuffer.data(kd_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xxxyyzz_yz = cbuffer.data(kd_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xxxyyzz_zz = cbuffer.data(kd_geom_01_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyyzz_xx, g_0_x_xxxyyzz_xy, g_0_x_xxxyyzz_xz, g_0_x_xxxyyzz_yy, g_0_x_xxxyyzz_yz, g_0_x_xxxyyzz_zz, g_0_x_xxxyzz_xx, g_0_x_xxxyzz_xxy, g_0_x_xxxyzz_xy, g_0_x_xxxyzz_xyy, g_0_x_xxxyzz_xyz, g_0_x_xxxyzz_xz, g_0_x_xxxyzz_yy, g_0_x_xxxyzz_yyy, g_0_x_xxxyzz_yyz, g_0_x_xxxyzz_yz, g_0_x_xxxyzz_yzz, g_0_x_xxxyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyyzz_xx[k] = -g_0_x_xxxyzz_xx[k] * ab_y + g_0_x_xxxyzz_xxy[k];

                g_0_x_xxxyyzz_xy[k] = -g_0_x_xxxyzz_xy[k] * ab_y + g_0_x_xxxyzz_xyy[k];

                g_0_x_xxxyyzz_xz[k] = -g_0_x_xxxyzz_xz[k] * ab_y + g_0_x_xxxyzz_xyz[k];

                g_0_x_xxxyyzz_yy[k] = -g_0_x_xxxyzz_yy[k] * ab_y + g_0_x_xxxyzz_yyy[k];

                g_0_x_xxxyyzz_yz[k] = -g_0_x_xxxyzz_yz[k] * ab_y + g_0_x_xxxyzz_yyz[k];

                g_0_x_xxxyyzz_zz[k] = -g_0_x_xxxyzz_zz[k] * ab_y + g_0_x_xxxyzz_yzz[k];
            }

            /// Set up 78-84 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyzzz_xx = cbuffer.data(kd_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xxxyzzz_xy = cbuffer.data(kd_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xxxyzzz_xz = cbuffer.data(kd_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xxxyzzz_yy = cbuffer.data(kd_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xxxyzzz_yz = cbuffer.data(kd_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xxxyzzz_zz = cbuffer.data(kd_geom_01_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyzzz_xx, g_0_x_xxxyzzz_xy, g_0_x_xxxyzzz_xz, g_0_x_xxxyzzz_yy, g_0_x_xxxyzzz_yz, g_0_x_xxxyzzz_zz, g_0_x_xxxzzz_xx, g_0_x_xxxzzz_xxy, g_0_x_xxxzzz_xy, g_0_x_xxxzzz_xyy, g_0_x_xxxzzz_xyz, g_0_x_xxxzzz_xz, g_0_x_xxxzzz_yy, g_0_x_xxxzzz_yyy, g_0_x_xxxzzz_yyz, g_0_x_xxxzzz_yz, g_0_x_xxxzzz_yzz, g_0_x_xxxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyzzz_xx[k] = -g_0_x_xxxzzz_xx[k] * ab_y + g_0_x_xxxzzz_xxy[k];

                g_0_x_xxxyzzz_xy[k] = -g_0_x_xxxzzz_xy[k] * ab_y + g_0_x_xxxzzz_xyy[k];

                g_0_x_xxxyzzz_xz[k] = -g_0_x_xxxzzz_xz[k] * ab_y + g_0_x_xxxzzz_xyz[k];

                g_0_x_xxxyzzz_yy[k] = -g_0_x_xxxzzz_yy[k] * ab_y + g_0_x_xxxzzz_yyy[k];

                g_0_x_xxxyzzz_yz[k] = -g_0_x_xxxzzz_yz[k] * ab_y + g_0_x_xxxzzz_yyz[k];

                g_0_x_xxxyzzz_zz[k] = -g_0_x_xxxzzz_zz[k] * ab_y + g_0_x_xxxzzz_yzz[k];
            }

            /// Set up 84-90 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxzzzz_xx = cbuffer.data(kd_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xxxzzzz_xy = cbuffer.data(kd_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xxxzzzz_xz = cbuffer.data(kd_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xxxzzzz_yy = cbuffer.data(kd_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xxxzzzz_yz = cbuffer.data(kd_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xxxzzzz_zz = cbuffer.data(kd_geom_01_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxzzz_xx, g_0_x_xxxzzz_xxz, g_0_x_xxxzzz_xy, g_0_x_xxxzzz_xyz, g_0_x_xxxzzz_xz, g_0_x_xxxzzz_xzz, g_0_x_xxxzzz_yy, g_0_x_xxxzzz_yyz, g_0_x_xxxzzz_yz, g_0_x_xxxzzz_yzz, g_0_x_xxxzzz_zz, g_0_x_xxxzzz_zzz, g_0_x_xxxzzzz_xx, g_0_x_xxxzzzz_xy, g_0_x_xxxzzzz_xz, g_0_x_xxxzzzz_yy, g_0_x_xxxzzzz_yz, g_0_x_xxxzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxzzzz_xx[k] = -g_0_x_xxxzzz_xx[k] * ab_z + g_0_x_xxxzzz_xxz[k];

                g_0_x_xxxzzzz_xy[k] = -g_0_x_xxxzzz_xy[k] * ab_z + g_0_x_xxxzzz_xyz[k];

                g_0_x_xxxzzzz_xz[k] = -g_0_x_xxxzzz_xz[k] * ab_z + g_0_x_xxxzzz_xzz[k];

                g_0_x_xxxzzzz_yy[k] = -g_0_x_xxxzzz_yy[k] * ab_z + g_0_x_xxxzzz_yyz[k];

                g_0_x_xxxzzzz_yz[k] = -g_0_x_xxxzzz_yz[k] * ab_z + g_0_x_xxxzzz_yzz[k];

                g_0_x_xxxzzzz_zz[k] = -g_0_x_xxxzzz_zz[k] * ab_z + g_0_x_xxxzzz_zzz[k];
            }

            /// Set up 90-96 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyyyy_xx = cbuffer.data(kd_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_xxyyyyy_xy = cbuffer.data(kd_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_xxyyyyy_xz = cbuffer.data(kd_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_xxyyyyy_yy = cbuffer.data(kd_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_xxyyyyy_yz = cbuffer.data(kd_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_xxyyyyy_zz = cbuffer.data(kd_geom_01_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyyy_xx, g_0_x_xxyyyy_xxy, g_0_x_xxyyyy_xy, g_0_x_xxyyyy_xyy, g_0_x_xxyyyy_xyz, g_0_x_xxyyyy_xz, g_0_x_xxyyyy_yy, g_0_x_xxyyyy_yyy, g_0_x_xxyyyy_yyz, g_0_x_xxyyyy_yz, g_0_x_xxyyyy_yzz, g_0_x_xxyyyy_zz, g_0_x_xxyyyyy_xx, g_0_x_xxyyyyy_xy, g_0_x_xxyyyyy_xz, g_0_x_xxyyyyy_yy, g_0_x_xxyyyyy_yz, g_0_x_xxyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyyyy_xx[k] = -g_0_x_xxyyyy_xx[k] * ab_y + g_0_x_xxyyyy_xxy[k];

                g_0_x_xxyyyyy_xy[k] = -g_0_x_xxyyyy_xy[k] * ab_y + g_0_x_xxyyyy_xyy[k];

                g_0_x_xxyyyyy_xz[k] = -g_0_x_xxyyyy_xz[k] * ab_y + g_0_x_xxyyyy_xyz[k];

                g_0_x_xxyyyyy_yy[k] = -g_0_x_xxyyyy_yy[k] * ab_y + g_0_x_xxyyyy_yyy[k];

                g_0_x_xxyyyyy_yz[k] = -g_0_x_xxyyyy_yz[k] * ab_y + g_0_x_xxyyyy_yyz[k];

                g_0_x_xxyyyyy_zz[k] = -g_0_x_xxyyyy_zz[k] * ab_y + g_0_x_xxyyyy_yzz[k];
            }

            /// Set up 96-102 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyyyz_xx = cbuffer.data(kd_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_xxyyyyz_xy = cbuffer.data(kd_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_xxyyyyz_xz = cbuffer.data(kd_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_xxyyyyz_yy = cbuffer.data(kd_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_xxyyyyz_yz = cbuffer.data(kd_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_xxyyyyz_zz = cbuffer.data(kd_geom_01_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyyyz_xx, g_0_x_xxyyyyz_xy, g_0_x_xxyyyyz_xz, g_0_x_xxyyyyz_yy, g_0_x_xxyyyyz_yz, g_0_x_xxyyyyz_zz, g_0_x_xxyyyz_xx, g_0_x_xxyyyz_xxy, g_0_x_xxyyyz_xy, g_0_x_xxyyyz_xyy, g_0_x_xxyyyz_xyz, g_0_x_xxyyyz_xz, g_0_x_xxyyyz_yy, g_0_x_xxyyyz_yyy, g_0_x_xxyyyz_yyz, g_0_x_xxyyyz_yz, g_0_x_xxyyyz_yzz, g_0_x_xxyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyyyz_xx[k] = -g_0_x_xxyyyz_xx[k] * ab_y + g_0_x_xxyyyz_xxy[k];

                g_0_x_xxyyyyz_xy[k] = -g_0_x_xxyyyz_xy[k] * ab_y + g_0_x_xxyyyz_xyy[k];

                g_0_x_xxyyyyz_xz[k] = -g_0_x_xxyyyz_xz[k] * ab_y + g_0_x_xxyyyz_xyz[k];

                g_0_x_xxyyyyz_yy[k] = -g_0_x_xxyyyz_yy[k] * ab_y + g_0_x_xxyyyz_yyy[k];

                g_0_x_xxyyyyz_yz[k] = -g_0_x_xxyyyz_yz[k] * ab_y + g_0_x_xxyyyz_yyz[k];

                g_0_x_xxyyyyz_zz[k] = -g_0_x_xxyyyz_zz[k] * ab_y + g_0_x_xxyyyz_yzz[k];
            }

            /// Set up 102-108 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyyzz_xx = cbuffer.data(kd_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_xxyyyzz_xy = cbuffer.data(kd_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_xxyyyzz_xz = cbuffer.data(kd_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_xxyyyzz_yy = cbuffer.data(kd_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_xxyyyzz_yz = cbuffer.data(kd_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_xxyyyzz_zz = cbuffer.data(kd_geom_01_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyyzz_xx, g_0_x_xxyyyzz_xy, g_0_x_xxyyyzz_xz, g_0_x_xxyyyzz_yy, g_0_x_xxyyyzz_yz, g_0_x_xxyyyzz_zz, g_0_x_xxyyzz_xx, g_0_x_xxyyzz_xxy, g_0_x_xxyyzz_xy, g_0_x_xxyyzz_xyy, g_0_x_xxyyzz_xyz, g_0_x_xxyyzz_xz, g_0_x_xxyyzz_yy, g_0_x_xxyyzz_yyy, g_0_x_xxyyzz_yyz, g_0_x_xxyyzz_yz, g_0_x_xxyyzz_yzz, g_0_x_xxyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyyzz_xx[k] = -g_0_x_xxyyzz_xx[k] * ab_y + g_0_x_xxyyzz_xxy[k];

                g_0_x_xxyyyzz_xy[k] = -g_0_x_xxyyzz_xy[k] * ab_y + g_0_x_xxyyzz_xyy[k];

                g_0_x_xxyyyzz_xz[k] = -g_0_x_xxyyzz_xz[k] * ab_y + g_0_x_xxyyzz_xyz[k];

                g_0_x_xxyyyzz_yy[k] = -g_0_x_xxyyzz_yy[k] * ab_y + g_0_x_xxyyzz_yyy[k];

                g_0_x_xxyyyzz_yz[k] = -g_0_x_xxyyzz_yz[k] * ab_y + g_0_x_xxyyzz_yyz[k];

                g_0_x_xxyyyzz_zz[k] = -g_0_x_xxyyzz_zz[k] * ab_y + g_0_x_xxyyzz_yzz[k];
            }

            /// Set up 108-114 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyzzz_xx = cbuffer.data(kd_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_xxyyzzz_xy = cbuffer.data(kd_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_xxyyzzz_xz = cbuffer.data(kd_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_xxyyzzz_yy = cbuffer.data(kd_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_xxyyzzz_yz = cbuffer.data(kd_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_xxyyzzz_zz = cbuffer.data(kd_geom_01_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyzzz_xx, g_0_x_xxyyzzz_xy, g_0_x_xxyyzzz_xz, g_0_x_xxyyzzz_yy, g_0_x_xxyyzzz_yz, g_0_x_xxyyzzz_zz, g_0_x_xxyzzz_xx, g_0_x_xxyzzz_xxy, g_0_x_xxyzzz_xy, g_0_x_xxyzzz_xyy, g_0_x_xxyzzz_xyz, g_0_x_xxyzzz_xz, g_0_x_xxyzzz_yy, g_0_x_xxyzzz_yyy, g_0_x_xxyzzz_yyz, g_0_x_xxyzzz_yz, g_0_x_xxyzzz_yzz, g_0_x_xxyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyzzz_xx[k] = -g_0_x_xxyzzz_xx[k] * ab_y + g_0_x_xxyzzz_xxy[k];

                g_0_x_xxyyzzz_xy[k] = -g_0_x_xxyzzz_xy[k] * ab_y + g_0_x_xxyzzz_xyy[k];

                g_0_x_xxyyzzz_xz[k] = -g_0_x_xxyzzz_xz[k] * ab_y + g_0_x_xxyzzz_xyz[k];

                g_0_x_xxyyzzz_yy[k] = -g_0_x_xxyzzz_yy[k] * ab_y + g_0_x_xxyzzz_yyy[k];

                g_0_x_xxyyzzz_yz[k] = -g_0_x_xxyzzz_yz[k] * ab_y + g_0_x_xxyzzz_yyz[k];

                g_0_x_xxyyzzz_zz[k] = -g_0_x_xxyzzz_zz[k] * ab_y + g_0_x_xxyzzz_yzz[k];
            }

            /// Set up 114-120 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyzzzz_xx = cbuffer.data(kd_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_xxyzzzz_xy = cbuffer.data(kd_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_xxyzzzz_xz = cbuffer.data(kd_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_xxyzzzz_yy = cbuffer.data(kd_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_xxyzzzz_yz = cbuffer.data(kd_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_xxyzzzz_zz = cbuffer.data(kd_geom_01_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyzzzz_xx, g_0_x_xxyzzzz_xy, g_0_x_xxyzzzz_xz, g_0_x_xxyzzzz_yy, g_0_x_xxyzzzz_yz, g_0_x_xxyzzzz_zz, g_0_x_xxzzzz_xx, g_0_x_xxzzzz_xxy, g_0_x_xxzzzz_xy, g_0_x_xxzzzz_xyy, g_0_x_xxzzzz_xyz, g_0_x_xxzzzz_xz, g_0_x_xxzzzz_yy, g_0_x_xxzzzz_yyy, g_0_x_xxzzzz_yyz, g_0_x_xxzzzz_yz, g_0_x_xxzzzz_yzz, g_0_x_xxzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyzzzz_xx[k] = -g_0_x_xxzzzz_xx[k] * ab_y + g_0_x_xxzzzz_xxy[k];

                g_0_x_xxyzzzz_xy[k] = -g_0_x_xxzzzz_xy[k] * ab_y + g_0_x_xxzzzz_xyy[k];

                g_0_x_xxyzzzz_xz[k] = -g_0_x_xxzzzz_xz[k] * ab_y + g_0_x_xxzzzz_xyz[k];

                g_0_x_xxyzzzz_yy[k] = -g_0_x_xxzzzz_yy[k] * ab_y + g_0_x_xxzzzz_yyy[k];

                g_0_x_xxyzzzz_yz[k] = -g_0_x_xxzzzz_yz[k] * ab_y + g_0_x_xxzzzz_yyz[k];

                g_0_x_xxyzzzz_zz[k] = -g_0_x_xxzzzz_zz[k] * ab_y + g_0_x_xxzzzz_yzz[k];
            }

            /// Set up 120-126 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxzzzzz_xx = cbuffer.data(kd_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_xxzzzzz_xy = cbuffer.data(kd_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_xxzzzzz_xz = cbuffer.data(kd_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_xxzzzzz_yy = cbuffer.data(kd_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_xxzzzzz_yz = cbuffer.data(kd_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_xxzzzzz_zz = cbuffer.data(kd_geom_01_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxzzzz_xx, g_0_x_xxzzzz_xxz, g_0_x_xxzzzz_xy, g_0_x_xxzzzz_xyz, g_0_x_xxzzzz_xz, g_0_x_xxzzzz_xzz, g_0_x_xxzzzz_yy, g_0_x_xxzzzz_yyz, g_0_x_xxzzzz_yz, g_0_x_xxzzzz_yzz, g_0_x_xxzzzz_zz, g_0_x_xxzzzz_zzz, g_0_x_xxzzzzz_xx, g_0_x_xxzzzzz_xy, g_0_x_xxzzzzz_xz, g_0_x_xxzzzzz_yy, g_0_x_xxzzzzz_yz, g_0_x_xxzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxzzzzz_xx[k] = -g_0_x_xxzzzz_xx[k] * ab_z + g_0_x_xxzzzz_xxz[k];

                g_0_x_xxzzzzz_xy[k] = -g_0_x_xxzzzz_xy[k] * ab_z + g_0_x_xxzzzz_xyz[k];

                g_0_x_xxzzzzz_xz[k] = -g_0_x_xxzzzz_xz[k] * ab_z + g_0_x_xxzzzz_xzz[k];

                g_0_x_xxzzzzz_yy[k] = -g_0_x_xxzzzz_yy[k] * ab_z + g_0_x_xxzzzz_yyz[k];

                g_0_x_xxzzzzz_yz[k] = -g_0_x_xxzzzz_yz[k] * ab_z + g_0_x_xxzzzz_yzz[k];

                g_0_x_xxzzzzz_zz[k] = -g_0_x_xxzzzz_zz[k] * ab_z + g_0_x_xxzzzz_zzz[k];
            }

            /// Set up 126-132 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyyyy_xx = cbuffer.data(kd_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_xyyyyyy_xy = cbuffer.data(kd_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_xyyyyyy_xz = cbuffer.data(kd_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_xyyyyyy_yy = cbuffer.data(kd_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_xyyyyyy_yz = cbuffer.data(kd_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_xyyyyyy_zz = cbuffer.data(kd_geom_01_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyyy_xx, g_0_x_xyyyyy_xxy, g_0_x_xyyyyy_xy, g_0_x_xyyyyy_xyy, g_0_x_xyyyyy_xyz, g_0_x_xyyyyy_xz, g_0_x_xyyyyy_yy, g_0_x_xyyyyy_yyy, g_0_x_xyyyyy_yyz, g_0_x_xyyyyy_yz, g_0_x_xyyyyy_yzz, g_0_x_xyyyyy_zz, g_0_x_xyyyyyy_xx, g_0_x_xyyyyyy_xy, g_0_x_xyyyyyy_xz, g_0_x_xyyyyyy_yy, g_0_x_xyyyyyy_yz, g_0_x_xyyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyyyy_xx[k] = -g_0_x_xyyyyy_xx[k] * ab_y + g_0_x_xyyyyy_xxy[k];

                g_0_x_xyyyyyy_xy[k] = -g_0_x_xyyyyy_xy[k] * ab_y + g_0_x_xyyyyy_xyy[k];

                g_0_x_xyyyyyy_xz[k] = -g_0_x_xyyyyy_xz[k] * ab_y + g_0_x_xyyyyy_xyz[k];

                g_0_x_xyyyyyy_yy[k] = -g_0_x_xyyyyy_yy[k] * ab_y + g_0_x_xyyyyy_yyy[k];

                g_0_x_xyyyyyy_yz[k] = -g_0_x_xyyyyy_yz[k] * ab_y + g_0_x_xyyyyy_yyz[k];

                g_0_x_xyyyyyy_zz[k] = -g_0_x_xyyyyy_zz[k] * ab_y + g_0_x_xyyyyy_yzz[k];
            }

            /// Set up 132-138 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyyyz_xx = cbuffer.data(kd_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_xyyyyyz_xy = cbuffer.data(kd_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_xyyyyyz_xz = cbuffer.data(kd_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_x_xyyyyyz_yy = cbuffer.data(kd_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_xyyyyyz_yz = cbuffer.data(kd_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_xyyyyyz_zz = cbuffer.data(kd_geom_01_off + 137 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyyyz_xx, g_0_x_xyyyyyz_xy, g_0_x_xyyyyyz_xz, g_0_x_xyyyyyz_yy, g_0_x_xyyyyyz_yz, g_0_x_xyyyyyz_zz, g_0_x_xyyyyz_xx, g_0_x_xyyyyz_xxy, g_0_x_xyyyyz_xy, g_0_x_xyyyyz_xyy, g_0_x_xyyyyz_xyz, g_0_x_xyyyyz_xz, g_0_x_xyyyyz_yy, g_0_x_xyyyyz_yyy, g_0_x_xyyyyz_yyz, g_0_x_xyyyyz_yz, g_0_x_xyyyyz_yzz, g_0_x_xyyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyyyz_xx[k] = -g_0_x_xyyyyz_xx[k] * ab_y + g_0_x_xyyyyz_xxy[k];

                g_0_x_xyyyyyz_xy[k] = -g_0_x_xyyyyz_xy[k] * ab_y + g_0_x_xyyyyz_xyy[k];

                g_0_x_xyyyyyz_xz[k] = -g_0_x_xyyyyz_xz[k] * ab_y + g_0_x_xyyyyz_xyz[k];

                g_0_x_xyyyyyz_yy[k] = -g_0_x_xyyyyz_yy[k] * ab_y + g_0_x_xyyyyz_yyy[k];

                g_0_x_xyyyyyz_yz[k] = -g_0_x_xyyyyz_yz[k] * ab_y + g_0_x_xyyyyz_yyz[k];

                g_0_x_xyyyyyz_zz[k] = -g_0_x_xyyyyz_zz[k] * ab_y + g_0_x_xyyyyz_yzz[k];
            }

            /// Set up 138-144 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyyzz_xx = cbuffer.data(kd_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_xyyyyzz_xy = cbuffer.data(kd_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_x_xyyyyzz_xz = cbuffer.data(kd_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_xyyyyzz_yy = cbuffer.data(kd_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_xyyyyzz_yz = cbuffer.data(kd_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_xyyyyzz_zz = cbuffer.data(kd_geom_01_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyyzz_xx, g_0_x_xyyyyzz_xy, g_0_x_xyyyyzz_xz, g_0_x_xyyyyzz_yy, g_0_x_xyyyyzz_yz, g_0_x_xyyyyzz_zz, g_0_x_xyyyzz_xx, g_0_x_xyyyzz_xxy, g_0_x_xyyyzz_xy, g_0_x_xyyyzz_xyy, g_0_x_xyyyzz_xyz, g_0_x_xyyyzz_xz, g_0_x_xyyyzz_yy, g_0_x_xyyyzz_yyy, g_0_x_xyyyzz_yyz, g_0_x_xyyyzz_yz, g_0_x_xyyyzz_yzz, g_0_x_xyyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyyzz_xx[k] = -g_0_x_xyyyzz_xx[k] * ab_y + g_0_x_xyyyzz_xxy[k];

                g_0_x_xyyyyzz_xy[k] = -g_0_x_xyyyzz_xy[k] * ab_y + g_0_x_xyyyzz_xyy[k];

                g_0_x_xyyyyzz_xz[k] = -g_0_x_xyyyzz_xz[k] * ab_y + g_0_x_xyyyzz_xyz[k];

                g_0_x_xyyyyzz_yy[k] = -g_0_x_xyyyzz_yy[k] * ab_y + g_0_x_xyyyzz_yyy[k];

                g_0_x_xyyyyzz_yz[k] = -g_0_x_xyyyzz_yz[k] * ab_y + g_0_x_xyyyzz_yyz[k];

                g_0_x_xyyyyzz_zz[k] = -g_0_x_xyyyzz_zz[k] * ab_y + g_0_x_xyyyzz_yzz[k];
            }

            /// Set up 144-150 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyzzz_xx = cbuffer.data(kd_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_xyyyzzz_xy = cbuffer.data(kd_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_xyyyzzz_xz = cbuffer.data(kd_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_xyyyzzz_yy = cbuffer.data(kd_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_xyyyzzz_yz = cbuffer.data(kd_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_xyyyzzz_zz = cbuffer.data(kd_geom_01_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyzzz_xx, g_0_x_xyyyzzz_xy, g_0_x_xyyyzzz_xz, g_0_x_xyyyzzz_yy, g_0_x_xyyyzzz_yz, g_0_x_xyyyzzz_zz, g_0_x_xyyzzz_xx, g_0_x_xyyzzz_xxy, g_0_x_xyyzzz_xy, g_0_x_xyyzzz_xyy, g_0_x_xyyzzz_xyz, g_0_x_xyyzzz_xz, g_0_x_xyyzzz_yy, g_0_x_xyyzzz_yyy, g_0_x_xyyzzz_yyz, g_0_x_xyyzzz_yz, g_0_x_xyyzzz_yzz, g_0_x_xyyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyzzz_xx[k] = -g_0_x_xyyzzz_xx[k] * ab_y + g_0_x_xyyzzz_xxy[k];

                g_0_x_xyyyzzz_xy[k] = -g_0_x_xyyzzz_xy[k] * ab_y + g_0_x_xyyzzz_xyy[k];

                g_0_x_xyyyzzz_xz[k] = -g_0_x_xyyzzz_xz[k] * ab_y + g_0_x_xyyzzz_xyz[k];

                g_0_x_xyyyzzz_yy[k] = -g_0_x_xyyzzz_yy[k] * ab_y + g_0_x_xyyzzz_yyy[k];

                g_0_x_xyyyzzz_yz[k] = -g_0_x_xyyzzz_yz[k] * ab_y + g_0_x_xyyzzz_yyz[k];

                g_0_x_xyyyzzz_zz[k] = -g_0_x_xyyzzz_zz[k] * ab_y + g_0_x_xyyzzz_yzz[k];
            }

            /// Set up 150-156 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyzzzz_xx = cbuffer.data(kd_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_x_xyyzzzz_xy = cbuffer.data(kd_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_x_xyyzzzz_xz = cbuffer.data(kd_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_x_xyyzzzz_yy = cbuffer.data(kd_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_x_xyyzzzz_yz = cbuffer.data(kd_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_x_xyyzzzz_zz = cbuffer.data(kd_geom_01_off + 155 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyzzzz_xx, g_0_x_xyyzzzz_xy, g_0_x_xyyzzzz_xz, g_0_x_xyyzzzz_yy, g_0_x_xyyzzzz_yz, g_0_x_xyyzzzz_zz, g_0_x_xyzzzz_xx, g_0_x_xyzzzz_xxy, g_0_x_xyzzzz_xy, g_0_x_xyzzzz_xyy, g_0_x_xyzzzz_xyz, g_0_x_xyzzzz_xz, g_0_x_xyzzzz_yy, g_0_x_xyzzzz_yyy, g_0_x_xyzzzz_yyz, g_0_x_xyzzzz_yz, g_0_x_xyzzzz_yzz, g_0_x_xyzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyzzzz_xx[k] = -g_0_x_xyzzzz_xx[k] * ab_y + g_0_x_xyzzzz_xxy[k];

                g_0_x_xyyzzzz_xy[k] = -g_0_x_xyzzzz_xy[k] * ab_y + g_0_x_xyzzzz_xyy[k];

                g_0_x_xyyzzzz_xz[k] = -g_0_x_xyzzzz_xz[k] * ab_y + g_0_x_xyzzzz_xyz[k];

                g_0_x_xyyzzzz_yy[k] = -g_0_x_xyzzzz_yy[k] * ab_y + g_0_x_xyzzzz_yyy[k];

                g_0_x_xyyzzzz_yz[k] = -g_0_x_xyzzzz_yz[k] * ab_y + g_0_x_xyzzzz_yyz[k];

                g_0_x_xyyzzzz_zz[k] = -g_0_x_xyzzzz_zz[k] * ab_y + g_0_x_xyzzzz_yzz[k];
            }

            /// Set up 156-162 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyzzzzz_xx = cbuffer.data(kd_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_x_xyzzzzz_xy = cbuffer.data(kd_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_x_xyzzzzz_xz = cbuffer.data(kd_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_x_xyzzzzz_yy = cbuffer.data(kd_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_x_xyzzzzz_yz = cbuffer.data(kd_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_x_xyzzzzz_zz = cbuffer.data(kd_geom_01_off + 161 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyzzzzz_xx, g_0_x_xyzzzzz_xy, g_0_x_xyzzzzz_xz, g_0_x_xyzzzzz_yy, g_0_x_xyzzzzz_yz, g_0_x_xyzzzzz_zz, g_0_x_xzzzzz_xx, g_0_x_xzzzzz_xxy, g_0_x_xzzzzz_xy, g_0_x_xzzzzz_xyy, g_0_x_xzzzzz_xyz, g_0_x_xzzzzz_xz, g_0_x_xzzzzz_yy, g_0_x_xzzzzz_yyy, g_0_x_xzzzzz_yyz, g_0_x_xzzzzz_yz, g_0_x_xzzzzz_yzz, g_0_x_xzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyzzzzz_xx[k] = -g_0_x_xzzzzz_xx[k] * ab_y + g_0_x_xzzzzz_xxy[k];

                g_0_x_xyzzzzz_xy[k] = -g_0_x_xzzzzz_xy[k] * ab_y + g_0_x_xzzzzz_xyy[k];

                g_0_x_xyzzzzz_xz[k] = -g_0_x_xzzzzz_xz[k] * ab_y + g_0_x_xzzzzz_xyz[k];

                g_0_x_xyzzzzz_yy[k] = -g_0_x_xzzzzz_yy[k] * ab_y + g_0_x_xzzzzz_yyy[k];

                g_0_x_xyzzzzz_yz[k] = -g_0_x_xzzzzz_yz[k] * ab_y + g_0_x_xzzzzz_yyz[k];

                g_0_x_xyzzzzz_zz[k] = -g_0_x_xzzzzz_zz[k] * ab_y + g_0_x_xzzzzz_yzz[k];
            }

            /// Set up 162-168 components of targeted buffer : cbuffer.data(

            auto g_0_x_xzzzzzz_xx = cbuffer.data(kd_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_x_xzzzzzz_xy = cbuffer.data(kd_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_x_xzzzzzz_xz = cbuffer.data(kd_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_x_xzzzzzz_yy = cbuffer.data(kd_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_x_xzzzzzz_yz = cbuffer.data(kd_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_x_xzzzzzz_zz = cbuffer.data(kd_geom_01_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xzzzzz_xx, g_0_x_xzzzzz_xxz, g_0_x_xzzzzz_xy, g_0_x_xzzzzz_xyz, g_0_x_xzzzzz_xz, g_0_x_xzzzzz_xzz, g_0_x_xzzzzz_yy, g_0_x_xzzzzz_yyz, g_0_x_xzzzzz_yz, g_0_x_xzzzzz_yzz, g_0_x_xzzzzz_zz, g_0_x_xzzzzz_zzz, g_0_x_xzzzzzz_xx, g_0_x_xzzzzzz_xy, g_0_x_xzzzzzz_xz, g_0_x_xzzzzzz_yy, g_0_x_xzzzzzz_yz, g_0_x_xzzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xzzzzzz_xx[k] = -g_0_x_xzzzzz_xx[k] * ab_z + g_0_x_xzzzzz_xxz[k];

                g_0_x_xzzzzzz_xy[k] = -g_0_x_xzzzzz_xy[k] * ab_z + g_0_x_xzzzzz_xyz[k];

                g_0_x_xzzzzzz_xz[k] = -g_0_x_xzzzzz_xz[k] * ab_z + g_0_x_xzzzzz_xzz[k];

                g_0_x_xzzzzzz_yy[k] = -g_0_x_xzzzzz_yy[k] * ab_z + g_0_x_xzzzzz_yyz[k];

                g_0_x_xzzzzzz_yz[k] = -g_0_x_xzzzzz_yz[k] * ab_z + g_0_x_xzzzzz_yzz[k];

                g_0_x_xzzzzzz_zz[k] = -g_0_x_xzzzzz_zz[k] * ab_z + g_0_x_xzzzzz_zzz[k];
            }

            /// Set up 168-174 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyyyy_xx = cbuffer.data(kd_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_x_yyyyyyy_xy = cbuffer.data(kd_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_x_yyyyyyy_xz = cbuffer.data(kd_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_x_yyyyyyy_yy = cbuffer.data(kd_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_x_yyyyyyy_yz = cbuffer.data(kd_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_x_yyyyyyy_zz = cbuffer.data(kd_geom_01_off + 173 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyyy_xx, g_0_x_yyyyyy_xxy, g_0_x_yyyyyy_xy, g_0_x_yyyyyy_xyy, g_0_x_yyyyyy_xyz, g_0_x_yyyyyy_xz, g_0_x_yyyyyy_yy, g_0_x_yyyyyy_yyy, g_0_x_yyyyyy_yyz, g_0_x_yyyyyy_yz, g_0_x_yyyyyy_yzz, g_0_x_yyyyyy_zz, g_0_x_yyyyyyy_xx, g_0_x_yyyyyyy_xy, g_0_x_yyyyyyy_xz, g_0_x_yyyyyyy_yy, g_0_x_yyyyyyy_yz, g_0_x_yyyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyyyy_xx[k] = -g_0_x_yyyyyy_xx[k] * ab_y + g_0_x_yyyyyy_xxy[k];

                g_0_x_yyyyyyy_xy[k] = -g_0_x_yyyyyy_xy[k] * ab_y + g_0_x_yyyyyy_xyy[k];

                g_0_x_yyyyyyy_xz[k] = -g_0_x_yyyyyy_xz[k] * ab_y + g_0_x_yyyyyy_xyz[k];

                g_0_x_yyyyyyy_yy[k] = -g_0_x_yyyyyy_yy[k] * ab_y + g_0_x_yyyyyy_yyy[k];

                g_0_x_yyyyyyy_yz[k] = -g_0_x_yyyyyy_yz[k] * ab_y + g_0_x_yyyyyy_yyz[k];

                g_0_x_yyyyyyy_zz[k] = -g_0_x_yyyyyy_zz[k] * ab_y + g_0_x_yyyyyy_yzz[k];
            }

            /// Set up 174-180 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyyyz_xx = cbuffer.data(kd_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_x_yyyyyyz_xy = cbuffer.data(kd_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_x_yyyyyyz_xz = cbuffer.data(kd_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_x_yyyyyyz_yy = cbuffer.data(kd_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_x_yyyyyyz_yz = cbuffer.data(kd_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_x_yyyyyyz_zz = cbuffer.data(kd_geom_01_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyyyz_xx, g_0_x_yyyyyyz_xy, g_0_x_yyyyyyz_xz, g_0_x_yyyyyyz_yy, g_0_x_yyyyyyz_yz, g_0_x_yyyyyyz_zz, g_0_x_yyyyyz_xx, g_0_x_yyyyyz_xxy, g_0_x_yyyyyz_xy, g_0_x_yyyyyz_xyy, g_0_x_yyyyyz_xyz, g_0_x_yyyyyz_xz, g_0_x_yyyyyz_yy, g_0_x_yyyyyz_yyy, g_0_x_yyyyyz_yyz, g_0_x_yyyyyz_yz, g_0_x_yyyyyz_yzz, g_0_x_yyyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyyyz_xx[k] = -g_0_x_yyyyyz_xx[k] * ab_y + g_0_x_yyyyyz_xxy[k];

                g_0_x_yyyyyyz_xy[k] = -g_0_x_yyyyyz_xy[k] * ab_y + g_0_x_yyyyyz_xyy[k];

                g_0_x_yyyyyyz_xz[k] = -g_0_x_yyyyyz_xz[k] * ab_y + g_0_x_yyyyyz_xyz[k];

                g_0_x_yyyyyyz_yy[k] = -g_0_x_yyyyyz_yy[k] * ab_y + g_0_x_yyyyyz_yyy[k];

                g_0_x_yyyyyyz_yz[k] = -g_0_x_yyyyyz_yz[k] * ab_y + g_0_x_yyyyyz_yyz[k];

                g_0_x_yyyyyyz_zz[k] = -g_0_x_yyyyyz_zz[k] * ab_y + g_0_x_yyyyyz_yzz[k];
            }

            /// Set up 180-186 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyyzz_xx = cbuffer.data(kd_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_x_yyyyyzz_xy = cbuffer.data(kd_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_x_yyyyyzz_xz = cbuffer.data(kd_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_x_yyyyyzz_yy = cbuffer.data(kd_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_x_yyyyyzz_yz = cbuffer.data(kd_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_x_yyyyyzz_zz = cbuffer.data(kd_geom_01_off + 185 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyyzz_xx, g_0_x_yyyyyzz_xy, g_0_x_yyyyyzz_xz, g_0_x_yyyyyzz_yy, g_0_x_yyyyyzz_yz, g_0_x_yyyyyzz_zz, g_0_x_yyyyzz_xx, g_0_x_yyyyzz_xxy, g_0_x_yyyyzz_xy, g_0_x_yyyyzz_xyy, g_0_x_yyyyzz_xyz, g_0_x_yyyyzz_xz, g_0_x_yyyyzz_yy, g_0_x_yyyyzz_yyy, g_0_x_yyyyzz_yyz, g_0_x_yyyyzz_yz, g_0_x_yyyyzz_yzz, g_0_x_yyyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyyzz_xx[k] = -g_0_x_yyyyzz_xx[k] * ab_y + g_0_x_yyyyzz_xxy[k];

                g_0_x_yyyyyzz_xy[k] = -g_0_x_yyyyzz_xy[k] * ab_y + g_0_x_yyyyzz_xyy[k];

                g_0_x_yyyyyzz_xz[k] = -g_0_x_yyyyzz_xz[k] * ab_y + g_0_x_yyyyzz_xyz[k];

                g_0_x_yyyyyzz_yy[k] = -g_0_x_yyyyzz_yy[k] * ab_y + g_0_x_yyyyzz_yyy[k];

                g_0_x_yyyyyzz_yz[k] = -g_0_x_yyyyzz_yz[k] * ab_y + g_0_x_yyyyzz_yyz[k];

                g_0_x_yyyyyzz_zz[k] = -g_0_x_yyyyzz_zz[k] * ab_y + g_0_x_yyyyzz_yzz[k];
            }

            /// Set up 186-192 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyzzz_xx = cbuffer.data(kd_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_x_yyyyzzz_xy = cbuffer.data(kd_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_x_yyyyzzz_xz = cbuffer.data(kd_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_x_yyyyzzz_yy = cbuffer.data(kd_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_x_yyyyzzz_yz = cbuffer.data(kd_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_x_yyyyzzz_zz = cbuffer.data(kd_geom_01_off + 191 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyzzz_xx, g_0_x_yyyyzzz_xy, g_0_x_yyyyzzz_xz, g_0_x_yyyyzzz_yy, g_0_x_yyyyzzz_yz, g_0_x_yyyyzzz_zz, g_0_x_yyyzzz_xx, g_0_x_yyyzzz_xxy, g_0_x_yyyzzz_xy, g_0_x_yyyzzz_xyy, g_0_x_yyyzzz_xyz, g_0_x_yyyzzz_xz, g_0_x_yyyzzz_yy, g_0_x_yyyzzz_yyy, g_0_x_yyyzzz_yyz, g_0_x_yyyzzz_yz, g_0_x_yyyzzz_yzz, g_0_x_yyyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyzzz_xx[k] = -g_0_x_yyyzzz_xx[k] * ab_y + g_0_x_yyyzzz_xxy[k];

                g_0_x_yyyyzzz_xy[k] = -g_0_x_yyyzzz_xy[k] * ab_y + g_0_x_yyyzzz_xyy[k];

                g_0_x_yyyyzzz_xz[k] = -g_0_x_yyyzzz_xz[k] * ab_y + g_0_x_yyyzzz_xyz[k];

                g_0_x_yyyyzzz_yy[k] = -g_0_x_yyyzzz_yy[k] * ab_y + g_0_x_yyyzzz_yyy[k];

                g_0_x_yyyyzzz_yz[k] = -g_0_x_yyyzzz_yz[k] * ab_y + g_0_x_yyyzzz_yyz[k];

                g_0_x_yyyyzzz_zz[k] = -g_0_x_yyyzzz_zz[k] * ab_y + g_0_x_yyyzzz_yzz[k];
            }

            /// Set up 192-198 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyzzzz_xx = cbuffer.data(kd_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_x_yyyzzzz_xy = cbuffer.data(kd_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_x_yyyzzzz_xz = cbuffer.data(kd_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_x_yyyzzzz_yy = cbuffer.data(kd_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_x_yyyzzzz_yz = cbuffer.data(kd_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_x_yyyzzzz_zz = cbuffer.data(kd_geom_01_off + 197 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyzzzz_xx, g_0_x_yyyzzzz_xy, g_0_x_yyyzzzz_xz, g_0_x_yyyzzzz_yy, g_0_x_yyyzzzz_yz, g_0_x_yyyzzzz_zz, g_0_x_yyzzzz_xx, g_0_x_yyzzzz_xxy, g_0_x_yyzzzz_xy, g_0_x_yyzzzz_xyy, g_0_x_yyzzzz_xyz, g_0_x_yyzzzz_xz, g_0_x_yyzzzz_yy, g_0_x_yyzzzz_yyy, g_0_x_yyzzzz_yyz, g_0_x_yyzzzz_yz, g_0_x_yyzzzz_yzz, g_0_x_yyzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyzzzz_xx[k] = -g_0_x_yyzzzz_xx[k] * ab_y + g_0_x_yyzzzz_xxy[k];

                g_0_x_yyyzzzz_xy[k] = -g_0_x_yyzzzz_xy[k] * ab_y + g_0_x_yyzzzz_xyy[k];

                g_0_x_yyyzzzz_xz[k] = -g_0_x_yyzzzz_xz[k] * ab_y + g_0_x_yyzzzz_xyz[k];

                g_0_x_yyyzzzz_yy[k] = -g_0_x_yyzzzz_yy[k] * ab_y + g_0_x_yyzzzz_yyy[k];

                g_0_x_yyyzzzz_yz[k] = -g_0_x_yyzzzz_yz[k] * ab_y + g_0_x_yyzzzz_yyz[k];

                g_0_x_yyyzzzz_zz[k] = -g_0_x_yyzzzz_zz[k] * ab_y + g_0_x_yyzzzz_yzz[k];
            }

            /// Set up 198-204 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyzzzzz_xx = cbuffer.data(kd_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_x_yyzzzzz_xy = cbuffer.data(kd_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_x_yyzzzzz_xz = cbuffer.data(kd_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_x_yyzzzzz_yy = cbuffer.data(kd_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_x_yyzzzzz_yz = cbuffer.data(kd_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_x_yyzzzzz_zz = cbuffer.data(kd_geom_01_off + 203 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyzzzzz_xx, g_0_x_yyzzzzz_xy, g_0_x_yyzzzzz_xz, g_0_x_yyzzzzz_yy, g_0_x_yyzzzzz_yz, g_0_x_yyzzzzz_zz, g_0_x_yzzzzz_xx, g_0_x_yzzzzz_xxy, g_0_x_yzzzzz_xy, g_0_x_yzzzzz_xyy, g_0_x_yzzzzz_xyz, g_0_x_yzzzzz_xz, g_0_x_yzzzzz_yy, g_0_x_yzzzzz_yyy, g_0_x_yzzzzz_yyz, g_0_x_yzzzzz_yz, g_0_x_yzzzzz_yzz, g_0_x_yzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyzzzzz_xx[k] = -g_0_x_yzzzzz_xx[k] * ab_y + g_0_x_yzzzzz_xxy[k];

                g_0_x_yyzzzzz_xy[k] = -g_0_x_yzzzzz_xy[k] * ab_y + g_0_x_yzzzzz_xyy[k];

                g_0_x_yyzzzzz_xz[k] = -g_0_x_yzzzzz_xz[k] * ab_y + g_0_x_yzzzzz_xyz[k];

                g_0_x_yyzzzzz_yy[k] = -g_0_x_yzzzzz_yy[k] * ab_y + g_0_x_yzzzzz_yyy[k];

                g_0_x_yyzzzzz_yz[k] = -g_0_x_yzzzzz_yz[k] * ab_y + g_0_x_yzzzzz_yyz[k];

                g_0_x_yyzzzzz_zz[k] = -g_0_x_yzzzzz_zz[k] * ab_y + g_0_x_yzzzzz_yzz[k];
            }

            /// Set up 204-210 components of targeted buffer : cbuffer.data(

            auto g_0_x_yzzzzzz_xx = cbuffer.data(kd_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_x_yzzzzzz_xy = cbuffer.data(kd_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_x_yzzzzzz_xz = cbuffer.data(kd_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_x_yzzzzzz_yy = cbuffer.data(kd_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_x_yzzzzzz_yz = cbuffer.data(kd_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_x_yzzzzzz_zz = cbuffer.data(kd_geom_01_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yzzzzzz_xx, g_0_x_yzzzzzz_xy, g_0_x_yzzzzzz_xz, g_0_x_yzzzzzz_yy, g_0_x_yzzzzzz_yz, g_0_x_yzzzzzz_zz, g_0_x_zzzzzz_xx, g_0_x_zzzzzz_xxy, g_0_x_zzzzzz_xy, g_0_x_zzzzzz_xyy, g_0_x_zzzzzz_xyz, g_0_x_zzzzzz_xz, g_0_x_zzzzzz_yy, g_0_x_zzzzzz_yyy, g_0_x_zzzzzz_yyz, g_0_x_zzzzzz_yz, g_0_x_zzzzzz_yzz, g_0_x_zzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yzzzzzz_xx[k] = -g_0_x_zzzzzz_xx[k] * ab_y + g_0_x_zzzzzz_xxy[k];

                g_0_x_yzzzzzz_xy[k] = -g_0_x_zzzzzz_xy[k] * ab_y + g_0_x_zzzzzz_xyy[k];

                g_0_x_yzzzzzz_xz[k] = -g_0_x_zzzzzz_xz[k] * ab_y + g_0_x_zzzzzz_xyz[k];

                g_0_x_yzzzzzz_yy[k] = -g_0_x_zzzzzz_yy[k] * ab_y + g_0_x_zzzzzz_yyy[k];

                g_0_x_yzzzzzz_yz[k] = -g_0_x_zzzzzz_yz[k] * ab_y + g_0_x_zzzzzz_yyz[k];

                g_0_x_yzzzzzz_zz[k] = -g_0_x_zzzzzz_zz[k] * ab_y + g_0_x_zzzzzz_yzz[k];
            }

            /// Set up 210-216 components of targeted buffer : cbuffer.data(

            auto g_0_x_zzzzzzz_xx = cbuffer.data(kd_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_x_zzzzzzz_xy = cbuffer.data(kd_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_x_zzzzzzz_xz = cbuffer.data(kd_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_x_zzzzzzz_yy = cbuffer.data(kd_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_x_zzzzzzz_yz = cbuffer.data(kd_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_x_zzzzzzz_zz = cbuffer.data(kd_geom_01_off + 215 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_zzzzzz_xx, g_0_x_zzzzzz_xxz, g_0_x_zzzzzz_xy, g_0_x_zzzzzz_xyz, g_0_x_zzzzzz_xz, g_0_x_zzzzzz_xzz, g_0_x_zzzzzz_yy, g_0_x_zzzzzz_yyz, g_0_x_zzzzzz_yz, g_0_x_zzzzzz_yzz, g_0_x_zzzzzz_zz, g_0_x_zzzzzz_zzz, g_0_x_zzzzzzz_xx, g_0_x_zzzzzzz_xy, g_0_x_zzzzzzz_xz, g_0_x_zzzzzzz_yy, g_0_x_zzzzzzz_yz, g_0_x_zzzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_zzzzzzz_xx[k] = -g_0_x_zzzzzz_xx[k] * ab_z + g_0_x_zzzzzz_xxz[k];

                g_0_x_zzzzzzz_xy[k] = -g_0_x_zzzzzz_xy[k] * ab_z + g_0_x_zzzzzz_xyz[k];

                g_0_x_zzzzzzz_xz[k] = -g_0_x_zzzzzz_xz[k] * ab_z + g_0_x_zzzzzz_xzz[k];

                g_0_x_zzzzzzz_yy[k] = -g_0_x_zzzzzz_yy[k] * ab_z + g_0_x_zzzzzz_yyz[k];

                g_0_x_zzzzzzz_yz[k] = -g_0_x_zzzzzz_yz[k] * ab_z + g_0_x_zzzzzz_yzz[k];

                g_0_x_zzzzzzz_zz[k] = -g_0_x_zzzzzz_zz[k] * ab_z + g_0_x_zzzzzz_zzz[k];
            }

            /// Set up 216-222 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxxx_xx = cbuffer.data(kd_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_y_xxxxxxx_xy = cbuffer.data(kd_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_y_xxxxxxx_xz = cbuffer.data(kd_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_y_xxxxxxx_yy = cbuffer.data(kd_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_y_xxxxxxx_yz = cbuffer.data(kd_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_y_xxxxxxx_zz = cbuffer.data(kd_geom_01_off + 221 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxx_xx, g_0_y_xxxxxx_xxx, g_0_y_xxxxxx_xxy, g_0_y_xxxxxx_xxz, g_0_y_xxxxxx_xy, g_0_y_xxxxxx_xyy, g_0_y_xxxxxx_xyz, g_0_y_xxxxxx_xz, g_0_y_xxxxxx_xzz, g_0_y_xxxxxx_yy, g_0_y_xxxxxx_yz, g_0_y_xxxxxx_zz, g_0_y_xxxxxxx_xx, g_0_y_xxxxxxx_xy, g_0_y_xxxxxxx_xz, g_0_y_xxxxxxx_yy, g_0_y_xxxxxxx_yz, g_0_y_xxxxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxxx_xx[k] = -g_0_y_xxxxxx_xx[k] * ab_x + g_0_y_xxxxxx_xxx[k];

                g_0_y_xxxxxxx_xy[k] = -g_0_y_xxxxxx_xy[k] * ab_x + g_0_y_xxxxxx_xxy[k];

                g_0_y_xxxxxxx_xz[k] = -g_0_y_xxxxxx_xz[k] * ab_x + g_0_y_xxxxxx_xxz[k];

                g_0_y_xxxxxxx_yy[k] = -g_0_y_xxxxxx_yy[k] * ab_x + g_0_y_xxxxxx_xyy[k];

                g_0_y_xxxxxxx_yz[k] = -g_0_y_xxxxxx_yz[k] * ab_x + g_0_y_xxxxxx_xyz[k];

                g_0_y_xxxxxxx_zz[k] = -g_0_y_xxxxxx_zz[k] * ab_x + g_0_y_xxxxxx_xzz[k];
            }

            /// Set up 222-228 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxxy_xx = cbuffer.data(kd_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_y_xxxxxxy_xy = cbuffer.data(kd_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_y_xxxxxxy_xz = cbuffer.data(kd_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_y_xxxxxxy_yy = cbuffer.data(kd_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_y_xxxxxxy_yz = cbuffer.data(kd_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_y_xxxxxxy_zz = cbuffer.data(kd_geom_01_off + 227 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxxy_xx, g_0_y_xxxxxxy_xy, g_0_y_xxxxxxy_xz, g_0_y_xxxxxxy_yy, g_0_y_xxxxxxy_yz, g_0_y_xxxxxxy_zz, g_0_y_xxxxxy_xx, g_0_y_xxxxxy_xxx, g_0_y_xxxxxy_xxy, g_0_y_xxxxxy_xxz, g_0_y_xxxxxy_xy, g_0_y_xxxxxy_xyy, g_0_y_xxxxxy_xyz, g_0_y_xxxxxy_xz, g_0_y_xxxxxy_xzz, g_0_y_xxxxxy_yy, g_0_y_xxxxxy_yz, g_0_y_xxxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxxy_xx[k] = -g_0_y_xxxxxy_xx[k] * ab_x + g_0_y_xxxxxy_xxx[k];

                g_0_y_xxxxxxy_xy[k] = -g_0_y_xxxxxy_xy[k] * ab_x + g_0_y_xxxxxy_xxy[k];

                g_0_y_xxxxxxy_xz[k] = -g_0_y_xxxxxy_xz[k] * ab_x + g_0_y_xxxxxy_xxz[k];

                g_0_y_xxxxxxy_yy[k] = -g_0_y_xxxxxy_yy[k] * ab_x + g_0_y_xxxxxy_xyy[k];

                g_0_y_xxxxxxy_yz[k] = -g_0_y_xxxxxy_yz[k] * ab_x + g_0_y_xxxxxy_xyz[k];

                g_0_y_xxxxxxy_zz[k] = -g_0_y_xxxxxy_zz[k] * ab_x + g_0_y_xxxxxy_xzz[k];
            }

            /// Set up 228-234 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxxz_xx = cbuffer.data(kd_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_y_xxxxxxz_xy = cbuffer.data(kd_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_y_xxxxxxz_xz = cbuffer.data(kd_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_y_xxxxxxz_yy = cbuffer.data(kd_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_y_xxxxxxz_yz = cbuffer.data(kd_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_y_xxxxxxz_zz = cbuffer.data(kd_geom_01_off + 233 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxxz_xx, g_0_y_xxxxxxz_xy, g_0_y_xxxxxxz_xz, g_0_y_xxxxxxz_yy, g_0_y_xxxxxxz_yz, g_0_y_xxxxxxz_zz, g_0_y_xxxxxz_xx, g_0_y_xxxxxz_xxx, g_0_y_xxxxxz_xxy, g_0_y_xxxxxz_xxz, g_0_y_xxxxxz_xy, g_0_y_xxxxxz_xyy, g_0_y_xxxxxz_xyz, g_0_y_xxxxxz_xz, g_0_y_xxxxxz_xzz, g_0_y_xxxxxz_yy, g_0_y_xxxxxz_yz, g_0_y_xxxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxxz_xx[k] = -g_0_y_xxxxxz_xx[k] * ab_x + g_0_y_xxxxxz_xxx[k];

                g_0_y_xxxxxxz_xy[k] = -g_0_y_xxxxxz_xy[k] * ab_x + g_0_y_xxxxxz_xxy[k];

                g_0_y_xxxxxxz_xz[k] = -g_0_y_xxxxxz_xz[k] * ab_x + g_0_y_xxxxxz_xxz[k];

                g_0_y_xxxxxxz_yy[k] = -g_0_y_xxxxxz_yy[k] * ab_x + g_0_y_xxxxxz_xyy[k];

                g_0_y_xxxxxxz_yz[k] = -g_0_y_xxxxxz_yz[k] * ab_x + g_0_y_xxxxxz_xyz[k];

                g_0_y_xxxxxxz_zz[k] = -g_0_y_xxxxxz_zz[k] * ab_x + g_0_y_xxxxxz_xzz[k];
            }

            /// Set up 234-240 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxyy_xx = cbuffer.data(kd_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_y_xxxxxyy_xy = cbuffer.data(kd_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_y_xxxxxyy_xz = cbuffer.data(kd_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_y_xxxxxyy_yy = cbuffer.data(kd_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_y_xxxxxyy_yz = cbuffer.data(kd_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_y_xxxxxyy_zz = cbuffer.data(kd_geom_01_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxyy_xx, g_0_y_xxxxxyy_xy, g_0_y_xxxxxyy_xz, g_0_y_xxxxxyy_yy, g_0_y_xxxxxyy_yz, g_0_y_xxxxxyy_zz, g_0_y_xxxxyy_xx, g_0_y_xxxxyy_xxx, g_0_y_xxxxyy_xxy, g_0_y_xxxxyy_xxz, g_0_y_xxxxyy_xy, g_0_y_xxxxyy_xyy, g_0_y_xxxxyy_xyz, g_0_y_xxxxyy_xz, g_0_y_xxxxyy_xzz, g_0_y_xxxxyy_yy, g_0_y_xxxxyy_yz, g_0_y_xxxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxyy_xx[k] = -g_0_y_xxxxyy_xx[k] * ab_x + g_0_y_xxxxyy_xxx[k];

                g_0_y_xxxxxyy_xy[k] = -g_0_y_xxxxyy_xy[k] * ab_x + g_0_y_xxxxyy_xxy[k];

                g_0_y_xxxxxyy_xz[k] = -g_0_y_xxxxyy_xz[k] * ab_x + g_0_y_xxxxyy_xxz[k];

                g_0_y_xxxxxyy_yy[k] = -g_0_y_xxxxyy_yy[k] * ab_x + g_0_y_xxxxyy_xyy[k];

                g_0_y_xxxxxyy_yz[k] = -g_0_y_xxxxyy_yz[k] * ab_x + g_0_y_xxxxyy_xyz[k];

                g_0_y_xxxxxyy_zz[k] = -g_0_y_xxxxyy_zz[k] * ab_x + g_0_y_xxxxyy_xzz[k];
            }

            /// Set up 240-246 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxyz_xx = cbuffer.data(kd_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_y_xxxxxyz_xy = cbuffer.data(kd_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_y_xxxxxyz_xz = cbuffer.data(kd_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_y_xxxxxyz_yy = cbuffer.data(kd_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_y_xxxxxyz_yz = cbuffer.data(kd_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_y_xxxxxyz_zz = cbuffer.data(kd_geom_01_off + 245 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxyz_xx, g_0_y_xxxxxyz_xy, g_0_y_xxxxxyz_xz, g_0_y_xxxxxyz_yy, g_0_y_xxxxxyz_yz, g_0_y_xxxxxyz_zz, g_0_y_xxxxyz_xx, g_0_y_xxxxyz_xxx, g_0_y_xxxxyz_xxy, g_0_y_xxxxyz_xxz, g_0_y_xxxxyz_xy, g_0_y_xxxxyz_xyy, g_0_y_xxxxyz_xyz, g_0_y_xxxxyz_xz, g_0_y_xxxxyz_xzz, g_0_y_xxxxyz_yy, g_0_y_xxxxyz_yz, g_0_y_xxxxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxyz_xx[k] = -g_0_y_xxxxyz_xx[k] * ab_x + g_0_y_xxxxyz_xxx[k];

                g_0_y_xxxxxyz_xy[k] = -g_0_y_xxxxyz_xy[k] * ab_x + g_0_y_xxxxyz_xxy[k];

                g_0_y_xxxxxyz_xz[k] = -g_0_y_xxxxyz_xz[k] * ab_x + g_0_y_xxxxyz_xxz[k];

                g_0_y_xxxxxyz_yy[k] = -g_0_y_xxxxyz_yy[k] * ab_x + g_0_y_xxxxyz_xyy[k];

                g_0_y_xxxxxyz_yz[k] = -g_0_y_xxxxyz_yz[k] * ab_x + g_0_y_xxxxyz_xyz[k];

                g_0_y_xxxxxyz_zz[k] = -g_0_y_xxxxyz_zz[k] * ab_x + g_0_y_xxxxyz_xzz[k];
            }

            /// Set up 246-252 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxzz_xx = cbuffer.data(kd_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_y_xxxxxzz_xy = cbuffer.data(kd_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_y_xxxxxzz_xz = cbuffer.data(kd_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_y_xxxxxzz_yy = cbuffer.data(kd_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_y_xxxxxzz_yz = cbuffer.data(kd_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_y_xxxxxzz_zz = cbuffer.data(kd_geom_01_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxzz_xx, g_0_y_xxxxxzz_xy, g_0_y_xxxxxzz_xz, g_0_y_xxxxxzz_yy, g_0_y_xxxxxzz_yz, g_0_y_xxxxxzz_zz, g_0_y_xxxxzz_xx, g_0_y_xxxxzz_xxx, g_0_y_xxxxzz_xxy, g_0_y_xxxxzz_xxz, g_0_y_xxxxzz_xy, g_0_y_xxxxzz_xyy, g_0_y_xxxxzz_xyz, g_0_y_xxxxzz_xz, g_0_y_xxxxzz_xzz, g_0_y_xxxxzz_yy, g_0_y_xxxxzz_yz, g_0_y_xxxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxzz_xx[k] = -g_0_y_xxxxzz_xx[k] * ab_x + g_0_y_xxxxzz_xxx[k];

                g_0_y_xxxxxzz_xy[k] = -g_0_y_xxxxzz_xy[k] * ab_x + g_0_y_xxxxzz_xxy[k];

                g_0_y_xxxxxzz_xz[k] = -g_0_y_xxxxzz_xz[k] * ab_x + g_0_y_xxxxzz_xxz[k];

                g_0_y_xxxxxzz_yy[k] = -g_0_y_xxxxzz_yy[k] * ab_x + g_0_y_xxxxzz_xyy[k];

                g_0_y_xxxxxzz_yz[k] = -g_0_y_xxxxzz_yz[k] * ab_x + g_0_y_xxxxzz_xyz[k];

                g_0_y_xxxxxzz_zz[k] = -g_0_y_xxxxzz_zz[k] * ab_x + g_0_y_xxxxzz_xzz[k];
            }

            /// Set up 252-258 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxyyy_xx = cbuffer.data(kd_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_y_xxxxyyy_xy = cbuffer.data(kd_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_y_xxxxyyy_xz = cbuffer.data(kd_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_y_xxxxyyy_yy = cbuffer.data(kd_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_y_xxxxyyy_yz = cbuffer.data(kd_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_y_xxxxyyy_zz = cbuffer.data(kd_geom_01_off + 257 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxyyy_xx, g_0_y_xxxxyyy_xy, g_0_y_xxxxyyy_xz, g_0_y_xxxxyyy_yy, g_0_y_xxxxyyy_yz, g_0_y_xxxxyyy_zz, g_0_y_xxxyyy_xx, g_0_y_xxxyyy_xxx, g_0_y_xxxyyy_xxy, g_0_y_xxxyyy_xxz, g_0_y_xxxyyy_xy, g_0_y_xxxyyy_xyy, g_0_y_xxxyyy_xyz, g_0_y_xxxyyy_xz, g_0_y_xxxyyy_xzz, g_0_y_xxxyyy_yy, g_0_y_xxxyyy_yz, g_0_y_xxxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxyyy_xx[k] = -g_0_y_xxxyyy_xx[k] * ab_x + g_0_y_xxxyyy_xxx[k];

                g_0_y_xxxxyyy_xy[k] = -g_0_y_xxxyyy_xy[k] * ab_x + g_0_y_xxxyyy_xxy[k];

                g_0_y_xxxxyyy_xz[k] = -g_0_y_xxxyyy_xz[k] * ab_x + g_0_y_xxxyyy_xxz[k];

                g_0_y_xxxxyyy_yy[k] = -g_0_y_xxxyyy_yy[k] * ab_x + g_0_y_xxxyyy_xyy[k];

                g_0_y_xxxxyyy_yz[k] = -g_0_y_xxxyyy_yz[k] * ab_x + g_0_y_xxxyyy_xyz[k];

                g_0_y_xxxxyyy_zz[k] = -g_0_y_xxxyyy_zz[k] * ab_x + g_0_y_xxxyyy_xzz[k];
            }

            /// Set up 258-264 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxyyz_xx = cbuffer.data(kd_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_y_xxxxyyz_xy = cbuffer.data(kd_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_y_xxxxyyz_xz = cbuffer.data(kd_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_y_xxxxyyz_yy = cbuffer.data(kd_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_y_xxxxyyz_yz = cbuffer.data(kd_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_y_xxxxyyz_zz = cbuffer.data(kd_geom_01_off + 263 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxyyz_xx, g_0_y_xxxxyyz_xy, g_0_y_xxxxyyz_xz, g_0_y_xxxxyyz_yy, g_0_y_xxxxyyz_yz, g_0_y_xxxxyyz_zz, g_0_y_xxxyyz_xx, g_0_y_xxxyyz_xxx, g_0_y_xxxyyz_xxy, g_0_y_xxxyyz_xxz, g_0_y_xxxyyz_xy, g_0_y_xxxyyz_xyy, g_0_y_xxxyyz_xyz, g_0_y_xxxyyz_xz, g_0_y_xxxyyz_xzz, g_0_y_xxxyyz_yy, g_0_y_xxxyyz_yz, g_0_y_xxxyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxyyz_xx[k] = -g_0_y_xxxyyz_xx[k] * ab_x + g_0_y_xxxyyz_xxx[k];

                g_0_y_xxxxyyz_xy[k] = -g_0_y_xxxyyz_xy[k] * ab_x + g_0_y_xxxyyz_xxy[k];

                g_0_y_xxxxyyz_xz[k] = -g_0_y_xxxyyz_xz[k] * ab_x + g_0_y_xxxyyz_xxz[k];

                g_0_y_xxxxyyz_yy[k] = -g_0_y_xxxyyz_yy[k] * ab_x + g_0_y_xxxyyz_xyy[k];

                g_0_y_xxxxyyz_yz[k] = -g_0_y_xxxyyz_yz[k] * ab_x + g_0_y_xxxyyz_xyz[k];

                g_0_y_xxxxyyz_zz[k] = -g_0_y_xxxyyz_zz[k] * ab_x + g_0_y_xxxyyz_xzz[k];
            }

            /// Set up 264-270 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxyzz_xx = cbuffer.data(kd_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_y_xxxxyzz_xy = cbuffer.data(kd_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_y_xxxxyzz_xz = cbuffer.data(kd_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_y_xxxxyzz_yy = cbuffer.data(kd_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_y_xxxxyzz_yz = cbuffer.data(kd_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_y_xxxxyzz_zz = cbuffer.data(kd_geom_01_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxyzz_xx, g_0_y_xxxxyzz_xy, g_0_y_xxxxyzz_xz, g_0_y_xxxxyzz_yy, g_0_y_xxxxyzz_yz, g_0_y_xxxxyzz_zz, g_0_y_xxxyzz_xx, g_0_y_xxxyzz_xxx, g_0_y_xxxyzz_xxy, g_0_y_xxxyzz_xxz, g_0_y_xxxyzz_xy, g_0_y_xxxyzz_xyy, g_0_y_xxxyzz_xyz, g_0_y_xxxyzz_xz, g_0_y_xxxyzz_xzz, g_0_y_xxxyzz_yy, g_0_y_xxxyzz_yz, g_0_y_xxxyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxyzz_xx[k] = -g_0_y_xxxyzz_xx[k] * ab_x + g_0_y_xxxyzz_xxx[k];

                g_0_y_xxxxyzz_xy[k] = -g_0_y_xxxyzz_xy[k] * ab_x + g_0_y_xxxyzz_xxy[k];

                g_0_y_xxxxyzz_xz[k] = -g_0_y_xxxyzz_xz[k] * ab_x + g_0_y_xxxyzz_xxz[k];

                g_0_y_xxxxyzz_yy[k] = -g_0_y_xxxyzz_yy[k] * ab_x + g_0_y_xxxyzz_xyy[k];

                g_0_y_xxxxyzz_yz[k] = -g_0_y_xxxyzz_yz[k] * ab_x + g_0_y_xxxyzz_xyz[k];

                g_0_y_xxxxyzz_zz[k] = -g_0_y_xxxyzz_zz[k] * ab_x + g_0_y_xxxyzz_xzz[k];
            }

            /// Set up 270-276 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxzzz_xx = cbuffer.data(kd_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_y_xxxxzzz_xy = cbuffer.data(kd_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_y_xxxxzzz_xz = cbuffer.data(kd_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_y_xxxxzzz_yy = cbuffer.data(kd_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_y_xxxxzzz_yz = cbuffer.data(kd_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_y_xxxxzzz_zz = cbuffer.data(kd_geom_01_off + 275 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxzzz_xx, g_0_y_xxxxzzz_xy, g_0_y_xxxxzzz_xz, g_0_y_xxxxzzz_yy, g_0_y_xxxxzzz_yz, g_0_y_xxxxzzz_zz, g_0_y_xxxzzz_xx, g_0_y_xxxzzz_xxx, g_0_y_xxxzzz_xxy, g_0_y_xxxzzz_xxz, g_0_y_xxxzzz_xy, g_0_y_xxxzzz_xyy, g_0_y_xxxzzz_xyz, g_0_y_xxxzzz_xz, g_0_y_xxxzzz_xzz, g_0_y_xxxzzz_yy, g_0_y_xxxzzz_yz, g_0_y_xxxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxzzz_xx[k] = -g_0_y_xxxzzz_xx[k] * ab_x + g_0_y_xxxzzz_xxx[k];

                g_0_y_xxxxzzz_xy[k] = -g_0_y_xxxzzz_xy[k] * ab_x + g_0_y_xxxzzz_xxy[k];

                g_0_y_xxxxzzz_xz[k] = -g_0_y_xxxzzz_xz[k] * ab_x + g_0_y_xxxzzz_xxz[k];

                g_0_y_xxxxzzz_yy[k] = -g_0_y_xxxzzz_yy[k] * ab_x + g_0_y_xxxzzz_xyy[k];

                g_0_y_xxxxzzz_yz[k] = -g_0_y_xxxzzz_yz[k] * ab_x + g_0_y_xxxzzz_xyz[k];

                g_0_y_xxxxzzz_zz[k] = -g_0_y_xxxzzz_zz[k] * ab_x + g_0_y_xxxzzz_xzz[k];
            }

            /// Set up 276-282 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyyyy_xx = cbuffer.data(kd_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_y_xxxyyyy_xy = cbuffer.data(kd_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_y_xxxyyyy_xz = cbuffer.data(kd_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_y_xxxyyyy_yy = cbuffer.data(kd_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_y_xxxyyyy_yz = cbuffer.data(kd_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_y_xxxyyyy_zz = cbuffer.data(kd_geom_01_off + 281 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyyyy_xx, g_0_y_xxxyyyy_xy, g_0_y_xxxyyyy_xz, g_0_y_xxxyyyy_yy, g_0_y_xxxyyyy_yz, g_0_y_xxxyyyy_zz, g_0_y_xxyyyy_xx, g_0_y_xxyyyy_xxx, g_0_y_xxyyyy_xxy, g_0_y_xxyyyy_xxz, g_0_y_xxyyyy_xy, g_0_y_xxyyyy_xyy, g_0_y_xxyyyy_xyz, g_0_y_xxyyyy_xz, g_0_y_xxyyyy_xzz, g_0_y_xxyyyy_yy, g_0_y_xxyyyy_yz, g_0_y_xxyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyyyy_xx[k] = -g_0_y_xxyyyy_xx[k] * ab_x + g_0_y_xxyyyy_xxx[k];

                g_0_y_xxxyyyy_xy[k] = -g_0_y_xxyyyy_xy[k] * ab_x + g_0_y_xxyyyy_xxy[k];

                g_0_y_xxxyyyy_xz[k] = -g_0_y_xxyyyy_xz[k] * ab_x + g_0_y_xxyyyy_xxz[k];

                g_0_y_xxxyyyy_yy[k] = -g_0_y_xxyyyy_yy[k] * ab_x + g_0_y_xxyyyy_xyy[k];

                g_0_y_xxxyyyy_yz[k] = -g_0_y_xxyyyy_yz[k] * ab_x + g_0_y_xxyyyy_xyz[k];

                g_0_y_xxxyyyy_zz[k] = -g_0_y_xxyyyy_zz[k] * ab_x + g_0_y_xxyyyy_xzz[k];
            }

            /// Set up 282-288 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyyyz_xx = cbuffer.data(kd_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_y_xxxyyyz_xy = cbuffer.data(kd_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_y_xxxyyyz_xz = cbuffer.data(kd_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_y_xxxyyyz_yy = cbuffer.data(kd_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_y_xxxyyyz_yz = cbuffer.data(kd_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_y_xxxyyyz_zz = cbuffer.data(kd_geom_01_off + 287 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyyyz_xx, g_0_y_xxxyyyz_xy, g_0_y_xxxyyyz_xz, g_0_y_xxxyyyz_yy, g_0_y_xxxyyyz_yz, g_0_y_xxxyyyz_zz, g_0_y_xxyyyz_xx, g_0_y_xxyyyz_xxx, g_0_y_xxyyyz_xxy, g_0_y_xxyyyz_xxz, g_0_y_xxyyyz_xy, g_0_y_xxyyyz_xyy, g_0_y_xxyyyz_xyz, g_0_y_xxyyyz_xz, g_0_y_xxyyyz_xzz, g_0_y_xxyyyz_yy, g_0_y_xxyyyz_yz, g_0_y_xxyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyyyz_xx[k] = -g_0_y_xxyyyz_xx[k] * ab_x + g_0_y_xxyyyz_xxx[k];

                g_0_y_xxxyyyz_xy[k] = -g_0_y_xxyyyz_xy[k] * ab_x + g_0_y_xxyyyz_xxy[k];

                g_0_y_xxxyyyz_xz[k] = -g_0_y_xxyyyz_xz[k] * ab_x + g_0_y_xxyyyz_xxz[k];

                g_0_y_xxxyyyz_yy[k] = -g_0_y_xxyyyz_yy[k] * ab_x + g_0_y_xxyyyz_xyy[k];

                g_0_y_xxxyyyz_yz[k] = -g_0_y_xxyyyz_yz[k] * ab_x + g_0_y_xxyyyz_xyz[k];

                g_0_y_xxxyyyz_zz[k] = -g_0_y_xxyyyz_zz[k] * ab_x + g_0_y_xxyyyz_xzz[k];
            }

            /// Set up 288-294 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyyzz_xx = cbuffer.data(kd_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_y_xxxyyzz_xy = cbuffer.data(kd_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_y_xxxyyzz_xz = cbuffer.data(kd_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_y_xxxyyzz_yy = cbuffer.data(kd_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_y_xxxyyzz_yz = cbuffer.data(kd_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_y_xxxyyzz_zz = cbuffer.data(kd_geom_01_off + 293 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyyzz_xx, g_0_y_xxxyyzz_xy, g_0_y_xxxyyzz_xz, g_0_y_xxxyyzz_yy, g_0_y_xxxyyzz_yz, g_0_y_xxxyyzz_zz, g_0_y_xxyyzz_xx, g_0_y_xxyyzz_xxx, g_0_y_xxyyzz_xxy, g_0_y_xxyyzz_xxz, g_0_y_xxyyzz_xy, g_0_y_xxyyzz_xyy, g_0_y_xxyyzz_xyz, g_0_y_xxyyzz_xz, g_0_y_xxyyzz_xzz, g_0_y_xxyyzz_yy, g_0_y_xxyyzz_yz, g_0_y_xxyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyyzz_xx[k] = -g_0_y_xxyyzz_xx[k] * ab_x + g_0_y_xxyyzz_xxx[k];

                g_0_y_xxxyyzz_xy[k] = -g_0_y_xxyyzz_xy[k] * ab_x + g_0_y_xxyyzz_xxy[k];

                g_0_y_xxxyyzz_xz[k] = -g_0_y_xxyyzz_xz[k] * ab_x + g_0_y_xxyyzz_xxz[k];

                g_0_y_xxxyyzz_yy[k] = -g_0_y_xxyyzz_yy[k] * ab_x + g_0_y_xxyyzz_xyy[k];

                g_0_y_xxxyyzz_yz[k] = -g_0_y_xxyyzz_yz[k] * ab_x + g_0_y_xxyyzz_xyz[k];

                g_0_y_xxxyyzz_zz[k] = -g_0_y_xxyyzz_zz[k] * ab_x + g_0_y_xxyyzz_xzz[k];
            }

            /// Set up 294-300 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyzzz_xx = cbuffer.data(kd_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_y_xxxyzzz_xy = cbuffer.data(kd_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_y_xxxyzzz_xz = cbuffer.data(kd_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_y_xxxyzzz_yy = cbuffer.data(kd_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_y_xxxyzzz_yz = cbuffer.data(kd_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_y_xxxyzzz_zz = cbuffer.data(kd_geom_01_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyzzz_xx, g_0_y_xxxyzzz_xy, g_0_y_xxxyzzz_xz, g_0_y_xxxyzzz_yy, g_0_y_xxxyzzz_yz, g_0_y_xxxyzzz_zz, g_0_y_xxyzzz_xx, g_0_y_xxyzzz_xxx, g_0_y_xxyzzz_xxy, g_0_y_xxyzzz_xxz, g_0_y_xxyzzz_xy, g_0_y_xxyzzz_xyy, g_0_y_xxyzzz_xyz, g_0_y_xxyzzz_xz, g_0_y_xxyzzz_xzz, g_0_y_xxyzzz_yy, g_0_y_xxyzzz_yz, g_0_y_xxyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyzzz_xx[k] = -g_0_y_xxyzzz_xx[k] * ab_x + g_0_y_xxyzzz_xxx[k];

                g_0_y_xxxyzzz_xy[k] = -g_0_y_xxyzzz_xy[k] * ab_x + g_0_y_xxyzzz_xxy[k];

                g_0_y_xxxyzzz_xz[k] = -g_0_y_xxyzzz_xz[k] * ab_x + g_0_y_xxyzzz_xxz[k];

                g_0_y_xxxyzzz_yy[k] = -g_0_y_xxyzzz_yy[k] * ab_x + g_0_y_xxyzzz_xyy[k];

                g_0_y_xxxyzzz_yz[k] = -g_0_y_xxyzzz_yz[k] * ab_x + g_0_y_xxyzzz_xyz[k];

                g_0_y_xxxyzzz_zz[k] = -g_0_y_xxyzzz_zz[k] * ab_x + g_0_y_xxyzzz_xzz[k];
            }

            /// Set up 300-306 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxzzzz_xx = cbuffer.data(kd_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_y_xxxzzzz_xy = cbuffer.data(kd_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_y_xxxzzzz_xz = cbuffer.data(kd_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_y_xxxzzzz_yy = cbuffer.data(kd_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_y_xxxzzzz_yz = cbuffer.data(kd_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_y_xxxzzzz_zz = cbuffer.data(kd_geom_01_off + 305 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxzzzz_xx, g_0_y_xxxzzzz_xy, g_0_y_xxxzzzz_xz, g_0_y_xxxzzzz_yy, g_0_y_xxxzzzz_yz, g_0_y_xxxzzzz_zz, g_0_y_xxzzzz_xx, g_0_y_xxzzzz_xxx, g_0_y_xxzzzz_xxy, g_0_y_xxzzzz_xxz, g_0_y_xxzzzz_xy, g_0_y_xxzzzz_xyy, g_0_y_xxzzzz_xyz, g_0_y_xxzzzz_xz, g_0_y_xxzzzz_xzz, g_0_y_xxzzzz_yy, g_0_y_xxzzzz_yz, g_0_y_xxzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxzzzz_xx[k] = -g_0_y_xxzzzz_xx[k] * ab_x + g_0_y_xxzzzz_xxx[k];

                g_0_y_xxxzzzz_xy[k] = -g_0_y_xxzzzz_xy[k] * ab_x + g_0_y_xxzzzz_xxy[k];

                g_0_y_xxxzzzz_xz[k] = -g_0_y_xxzzzz_xz[k] * ab_x + g_0_y_xxzzzz_xxz[k];

                g_0_y_xxxzzzz_yy[k] = -g_0_y_xxzzzz_yy[k] * ab_x + g_0_y_xxzzzz_xyy[k];

                g_0_y_xxxzzzz_yz[k] = -g_0_y_xxzzzz_yz[k] * ab_x + g_0_y_xxzzzz_xyz[k];

                g_0_y_xxxzzzz_zz[k] = -g_0_y_xxzzzz_zz[k] * ab_x + g_0_y_xxzzzz_xzz[k];
            }

            /// Set up 306-312 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyyyy_xx = cbuffer.data(kd_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_y_xxyyyyy_xy = cbuffer.data(kd_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_y_xxyyyyy_xz = cbuffer.data(kd_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_y_xxyyyyy_yy = cbuffer.data(kd_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_y_xxyyyyy_yz = cbuffer.data(kd_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_y_xxyyyyy_zz = cbuffer.data(kd_geom_01_off + 311 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyyyy_xx, g_0_y_xxyyyyy_xy, g_0_y_xxyyyyy_xz, g_0_y_xxyyyyy_yy, g_0_y_xxyyyyy_yz, g_0_y_xxyyyyy_zz, g_0_y_xyyyyy_xx, g_0_y_xyyyyy_xxx, g_0_y_xyyyyy_xxy, g_0_y_xyyyyy_xxz, g_0_y_xyyyyy_xy, g_0_y_xyyyyy_xyy, g_0_y_xyyyyy_xyz, g_0_y_xyyyyy_xz, g_0_y_xyyyyy_xzz, g_0_y_xyyyyy_yy, g_0_y_xyyyyy_yz, g_0_y_xyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyyyy_xx[k] = -g_0_y_xyyyyy_xx[k] * ab_x + g_0_y_xyyyyy_xxx[k];

                g_0_y_xxyyyyy_xy[k] = -g_0_y_xyyyyy_xy[k] * ab_x + g_0_y_xyyyyy_xxy[k];

                g_0_y_xxyyyyy_xz[k] = -g_0_y_xyyyyy_xz[k] * ab_x + g_0_y_xyyyyy_xxz[k];

                g_0_y_xxyyyyy_yy[k] = -g_0_y_xyyyyy_yy[k] * ab_x + g_0_y_xyyyyy_xyy[k];

                g_0_y_xxyyyyy_yz[k] = -g_0_y_xyyyyy_yz[k] * ab_x + g_0_y_xyyyyy_xyz[k];

                g_0_y_xxyyyyy_zz[k] = -g_0_y_xyyyyy_zz[k] * ab_x + g_0_y_xyyyyy_xzz[k];
            }

            /// Set up 312-318 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyyyz_xx = cbuffer.data(kd_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_y_xxyyyyz_xy = cbuffer.data(kd_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_y_xxyyyyz_xz = cbuffer.data(kd_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_y_xxyyyyz_yy = cbuffer.data(kd_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_y_xxyyyyz_yz = cbuffer.data(kd_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_y_xxyyyyz_zz = cbuffer.data(kd_geom_01_off + 317 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyyyz_xx, g_0_y_xxyyyyz_xy, g_0_y_xxyyyyz_xz, g_0_y_xxyyyyz_yy, g_0_y_xxyyyyz_yz, g_0_y_xxyyyyz_zz, g_0_y_xyyyyz_xx, g_0_y_xyyyyz_xxx, g_0_y_xyyyyz_xxy, g_0_y_xyyyyz_xxz, g_0_y_xyyyyz_xy, g_0_y_xyyyyz_xyy, g_0_y_xyyyyz_xyz, g_0_y_xyyyyz_xz, g_0_y_xyyyyz_xzz, g_0_y_xyyyyz_yy, g_0_y_xyyyyz_yz, g_0_y_xyyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyyyz_xx[k] = -g_0_y_xyyyyz_xx[k] * ab_x + g_0_y_xyyyyz_xxx[k];

                g_0_y_xxyyyyz_xy[k] = -g_0_y_xyyyyz_xy[k] * ab_x + g_0_y_xyyyyz_xxy[k];

                g_0_y_xxyyyyz_xz[k] = -g_0_y_xyyyyz_xz[k] * ab_x + g_0_y_xyyyyz_xxz[k];

                g_0_y_xxyyyyz_yy[k] = -g_0_y_xyyyyz_yy[k] * ab_x + g_0_y_xyyyyz_xyy[k];

                g_0_y_xxyyyyz_yz[k] = -g_0_y_xyyyyz_yz[k] * ab_x + g_0_y_xyyyyz_xyz[k];

                g_0_y_xxyyyyz_zz[k] = -g_0_y_xyyyyz_zz[k] * ab_x + g_0_y_xyyyyz_xzz[k];
            }

            /// Set up 318-324 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyyzz_xx = cbuffer.data(kd_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_y_xxyyyzz_xy = cbuffer.data(kd_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_y_xxyyyzz_xz = cbuffer.data(kd_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_y_xxyyyzz_yy = cbuffer.data(kd_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_y_xxyyyzz_yz = cbuffer.data(kd_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_y_xxyyyzz_zz = cbuffer.data(kd_geom_01_off + 323 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyyzz_xx, g_0_y_xxyyyzz_xy, g_0_y_xxyyyzz_xz, g_0_y_xxyyyzz_yy, g_0_y_xxyyyzz_yz, g_0_y_xxyyyzz_zz, g_0_y_xyyyzz_xx, g_0_y_xyyyzz_xxx, g_0_y_xyyyzz_xxy, g_0_y_xyyyzz_xxz, g_0_y_xyyyzz_xy, g_0_y_xyyyzz_xyy, g_0_y_xyyyzz_xyz, g_0_y_xyyyzz_xz, g_0_y_xyyyzz_xzz, g_0_y_xyyyzz_yy, g_0_y_xyyyzz_yz, g_0_y_xyyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyyzz_xx[k] = -g_0_y_xyyyzz_xx[k] * ab_x + g_0_y_xyyyzz_xxx[k];

                g_0_y_xxyyyzz_xy[k] = -g_0_y_xyyyzz_xy[k] * ab_x + g_0_y_xyyyzz_xxy[k];

                g_0_y_xxyyyzz_xz[k] = -g_0_y_xyyyzz_xz[k] * ab_x + g_0_y_xyyyzz_xxz[k];

                g_0_y_xxyyyzz_yy[k] = -g_0_y_xyyyzz_yy[k] * ab_x + g_0_y_xyyyzz_xyy[k];

                g_0_y_xxyyyzz_yz[k] = -g_0_y_xyyyzz_yz[k] * ab_x + g_0_y_xyyyzz_xyz[k];

                g_0_y_xxyyyzz_zz[k] = -g_0_y_xyyyzz_zz[k] * ab_x + g_0_y_xyyyzz_xzz[k];
            }

            /// Set up 324-330 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyzzz_xx = cbuffer.data(kd_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_y_xxyyzzz_xy = cbuffer.data(kd_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_y_xxyyzzz_xz = cbuffer.data(kd_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_y_xxyyzzz_yy = cbuffer.data(kd_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_y_xxyyzzz_yz = cbuffer.data(kd_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_y_xxyyzzz_zz = cbuffer.data(kd_geom_01_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyzzz_xx, g_0_y_xxyyzzz_xy, g_0_y_xxyyzzz_xz, g_0_y_xxyyzzz_yy, g_0_y_xxyyzzz_yz, g_0_y_xxyyzzz_zz, g_0_y_xyyzzz_xx, g_0_y_xyyzzz_xxx, g_0_y_xyyzzz_xxy, g_0_y_xyyzzz_xxz, g_0_y_xyyzzz_xy, g_0_y_xyyzzz_xyy, g_0_y_xyyzzz_xyz, g_0_y_xyyzzz_xz, g_0_y_xyyzzz_xzz, g_0_y_xyyzzz_yy, g_0_y_xyyzzz_yz, g_0_y_xyyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyzzz_xx[k] = -g_0_y_xyyzzz_xx[k] * ab_x + g_0_y_xyyzzz_xxx[k];

                g_0_y_xxyyzzz_xy[k] = -g_0_y_xyyzzz_xy[k] * ab_x + g_0_y_xyyzzz_xxy[k];

                g_0_y_xxyyzzz_xz[k] = -g_0_y_xyyzzz_xz[k] * ab_x + g_0_y_xyyzzz_xxz[k];

                g_0_y_xxyyzzz_yy[k] = -g_0_y_xyyzzz_yy[k] * ab_x + g_0_y_xyyzzz_xyy[k];

                g_0_y_xxyyzzz_yz[k] = -g_0_y_xyyzzz_yz[k] * ab_x + g_0_y_xyyzzz_xyz[k];

                g_0_y_xxyyzzz_zz[k] = -g_0_y_xyyzzz_zz[k] * ab_x + g_0_y_xyyzzz_xzz[k];
            }

            /// Set up 330-336 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyzzzz_xx = cbuffer.data(kd_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_y_xxyzzzz_xy = cbuffer.data(kd_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_y_xxyzzzz_xz = cbuffer.data(kd_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_y_xxyzzzz_yy = cbuffer.data(kd_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_y_xxyzzzz_yz = cbuffer.data(kd_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_y_xxyzzzz_zz = cbuffer.data(kd_geom_01_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyzzzz_xx, g_0_y_xxyzzzz_xy, g_0_y_xxyzzzz_xz, g_0_y_xxyzzzz_yy, g_0_y_xxyzzzz_yz, g_0_y_xxyzzzz_zz, g_0_y_xyzzzz_xx, g_0_y_xyzzzz_xxx, g_0_y_xyzzzz_xxy, g_0_y_xyzzzz_xxz, g_0_y_xyzzzz_xy, g_0_y_xyzzzz_xyy, g_0_y_xyzzzz_xyz, g_0_y_xyzzzz_xz, g_0_y_xyzzzz_xzz, g_0_y_xyzzzz_yy, g_0_y_xyzzzz_yz, g_0_y_xyzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyzzzz_xx[k] = -g_0_y_xyzzzz_xx[k] * ab_x + g_0_y_xyzzzz_xxx[k];

                g_0_y_xxyzzzz_xy[k] = -g_0_y_xyzzzz_xy[k] * ab_x + g_0_y_xyzzzz_xxy[k];

                g_0_y_xxyzzzz_xz[k] = -g_0_y_xyzzzz_xz[k] * ab_x + g_0_y_xyzzzz_xxz[k];

                g_0_y_xxyzzzz_yy[k] = -g_0_y_xyzzzz_yy[k] * ab_x + g_0_y_xyzzzz_xyy[k];

                g_0_y_xxyzzzz_yz[k] = -g_0_y_xyzzzz_yz[k] * ab_x + g_0_y_xyzzzz_xyz[k];

                g_0_y_xxyzzzz_zz[k] = -g_0_y_xyzzzz_zz[k] * ab_x + g_0_y_xyzzzz_xzz[k];
            }

            /// Set up 336-342 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxzzzzz_xx = cbuffer.data(kd_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_y_xxzzzzz_xy = cbuffer.data(kd_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_y_xxzzzzz_xz = cbuffer.data(kd_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_y_xxzzzzz_yy = cbuffer.data(kd_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_y_xxzzzzz_yz = cbuffer.data(kd_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_y_xxzzzzz_zz = cbuffer.data(kd_geom_01_off + 341 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxzzzzz_xx, g_0_y_xxzzzzz_xy, g_0_y_xxzzzzz_xz, g_0_y_xxzzzzz_yy, g_0_y_xxzzzzz_yz, g_0_y_xxzzzzz_zz, g_0_y_xzzzzz_xx, g_0_y_xzzzzz_xxx, g_0_y_xzzzzz_xxy, g_0_y_xzzzzz_xxz, g_0_y_xzzzzz_xy, g_0_y_xzzzzz_xyy, g_0_y_xzzzzz_xyz, g_0_y_xzzzzz_xz, g_0_y_xzzzzz_xzz, g_0_y_xzzzzz_yy, g_0_y_xzzzzz_yz, g_0_y_xzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxzzzzz_xx[k] = -g_0_y_xzzzzz_xx[k] * ab_x + g_0_y_xzzzzz_xxx[k];

                g_0_y_xxzzzzz_xy[k] = -g_0_y_xzzzzz_xy[k] * ab_x + g_0_y_xzzzzz_xxy[k];

                g_0_y_xxzzzzz_xz[k] = -g_0_y_xzzzzz_xz[k] * ab_x + g_0_y_xzzzzz_xxz[k];

                g_0_y_xxzzzzz_yy[k] = -g_0_y_xzzzzz_yy[k] * ab_x + g_0_y_xzzzzz_xyy[k];

                g_0_y_xxzzzzz_yz[k] = -g_0_y_xzzzzz_yz[k] * ab_x + g_0_y_xzzzzz_xyz[k];

                g_0_y_xxzzzzz_zz[k] = -g_0_y_xzzzzz_zz[k] * ab_x + g_0_y_xzzzzz_xzz[k];
            }

            /// Set up 342-348 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyyyy_xx = cbuffer.data(kd_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_y_xyyyyyy_xy = cbuffer.data(kd_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_y_xyyyyyy_xz = cbuffer.data(kd_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_y_xyyyyyy_yy = cbuffer.data(kd_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_y_xyyyyyy_yz = cbuffer.data(kd_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_y_xyyyyyy_zz = cbuffer.data(kd_geom_01_off + 347 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyyyy_xx, g_0_y_xyyyyyy_xy, g_0_y_xyyyyyy_xz, g_0_y_xyyyyyy_yy, g_0_y_xyyyyyy_yz, g_0_y_xyyyyyy_zz, g_0_y_yyyyyy_xx, g_0_y_yyyyyy_xxx, g_0_y_yyyyyy_xxy, g_0_y_yyyyyy_xxz, g_0_y_yyyyyy_xy, g_0_y_yyyyyy_xyy, g_0_y_yyyyyy_xyz, g_0_y_yyyyyy_xz, g_0_y_yyyyyy_xzz, g_0_y_yyyyyy_yy, g_0_y_yyyyyy_yz, g_0_y_yyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyyyy_xx[k] = -g_0_y_yyyyyy_xx[k] * ab_x + g_0_y_yyyyyy_xxx[k];

                g_0_y_xyyyyyy_xy[k] = -g_0_y_yyyyyy_xy[k] * ab_x + g_0_y_yyyyyy_xxy[k];

                g_0_y_xyyyyyy_xz[k] = -g_0_y_yyyyyy_xz[k] * ab_x + g_0_y_yyyyyy_xxz[k];

                g_0_y_xyyyyyy_yy[k] = -g_0_y_yyyyyy_yy[k] * ab_x + g_0_y_yyyyyy_xyy[k];

                g_0_y_xyyyyyy_yz[k] = -g_0_y_yyyyyy_yz[k] * ab_x + g_0_y_yyyyyy_xyz[k];

                g_0_y_xyyyyyy_zz[k] = -g_0_y_yyyyyy_zz[k] * ab_x + g_0_y_yyyyyy_xzz[k];
            }

            /// Set up 348-354 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyyyz_xx = cbuffer.data(kd_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_y_xyyyyyz_xy = cbuffer.data(kd_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_y_xyyyyyz_xz = cbuffer.data(kd_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_y_xyyyyyz_yy = cbuffer.data(kd_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_y_xyyyyyz_yz = cbuffer.data(kd_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_y_xyyyyyz_zz = cbuffer.data(kd_geom_01_off + 353 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyyyz_xx, g_0_y_xyyyyyz_xy, g_0_y_xyyyyyz_xz, g_0_y_xyyyyyz_yy, g_0_y_xyyyyyz_yz, g_0_y_xyyyyyz_zz, g_0_y_yyyyyz_xx, g_0_y_yyyyyz_xxx, g_0_y_yyyyyz_xxy, g_0_y_yyyyyz_xxz, g_0_y_yyyyyz_xy, g_0_y_yyyyyz_xyy, g_0_y_yyyyyz_xyz, g_0_y_yyyyyz_xz, g_0_y_yyyyyz_xzz, g_0_y_yyyyyz_yy, g_0_y_yyyyyz_yz, g_0_y_yyyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyyyz_xx[k] = -g_0_y_yyyyyz_xx[k] * ab_x + g_0_y_yyyyyz_xxx[k];

                g_0_y_xyyyyyz_xy[k] = -g_0_y_yyyyyz_xy[k] * ab_x + g_0_y_yyyyyz_xxy[k];

                g_0_y_xyyyyyz_xz[k] = -g_0_y_yyyyyz_xz[k] * ab_x + g_0_y_yyyyyz_xxz[k];

                g_0_y_xyyyyyz_yy[k] = -g_0_y_yyyyyz_yy[k] * ab_x + g_0_y_yyyyyz_xyy[k];

                g_0_y_xyyyyyz_yz[k] = -g_0_y_yyyyyz_yz[k] * ab_x + g_0_y_yyyyyz_xyz[k];

                g_0_y_xyyyyyz_zz[k] = -g_0_y_yyyyyz_zz[k] * ab_x + g_0_y_yyyyyz_xzz[k];
            }

            /// Set up 354-360 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyyzz_xx = cbuffer.data(kd_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_y_xyyyyzz_xy = cbuffer.data(kd_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_y_xyyyyzz_xz = cbuffer.data(kd_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_y_xyyyyzz_yy = cbuffer.data(kd_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_y_xyyyyzz_yz = cbuffer.data(kd_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_y_xyyyyzz_zz = cbuffer.data(kd_geom_01_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyyzz_xx, g_0_y_xyyyyzz_xy, g_0_y_xyyyyzz_xz, g_0_y_xyyyyzz_yy, g_0_y_xyyyyzz_yz, g_0_y_xyyyyzz_zz, g_0_y_yyyyzz_xx, g_0_y_yyyyzz_xxx, g_0_y_yyyyzz_xxy, g_0_y_yyyyzz_xxz, g_0_y_yyyyzz_xy, g_0_y_yyyyzz_xyy, g_0_y_yyyyzz_xyz, g_0_y_yyyyzz_xz, g_0_y_yyyyzz_xzz, g_0_y_yyyyzz_yy, g_0_y_yyyyzz_yz, g_0_y_yyyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyyzz_xx[k] = -g_0_y_yyyyzz_xx[k] * ab_x + g_0_y_yyyyzz_xxx[k];

                g_0_y_xyyyyzz_xy[k] = -g_0_y_yyyyzz_xy[k] * ab_x + g_0_y_yyyyzz_xxy[k];

                g_0_y_xyyyyzz_xz[k] = -g_0_y_yyyyzz_xz[k] * ab_x + g_0_y_yyyyzz_xxz[k];

                g_0_y_xyyyyzz_yy[k] = -g_0_y_yyyyzz_yy[k] * ab_x + g_0_y_yyyyzz_xyy[k];

                g_0_y_xyyyyzz_yz[k] = -g_0_y_yyyyzz_yz[k] * ab_x + g_0_y_yyyyzz_xyz[k];

                g_0_y_xyyyyzz_zz[k] = -g_0_y_yyyyzz_zz[k] * ab_x + g_0_y_yyyyzz_xzz[k];
            }

            /// Set up 360-366 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyzzz_xx = cbuffer.data(kd_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_y_xyyyzzz_xy = cbuffer.data(kd_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_y_xyyyzzz_xz = cbuffer.data(kd_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_y_xyyyzzz_yy = cbuffer.data(kd_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_y_xyyyzzz_yz = cbuffer.data(kd_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_y_xyyyzzz_zz = cbuffer.data(kd_geom_01_off + 365 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyzzz_xx, g_0_y_xyyyzzz_xy, g_0_y_xyyyzzz_xz, g_0_y_xyyyzzz_yy, g_0_y_xyyyzzz_yz, g_0_y_xyyyzzz_zz, g_0_y_yyyzzz_xx, g_0_y_yyyzzz_xxx, g_0_y_yyyzzz_xxy, g_0_y_yyyzzz_xxz, g_0_y_yyyzzz_xy, g_0_y_yyyzzz_xyy, g_0_y_yyyzzz_xyz, g_0_y_yyyzzz_xz, g_0_y_yyyzzz_xzz, g_0_y_yyyzzz_yy, g_0_y_yyyzzz_yz, g_0_y_yyyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyzzz_xx[k] = -g_0_y_yyyzzz_xx[k] * ab_x + g_0_y_yyyzzz_xxx[k];

                g_0_y_xyyyzzz_xy[k] = -g_0_y_yyyzzz_xy[k] * ab_x + g_0_y_yyyzzz_xxy[k];

                g_0_y_xyyyzzz_xz[k] = -g_0_y_yyyzzz_xz[k] * ab_x + g_0_y_yyyzzz_xxz[k];

                g_0_y_xyyyzzz_yy[k] = -g_0_y_yyyzzz_yy[k] * ab_x + g_0_y_yyyzzz_xyy[k];

                g_0_y_xyyyzzz_yz[k] = -g_0_y_yyyzzz_yz[k] * ab_x + g_0_y_yyyzzz_xyz[k];

                g_0_y_xyyyzzz_zz[k] = -g_0_y_yyyzzz_zz[k] * ab_x + g_0_y_yyyzzz_xzz[k];
            }

            /// Set up 366-372 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyzzzz_xx = cbuffer.data(kd_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_y_xyyzzzz_xy = cbuffer.data(kd_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_y_xyyzzzz_xz = cbuffer.data(kd_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_y_xyyzzzz_yy = cbuffer.data(kd_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_y_xyyzzzz_yz = cbuffer.data(kd_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_y_xyyzzzz_zz = cbuffer.data(kd_geom_01_off + 371 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyzzzz_xx, g_0_y_xyyzzzz_xy, g_0_y_xyyzzzz_xz, g_0_y_xyyzzzz_yy, g_0_y_xyyzzzz_yz, g_0_y_xyyzzzz_zz, g_0_y_yyzzzz_xx, g_0_y_yyzzzz_xxx, g_0_y_yyzzzz_xxy, g_0_y_yyzzzz_xxz, g_0_y_yyzzzz_xy, g_0_y_yyzzzz_xyy, g_0_y_yyzzzz_xyz, g_0_y_yyzzzz_xz, g_0_y_yyzzzz_xzz, g_0_y_yyzzzz_yy, g_0_y_yyzzzz_yz, g_0_y_yyzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyzzzz_xx[k] = -g_0_y_yyzzzz_xx[k] * ab_x + g_0_y_yyzzzz_xxx[k];

                g_0_y_xyyzzzz_xy[k] = -g_0_y_yyzzzz_xy[k] * ab_x + g_0_y_yyzzzz_xxy[k];

                g_0_y_xyyzzzz_xz[k] = -g_0_y_yyzzzz_xz[k] * ab_x + g_0_y_yyzzzz_xxz[k];

                g_0_y_xyyzzzz_yy[k] = -g_0_y_yyzzzz_yy[k] * ab_x + g_0_y_yyzzzz_xyy[k];

                g_0_y_xyyzzzz_yz[k] = -g_0_y_yyzzzz_yz[k] * ab_x + g_0_y_yyzzzz_xyz[k];

                g_0_y_xyyzzzz_zz[k] = -g_0_y_yyzzzz_zz[k] * ab_x + g_0_y_yyzzzz_xzz[k];
            }

            /// Set up 372-378 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyzzzzz_xx = cbuffer.data(kd_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_y_xyzzzzz_xy = cbuffer.data(kd_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_y_xyzzzzz_xz = cbuffer.data(kd_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_y_xyzzzzz_yy = cbuffer.data(kd_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_y_xyzzzzz_yz = cbuffer.data(kd_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_y_xyzzzzz_zz = cbuffer.data(kd_geom_01_off + 377 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyzzzzz_xx, g_0_y_xyzzzzz_xy, g_0_y_xyzzzzz_xz, g_0_y_xyzzzzz_yy, g_0_y_xyzzzzz_yz, g_0_y_xyzzzzz_zz, g_0_y_yzzzzz_xx, g_0_y_yzzzzz_xxx, g_0_y_yzzzzz_xxy, g_0_y_yzzzzz_xxz, g_0_y_yzzzzz_xy, g_0_y_yzzzzz_xyy, g_0_y_yzzzzz_xyz, g_0_y_yzzzzz_xz, g_0_y_yzzzzz_xzz, g_0_y_yzzzzz_yy, g_0_y_yzzzzz_yz, g_0_y_yzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyzzzzz_xx[k] = -g_0_y_yzzzzz_xx[k] * ab_x + g_0_y_yzzzzz_xxx[k];

                g_0_y_xyzzzzz_xy[k] = -g_0_y_yzzzzz_xy[k] * ab_x + g_0_y_yzzzzz_xxy[k];

                g_0_y_xyzzzzz_xz[k] = -g_0_y_yzzzzz_xz[k] * ab_x + g_0_y_yzzzzz_xxz[k];

                g_0_y_xyzzzzz_yy[k] = -g_0_y_yzzzzz_yy[k] * ab_x + g_0_y_yzzzzz_xyy[k];

                g_0_y_xyzzzzz_yz[k] = -g_0_y_yzzzzz_yz[k] * ab_x + g_0_y_yzzzzz_xyz[k];

                g_0_y_xyzzzzz_zz[k] = -g_0_y_yzzzzz_zz[k] * ab_x + g_0_y_yzzzzz_xzz[k];
            }

            /// Set up 378-384 components of targeted buffer : cbuffer.data(

            auto g_0_y_xzzzzzz_xx = cbuffer.data(kd_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_y_xzzzzzz_xy = cbuffer.data(kd_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_y_xzzzzzz_xz = cbuffer.data(kd_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_y_xzzzzzz_yy = cbuffer.data(kd_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_y_xzzzzzz_yz = cbuffer.data(kd_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_y_xzzzzzz_zz = cbuffer.data(kd_geom_01_off + 383 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xzzzzzz_xx, g_0_y_xzzzzzz_xy, g_0_y_xzzzzzz_xz, g_0_y_xzzzzzz_yy, g_0_y_xzzzzzz_yz, g_0_y_xzzzzzz_zz, g_0_y_zzzzzz_xx, g_0_y_zzzzzz_xxx, g_0_y_zzzzzz_xxy, g_0_y_zzzzzz_xxz, g_0_y_zzzzzz_xy, g_0_y_zzzzzz_xyy, g_0_y_zzzzzz_xyz, g_0_y_zzzzzz_xz, g_0_y_zzzzzz_xzz, g_0_y_zzzzzz_yy, g_0_y_zzzzzz_yz, g_0_y_zzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xzzzzzz_xx[k] = -g_0_y_zzzzzz_xx[k] * ab_x + g_0_y_zzzzzz_xxx[k];

                g_0_y_xzzzzzz_xy[k] = -g_0_y_zzzzzz_xy[k] * ab_x + g_0_y_zzzzzz_xxy[k];

                g_0_y_xzzzzzz_xz[k] = -g_0_y_zzzzzz_xz[k] * ab_x + g_0_y_zzzzzz_xxz[k];

                g_0_y_xzzzzzz_yy[k] = -g_0_y_zzzzzz_yy[k] * ab_x + g_0_y_zzzzzz_xyy[k];

                g_0_y_xzzzzzz_yz[k] = -g_0_y_zzzzzz_yz[k] * ab_x + g_0_y_zzzzzz_xyz[k];

                g_0_y_xzzzzzz_zz[k] = -g_0_y_zzzzzz_zz[k] * ab_x + g_0_y_zzzzzz_xzz[k];
            }

            /// Set up 384-390 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyyyy_xx = cbuffer.data(kd_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_y_yyyyyyy_xy = cbuffer.data(kd_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_y_yyyyyyy_xz = cbuffer.data(kd_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_y_yyyyyyy_yy = cbuffer.data(kd_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_y_yyyyyyy_yz = cbuffer.data(kd_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_y_yyyyyyy_zz = cbuffer.data(kd_geom_01_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyyyy_xx, g_0_y_yyyyyy_xxy, g_0_y_yyyyyy_xy, g_0_y_yyyyyy_xyy, g_0_y_yyyyyy_xyz, g_0_y_yyyyyy_xz, g_0_y_yyyyyy_yy, g_0_y_yyyyyy_yyy, g_0_y_yyyyyy_yyz, g_0_y_yyyyyy_yz, g_0_y_yyyyyy_yzz, g_0_y_yyyyyy_zz, g_0_y_yyyyyyy_xx, g_0_y_yyyyyyy_xy, g_0_y_yyyyyyy_xz, g_0_y_yyyyyyy_yy, g_0_y_yyyyyyy_yz, g_0_y_yyyyyyy_zz, g_yyyyyy_xx, g_yyyyyy_xy, g_yyyyyy_xz, g_yyyyyy_yy, g_yyyyyy_yz, g_yyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyyyy_xx[k] = g_yyyyyy_xx[k] - g_0_y_yyyyyy_xx[k] * ab_y + g_0_y_yyyyyy_xxy[k];

                g_0_y_yyyyyyy_xy[k] = g_yyyyyy_xy[k] - g_0_y_yyyyyy_xy[k] * ab_y + g_0_y_yyyyyy_xyy[k];

                g_0_y_yyyyyyy_xz[k] = g_yyyyyy_xz[k] - g_0_y_yyyyyy_xz[k] * ab_y + g_0_y_yyyyyy_xyz[k];

                g_0_y_yyyyyyy_yy[k] = g_yyyyyy_yy[k] - g_0_y_yyyyyy_yy[k] * ab_y + g_0_y_yyyyyy_yyy[k];

                g_0_y_yyyyyyy_yz[k] = g_yyyyyy_yz[k] - g_0_y_yyyyyy_yz[k] * ab_y + g_0_y_yyyyyy_yyz[k];

                g_0_y_yyyyyyy_zz[k] = g_yyyyyy_zz[k] - g_0_y_yyyyyy_zz[k] * ab_y + g_0_y_yyyyyy_yzz[k];
            }

            /// Set up 390-396 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyyyz_xx = cbuffer.data(kd_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_y_yyyyyyz_xy = cbuffer.data(kd_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_y_yyyyyyz_xz = cbuffer.data(kd_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_y_yyyyyyz_yy = cbuffer.data(kd_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_y_yyyyyyz_yz = cbuffer.data(kd_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_y_yyyyyyz_zz = cbuffer.data(kd_geom_01_off + 395 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyyyy_xx, g_0_y_yyyyyy_xxz, g_0_y_yyyyyy_xy, g_0_y_yyyyyy_xyz, g_0_y_yyyyyy_xz, g_0_y_yyyyyy_xzz, g_0_y_yyyyyy_yy, g_0_y_yyyyyy_yyz, g_0_y_yyyyyy_yz, g_0_y_yyyyyy_yzz, g_0_y_yyyyyy_zz, g_0_y_yyyyyy_zzz, g_0_y_yyyyyyz_xx, g_0_y_yyyyyyz_xy, g_0_y_yyyyyyz_xz, g_0_y_yyyyyyz_yy, g_0_y_yyyyyyz_yz, g_0_y_yyyyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyyyz_xx[k] = -g_0_y_yyyyyy_xx[k] * ab_z + g_0_y_yyyyyy_xxz[k];

                g_0_y_yyyyyyz_xy[k] = -g_0_y_yyyyyy_xy[k] * ab_z + g_0_y_yyyyyy_xyz[k];

                g_0_y_yyyyyyz_xz[k] = -g_0_y_yyyyyy_xz[k] * ab_z + g_0_y_yyyyyy_xzz[k];

                g_0_y_yyyyyyz_yy[k] = -g_0_y_yyyyyy_yy[k] * ab_z + g_0_y_yyyyyy_yyz[k];

                g_0_y_yyyyyyz_yz[k] = -g_0_y_yyyyyy_yz[k] * ab_z + g_0_y_yyyyyy_yzz[k];

                g_0_y_yyyyyyz_zz[k] = -g_0_y_yyyyyy_zz[k] * ab_z + g_0_y_yyyyyy_zzz[k];
            }

            /// Set up 396-402 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyyzz_xx = cbuffer.data(kd_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_y_yyyyyzz_xy = cbuffer.data(kd_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_y_yyyyyzz_xz = cbuffer.data(kd_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_y_yyyyyzz_yy = cbuffer.data(kd_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_y_yyyyyzz_yz = cbuffer.data(kd_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_y_yyyyyzz_zz = cbuffer.data(kd_geom_01_off + 401 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyyyz_xx, g_0_y_yyyyyz_xxz, g_0_y_yyyyyz_xy, g_0_y_yyyyyz_xyz, g_0_y_yyyyyz_xz, g_0_y_yyyyyz_xzz, g_0_y_yyyyyz_yy, g_0_y_yyyyyz_yyz, g_0_y_yyyyyz_yz, g_0_y_yyyyyz_yzz, g_0_y_yyyyyz_zz, g_0_y_yyyyyz_zzz, g_0_y_yyyyyzz_xx, g_0_y_yyyyyzz_xy, g_0_y_yyyyyzz_xz, g_0_y_yyyyyzz_yy, g_0_y_yyyyyzz_yz, g_0_y_yyyyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyyzz_xx[k] = -g_0_y_yyyyyz_xx[k] * ab_z + g_0_y_yyyyyz_xxz[k];

                g_0_y_yyyyyzz_xy[k] = -g_0_y_yyyyyz_xy[k] * ab_z + g_0_y_yyyyyz_xyz[k];

                g_0_y_yyyyyzz_xz[k] = -g_0_y_yyyyyz_xz[k] * ab_z + g_0_y_yyyyyz_xzz[k];

                g_0_y_yyyyyzz_yy[k] = -g_0_y_yyyyyz_yy[k] * ab_z + g_0_y_yyyyyz_yyz[k];

                g_0_y_yyyyyzz_yz[k] = -g_0_y_yyyyyz_yz[k] * ab_z + g_0_y_yyyyyz_yzz[k];

                g_0_y_yyyyyzz_zz[k] = -g_0_y_yyyyyz_zz[k] * ab_z + g_0_y_yyyyyz_zzz[k];
            }

            /// Set up 402-408 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyzzz_xx = cbuffer.data(kd_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_y_yyyyzzz_xy = cbuffer.data(kd_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_y_yyyyzzz_xz = cbuffer.data(kd_geom_01_off + 404 * ccomps * dcomps);

            auto g_0_y_yyyyzzz_yy = cbuffer.data(kd_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_y_yyyyzzz_yz = cbuffer.data(kd_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_y_yyyyzzz_zz = cbuffer.data(kd_geom_01_off + 407 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyyzz_xx, g_0_y_yyyyzz_xxz, g_0_y_yyyyzz_xy, g_0_y_yyyyzz_xyz, g_0_y_yyyyzz_xz, g_0_y_yyyyzz_xzz, g_0_y_yyyyzz_yy, g_0_y_yyyyzz_yyz, g_0_y_yyyyzz_yz, g_0_y_yyyyzz_yzz, g_0_y_yyyyzz_zz, g_0_y_yyyyzz_zzz, g_0_y_yyyyzzz_xx, g_0_y_yyyyzzz_xy, g_0_y_yyyyzzz_xz, g_0_y_yyyyzzz_yy, g_0_y_yyyyzzz_yz, g_0_y_yyyyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyzzz_xx[k] = -g_0_y_yyyyzz_xx[k] * ab_z + g_0_y_yyyyzz_xxz[k];

                g_0_y_yyyyzzz_xy[k] = -g_0_y_yyyyzz_xy[k] * ab_z + g_0_y_yyyyzz_xyz[k];

                g_0_y_yyyyzzz_xz[k] = -g_0_y_yyyyzz_xz[k] * ab_z + g_0_y_yyyyzz_xzz[k];

                g_0_y_yyyyzzz_yy[k] = -g_0_y_yyyyzz_yy[k] * ab_z + g_0_y_yyyyzz_yyz[k];

                g_0_y_yyyyzzz_yz[k] = -g_0_y_yyyyzz_yz[k] * ab_z + g_0_y_yyyyzz_yzz[k];

                g_0_y_yyyyzzz_zz[k] = -g_0_y_yyyyzz_zz[k] * ab_z + g_0_y_yyyyzz_zzz[k];
            }

            /// Set up 408-414 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyzzzz_xx = cbuffer.data(kd_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_y_yyyzzzz_xy = cbuffer.data(kd_geom_01_off + 409 * ccomps * dcomps);

            auto g_0_y_yyyzzzz_xz = cbuffer.data(kd_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_y_yyyzzzz_yy = cbuffer.data(kd_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_y_yyyzzzz_yz = cbuffer.data(kd_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_y_yyyzzzz_zz = cbuffer.data(kd_geom_01_off + 413 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyzzz_xx, g_0_y_yyyzzz_xxz, g_0_y_yyyzzz_xy, g_0_y_yyyzzz_xyz, g_0_y_yyyzzz_xz, g_0_y_yyyzzz_xzz, g_0_y_yyyzzz_yy, g_0_y_yyyzzz_yyz, g_0_y_yyyzzz_yz, g_0_y_yyyzzz_yzz, g_0_y_yyyzzz_zz, g_0_y_yyyzzz_zzz, g_0_y_yyyzzzz_xx, g_0_y_yyyzzzz_xy, g_0_y_yyyzzzz_xz, g_0_y_yyyzzzz_yy, g_0_y_yyyzzzz_yz, g_0_y_yyyzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyzzzz_xx[k] = -g_0_y_yyyzzz_xx[k] * ab_z + g_0_y_yyyzzz_xxz[k];

                g_0_y_yyyzzzz_xy[k] = -g_0_y_yyyzzz_xy[k] * ab_z + g_0_y_yyyzzz_xyz[k];

                g_0_y_yyyzzzz_xz[k] = -g_0_y_yyyzzz_xz[k] * ab_z + g_0_y_yyyzzz_xzz[k];

                g_0_y_yyyzzzz_yy[k] = -g_0_y_yyyzzz_yy[k] * ab_z + g_0_y_yyyzzz_yyz[k];

                g_0_y_yyyzzzz_yz[k] = -g_0_y_yyyzzz_yz[k] * ab_z + g_0_y_yyyzzz_yzz[k];

                g_0_y_yyyzzzz_zz[k] = -g_0_y_yyyzzz_zz[k] * ab_z + g_0_y_yyyzzz_zzz[k];
            }

            /// Set up 414-420 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyzzzzz_xx = cbuffer.data(kd_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_y_yyzzzzz_xy = cbuffer.data(kd_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_y_yyzzzzz_xz = cbuffer.data(kd_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_y_yyzzzzz_yy = cbuffer.data(kd_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_y_yyzzzzz_yz = cbuffer.data(kd_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_y_yyzzzzz_zz = cbuffer.data(kd_geom_01_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyzzzz_xx, g_0_y_yyzzzz_xxz, g_0_y_yyzzzz_xy, g_0_y_yyzzzz_xyz, g_0_y_yyzzzz_xz, g_0_y_yyzzzz_xzz, g_0_y_yyzzzz_yy, g_0_y_yyzzzz_yyz, g_0_y_yyzzzz_yz, g_0_y_yyzzzz_yzz, g_0_y_yyzzzz_zz, g_0_y_yyzzzz_zzz, g_0_y_yyzzzzz_xx, g_0_y_yyzzzzz_xy, g_0_y_yyzzzzz_xz, g_0_y_yyzzzzz_yy, g_0_y_yyzzzzz_yz, g_0_y_yyzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyzzzzz_xx[k] = -g_0_y_yyzzzz_xx[k] * ab_z + g_0_y_yyzzzz_xxz[k];

                g_0_y_yyzzzzz_xy[k] = -g_0_y_yyzzzz_xy[k] * ab_z + g_0_y_yyzzzz_xyz[k];

                g_0_y_yyzzzzz_xz[k] = -g_0_y_yyzzzz_xz[k] * ab_z + g_0_y_yyzzzz_xzz[k];

                g_0_y_yyzzzzz_yy[k] = -g_0_y_yyzzzz_yy[k] * ab_z + g_0_y_yyzzzz_yyz[k];

                g_0_y_yyzzzzz_yz[k] = -g_0_y_yyzzzz_yz[k] * ab_z + g_0_y_yyzzzz_yzz[k];

                g_0_y_yyzzzzz_zz[k] = -g_0_y_yyzzzz_zz[k] * ab_z + g_0_y_yyzzzz_zzz[k];
            }

            /// Set up 420-426 components of targeted buffer : cbuffer.data(

            auto g_0_y_yzzzzzz_xx = cbuffer.data(kd_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_y_yzzzzzz_xy = cbuffer.data(kd_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_y_yzzzzzz_xz = cbuffer.data(kd_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_y_yzzzzzz_yy = cbuffer.data(kd_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_y_yzzzzzz_yz = cbuffer.data(kd_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_y_yzzzzzz_zz = cbuffer.data(kd_geom_01_off + 425 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yzzzzz_xx, g_0_y_yzzzzz_xxz, g_0_y_yzzzzz_xy, g_0_y_yzzzzz_xyz, g_0_y_yzzzzz_xz, g_0_y_yzzzzz_xzz, g_0_y_yzzzzz_yy, g_0_y_yzzzzz_yyz, g_0_y_yzzzzz_yz, g_0_y_yzzzzz_yzz, g_0_y_yzzzzz_zz, g_0_y_yzzzzz_zzz, g_0_y_yzzzzzz_xx, g_0_y_yzzzzzz_xy, g_0_y_yzzzzzz_xz, g_0_y_yzzzzzz_yy, g_0_y_yzzzzzz_yz, g_0_y_yzzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yzzzzzz_xx[k] = -g_0_y_yzzzzz_xx[k] * ab_z + g_0_y_yzzzzz_xxz[k];

                g_0_y_yzzzzzz_xy[k] = -g_0_y_yzzzzz_xy[k] * ab_z + g_0_y_yzzzzz_xyz[k];

                g_0_y_yzzzzzz_xz[k] = -g_0_y_yzzzzz_xz[k] * ab_z + g_0_y_yzzzzz_xzz[k];

                g_0_y_yzzzzzz_yy[k] = -g_0_y_yzzzzz_yy[k] * ab_z + g_0_y_yzzzzz_yyz[k];

                g_0_y_yzzzzzz_yz[k] = -g_0_y_yzzzzz_yz[k] * ab_z + g_0_y_yzzzzz_yzz[k];

                g_0_y_yzzzzzz_zz[k] = -g_0_y_yzzzzz_zz[k] * ab_z + g_0_y_yzzzzz_zzz[k];
            }

            /// Set up 426-432 components of targeted buffer : cbuffer.data(

            auto g_0_y_zzzzzzz_xx = cbuffer.data(kd_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_y_zzzzzzz_xy = cbuffer.data(kd_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_y_zzzzzzz_xz = cbuffer.data(kd_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_y_zzzzzzz_yy = cbuffer.data(kd_geom_01_off + 429 * ccomps * dcomps);

            auto g_0_y_zzzzzzz_yz = cbuffer.data(kd_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_y_zzzzzzz_zz = cbuffer.data(kd_geom_01_off + 431 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_zzzzzz_xx, g_0_y_zzzzzz_xxz, g_0_y_zzzzzz_xy, g_0_y_zzzzzz_xyz, g_0_y_zzzzzz_xz, g_0_y_zzzzzz_xzz, g_0_y_zzzzzz_yy, g_0_y_zzzzzz_yyz, g_0_y_zzzzzz_yz, g_0_y_zzzzzz_yzz, g_0_y_zzzzzz_zz, g_0_y_zzzzzz_zzz, g_0_y_zzzzzzz_xx, g_0_y_zzzzzzz_xy, g_0_y_zzzzzzz_xz, g_0_y_zzzzzzz_yy, g_0_y_zzzzzzz_yz, g_0_y_zzzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_zzzzzzz_xx[k] = -g_0_y_zzzzzz_xx[k] * ab_z + g_0_y_zzzzzz_xxz[k];

                g_0_y_zzzzzzz_xy[k] = -g_0_y_zzzzzz_xy[k] * ab_z + g_0_y_zzzzzz_xyz[k];

                g_0_y_zzzzzzz_xz[k] = -g_0_y_zzzzzz_xz[k] * ab_z + g_0_y_zzzzzz_xzz[k];

                g_0_y_zzzzzzz_yy[k] = -g_0_y_zzzzzz_yy[k] * ab_z + g_0_y_zzzzzz_yyz[k];

                g_0_y_zzzzzzz_yz[k] = -g_0_y_zzzzzz_yz[k] * ab_z + g_0_y_zzzzzz_yzz[k];

                g_0_y_zzzzzzz_zz[k] = -g_0_y_zzzzzz_zz[k] * ab_z + g_0_y_zzzzzz_zzz[k];
            }

            /// Set up 432-438 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxxx_xx = cbuffer.data(kd_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_z_xxxxxxx_xy = cbuffer.data(kd_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_z_xxxxxxx_xz = cbuffer.data(kd_geom_01_off + 434 * ccomps * dcomps);

            auto g_0_z_xxxxxxx_yy = cbuffer.data(kd_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_z_xxxxxxx_yz = cbuffer.data(kd_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_z_xxxxxxx_zz = cbuffer.data(kd_geom_01_off + 437 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxx_xx, g_0_z_xxxxxx_xxx, g_0_z_xxxxxx_xxy, g_0_z_xxxxxx_xxz, g_0_z_xxxxxx_xy, g_0_z_xxxxxx_xyy, g_0_z_xxxxxx_xyz, g_0_z_xxxxxx_xz, g_0_z_xxxxxx_xzz, g_0_z_xxxxxx_yy, g_0_z_xxxxxx_yz, g_0_z_xxxxxx_zz, g_0_z_xxxxxxx_xx, g_0_z_xxxxxxx_xy, g_0_z_xxxxxxx_xz, g_0_z_xxxxxxx_yy, g_0_z_xxxxxxx_yz, g_0_z_xxxxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxxx_xx[k] = -g_0_z_xxxxxx_xx[k] * ab_x + g_0_z_xxxxxx_xxx[k];

                g_0_z_xxxxxxx_xy[k] = -g_0_z_xxxxxx_xy[k] * ab_x + g_0_z_xxxxxx_xxy[k];

                g_0_z_xxxxxxx_xz[k] = -g_0_z_xxxxxx_xz[k] * ab_x + g_0_z_xxxxxx_xxz[k];

                g_0_z_xxxxxxx_yy[k] = -g_0_z_xxxxxx_yy[k] * ab_x + g_0_z_xxxxxx_xyy[k];

                g_0_z_xxxxxxx_yz[k] = -g_0_z_xxxxxx_yz[k] * ab_x + g_0_z_xxxxxx_xyz[k];

                g_0_z_xxxxxxx_zz[k] = -g_0_z_xxxxxx_zz[k] * ab_x + g_0_z_xxxxxx_xzz[k];
            }

            /// Set up 438-444 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxxy_xx = cbuffer.data(kd_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_z_xxxxxxy_xy = cbuffer.data(kd_geom_01_off + 439 * ccomps * dcomps);

            auto g_0_z_xxxxxxy_xz = cbuffer.data(kd_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_z_xxxxxxy_yy = cbuffer.data(kd_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_z_xxxxxxy_yz = cbuffer.data(kd_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_z_xxxxxxy_zz = cbuffer.data(kd_geom_01_off + 443 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxxy_xx, g_0_z_xxxxxxy_xy, g_0_z_xxxxxxy_xz, g_0_z_xxxxxxy_yy, g_0_z_xxxxxxy_yz, g_0_z_xxxxxxy_zz, g_0_z_xxxxxy_xx, g_0_z_xxxxxy_xxx, g_0_z_xxxxxy_xxy, g_0_z_xxxxxy_xxz, g_0_z_xxxxxy_xy, g_0_z_xxxxxy_xyy, g_0_z_xxxxxy_xyz, g_0_z_xxxxxy_xz, g_0_z_xxxxxy_xzz, g_0_z_xxxxxy_yy, g_0_z_xxxxxy_yz, g_0_z_xxxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxxy_xx[k] = -g_0_z_xxxxxy_xx[k] * ab_x + g_0_z_xxxxxy_xxx[k];

                g_0_z_xxxxxxy_xy[k] = -g_0_z_xxxxxy_xy[k] * ab_x + g_0_z_xxxxxy_xxy[k];

                g_0_z_xxxxxxy_xz[k] = -g_0_z_xxxxxy_xz[k] * ab_x + g_0_z_xxxxxy_xxz[k];

                g_0_z_xxxxxxy_yy[k] = -g_0_z_xxxxxy_yy[k] * ab_x + g_0_z_xxxxxy_xyy[k];

                g_0_z_xxxxxxy_yz[k] = -g_0_z_xxxxxy_yz[k] * ab_x + g_0_z_xxxxxy_xyz[k];

                g_0_z_xxxxxxy_zz[k] = -g_0_z_xxxxxy_zz[k] * ab_x + g_0_z_xxxxxy_xzz[k];
            }

            /// Set up 444-450 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxxz_xx = cbuffer.data(kd_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_z_xxxxxxz_xy = cbuffer.data(kd_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_z_xxxxxxz_xz = cbuffer.data(kd_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_z_xxxxxxz_yy = cbuffer.data(kd_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_z_xxxxxxz_yz = cbuffer.data(kd_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_z_xxxxxxz_zz = cbuffer.data(kd_geom_01_off + 449 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxxz_xx, g_0_z_xxxxxxz_xy, g_0_z_xxxxxxz_xz, g_0_z_xxxxxxz_yy, g_0_z_xxxxxxz_yz, g_0_z_xxxxxxz_zz, g_0_z_xxxxxz_xx, g_0_z_xxxxxz_xxx, g_0_z_xxxxxz_xxy, g_0_z_xxxxxz_xxz, g_0_z_xxxxxz_xy, g_0_z_xxxxxz_xyy, g_0_z_xxxxxz_xyz, g_0_z_xxxxxz_xz, g_0_z_xxxxxz_xzz, g_0_z_xxxxxz_yy, g_0_z_xxxxxz_yz, g_0_z_xxxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxxz_xx[k] = -g_0_z_xxxxxz_xx[k] * ab_x + g_0_z_xxxxxz_xxx[k];

                g_0_z_xxxxxxz_xy[k] = -g_0_z_xxxxxz_xy[k] * ab_x + g_0_z_xxxxxz_xxy[k];

                g_0_z_xxxxxxz_xz[k] = -g_0_z_xxxxxz_xz[k] * ab_x + g_0_z_xxxxxz_xxz[k];

                g_0_z_xxxxxxz_yy[k] = -g_0_z_xxxxxz_yy[k] * ab_x + g_0_z_xxxxxz_xyy[k];

                g_0_z_xxxxxxz_yz[k] = -g_0_z_xxxxxz_yz[k] * ab_x + g_0_z_xxxxxz_xyz[k];

                g_0_z_xxxxxxz_zz[k] = -g_0_z_xxxxxz_zz[k] * ab_x + g_0_z_xxxxxz_xzz[k];
            }

            /// Set up 450-456 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxyy_xx = cbuffer.data(kd_geom_01_off + 450 * ccomps * dcomps);

            auto g_0_z_xxxxxyy_xy = cbuffer.data(kd_geom_01_off + 451 * ccomps * dcomps);

            auto g_0_z_xxxxxyy_xz = cbuffer.data(kd_geom_01_off + 452 * ccomps * dcomps);

            auto g_0_z_xxxxxyy_yy = cbuffer.data(kd_geom_01_off + 453 * ccomps * dcomps);

            auto g_0_z_xxxxxyy_yz = cbuffer.data(kd_geom_01_off + 454 * ccomps * dcomps);

            auto g_0_z_xxxxxyy_zz = cbuffer.data(kd_geom_01_off + 455 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxyy_xx, g_0_z_xxxxxyy_xy, g_0_z_xxxxxyy_xz, g_0_z_xxxxxyy_yy, g_0_z_xxxxxyy_yz, g_0_z_xxxxxyy_zz, g_0_z_xxxxyy_xx, g_0_z_xxxxyy_xxx, g_0_z_xxxxyy_xxy, g_0_z_xxxxyy_xxz, g_0_z_xxxxyy_xy, g_0_z_xxxxyy_xyy, g_0_z_xxxxyy_xyz, g_0_z_xxxxyy_xz, g_0_z_xxxxyy_xzz, g_0_z_xxxxyy_yy, g_0_z_xxxxyy_yz, g_0_z_xxxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxyy_xx[k] = -g_0_z_xxxxyy_xx[k] * ab_x + g_0_z_xxxxyy_xxx[k];

                g_0_z_xxxxxyy_xy[k] = -g_0_z_xxxxyy_xy[k] * ab_x + g_0_z_xxxxyy_xxy[k];

                g_0_z_xxxxxyy_xz[k] = -g_0_z_xxxxyy_xz[k] * ab_x + g_0_z_xxxxyy_xxz[k];

                g_0_z_xxxxxyy_yy[k] = -g_0_z_xxxxyy_yy[k] * ab_x + g_0_z_xxxxyy_xyy[k];

                g_0_z_xxxxxyy_yz[k] = -g_0_z_xxxxyy_yz[k] * ab_x + g_0_z_xxxxyy_xyz[k];

                g_0_z_xxxxxyy_zz[k] = -g_0_z_xxxxyy_zz[k] * ab_x + g_0_z_xxxxyy_xzz[k];
            }

            /// Set up 456-462 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxyz_xx = cbuffer.data(kd_geom_01_off + 456 * ccomps * dcomps);

            auto g_0_z_xxxxxyz_xy = cbuffer.data(kd_geom_01_off + 457 * ccomps * dcomps);

            auto g_0_z_xxxxxyz_xz = cbuffer.data(kd_geom_01_off + 458 * ccomps * dcomps);

            auto g_0_z_xxxxxyz_yy = cbuffer.data(kd_geom_01_off + 459 * ccomps * dcomps);

            auto g_0_z_xxxxxyz_yz = cbuffer.data(kd_geom_01_off + 460 * ccomps * dcomps);

            auto g_0_z_xxxxxyz_zz = cbuffer.data(kd_geom_01_off + 461 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxyz_xx, g_0_z_xxxxxyz_xy, g_0_z_xxxxxyz_xz, g_0_z_xxxxxyz_yy, g_0_z_xxxxxyz_yz, g_0_z_xxxxxyz_zz, g_0_z_xxxxyz_xx, g_0_z_xxxxyz_xxx, g_0_z_xxxxyz_xxy, g_0_z_xxxxyz_xxz, g_0_z_xxxxyz_xy, g_0_z_xxxxyz_xyy, g_0_z_xxxxyz_xyz, g_0_z_xxxxyz_xz, g_0_z_xxxxyz_xzz, g_0_z_xxxxyz_yy, g_0_z_xxxxyz_yz, g_0_z_xxxxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxyz_xx[k] = -g_0_z_xxxxyz_xx[k] * ab_x + g_0_z_xxxxyz_xxx[k];

                g_0_z_xxxxxyz_xy[k] = -g_0_z_xxxxyz_xy[k] * ab_x + g_0_z_xxxxyz_xxy[k];

                g_0_z_xxxxxyz_xz[k] = -g_0_z_xxxxyz_xz[k] * ab_x + g_0_z_xxxxyz_xxz[k];

                g_0_z_xxxxxyz_yy[k] = -g_0_z_xxxxyz_yy[k] * ab_x + g_0_z_xxxxyz_xyy[k];

                g_0_z_xxxxxyz_yz[k] = -g_0_z_xxxxyz_yz[k] * ab_x + g_0_z_xxxxyz_xyz[k];

                g_0_z_xxxxxyz_zz[k] = -g_0_z_xxxxyz_zz[k] * ab_x + g_0_z_xxxxyz_xzz[k];
            }

            /// Set up 462-468 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxzz_xx = cbuffer.data(kd_geom_01_off + 462 * ccomps * dcomps);

            auto g_0_z_xxxxxzz_xy = cbuffer.data(kd_geom_01_off + 463 * ccomps * dcomps);

            auto g_0_z_xxxxxzz_xz = cbuffer.data(kd_geom_01_off + 464 * ccomps * dcomps);

            auto g_0_z_xxxxxzz_yy = cbuffer.data(kd_geom_01_off + 465 * ccomps * dcomps);

            auto g_0_z_xxxxxzz_yz = cbuffer.data(kd_geom_01_off + 466 * ccomps * dcomps);

            auto g_0_z_xxxxxzz_zz = cbuffer.data(kd_geom_01_off + 467 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxzz_xx, g_0_z_xxxxxzz_xy, g_0_z_xxxxxzz_xz, g_0_z_xxxxxzz_yy, g_0_z_xxxxxzz_yz, g_0_z_xxxxxzz_zz, g_0_z_xxxxzz_xx, g_0_z_xxxxzz_xxx, g_0_z_xxxxzz_xxy, g_0_z_xxxxzz_xxz, g_0_z_xxxxzz_xy, g_0_z_xxxxzz_xyy, g_0_z_xxxxzz_xyz, g_0_z_xxxxzz_xz, g_0_z_xxxxzz_xzz, g_0_z_xxxxzz_yy, g_0_z_xxxxzz_yz, g_0_z_xxxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxzz_xx[k] = -g_0_z_xxxxzz_xx[k] * ab_x + g_0_z_xxxxzz_xxx[k];

                g_0_z_xxxxxzz_xy[k] = -g_0_z_xxxxzz_xy[k] * ab_x + g_0_z_xxxxzz_xxy[k];

                g_0_z_xxxxxzz_xz[k] = -g_0_z_xxxxzz_xz[k] * ab_x + g_0_z_xxxxzz_xxz[k];

                g_0_z_xxxxxzz_yy[k] = -g_0_z_xxxxzz_yy[k] * ab_x + g_0_z_xxxxzz_xyy[k];

                g_0_z_xxxxxzz_yz[k] = -g_0_z_xxxxzz_yz[k] * ab_x + g_0_z_xxxxzz_xyz[k];

                g_0_z_xxxxxzz_zz[k] = -g_0_z_xxxxzz_zz[k] * ab_x + g_0_z_xxxxzz_xzz[k];
            }

            /// Set up 468-474 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxyyy_xx = cbuffer.data(kd_geom_01_off + 468 * ccomps * dcomps);

            auto g_0_z_xxxxyyy_xy = cbuffer.data(kd_geom_01_off + 469 * ccomps * dcomps);

            auto g_0_z_xxxxyyy_xz = cbuffer.data(kd_geom_01_off + 470 * ccomps * dcomps);

            auto g_0_z_xxxxyyy_yy = cbuffer.data(kd_geom_01_off + 471 * ccomps * dcomps);

            auto g_0_z_xxxxyyy_yz = cbuffer.data(kd_geom_01_off + 472 * ccomps * dcomps);

            auto g_0_z_xxxxyyy_zz = cbuffer.data(kd_geom_01_off + 473 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxyyy_xx, g_0_z_xxxxyyy_xy, g_0_z_xxxxyyy_xz, g_0_z_xxxxyyy_yy, g_0_z_xxxxyyy_yz, g_0_z_xxxxyyy_zz, g_0_z_xxxyyy_xx, g_0_z_xxxyyy_xxx, g_0_z_xxxyyy_xxy, g_0_z_xxxyyy_xxz, g_0_z_xxxyyy_xy, g_0_z_xxxyyy_xyy, g_0_z_xxxyyy_xyz, g_0_z_xxxyyy_xz, g_0_z_xxxyyy_xzz, g_0_z_xxxyyy_yy, g_0_z_xxxyyy_yz, g_0_z_xxxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxyyy_xx[k] = -g_0_z_xxxyyy_xx[k] * ab_x + g_0_z_xxxyyy_xxx[k];

                g_0_z_xxxxyyy_xy[k] = -g_0_z_xxxyyy_xy[k] * ab_x + g_0_z_xxxyyy_xxy[k];

                g_0_z_xxxxyyy_xz[k] = -g_0_z_xxxyyy_xz[k] * ab_x + g_0_z_xxxyyy_xxz[k];

                g_0_z_xxxxyyy_yy[k] = -g_0_z_xxxyyy_yy[k] * ab_x + g_0_z_xxxyyy_xyy[k];

                g_0_z_xxxxyyy_yz[k] = -g_0_z_xxxyyy_yz[k] * ab_x + g_0_z_xxxyyy_xyz[k];

                g_0_z_xxxxyyy_zz[k] = -g_0_z_xxxyyy_zz[k] * ab_x + g_0_z_xxxyyy_xzz[k];
            }

            /// Set up 474-480 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxyyz_xx = cbuffer.data(kd_geom_01_off + 474 * ccomps * dcomps);

            auto g_0_z_xxxxyyz_xy = cbuffer.data(kd_geom_01_off + 475 * ccomps * dcomps);

            auto g_0_z_xxxxyyz_xz = cbuffer.data(kd_geom_01_off + 476 * ccomps * dcomps);

            auto g_0_z_xxxxyyz_yy = cbuffer.data(kd_geom_01_off + 477 * ccomps * dcomps);

            auto g_0_z_xxxxyyz_yz = cbuffer.data(kd_geom_01_off + 478 * ccomps * dcomps);

            auto g_0_z_xxxxyyz_zz = cbuffer.data(kd_geom_01_off + 479 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxyyz_xx, g_0_z_xxxxyyz_xy, g_0_z_xxxxyyz_xz, g_0_z_xxxxyyz_yy, g_0_z_xxxxyyz_yz, g_0_z_xxxxyyz_zz, g_0_z_xxxyyz_xx, g_0_z_xxxyyz_xxx, g_0_z_xxxyyz_xxy, g_0_z_xxxyyz_xxz, g_0_z_xxxyyz_xy, g_0_z_xxxyyz_xyy, g_0_z_xxxyyz_xyz, g_0_z_xxxyyz_xz, g_0_z_xxxyyz_xzz, g_0_z_xxxyyz_yy, g_0_z_xxxyyz_yz, g_0_z_xxxyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxyyz_xx[k] = -g_0_z_xxxyyz_xx[k] * ab_x + g_0_z_xxxyyz_xxx[k];

                g_0_z_xxxxyyz_xy[k] = -g_0_z_xxxyyz_xy[k] * ab_x + g_0_z_xxxyyz_xxy[k];

                g_0_z_xxxxyyz_xz[k] = -g_0_z_xxxyyz_xz[k] * ab_x + g_0_z_xxxyyz_xxz[k];

                g_0_z_xxxxyyz_yy[k] = -g_0_z_xxxyyz_yy[k] * ab_x + g_0_z_xxxyyz_xyy[k];

                g_0_z_xxxxyyz_yz[k] = -g_0_z_xxxyyz_yz[k] * ab_x + g_0_z_xxxyyz_xyz[k];

                g_0_z_xxxxyyz_zz[k] = -g_0_z_xxxyyz_zz[k] * ab_x + g_0_z_xxxyyz_xzz[k];
            }

            /// Set up 480-486 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxyzz_xx = cbuffer.data(kd_geom_01_off + 480 * ccomps * dcomps);

            auto g_0_z_xxxxyzz_xy = cbuffer.data(kd_geom_01_off + 481 * ccomps * dcomps);

            auto g_0_z_xxxxyzz_xz = cbuffer.data(kd_geom_01_off + 482 * ccomps * dcomps);

            auto g_0_z_xxxxyzz_yy = cbuffer.data(kd_geom_01_off + 483 * ccomps * dcomps);

            auto g_0_z_xxxxyzz_yz = cbuffer.data(kd_geom_01_off + 484 * ccomps * dcomps);

            auto g_0_z_xxxxyzz_zz = cbuffer.data(kd_geom_01_off + 485 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxyzz_xx, g_0_z_xxxxyzz_xy, g_0_z_xxxxyzz_xz, g_0_z_xxxxyzz_yy, g_0_z_xxxxyzz_yz, g_0_z_xxxxyzz_zz, g_0_z_xxxyzz_xx, g_0_z_xxxyzz_xxx, g_0_z_xxxyzz_xxy, g_0_z_xxxyzz_xxz, g_0_z_xxxyzz_xy, g_0_z_xxxyzz_xyy, g_0_z_xxxyzz_xyz, g_0_z_xxxyzz_xz, g_0_z_xxxyzz_xzz, g_0_z_xxxyzz_yy, g_0_z_xxxyzz_yz, g_0_z_xxxyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxyzz_xx[k] = -g_0_z_xxxyzz_xx[k] * ab_x + g_0_z_xxxyzz_xxx[k];

                g_0_z_xxxxyzz_xy[k] = -g_0_z_xxxyzz_xy[k] * ab_x + g_0_z_xxxyzz_xxy[k];

                g_0_z_xxxxyzz_xz[k] = -g_0_z_xxxyzz_xz[k] * ab_x + g_0_z_xxxyzz_xxz[k];

                g_0_z_xxxxyzz_yy[k] = -g_0_z_xxxyzz_yy[k] * ab_x + g_0_z_xxxyzz_xyy[k];

                g_0_z_xxxxyzz_yz[k] = -g_0_z_xxxyzz_yz[k] * ab_x + g_0_z_xxxyzz_xyz[k];

                g_0_z_xxxxyzz_zz[k] = -g_0_z_xxxyzz_zz[k] * ab_x + g_0_z_xxxyzz_xzz[k];
            }

            /// Set up 486-492 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxzzz_xx = cbuffer.data(kd_geom_01_off + 486 * ccomps * dcomps);

            auto g_0_z_xxxxzzz_xy = cbuffer.data(kd_geom_01_off + 487 * ccomps * dcomps);

            auto g_0_z_xxxxzzz_xz = cbuffer.data(kd_geom_01_off + 488 * ccomps * dcomps);

            auto g_0_z_xxxxzzz_yy = cbuffer.data(kd_geom_01_off + 489 * ccomps * dcomps);

            auto g_0_z_xxxxzzz_yz = cbuffer.data(kd_geom_01_off + 490 * ccomps * dcomps);

            auto g_0_z_xxxxzzz_zz = cbuffer.data(kd_geom_01_off + 491 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxzzz_xx, g_0_z_xxxxzzz_xy, g_0_z_xxxxzzz_xz, g_0_z_xxxxzzz_yy, g_0_z_xxxxzzz_yz, g_0_z_xxxxzzz_zz, g_0_z_xxxzzz_xx, g_0_z_xxxzzz_xxx, g_0_z_xxxzzz_xxy, g_0_z_xxxzzz_xxz, g_0_z_xxxzzz_xy, g_0_z_xxxzzz_xyy, g_0_z_xxxzzz_xyz, g_0_z_xxxzzz_xz, g_0_z_xxxzzz_xzz, g_0_z_xxxzzz_yy, g_0_z_xxxzzz_yz, g_0_z_xxxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxzzz_xx[k] = -g_0_z_xxxzzz_xx[k] * ab_x + g_0_z_xxxzzz_xxx[k];

                g_0_z_xxxxzzz_xy[k] = -g_0_z_xxxzzz_xy[k] * ab_x + g_0_z_xxxzzz_xxy[k];

                g_0_z_xxxxzzz_xz[k] = -g_0_z_xxxzzz_xz[k] * ab_x + g_0_z_xxxzzz_xxz[k];

                g_0_z_xxxxzzz_yy[k] = -g_0_z_xxxzzz_yy[k] * ab_x + g_0_z_xxxzzz_xyy[k];

                g_0_z_xxxxzzz_yz[k] = -g_0_z_xxxzzz_yz[k] * ab_x + g_0_z_xxxzzz_xyz[k];

                g_0_z_xxxxzzz_zz[k] = -g_0_z_xxxzzz_zz[k] * ab_x + g_0_z_xxxzzz_xzz[k];
            }

            /// Set up 492-498 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyyyy_xx = cbuffer.data(kd_geom_01_off + 492 * ccomps * dcomps);

            auto g_0_z_xxxyyyy_xy = cbuffer.data(kd_geom_01_off + 493 * ccomps * dcomps);

            auto g_0_z_xxxyyyy_xz = cbuffer.data(kd_geom_01_off + 494 * ccomps * dcomps);

            auto g_0_z_xxxyyyy_yy = cbuffer.data(kd_geom_01_off + 495 * ccomps * dcomps);

            auto g_0_z_xxxyyyy_yz = cbuffer.data(kd_geom_01_off + 496 * ccomps * dcomps);

            auto g_0_z_xxxyyyy_zz = cbuffer.data(kd_geom_01_off + 497 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyyyy_xx, g_0_z_xxxyyyy_xy, g_0_z_xxxyyyy_xz, g_0_z_xxxyyyy_yy, g_0_z_xxxyyyy_yz, g_0_z_xxxyyyy_zz, g_0_z_xxyyyy_xx, g_0_z_xxyyyy_xxx, g_0_z_xxyyyy_xxy, g_0_z_xxyyyy_xxz, g_0_z_xxyyyy_xy, g_0_z_xxyyyy_xyy, g_0_z_xxyyyy_xyz, g_0_z_xxyyyy_xz, g_0_z_xxyyyy_xzz, g_0_z_xxyyyy_yy, g_0_z_xxyyyy_yz, g_0_z_xxyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyyyy_xx[k] = -g_0_z_xxyyyy_xx[k] * ab_x + g_0_z_xxyyyy_xxx[k];

                g_0_z_xxxyyyy_xy[k] = -g_0_z_xxyyyy_xy[k] * ab_x + g_0_z_xxyyyy_xxy[k];

                g_0_z_xxxyyyy_xz[k] = -g_0_z_xxyyyy_xz[k] * ab_x + g_0_z_xxyyyy_xxz[k];

                g_0_z_xxxyyyy_yy[k] = -g_0_z_xxyyyy_yy[k] * ab_x + g_0_z_xxyyyy_xyy[k];

                g_0_z_xxxyyyy_yz[k] = -g_0_z_xxyyyy_yz[k] * ab_x + g_0_z_xxyyyy_xyz[k];

                g_0_z_xxxyyyy_zz[k] = -g_0_z_xxyyyy_zz[k] * ab_x + g_0_z_xxyyyy_xzz[k];
            }

            /// Set up 498-504 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyyyz_xx = cbuffer.data(kd_geom_01_off + 498 * ccomps * dcomps);

            auto g_0_z_xxxyyyz_xy = cbuffer.data(kd_geom_01_off + 499 * ccomps * dcomps);

            auto g_0_z_xxxyyyz_xz = cbuffer.data(kd_geom_01_off + 500 * ccomps * dcomps);

            auto g_0_z_xxxyyyz_yy = cbuffer.data(kd_geom_01_off + 501 * ccomps * dcomps);

            auto g_0_z_xxxyyyz_yz = cbuffer.data(kd_geom_01_off + 502 * ccomps * dcomps);

            auto g_0_z_xxxyyyz_zz = cbuffer.data(kd_geom_01_off + 503 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyyyz_xx, g_0_z_xxxyyyz_xy, g_0_z_xxxyyyz_xz, g_0_z_xxxyyyz_yy, g_0_z_xxxyyyz_yz, g_0_z_xxxyyyz_zz, g_0_z_xxyyyz_xx, g_0_z_xxyyyz_xxx, g_0_z_xxyyyz_xxy, g_0_z_xxyyyz_xxz, g_0_z_xxyyyz_xy, g_0_z_xxyyyz_xyy, g_0_z_xxyyyz_xyz, g_0_z_xxyyyz_xz, g_0_z_xxyyyz_xzz, g_0_z_xxyyyz_yy, g_0_z_xxyyyz_yz, g_0_z_xxyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyyyz_xx[k] = -g_0_z_xxyyyz_xx[k] * ab_x + g_0_z_xxyyyz_xxx[k];

                g_0_z_xxxyyyz_xy[k] = -g_0_z_xxyyyz_xy[k] * ab_x + g_0_z_xxyyyz_xxy[k];

                g_0_z_xxxyyyz_xz[k] = -g_0_z_xxyyyz_xz[k] * ab_x + g_0_z_xxyyyz_xxz[k];

                g_0_z_xxxyyyz_yy[k] = -g_0_z_xxyyyz_yy[k] * ab_x + g_0_z_xxyyyz_xyy[k];

                g_0_z_xxxyyyz_yz[k] = -g_0_z_xxyyyz_yz[k] * ab_x + g_0_z_xxyyyz_xyz[k];

                g_0_z_xxxyyyz_zz[k] = -g_0_z_xxyyyz_zz[k] * ab_x + g_0_z_xxyyyz_xzz[k];
            }

            /// Set up 504-510 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyyzz_xx = cbuffer.data(kd_geom_01_off + 504 * ccomps * dcomps);

            auto g_0_z_xxxyyzz_xy = cbuffer.data(kd_geom_01_off + 505 * ccomps * dcomps);

            auto g_0_z_xxxyyzz_xz = cbuffer.data(kd_geom_01_off + 506 * ccomps * dcomps);

            auto g_0_z_xxxyyzz_yy = cbuffer.data(kd_geom_01_off + 507 * ccomps * dcomps);

            auto g_0_z_xxxyyzz_yz = cbuffer.data(kd_geom_01_off + 508 * ccomps * dcomps);

            auto g_0_z_xxxyyzz_zz = cbuffer.data(kd_geom_01_off + 509 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyyzz_xx, g_0_z_xxxyyzz_xy, g_0_z_xxxyyzz_xz, g_0_z_xxxyyzz_yy, g_0_z_xxxyyzz_yz, g_0_z_xxxyyzz_zz, g_0_z_xxyyzz_xx, g_0_z_xxyyzz_xxx, g_0_z_xxyyzz_xxy, g_0_z_xxyyzz_xxz, g_0_z_xxyyzz_xy, g_0_z_xxyyzz_xyy, g_0_z_xxyyzz_xyz, g_0_z_xxyyzz_xz, g_0_z_xxyyzz_xzz, g_0_z_xxyyzz_yy, g_0_z_xxyyzz_yz, g_0_z_xxyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyyzz_xx[k] = -g_0_z_xxyyzz_xx[k] * ab_x + g_0_z_xxyyzz_xxx[k];

                g_0_z_xxxyyzz_xy[k] = -g_0_z_xxyyzz_xy[k] * ab_x + g_0_z_xxyyzz_xxy[k];

                g_0_z_xxxyyzz_xz[k] = -g_0_z_xxyyzz_xz[k] * ab_x + g_0_z_xxyyzz_xxz[k];

                g_0_z_xxxyyzz_yy[k] = -g_0_z_xxyyzz_yy[k] * ab_x + g_0_z_xxyyzz_xyy[k];

                g_0_z_xxxyyzz_yz[k] = -g_0_z_xxyyzz_yz[k] * ab_x + g_0_z_xxyyzz_xyz[k];

                g_0_z_xxxyyzz_zz[k] = -g_0_z_xxyyzz_zz[k] * ab_x + g_0_z_xxyyzz_xzz[k];
            }

            /// Set up 510-516 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyzzz_xx = cbuffer.data(kd_geom_01_off + 510 * ccomps * dcomps);

            auto g_0_z_xxxyzzz_xy = cbuffer.data(kd_geom_01_off + 511 * ccomps * dcomps);

            auto g_0_z_xxxyzzz_xz = cbuffer.data(kd_geom_01_off + 512 * ccomps * dcomps);

            auto g_0_z_xxxyzzz_yy = cbuffer.data(kd_geom_01_off + 513 * ccomps * dcomps);

            auto g_0_z_xxxyzzz_yz = cbuffer.data(kd_geom_01_off + 514 * ccomps * dcomps);

            auto g_0_z_xxxyzzz_zz = cbuffer.data(kd_geom_01_off + 515 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyzzz_xx, g_0_z_xxxyzzz_xy, g_0_z_xxxyzzz_xz, g_0_z_xxxyzzz_yy, g_0_z_xxxyzzz_yz, g_0_z_xxxyzzz_zz, g_0_z_xxyzzz_xx, g_0_z_xxyzzz_xxx, g_0_z_xxyzzz_xxy, g_0_z_xxyzzz_xxz, g_0_z_xxyzzz_xy, g_0_z_xxyzzz_xyy, g_0_z_xxyzzz_xyz, g_0_z_xxyzzz_xz, g_0_z_xxyzzz_xzz, g_0_z_xxyzzz_yy, g_0_z_xxyzzz_yz, g_0_z_xxyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyzzz_xx[k] = -g_0_z_xxyzzz_xx[k] * ab_x + g_0_z_xxyzzz_xxx[k];

                g_0_z_xxxyzzz_xy[k] = -g_0_z_xxyzzz_xy[k] * ab_x + g_0_z_xxyzzz_xxy[k];

                g_0_z_xxxyzzz_xz[k] = -g_0_z_xxyzzz_xz[k] * ab_x + g_0_z_xxyzzz_xxz[k];

                g_0_z_xxxyzzz_yy[k] = -g_0_z_xxyzzz_yy[k] * ab_x + g_0_z_xxyzzz_xyy[k];

                g_0_z_xxxyzzz_yz[k] = -g_0_z_xxyzzz_yz[k] * ab_x + g_0_z_xxyzzz_xyz[k];

                g_0_z_xxxyzzz_zz[k] = -g_0_z_xxyzzz_zz[k] * ab_x + g_0_z_xxyzzz_xzz[k];
            }

            /// Set up 516-522 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxzzzz_xx = cbuffer.data(kd_geom_01_off + 516 * ccomps * dcomps);

            auto g_0_z_xxxzzzz_xy = cbuffer.data(kd_geom_01_off + 517 * ccomps * dcomps);

            auto g_0_z_xxxzzzz_xz = cbuffer.data(kd_geom_01_off + 518 * ccomps * dcomps);

            auto g_0_z_xxxzzzz_yy = cbuffer.data(kd_geom_01_off + 519 * ccomps * dcomps);

            auto g_0_z_xxxzzzz_yz = cbuffer.data(kd_geom_01_off + 520 * ccomps * dcomps);

            auto g_0_z_xxxzzzz_zz = cbuffer.data(kd_geom_01_off + 521 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxzzzz_xx, g_0_z_xxxzzzz_xy, g_0_z_xxxzzzz_xz, g_0_z_xxxzzzz_yy, g_0_z_xxxzzzz_yz, g_0_z_xxxzzzz_zz, g_0_z_xxzzzz_xx, g_0_z_xxzzzz_xxx, g_0_z_xxzzzz_xxy, g_0_z_xxzzzz_xxz, g_0_z_xxzzzz_xy, g_0_z_xxzzzz_xyy, g_0_z_xxzzzz_xyz, g_0_z_xxzzzz_xz, g_0_z_xxzzzz_xzz, g_0_z_xxzzzz_yy, g_0_z_xxzzzz_yz, g_0_z_xxzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxzzzz_xx[k] = -g_0_z_xxzzzz_xx[k] * ab_x + g_0_z_xxzzzz_xxx[k];

                g_0_z_xxxzzzz_xy[k] = -g_0_z_xxzzzz_xy[k] * ab_x + g_0_z_xxzzzz_xxy[k];

                g_0_z_xxxzzzz_xz[k] = -g_0_z_xxzzzz_xz[k] * ab_x + g_0_z_xxzzzz_xxz[k];

                g_0_z_xxxzzzz_yy[k] = -g_0_z_xxzzzz_yy[k] * ab_x + g_0_z_xxzzzz_xyy[k];

                g_0_z_xxxzzzz_yz[k] = -g_0_z_xxzzzz_yz[k] * ab_x + g_0_z_xxzzzz_xyz[k];

                g_0_z_xxxzzzz_zz[k] = -g_0_z_xxzzzz_zz[k] * ab_x + g_0_z_xxzzzz_xzz[k];
            }

            /// Set up 522-528 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyyyy_xx = cbuffer.data(kd_geom_01_off + 522 * ccomps * dcomps);

            auto g_0_z_xxyyyyy_xy = cbuffer.data(kd_geom_01_off + 523 * ccomps * dcomps);

            auto g_0_z_xxyyyyy_xz = cbuffer.data(kd_geom_01_off + 524 * ccomps * dcomps);

            auto g_0_z_xxyyyyy_yy = cbuffer.data(kd_geom_01_off + 525 * ccomps * dcomps);

            auto g_0_z_xxyyyyy_yz = cbuffer.data(kd_geom_01_off + 526 * ccomps * dcomps);

            auto g_0_z_xxyyyyy_zz = cbuffer.data(kd_geom_01_off + 527 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyyyy_xx, g_0_z_xxyyyyy_xy, g_0_z_xxyyyyy_xz, g_0_z_xxyyyyy_yy, g_0_z_xxyyyyy_yz, g_0_z_xxyyyyy_zz, g_0_z_xyyyyy_xx, g_0_z_xyyyyy_xxx, g_0_z_xyyyyy_xxy, g_0_z_xyyyyy_xxz, g_0_z_xyyyyy_xy, g_0_z_xyyyyy_xyy, g_0_z_xyyyyy_xyz, g_0_z_xyyyyy_xz, g_0_z_xyyyyy_xzz, g_0_z_xyyyyy_yy, g_0_z_xyyyyy_yz, g_0_z_xyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyyyy_xx[k] = -g_0_z_xyyyyy_xx[k] * ab_x + g_0_z_xyyyyy_xxx[k];

                g_0_z_xxyyyyy_xy[k] = -g_0_z_xyyyyy_xy[k] * ab_x + g_0_z_xyyyyy_xxy[k];

                g_0_z_xxyyyyy_xz[k] = -g_0_z_xyyyyy_xz[k] * ab_x + g_0_z_xyyyyy_xxz[k];

                g_0_z_xxyyyyy_yy[k] = -g_0_z_xyyyyy_yy[k] * ab_x + g_0_z_xyyyyy_xyy[k];

                g_0_z_xxyyyyy_yz[k] = -g_0_z_xyyyyy_yz[k] * ab_x + g_0_z_xyyyyy_xyz[k];

                g_0_z_xxyyyyy_zz[k] = -g_0_z_xyyyyy_zz[k] * ab_x + g_0_z_xyyyyy_xzz[k];
            }

            /// Set up 528-534 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyyyz_xx = cbuffer.data(kd_geom_01_off + 528 * ccomps * dcomps);

            auto g_0_z_xxyyyyz_xy = cbuffer.data(kd_geom_01_off + 529 * ccomps * dcomps);

            auto g_0_z_xxyyyyz_xz = cbuffer.data(kd_geom_01_off + 530 * ccomps * dcomps);

            auto g_0_z_xxyyyyz_yy = cbuffer.data(kd_geom_01_off + 531 * ccomps * dcomps);

            auto g_0_z_xxyyyyz_yz = cbuffer.data(kd_geom_01_off + 532 * ccomps * dcomps);

            auto g_0_z_xxyyyyz_zz = cbuffer.data(kd_geom_01_off + 533 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyyyz_xx, g_0_z_xxyyyyz_xy, g_0_z_xxyyyyz_xz, g_0_z_xxyyyyz_yy, g_0_z_xxyyyyz_yz, g_0_z_xxyyyyz_zz, g_0_z_xyyyyz_xx, g_0_z_xyyyyz_xxx, g_0_z_xyyyyz_xxy, g_0_z_xyyyyz_xxz, g_0_z_xyyyyz_xy, g_0_z_xyyyyz_xyy, g_0_z_xyyyyz_xyz, g_0_z_xyyyyz_xz, g_0_z_xyyyyz_xzz, g_0_z_xyyyyz_yy, g_0_z_xyyyyz_yz, g_0_z_xyyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyyyz_xx[k] = -g_0_z_xyyyyz_xx[k] * ab_x + g_0_z_xyyyyz_xxx[k];

                g_0_z_xxyyyyz_xy[k] = -g_0_z_xyyyyz_xy[k] * ab_x + g_0_z_xyyyyz_xxy[k];

                g_0_z_xxyyyyz_xz[k] = -g_0_z_xyyyyz_xz[k] * ab_x + g_0_z_xyyyyz_xxz[k];

                g_0_z_xxyyyyz_yy[k] = -g_0_z_xyyyyz_yy[k] * ab_x + g_0_z_xyyyyz_xyy[k];

                g_0_z_xxyyyyz_yz[k] = -g_0_z_xyyyyz_yz[k] * ab_x + g_0_z_xyyyyz_xyz[k];

                g_0_z_xxyyyyz_zz[k] = -g_0_z_xyyyyz_zz[k] * ab_x + g_0_z_xyyyyz_xzz[k];
            }

            /// Set up 534-540 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyyzz_xx = cbuffer.data(kd_geom_01_off + 534 * ccomps * dcomps);

            auto g_0_z_xxyyyzz_xy = cbuffer.data(kd_geom_01_off + 535 * ccomps * dcomps);

            auto g_0_z_xxyyyzz_xz = cbuffer.data(kd_geom_01_off + 536 * ccomps * dcomps);

            auto g_0_z_xxyyyzz_yy = cbuffer.data(kd_geom_01_off + 537 * ccomps * dcomps);

            auto g_0_z_xxyyyzz_yz = cbuffer.data(kd_geom_01_off + 538 * ccomps * dcomps);

            auto g_0_z_xxyyyzz_zz = cbuffer.data(kd_geom_01_off + 539 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyyzz_xx, g_0_z_xxyyyzz_xy, g_0_z_xxyyyzz_xz, g_0_z_xxyyyzz_yy, g_0_z_xxyyyzz_yz, g_0_z_xxyyyzz_zz, g_0_z_xyyyzz_xx, g_0_z_xyyyzz_xxx, g_0_z_xyyyzz_xxy, g_0_z_xyyyzz_xxz, g_0_z_xyyyzz_xy, g_0_z_xyyyzz_xyy, g_0_z_xyyyzz_xyz, g_0_z_xyyyzz_xz, g_0_z_xyyyzz_xzz, g_0_z_xyyyzz_yy, g_0_z_xyyyzz_yz, g_0_z_xyyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyyzz_xx[k] = -g_0_z_xyyyzz_xx[k] * ab_x + g_0_z_xyyyzz_xxx[k];

                g_0_z_xxyyyzz_xy[k] = -g_0_z_xyyyzz_xy[k] * ab_x + g_0_z_xyyyzz_xxy[k];

                g_0_z_xxyyyzz_xz[k] = -g_0_z_xyyyzz_xz[k] * ab_x + g_0_z_xyyyzz_xxz[k];

                g_0_z_xxyyyzz_yy[k] = -g_0_z_xyyyzz_yy[k] * ab_x + g_0_z_xyyyzz_xyy[k];

                g_0_z_xxyyyzz_yz[k] = -g_0_z_xyyyzz_yz[k] * ab_x + g_0_z_xyyyzz_xyz[k];

                g_0_z_xxyyyzz_zz[k] = -g_0_z_xyyyzz_zz[k] * ab_x + g_0_z_xyyyzz_xzz[k];
            }

            /// Set up 540-546 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyzzz_xx = cbuffer.data(kd_geom_01_off + 540 * ccomps * dcomps);

            auto g_0_z_xxyyzzz_xy = cbuffer.data(kd_geom_01_off + 541 * ccomps * dcomps);

            auto g_0_z_xxyyzzz_xz = cbuffer.data(kd_geom_01_off + 542 * ccomps * dcomps);

            auto g_0_z_xxyyzzz_yy = cbuffer.data(kd_geom_01_off + 543 * ccomps * dcomps);

            auto g_0_z_xxyyzzz_yz = cbuffer.data(kd_geom_01_off + 544 * ccomps * dcomps);

            auto g_0_z_xxyyzzz_zz = cbuffer.data(kd_geom_01_off + 545 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyzzz_xx, g_0_z_xxyyzzz_xy, g_0_z_xxyyzzz_xz, g_0_z_xxyyzzz_yy, g_0_z_xxyyzzz_yz, g_0_z_xxyyzzz_zz, g_0_z_xyyzzz_xx, g_0_z_xyyzzz_xxx, g_0_z_xyyzzz_xxy, g_0_z_xyyzzz_xxz, g_0_z_xyyzzz_xy, g_0_z_xyyzzz_xyy, g_0_z_xyyzzz_xyz, g_0_z_xyyzzz_xz, g_0_z_xyyzzz_xzz, g_0_z_xyyzzz_yy, g_0_z_xyyzzz_yz, g_0_z_xyyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyzzz_xx[k] = -g_0_z_xyyzzz_xx[k] * ab_x + g_0_z_xyyzzz_xxx[k];

                g_0_z_xxyyzzz_xy[k] = -g_0_z_xyyzzz_xy[k] * ab_x + g_0_z_xyyzzz_xxy[k];

                g_0_z_xxyyzzz_xz[k] = -g_0_z_xyyzzz_xz[k] * ab_x + g_0_z_xyyzzz_xxz[k];

                g_0_z_xxyyzzz_yy[k] = -g_0_z_xyyzzz_yy[k] * ab_x + g_0_z_xyyzzz_xyy[k];

                g_0_z_xxyyzzz_yz[k] = -g_0_z_xyyzzz_yz[k] * ab_x + g_0_z_xyyzzz_xyz[k];

                g_0_z_xxyyzzz_zz[k] = -g_0_z_xyyzzz_zz[k] * ab_x + g_0_z_xyyzzz_xzz[k];
            }

            /// Set up 546-552 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyzzzz_xx = cbuffer.data(kd_geom_01_off + 546 * ccomps * dcomps);

            auto g_0_z_xxyzzzz_xy = cbuffer.data(kd_geom_01_off + 547 * ccomps * dcomps);

            auto g_0_z_xxyzzzz_xz = cbuffer.data(kd_geom_01_off + 548 * ccomps * dcomps);

            auto g_0_z_xxyzzzz_yy = cbuffer.data(kd_geom_01_off + 549 * ccomps * dcomps);

            auto g_0_z_xxyzzzz_yz = cbuffer.data(kd_geom_01_off + 550 * ccomps * dcomps);

            auto g_0_z_xxyzzzz_zz = cbuffer.data(kd_geom_01_off + 551 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyzzzz_xx, g_0_z_xxyzzzz_xy, g_0_z_xxyzzzz_xz, g_0_z_xxyzzzz_yy, g_0_z_xxyzzzz_yz, g_0_z_xxyzzzz_zz, g_0_z_xyzzzz_xx, g_0_z_xyzzzz_xxx, g_0_z_xyzzzz_xxy, g_0_z_xyzzzz_xxz, g_0_z_xyzzzz_xy, g_0_z_xyzzzz_xyy, g_0_z_xyzzzz_xyz, g_0_z_xyzzzz_xz, g_0_z_xyzzzz_xzz, g_0_z_xyzzzz_yy, g_0_z_xyzzzz_yz, g_0_z_xyzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyzzzz_xx[k] = -g_0_z_xyzzzz_xx[k] * ab_x + g_0_z_xyzzzz_xxx[k];

                g_0_z_xxyzzzz_xy[k] = -g_0_z_xyzzzz_xy[k] * ab_x + g_0_z_xyzzzz_xxy[k];

                g_0_z_xxyzzzz_xz[k] = -g_0_z_xyzzzz_xz[k] * ab_x + g_0_z_xyzzzz_xxz[k];

                g_0_z_xxyzzzz_yy[k] = -g_0_z_xyzzzz_yy[k] * ab_x + g_0_z_xyzzzz_xyy[k];

                g_0_z_xxyzzzz_yz[k] = -g_0_z_xyzzzz_yz[k] * ab_x + g_0_z_xyzzzz_xyz[k];

                g_0_z_xxyzzzz_zz[k] = -g_0_z_xyzzzz_zz[k] * ab_x + g_0_z_xyzzzz_xzz[k];
            }

            /// Set up 552-558 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxzzzzz_xx = cbuffer.data(kd_geom_01_off + 552 * ccomps * dcomps);

            auto g_0_z_xxzzzzz_xy = cbuffer.data(kd_geom_01_off + 553 * ccomps * dcomps);

            auto g_0_z_xxzzzzz_xz = cbuffer.data(kd_geom_01_off + 554 * ccomps * dcomps);

            auto g_0_z_xxzzzzz_yy = cbuffer.data(kd_geom_01_off + 555 * ccomps * dcomps);

            auto g_0_z_xxzzzzz_yz = cbuffer.data(kd_geom_01_off + 556 * ccomps * dcomps);

            auto g_0_z_xxzzzzz_zz = cbuffer.data(kd_geom_01_off + 557 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxzzzzz_xx, g_0_z_xxzzzzz_xy, g_0_z_xxzzzzz_xz, g_0_z_xxzzzzz_yy, g_0_z_xxzzzzz_yz, g_0_z_xxzzzzz_zz, g_0_z_xzzzzz_xx, g_0_z_xzzzzz_xxx, g_0_z_xzzzzz_xxy, g_0_z_xzzzzz_xxz, g_0_z_xzzzzz_xy, g_0_z_xzzzzz_xyy, g_0_z_xzzzzz_xyz, g_0_z_xzzzzz_xz, g_0_z_xzzzzz_xzz, g_0_z_xzzzzz_yy, g_0_z_xzzzzz_yz, g_0_z_xzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxzzzzz_xx[k] = -g_0_z_xzzzzz_xx[k] * ab_x + g_0_z_xzzzzz_xxx[k];

                g_0_z_xxzzzzz_xy[k] = -g_0_z_xzzzzz_xy[k] * ab_x + g_0_z_xzzzzz_xxy[k];

                g_0_z_xxzzzzz_xz[k] = -g_0_z_xzzzzz_xz[k] * ab_x + g_0_z_xzzzzz_xxz[k];

                g_0_z_xxzzzzz_yy[k] = -g_0_z_xzzzzz_yy[k] * ab_x + g_0_z_xzzzzz_xyy[k];

                g_0_z_xxzzzzz_yz[k] = -g_0_z_xzzzzz_yz[k] * ab_x + g_0_z_xzzzzz_xyz[k];

                g_0_z_xxzzzzz_zz[k] = -g_0_z_xzzzzz_zz[k] * ab_x + g_0_z_xzzzzz_xzz[k];
            }

            /// Set up 558-564 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyyyy_xx = cbuffer.data(kd_geom_01_off + 558 * ccomps * dcomps);

            auto g_0_z_xyyyyyy_xy = cbuffer.data(kd_geom_01_off + 559 * ccomps * dcomps);

            auto g_0_z_xyyyyyy_xz = cbuffer.data(kd_geom_01_off + 560 * ccomps * dcomps);

            auto g_0_z_xyyyyyy_yy = cbuffer.data(kd_geom_01_off + 561 * ccomps * dcomps);

            auto g_0_z_xyyyyyy_yz = cbuffer.data(kd_geom_01_off + 562 * ccomps * dcomps);

            auto g_0_z_xyyyyyy_zz = cbuffer.data(kd_geom_01_off + 563 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyyyy_xx, g_0_z_xyyyyyy_xy, g_0_z_xyyyyyy_xz, g_0_z_xyyyyyy_yy, g_0_z_xyyyyyy_yz, g_0_z_xyyyyyy_zz, g_0_z_yyyyyy_xx, g_0_z_yyyyyy_xxx, g_0_z_yyyyyy_xxy, g_0_z_yyyyyy_xxz, g_0_z_yyyyyy_xy, g_0_z_yyyyyy_xyy, g_0_z_yyyyyy_xyz, g_0_z_yyyyyy_xz, g_0_z_yyyyyy_xzz, g_0_z_yyyyyy_yy, g_0_z_yyyyyy_yz, g_0_z_yyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyyyy_xx[k] = -g_0_z_yyyyyy_xx[k] * ab_x + g_0_z_yyyyyy_xxx[k];

                g_0_z_xyyyyyy_xy[k] = -g_0_z_yyyyyy_xy[k] * ab_x + g_0_z_yyyyyy_xxy[k];

                g_0_z_xyyyyyy_xz[k] = -g_0_z_yyyyyy_xz[k] * ab_x + g_0_z_yyyyyy_xxz[k];

                g_0_z_xyyyyyy_yy[k] = -g_0_z_yyyyyy_yy[k] * ab_x + g_0_z_yyyyyy_xyy[k];

                g_0_z_xyyyyyy_yz[k] = -g_0_z_yyyyyy_yz[k] * ab_x + g_0_z_yyyyyy_xyz[k];

                g_0_z_xyyyyyy_zz[k] = -g_0_z_yyyyyy_zz[k] * ab_x + g_0_z_yyyyyy_xzz[k];
            }

            /// Set up 564-570 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyyyz_xx = cbuffer.data(kd_geom_01_off + 564 * ccomps * dcomps);

            auto g_0_z_xyyyyyz_xy = cbuffer.data(kd_geom_01_off + 565 * ccomps * dcomps);

            auto g_0_z_xyyyyyz_xz = cbuffer.data(kd_geom_01_off + 566 * ccomps * dcomps);

            auto g_0_z_xyyyyyz_yy = cbuffer.data(kd_geom_01_off + 567 * ccomps * dcomps);

            auto g_0_z_xyyyyyz_yz = cbuffer.data(kd_geom_01_off + 568 * ccomps * dcomps);

            auto g_0_z_xyyyyyz_zz = cbuffer.data(kd_geom_01_off + 569 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyyyz_xx, g_0_z_xyyyyyz_xy, g_0_z_xyyyyyz_xz, g_0_z_xyyyyyz_yy, g_0_z_xyyyyyz_yz, g_0_z_xyyyyyz_zz, g_0_z_yyyyyz_xx, g_0_z_yyyyyz_xxx, g_0_z_yyyyyz_xxy, g_0_z_yyyyyz_xxz, g_0_z_yyyyyz_xy, g_0_z_yyyyyz_xyy, g_0_z_yyyyyz_xyz, g_0_z_yyyyyz_xz, g_0_z_yyyyyz_xzz, g_0_z_yyyyyz_yy, g_0_z_yyyyyz_yz, g_0_z_yyyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyyyz_xx[k] = -g_0_z_yyyyyz_xx[k] * ab_x + g_0_z_yyyyyz_xxx[k];

                g_0_z_xyyyyyz_xy[k] = -g_0_z_yyyyyz_xy[k] * ab_x + g_0_z_yyyyyz_xxy[k];

                g_0_z_xyyyyyz_xz[k] = -g_0_z_yyyyyz_xz[k] * ab_x + g_0_z_yyyyyz_xxz[k];

                g_0_z_xyyyyyz_yy[k] = -g_0_z_yyyyyz_yy[k] * ab_x + g_0_z_yyyyyz_xyy[k];

                g_0_z_xyyyyyz_yz[k] = -g_0_z_yyyyyz_yz[k] * ab_x + g_0_z_yyyyyz_xyz[k];

                g_0_z_xyyyyyz_zz[k] = -g_0_z_yyyyyz_zz[k] * ab_x + g_0_z_yyyyyz_xzz[k];
            }

            /// Set up 570-576 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyyzz_xx = cbuffer.data(kd_geom_01_off + 570 * ccomps * dcomps);

            auto g_0_z_xyyyyzz_xy = cbuffer.data(kd_geom_01_off + 571 * ccomps * dcomps);

            auto g_0_z_xyyyyzz_xz = cbuffer.data(kd_geom_01_off + 572 * ccomps * dcomps);

            auto g_0_z_xyyyyzz_yy = cbuffer.data(kd_geom_01_off + 573 * ccomps * dcomps);

            auto g_0_z_xyyyyzz_yz = cbuffer.data(kd_geom_01_off + 574 * ccomps * dcomps);

            auto g_0_z_xyyyyzz_zz = cbuffer.data(kd_geom_01_off + 575 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyyzz_xx, g_0_z_xyyyyzz_xy, g_0_z_xyyyyzz_xz, g_0_z_xyyyyzz_yy, g_0_z_xyyyyzz_yz, g_0_z_xyyyyzz_zz, g_0_z_yyyyzz_xx, g_0_z_yyyyzz_xxx, g_0_z_yyyyzz_xxy, g_0_z_yyyyzz_xxz, g_0_z_yyyyzz_xy, g_0_z_yyyyzz_xyy, g_0_z_yyyyzz_xyz, g_0_z_yyyyzz_xz, g_0_z_yyyyzz_xzz, g_0_z_yyyyzz_yy, g_0_z_yyyyzz_yz, g_0_z_yyyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyyzz_xx[k] = -g_0_z_yyyyzz_xx[k] * ab_x + g_0_z_yyyyzz_xxx[k];

                g_0_z_xyyyyzz_xy[k] = -g_0_z_yyyyzz_xy[k] * ab_x + g_0_z_yyyyzz_xxy[k];

                g_0_z_xyyyyzz_xz[k] = -g_0_z_yyyyzz_xz[k] * ab_x + g_0_z_yyyyzz_xxz[k];

                g_0_z_xyyyyzz_yy[k] = -g_0_z_yyyyzz_yy[k] * ab_x + g_0_z_yyyyzz_xyy[k];

                g_0_z_xyyyyzz_yz[k] = -g_0_z_yyyyzz_yz[k] * ab_x + g_0_z_yyyyzz_xyz[k];

                g_0_z_xyyyyzz_zz[k] = -g_0_z_yyyyzz_zz[k] * ab_x + g_0_z_yyyyzz_xzz[k];
            }

            /// Set up 576-582 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyzzz_xx = cbuffer.data(kd_geom_01_off + 576 * ccomps * dcomps);

            auto g_0_z_xyyyzzz_xy = cbuffer.data(kd_geom_01_off + 577 * ccomps * dcomps);

            auto g_0_z_xyyyzzz_xz = cbuffer.data(kd_geom_01_off + 578 * ccomps * dcomps);

            auto g_0_z_xyyyzzz_yy = cbuffer.data(kd_geom_01_off + 579 * ccomps * dcomps);

            auto g_0_z_xyyyzzz_yz = cbuffer.data(kd_geom_01_off + 580 * ccomps * dcomps);

            auto g_0_z_xyyyzzz_zz = cbuffer.data(kd_geom_01_off + 581 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyzzz_xx, g_0_z_xyyyzzz_xy, g_0_z_xyyyzzz_xz, g_0_z_xyyyzzz_yy, g_0_z_xyyyzzz_yz, g_0_z_xyyyzzz_zz, g_0_z_yyyzzz_xx, g_0_z_yyyzzz_xxx, g_0_z_yyyzzz_xxy, g_0_z_yyyzzz_xxz, g_0_z_yyyzzz_xy, g_0_z_yyyzzz_xyy, g_0_z_yyyzzz_xyz, g_0_z_yyyzzz_xz, g_0_z_yyyzzz_xzz, g_0_z_yyyzzz_yy, g_0_z_yyyzzz_yz, g_0_z_yyyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyzzz_xx[k] = -g_0_z_yyyzzz_xx[k] * ab_x + g_0_z_yyyzzz_xxx[k];

                g_0_z_xyyyzzz_xy[k] = -g_0_z_yyyzzz_xy[k] * ab_x + g_0_z_yyyzzz_xxy[k];

                g_0_z_xyyyzzz_xz[k] = -g_0_z_yyyzzz_xz[k] * ab_x + g_0_z_yyyzzz_xxz[k];

                g_0_z_xyyyzzz_yy[k] = -g_0_z_yyyzzz_yy[k] * ab_x + g_0_z_yyyzzz_xyy[k];

                g_0_z_xyyyzzz_yz[k] = -g_0_z_yyyzzz_yz[k] * ab_x + g_0_z_yyyzzz_xyz[k];

                g_0_z_xyyyzzz_zz[k] = -g_0_z_yyyzzz_zz[k] * ab_x + g_0_z_yyyzzz_xzz[k];
            }

            /// Set up 582-588 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyzzzz_xx = cbuffer.data(kd_geom_01_off + 582 * ccomps * dcomps);

            auto g_0_z_xyyzzzz_xy = cbuffer.data(kd_geom_01_off + 583 * ccomps * dcomps);

            auto g_0_z_xyyzzzz_xz = cbuffer.data(kd_geom_01_off + 584 * ccomps * dcomps);

            auto g_0_z_xyyzzzz_yy = cbuffer.data(kd_geom_01_off + 585 * ccomps * dcomps);

            auto g_0_z_xyyzzzz_yz = cbuffer.data(kd_geom_01_off + 586 * ccomps * dcomps);

            auto g_0_z_xyyzzzz_zz = cbuffer.data(kd_geom_01_off + 587 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyzzzz_xx, g_0_z_xyyzzzz_xy, g_0_z_xyyzzzz_xz, g_0_z_xyyzzzz_yy, g_0_z_xyyzzzz_yz, g_0_z_xyyzzzz_zz, g_0_z_yyzzzz_xx, g_0_z_yyzzzz_xxx, g_0_z_yyzzzz_xxy, g_0_z_yyzzzz_xxz, g_0_z_yyzzzz_xy, g_0_z_yyzzzz_xyy, g_0_z_yyzzzz_xyz, g_0_z_yyzzzz_xz, g_0_z_yyzzzz_xzz, g_0_z_yyzzzz_yy, g_0_z_yyzzzz_yz, g_0_z_yyzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyzzzz_xx[k] = -g_0_z_yyzzzz_xx[k] * ab_x + g_0_z_yyzzzz_xxx[k];

                g_0_z_xyyzzzz_xy[k] = -g_0_z_yyzzzz_xy[k] * ab_x + g_0_z_yyzzzz_xxy[k];

                g_0_z_xyyzzzz_xz[k] = -g_0_z_yyzzzz_xz[k] * ab_x + g_0_z_yyzzzz_xxz[k];

                g_0_z_xyyzzzz_yy[k] = -g_0_z_yyzzzz_yy[k] * ab_x + g_0_z_yyzzzz_xyy[k];

                g_0_z_xyyzzzz_yz[k] = -g_0_z_yyzzzz_yz[k] * ab_x + g_0_z_yyzzzz_xyz[k];

                g_0_z_xyyzzzz_zz[k] = -g_0_z_yyzzzz_zz[k] * ab_x + g_0_z_yyzzzz_xzz[k];
            }

            /// Set up 588-594 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyzzzzz_xx = cbuffer.data(kd_geom_01_off + 588 * ccomps * dcomps);

            auto g_0_z_xyzzzzz_xy = cbuffer.data(kd_geom_01_off + 589 * ccomps * dcomps);

            auto g_0_z_xyzzzzz_xz = cbuffer.data(kd_geom_01_off + 590 * ccomps * dcomps);

            auto g_0_z_xyzzzzz_yy = cbuffer.data(kd_geom_01_off + 591 * ccomps * dcomps);

            auto g_0_z_xyzzzzz_yz = cbuffer.data(kd_geom_01_off + 592 * ccomps * dcomps);

            auto g_0_z_xyzzzzz_zz = cbuffer.data(kd_geom_01_off + 593 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyzzzzz_xx, g_0_z_xyzzzzz_xy, g_0_z_xyzzzzz_xz, g_0_z_xyzzzzz_yy, g_0_z_xyzzzzz_yz, g_0_z_xyzzzzz_zz, g_0_z_yzzzzz_xx, g_0_z_yzzzzz_xxx, g_0_z_yzzzzz_xxy, g_0_z_yzzzzz_xxz, g_0_z_yzzzzz_xy, g_0_z_yzzzzz_xyy, g_0_z_yzzzzz_xyz, g_0_z_yzzzzz_xz, g_0_z_yzzzzz_xzz, g_0_z_yzzzzz_yy, g_0_z_yzzzzz_yz, g_0_z_yzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyzzzzz_xx[k] = -g_0_z_yzzzzz_xx[k] * ab_x + g_0_z_yzzzzz_xxx[k];

                g_0_z_xyzzzzz_xy[k] = -g_0_z_yzzzzz_xy[k] * ab_x + g_0_z_yzzzzz_xxy[k];

                g_0_z_xyzzzzz_xz[k] = -g_0_z_yzzzzz_xz[k] * ab_x + g_0_z_yzzzzz_xxz[k];

                g_0_z_xyzzzzz_yy[k] = -g_0_z_yzzzzz_yy[k] * ab_x + g_0_z_yzzzzz_xyy[k];

                g_0_z_xyzzzzz_yz[k] = -g_0_z_yzzzzz_yz[k] * ab_x + g_0_z_yzzzzz_xyz[k];

                g_0_z_xyzzzzz_zz[k] = -g_0_z_yzzzzz_zz[k] * ab_x + g_0_z_yzzzzz_xzz[k];
            }

            /// Set up 594-600 components of targeted buffer : cbuffer.data(

            auto g_0_z_xzzzzzz_xx = cbuffer.data(kd_geom_01_off + 594 * ccomps * dcomps);

            auto g_0_z_xzzzzzz_xy = cbuffer.data(kd_geom_01_off + 595 * ccomps * dcomps);

            auto g_0_z_xzzzzzz_xz = cbuffer.data(kd_geom_01_off + 596 * ccomps * dcomps);

            auto g_0_z_xzzzzzz_yy = cbuffer.data(kd_geom_01_off + 597 * ccomps * dcomps);

            auto g_0_z_xzzzzzz_yz = cbuffer.data(kd_geom_01_off + 598 * ccomps * dcomps);

            auto g_0_z_xzzzzzz_zz = cbuffer.data(kd_geom_01_off + 599 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xzzzzzz_xx, g_0_z_xzzzzzz_xy, g_0_z_xzzzzzz_xz, g_0_z_xzzzzzz_yy, g_0_z_xzzzzzz_yz, g_0_z_xzzzzzz_zz, g_0_z_zzzzzz_xx, g_0_z_zzzzzz_xxx, g_0_z_zzzzzz_xxy, g_0_z_zzzzzz_xxz, g_0_z_zzzzzz_xy, g_0_z_zzzzzz_xyy, g_0_z_zzzzzz_xyz, g_0_z_zzzzzz_xz, g_0_z_zzzzzz_xzz, g_0_z_zzzzzz_yy, g_0_z_zzzzzz_yz, g_0_z_zzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xzzzzzz_xx[k] = -g_0_z_zzzzzz_xx[k] * ab_x + g_0_z_zzzzzz_xxx[k];

                g_0_z_xzzzzzz_xy[k] = -g_0_z_zzzzzz_xy[k] * ab_x + g_0_z_zzzzzz_xxy[k];

                g_0_z_xzzzzzz_xz[k] = -g_0_z_zzzzzz_xz[k] * ab_x + g_0_z_zzzzzz_xxz[k];

                g_0_z_xzzzzzz_yy[k] = -g_0_z_zzzzzz_yy[k] * ab_x + g_0_z_zzzzzz_xyy[k];

                g_0_z_xzzzzzz_yz[k] = -g_0_z_zzzzzz_yz[k] * ab_x + g_0_z_zzzzzz_xyz[k];

                g_0_z_xzzzzzz_zz[k] = -g_0_z_zzzzzz_zz[k] * ab_x + g_0_z_zzzzzz_xzz[k];
            }

            /// Set up 600-606 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyyyy_xx = cbuffer.data(kd_geom_01_off + 600 * ccomps * dcomps);

            auto g_0_z_yyyyyyy_xy = cbuffer.data(kd_geom_01_off + 601 * ccomps * dcomps);

            auto g_0_z_yyyyyyy_xz = cbuffer.data(kd_geom_01_off + 602 * ccomps * dcomps);

            auto g_0_z_yyyyyyy_yy = cbuffer.data(kd_geom_01_off + 603 * ccomps * dcomps);

            auto g_0_z_yyyyyyy_yz = cbuffer.data(kd_geom_01_off + 604 * ccomps * dcomps);

            auto g_0_z_yyyyyyy_zz = cbuffer.data(kd_geom_01_off + 605 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyyy_xx, g_0_z_yyyyyy_xxy, g_0_z_yyyyyy_xy, g_0_z_yyyyyy_xyy, g_0_z_yyyyyy_xyz, g_0_z_yyyyyy_xz, g_0_z_yyyyyy_yy, g_0_z_yyyyyy_yyy, g_0_z_yyyyyy_yyz, g_0_z_yyyyyy_yz, g_0_z_yyyyyy_yzz, g_0_z_yyyyyy_zz, g_0_z_yyyyyyy_xx, g_0_z_yyyyyyy_xy, g_0_z_yyyyyyy_xz, g_0_z_yyyyyyy_yy, g_0_z_yyyyyyy_yz, g_0_z_yyyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyyyy_xx[k] = -g_0_z_yyyyyy_xx[k] * ab_y + g_0_z_yyyyyy_xxy[k];

                g_0_z_yyyyyyy_xy[k] = -g_0_z_yyyyyy_xy[k] * ab_y + g_0_z_yyyyyy_xyy[k];

                g_0_z_yyyyyyy_xz[k] = -g_0_z_yyyyyy_xz[k] * ab_y + g_0_z_yyyyyy_xyz[k];

                g_0_z_yyyyyyy_yy[k] = -g_0_z_yyyyyy_yy[k] * ab_y + g_0_z_yyyyyy_yyy[k];

                g_0_z_yyyyyyy_yz[k] = -g_0_z_yyyyyy_yz[k] * ab_y + g_0_z_yyyyyy_yyz[k];

                g_0_z_yyyyyyy_zz[k] = -g_0_z_yyyyyy_zz[k] * ab_y + g_0_z_yyyyyy_yzz[k];
            }

            /// Set up 606-612 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyyyz_xx = cbuffer.data(kd_geom_01_off + 606 * ccomps * dcomps);

            auto g_0_z_yyyyyyz_xy = cbuffer.data(kd_geom_01_off + 607 * ccomps * dcomps);

            auto g_0_z_yyyyyyz_xz = cbuffer.data(kd_geom_01_off + 608 * ccomps * dcomps);

            auto g_0_z_yyyyyyz_yy = cbuffer.data(kd_geom_01_off + 609 * ccomps * dcomps);

            auto g_0_z_yyyyyyz_yz = cbuffer.data(kd_geom_01_off + 610 * ccomps * dcomps);

            auto g_0_z_yyyyyyz_zz = cbuffer.data(kd_geom_01_off + 611 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyyyz_xx, g_0_z_yyyyyyz_xy, g_0_z_yyyyyyz_xz, g_0_z_yyyyyyz_yy, g_0_z_yyyyyyz_yz, g_0_z_yyyyyyz_zz, g_0_z_yyyyyz_xx, g_0_z_yyyyyz_xxy, g_0_z_yyyyyz_xy, g_0_z_yyyyyz_xyy, g_0_z_yyyyyz_xyz, g_0_z_yyyyyz_xz, g_0_z_yyyyyz_yy, g_0_z_yyyyyz_yyy, g_0_z_yyyyyz_yyz, g_0_z_yyyyyz_yz, g_0_z_yyyyyz_yzz, g_0_z_yyyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyyyz_xx[k] = -g_0_z_yyyyyz_xx[k] * ab_y + g_0_z_yyyyyz_xxy[k];

                g_0_z_yyyyyyz_xy[k] = -g_0_z_yyyyyz_xy[k] * ab_y + g_0_z_yyyyyz_xyy[k];

                g_0_z_yyyyyyz_xz[k] = -g_0_z_yyyyyz_xz[k] * ab_y + g_0_z_yyyyyz_xyz[k];

                g_0_z_yyyyyyz_yy[k] = -g_0_z_yyyyyz_yy[k] * ab_y + g_0_z_yyyyyz_yyy[k];

                g_0_z_yyyyyyz_yz[k] = -g_0_z_yyyyyz_yz[k] * ab_y + g_0_z_yyyyyz_yyz[k];

                g_0_z_yyyyyyz_zz[k] = -g_0_z_yyyyyz_zz[k] * ab_y + g_0_z_yyyyyz_yzz[k];
            }

            /// Set up 612-618 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyyzz_xx = cbuffer.data(kd_geom_01_off + 612 * ccomps * dcomps);

            auto g_0_z_yyyyyzz_xy = cbuffer.data(kd_geom_01_off + 613 * ccomps * dcomps);

            auto g_0_z_yyyyyzz_xz = cbuffer.data(kd_geom_01_off + 614 * ccomps * dcomps);

            auto g_0_z_yyyyyzz_yy = cbuffer.data(kd_geom_01_off + 615 * ccomps * dcomps);

            auto g_0_z_yyyyyzz_yz = cbuffer.data(kd_geom_01_off + 616 * ccomps * dcomps);

            auto g_0_z_yyyyyzz_zz = cbuffer.data(kd_geom_01_off + 617 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyyzz_xx, g_0_z_yyyyyzz_xy, g_0_z_yyyyyzz_xz, g_0_z_yyyyyzz_yy, g_0_z_yyyyyzz_yz, g_0_z_yyyyyzz_zz, g_0_z_yyyyzz_xx, g_0_z_yyyyzz_xxy, g_0_z_yyyyzz_xy, g_0_z_yyyyzz_xyy, g_0_z_yyyyzz_xyz, g_0_z_yyyyzz_xz, g_0_z_yyyyzz_yy, g_0_z_yyyyzz_yyy, g_0_z_yyyyzz_yyz, g_0_z_yyyyzz_yz, g_0_z_yyyyzz_yzz, g_0_z_yyyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyyzz_xx[k] = -g_0_z_yyyyzz_xx[k] * ab_y + g_0_z_yyyyzz_xxy[k];

                g_0_z_yyyyyzz_xy[k] = -g_0_z_yyyyzz_xy[k] * ab_y + g_0_z_yyyyzz_xyy[k];

                g_0_z_yyyyyzz_xz[k] = -g_0_z_yyyyzz_xz[k] * ab_y + g_0_z_yyyyzz_xyz[k];

                g_0_z_yyyyyzz_yy[k] = -g_0_z_yyyyzz_yy[k] * ab_y + g_0_z_yyyyzz_yyy[k];

                g_0_z_yyyyyzz_yz[k] = -g_0_z_yyyyzz_yz[k] * ab_y + g_0_z_yyyyzz_yyz[k];

                g_0_z_yyyyyzz_zz[k] = -g_0_z_yyyyzz_zz[k] * ab_y + g_0_z_yyyyzz_yzz[k];
            }

            /// Set up 618-624 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyzzz_xx = cbuffer.data(kd_geom_01_off + 618 * ccomps * dcomps);

            auto g_0_z_yyyyzzz_xy = cbuffer.data(kd_geom_01_off + 619 * ccomps * dcomps);

            auto g_0_z_yyyyzzz_xz = cbuffer.data(kd_geom_01_off + 620 * ccomps * dcomps);

            auto g_0_z_yyyyzzz_yy = cbuffer.data(kd_geom_01_off + 621 * ccomps * dcomps);

            auto g_0_z_yyyyzzz_yz = cbuffer.data(kd_geom_01_off + 622 * ccomps * dcomps);

            auto g_0_z_yyyyzzz_zz = cbuffer.data(kd_geom_01_off + 623 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyzzz_xx, g_0_z_yyyyzzz_xy, g_0_z_yyyyzzz_xz, g_0_z_yyyyzzz_yy, g_0_z_yyyyzzz_yz, g_0_z_yyyyzzz_zz, g_0_z_yyyzzz_xx, g_0_z_yyyzzz_xxy, g_0_z_yyyzzz_xy, g_0_z_yyyzzz_xyy, g_0_z_yyyzzz_xyz, g_0_z_yyyzzz_xz, g_0_z_yyyzzz_yy, g_0_z_yyyzzz_yyy, g_0_z_yyyzzz_yyz, g_0_z_yyyzzz_yz, g_0_z_yyyzzz_yzz, g_0_z_yyyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyzzz_xx[k] = -g_0_z_yyyzzz_xx[k] * ab_y + g_0_z_yyyzzz_xxy[k];

                g_0_z_yyyyzzz_xy[k] = -g_0_z_yyyzzz_xy[k] * ab_y + g_0_z_yyyzzz_xyy[k];

                g_0_z_yyyyzzz_xz[k] = -g_0_z_yyyzzz_xz[k] * ab_y + g_0_z_yyyzzz_xyz[k];

                g_0_z_yyyyzzz_yy[k] = -g_0_z_yyyzzz_yy[k] * ab_y + g_0_z_yyyzzz_yyy[k];

                g_0_z_yyyyzzz_yz[k] = -g_0_z_yyyzzz_yz[k] * ab_y + g_0_z_yyyzzz_yyz[k];

                g_0_z_yyyyzzz_zz[k] = -g_0_z_yyyzzz_zz[k] * ab_y + g_0_z_yyyzzz_yzz[k];
            }

            /// Set up 624-630 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyzzzz_xx = cbuffer.data(kd_geom_01_off + 624 * ccomps * dcomps);

            auto g_0_z_yyyzzzz_xy = cbuffer.data(kd_geom_01_off + 625 * ccomps * dcomps);

            auto g_0_z_yyyzzzz_xz = cbuffer.data(kd_geom_01_off + 626 * ccomps * dcomps);

            auto g_0_z_yyyzzzz_yy = cbuffer.data(kd_geom_01_off + 627 * ccomps * dcomps);

            auto g_0_z_yyyzzzz_yz = cbuffer.data(kd_geom_01_off + 628 * ccomps * dcomps);

            auto g_0_z_yyyzzzz_zz = cbuffer.data(kd_geom_01_off + 629 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyzzzz_xx, g_0_z_yyyzzzz_xy, g_0_z_yyyzzzz_xz, g_0_z_yyyzzzz_yy, g_0_z_yyyzzzz_yz, g_0_z_yyyzzzz_zz, g_0_z_yyzzzz_xx, g_0_z_yyzzzz_xxy, g_0_z_yyzzzz_xy, g_0_z_yyzzzz_xyy, g_0_z_yyzzzz_xyz, g_0_z_yyzzzz_xz, g_0_z_yyzzzz_yy, g_0_z_yyzzzz_yyy, g_0_z_yyzzzz_yyz, g_0_z_yyzzzz_yz, g_0_z_yyzzzz_yzz, g_0_z_yyzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyzzzz_xx[k] = -g_0_z_yyzzzz_xx[k] * ab_y + g_0_z_yyzzzz_xxy[k];

                g_0_z_yyyzzzz_xy[k] = -g_0_z_yyzzzz_xy[k] * ab_y + g_0_z_yyzzzz_xyy[k];

                g_0_z_yyyzzzz_xz[k] = -g_0_z_yyzzzz_xz[k] * ab_y + g_0_z_yyzzzz_xyz[k];

                g_0_z_yyyzzzz_yy[k] = -g_0_z_yyzzzz_yy[k] * ab_y + g_0_z_yyzzzz_yyy[k];

                g_0_z_yyyzzzz_yz[k] = -g_0_z_yyzzzz_yz[k] * ab_y + g_0_z_yyzzzz_yyz[k];

                g_0_z_yyyzzzz_zz[k] = -g_0_z_yyzzzz_zz[k] * ab_y + g_0_z_yyzzzz_yzz[k];
            }

            /// Set up 630-636 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyzzzzz_xx = cbuffer.data(kd_geom_01_off + 630 * ccomps * dcomps);

            auto g_0_z_yyzzzzz_xy = cbuffer.data(kd_geom_01_off + 631 * ccomps * dcomps);

            auto g_0_z_yyzzzzz_xz = cbuffer.data(kd_geom_01_off + 632 * ccomps * dcomps);

            auto g_0_z_yyzzzzz_yy = cbuffer.data(kd_geom_01_off + 633 * ccomps * dcomps);

            auto g_0_z_yyzzzzz_yz = cbuffer.data(kd_geom_01_off + 634 * ccomps * dcomps);

            auto g_0_z_yyzzzzz_zz = cbuffer.data(kd_geom_01_off + 635 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyzzzzz_xx, g_0_z_yyzzzzz_xy, g_0_z_yyzzzzz_xz, g_0_z_yyzzzzz_yy, g_0_z_yyzzzzz_yz, g_0_z_yyzzzzz_zz, g_0_z_yzzzzz_xx, g_0_z_yzzzzz_xxy, g_0_z_yzzzzz_xy, g_0_z_yzzzzz_xyy, g_0_z_yzzzzz_xyz, g_0_z_yzzzzz_xz, g_0_z_yzzzzz_yy, g_0_z_yzzzzz_yyy, g_0_z_yzzzzz_yyz, g_0_z_yzzzzz_yz, g_0_z_yzzzzz_yzz, g_0_z_yzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyzzzzz_xx[k] = -g_0_z_yzzzzz_xx[k] * ab_y + g_0_z_yzzzzz_xxy[k];

                g_0_z_yyzzzzz_xy[k] = -g_0_z_yzzzzz_xy[k] * ab_y + g_0_z_yzzzzz_xyy[k];

                g_0_z_yyzzzzz_xz[k] = -g_0_z_yzzzzz_xz[k] * ab_y + g_0_z_yzzzzz_xyz[k];

                g_0_z_yyzzzzz_yy[k] = -g_0_z_yzzzzz_yy[k] * ab_y + g_0_z_yzzzzz_yyy[k];

                g_0_z_yyzzzzz_yz[k] = -g_0_z_yzzzzz_yz[k] * ab_y + g_0_z_yzzzzz_yyz[k];

                g_0_z_yyzzzzz_zz[k] = -g_0_z_yzzzzz_zz[k] * ab_y + g_0_z_yzzzzz_yzz[k];
            }

            /// Set up 636-642 components of targeted buffer : cbuffer.data(

            auto g_0_z_yzzzzzz_xx = cbuffer.data(kd_geom_01_off + 636 * ccomps * dcomps);

            auto g_0_z_yzzzzzz_xy = cbuffer.data(kd_geom_01_off + 637 * ccomps * dcomps);

            auto g_0_z_yzzzzzz_xz = cbuffer.data(kd_geom_01_off + 638 * ccomps * dcomps);

            auto g_0_z_yzzzzzz_yy = cbuffer.data(kd_geom_01_off + 639 * ccomps * dcomps);

            auto g_0_z_yzzzzzz_yz = cbuffer.data(kd_geom_01_off + 640 * ccomps * dcomps);

            auto g_0_z_yzzzzzz_zz = cbuffer.data(kd_geom_01_off + 641 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yzzzzzz_xx, g_0_z_yzzzzzz_xy, g_0_z_yzzzzzz_xz, g_0_z_yzzzzzz_yy, g_0_z_yzzzzzz_yz, g_0_z_yzzzzzz_zz, g_0_z_zzzzzz_xx, g_0_z_zzzzzz_xxy, g_0_z_zzzzzz_xy, g_0_z_zzzzzz_xyy, g_0_z_zzzzzz_xyz, g_0_z_zzzzzz_xz, g_0_z_zzzzzz_yy, g_0_z_zzzzzz_yyy, g_0_z_zzzzzz_yyz, g_0_z_zzzzzz_yz, g_0_z_zzzzzz_yzz, g_0_z_zzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yzzzzzz_xx[k] = -g_0_z_zzzzzz_xx[k] * ab_y + g_0_z_zzzzzz_xxy[k];

                g_0_z_yzzzzzz_xy[k] = -g_0_z_zzzzzz_xy[k] * ab_y + g_0_z_zzzzzz_xyy[k];

                g_0_z_yzzzzzz_xz[k] = -g_0_z_zzzzzz_xz[k] * ab_y + g_0_z_zzzzzz_xyz[k];

                g_0_z_yzzzzzz_yy[k] = -g_0_z_zzzzzz_yy[k] * ab_y + g_0_z_zzzzzz_yyy[k];

                g_0_z_yzzzzzz_yz[k] = -g_0_z_zzzzzz_yz[k] * ab_y + g_0_z_zzzzzz_yyz[k];

                g_0_z_yzzzzzz_zz[k] = -g_0_z_zzzzzz_zz[k] * ab_y + g_0_z_zzzzzz_yzz[k];
            }

            /// Set up 642-648 components of targeted buffer : cbuffer.data(

            auto g_0_z_zzzzzzz_xx = cbuffer.data(kd_geom_01_off + 642 * ccomps * dcomps);

            auto g_0_z_zzzzzzz_xy = cbuffer.data(kd_geom_01_off + 643 * ccomps * dcomps);

            auto g_0_z_zzzzzzz_xz = cbuffer.data(kd_geom_01_off + 644 * ccomps * dcomps);

            auto g_0_z_zzzzzzz_yy = cbuffer.data(kd_geom_01_off + 645 * ccomps * dcomps);

            auto g_0_z_zzzzzzz_yz = cbuffer.data(kd_geom_01_off + 646 * ccomps * dcomps);

            auto g_0_z_zzzzzzz_zz = cbuffer.data(kd_geom_01_off + 647 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zzzzzz_xx, g_0_z_zzzzzz_xxz, g_0_z_zzzzzz_xy, g_0_z_zzzzzz_xyz, g_0_z_zzzzzz_xz, g_0_z_zzzzzz_xzz, g_0_z_zzzzzz_yy, g_0_z_zzzzzz_yyz, g_0_z_zzzzzz_yz, g_0_z_zzzzzz_yzz, g_0_z_zzzzzz_zz, g_0_z_zzzzzz_zzz, g_0_z_zzzzzzz_xx, g_0_z_zzzzzzz_xy, g_0_z_zzzzzzz_xz, g_0_z_zzzzzzz_yy, g_0_z_zzzzzzz_yz, g_0_z_zzzzzzz_zz, g_zzzzzz_xx, g_zzzzzz_xy, g_zzzzzz_xz, g_zzzzzz_yy, g_zzzzzz_yz, g_zzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_zzzzzzz_xx[k] = g_zzzzzz_xx[k] - g_0_z_zzzzzz_xx[k] * ab_z + g_0_z_zzzzzz_xxz[k];

                g_0_z_zzzzzzz_xy[k] = g_zzzzzz_xy[k] - g_0_z_zzzzzz_xy[k] * ab_z + g_0_z_zzzzzz_xyz[k];

                g_0_z_zzzzzzz_xz[k] = g_zzzzzz_xz[k] - g_0_z_zzzzzz_xz[k] * ab_z + g_0_z_zzzzzz_xzz[k];

                g_0_z_zzzzzzz_yy[k] = g_zzzzzz_yy[k] - g_0_z_zzzzzz_yy[k] * ab_z + g_0_z_zzzzzz_yyz[k];

                g_0_z_zzzzzzz_yz[k] = g_zzzzzz_yz[k] - g_0_z_zzzzzz_yz[k] * ab_z + g_0_z_zzzzzz_yzz[k];

                g_0_z_zzzzzzz_zz[k] = g_zzzzzz_zz[k] - g_0_z_zzzzzz_zz[k] * ab_z + g_0_z_zzzzzz_zzz[k];
            }
        }
    }
}

} // erirec namespace

