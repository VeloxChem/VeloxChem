#include "ElectronRepulsionContrRecDIXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_hrr_electron_repulsion_dixx(CSimdArray<double>& contr_buffer_dixx,
                                     const CSimdArray<double>& contr_buffer_pixx,
                                     const CSimdArray<double>& contr_buffer_pkxx,
                                     const double ab_x,
                                     const double ab_y,
                                     const double ab_z,
                                     const int c_angmom,
                                     const int d_angmom) -> void
{
    const auto ndims = contr_buffer_dixx.number_of_columns();

    const auto ccomps = tensor::number_of_spherical_components(c_angmom);

    const auto dcomps = tensor::number_of_spherical_components(d_angmom);

    for (int i = 0; i < ccomps; i++)
    {
        for (int j = 0; j < dcomps; j++)
        {
            /// Set up components of auxilary buffer : contr_buffer_pixx

            const auto pi_off = i * dcomps + j;

            auto g_x_xxxxxx = contr_buffer_pixx[pi_off + 0 * ccomps * dcomps];

            auto g_x_xxxxxy = contr_buffer_pixx[pi_off + 1 * ccomps * dcomps];

            auto g_x_xxxxxz = contr_buffer_pixx[pi_off + 2 * ccomps * dcomps];

            auto g_x_xxxxyy = contr_buffer_pixx[pi_off + 3 * ccomps * dcomps];

            auto g_x_xxxxyz = contr_buffer_pixx[pi_off + 4 * ccomps * dcomps];

            auto g_x_xxxxzz = contr_buffer_pixx[pi_off + 5 * ccomps * dcomps];

            auto g_x_xxxyyy = contr_buffer_pixx[pi_off + 6 * ccomps * dcomps];

            auto g_x_xxxyyz = contr_buffer_pixx[pi_off + 7 * ccomps * dcomps];

            auto g_x_xxxyzz = contr_buffer_pixx[pi_off + 8 * ccomps * dcomps];

            auto g_x_xxxzzz = contr_buffer_pixx[pi_off + 9 * ccomps * dcomps];

            auto g_x_xxyyyy = contr_buffer_pixx[pi_off + 10 * ccomps * dcomps];

            auto g_x_xxyyyz = contr_buffer_pixx[pi_off + 11 * ccomps * dcomps];

            auto g_x_xxyyzz = contr_buffer_pixx[pi_off + 12 * ccomps * dcomps];

            auto g_x_xxyzzz = contr_buffer_pixx[pi_off + 13 * ccomps * dcomps];

            auto g_x_xxzzzz = contr_buffer_pixx[pi_off + 14 * ccomps * dcomps];

            auto g_x_xyyyyy = contr_buffer_pixx[pi_off + 15 * ccomps * dcomps];

            auto g_x_xyyyyz = contr_buffer_pixx[pi_off + 16 * ccomps * dcomps];

            auto g_x_xyyyzz = contr_buffer_pixx[pi_off + 17 * ccomps * dcomps];

            auto g_x_xyyzzz = contr_buffer_pixx[pi_off + 18 * ccomps * dcomps];

            auto g_x_xyzzzz = contr_buffer_pixx[pi_off + 19 * ccomps * dcomps];

            auto g_x_xzzzzz = contr_buffer_pixx[pi_off + 20 * ccomps * dcomps];

            auto g_x_yyyyyy = contr_buffer_pixx[pi_off + 21 * ccomps * dcomps];

            auto g_x_yyyyyz = contr_buffer_pixx[pi_off + 22 * ccomps * dcomps];

            auto g_x_yyyyzz = contr_buffer_pixx[pi_off + 23 * ccomps * dcomps];

            auto g_x_yyyzzz = contr_buffer_pixx[pi_off + 24 * ccomps * dcomps];

            auto g_x_yyzzzz = contr_buffer_pixx[pi_off + 25 * ccomps * dcomps];

            auto g_x_yzzzzz = contr_buffer_pixx[pi_off + 26 * ccomps * dcomps];

            auto g_x_zzzzzz = contr_buffer_pixx[pi_off + 27 * ccomps * dcomps];

            auto g_y_xxxxxx = contr_buffer_pixx[pi_off + 28 * ccomps * dcomps];

            auto g_y_xxxxxy = contr_buffer_pixx[pi_off + 29 * ccomps * dcomps];

            auto g_y_xxxxxz = contr_buffer_pixx[pi_off + 30 * ccomps * dcomps];

            auto g_y_xxxxyy = contr_buffer_pixx[pi_off + 31 * ccomps * dcomps];

            auto g_y_xxxxyz = contr_buffer_pixx[pi_off + 32 * ccomps * dcomps];

            auto g_y_xxxxzz = contr_buffer_pixx[pi_off + 33 * ccomps * dcomps];

            auto g_y_xxxyyy = contr_buffer_pixx[pi_off + 34 * ccomps * dcomps];

            auto g_y_xxxyyz = contr_buffer_pixx[pi_off + 35 * ccomps * dcomps];

            auto g_y_xxxyzz = contr_buffer_pixx[pi_off + 36 * ccomps * dcomps];

            auto g_y_xxxzzz = contr_buffer_pixx[pi_off + 37 * ccomps * dcomps];

            auto g_y_xxyyyy = contr_buffer_pixx[pi_off + 38 * ccomps * dcomps];

            auto g_y_xxyyyz = contr_buffer_pixx[pi_off + 39 * ccomps * dcomps];

            auto g_y_xxyyzz = contr_buffer_pixx[pi_off + 40 * ccomps * dcomps];

            auto g_y_xxyzzz = contr_buffer_pixx[pi_off + 41 * ccomps * dcomps];

            auto g_y_xxzzzz = contr_buffer_pixx[pi_off + 42 * ccomps * dcomps];

            auto g_y_xyyyyy = contr_buffer_pixx[pi_off + 43 * ccomps * dcomps];

            auto g_y_xyyyyz = contr_buffer_pixx[pi_off + 44 * ccomps * dcomps];

            auto g_y_xyyyzz = contr_buffer_pixx[pi_off + 45 * ccomps * dcomps];

            auto g_y_xyyzzz = contr_buffer_pixx[pi_off + 46 * ccomps * dcomps];

            auto g_y_xyzzzz = contr_buffer_pixx[pi_off + 47 * ccomps * dcomps];

            auto g_y_xzzzzz = contr_buffer_pixx[pi_off + 48 * ccomps * dcomps];

            auto g_y_yyyyyy = contr_buffer_pixx[pi_off + 49 * ccomps * dcomps];

            auto g_y_yyyyyz = contr_buffer_pixx[pi_off + 50 * ccomps * dcomps];

            auto g_y_yyyyzz = contr_buffer_pixx[pi_off + 51 * ccomps * dcomps];

            auto g_y_yyyzzz = contr_buffer_pixx[pi_off + 52 * ccomps * dcomps];

            auto g_y_yyzzzz = contr_buffer_pixx[pi_off + 53 * ccomps * dcomps];

            auto g_y_yzzzzz = contr_buffer_pixx[pi_off + 54 * ccomps * dcomps];

            auto g_y_zzzzzz = contr_buffer_pixx[pi_off + 55 * ccomps * dcomps];

            auto g_z_xxxxxx = contr_buffer_pixx[pi_off + 56 * ccomps * dcomps];

            auto g_z_xxxxxy = contr_buffer_pixx[pi_off + 57 * ccomps * dcomps];

            auto g_z_xxxxxz = contr_buffer_pixx[pi_off + 58 * ccomps * dcomps];

            auto g_z_xxxxyy = contr_buffer_pixx[pi_off + 59 * ccomps * dcomps];

            auto g_z_xxxxyz = contr_buffer_pixx[pi_off + 60 * ccomps * dcomps];

            auto g_z_xxxxzz = contr_buffer_pixx[pi_off + 61 * ccomps * dcomps];

            auto g_z_xxxyyy = contr_buffer_pixx[pi_off + 62 * ccomps * dcomps];

            auto g_z_xxxyyz = contr_buffer_pixx[pi_off + 63 * ccomps * dcomps];

            auto g_z_xxxyzz = contr_buffer_pixx[pi_off + 64 * ccomps * dcomps];

            auto g_z_xxxzzz = contr_buffer_pixx[pi_off + 65 * ccomps * dcomps];

            auto g_z_xxyyyy = contr_buffer_pixx[pi_off + 66 * ccomps * dcomps];

            auto g_z_xxyyyz = contr_buffer_pixx[pi_off + 67 * ccomps * dcomps];

            auto g_z_xxyyzz = contr_buffer_pixx[pi_off + 68 * ccomps * dcomps];

            auto g_z_xxyzzz = contr_buffer_pixx[pi_off + 69 * ccomps * dcomps];

            auto g_z_xxzzzz = contr_buffer_pixx[pi_off + 70 * ccomps * dcomps];

            auto g_z_xyyyyy = contr_buffer_pixx[pi_off + 71 * ccomps * dcomps];

            auto g_z_xyyyyz = contr_buffer_pixx[pi_off + 72 * ccomps * dcomps];

            auto g_z_xyyyzz = contr_buffer_pixx[pi_off + 73 * ccomps * dcomps];

            auto g_z_xyyzzz = contr_buffer_pixx[pi_off + 74 * ccomps * dcomps];

            auto g_z_xyzzzz = contr_buffer_pixx[pi_off + 75 * ccomps * dcomps];

            auto g_z_xzzzzz = contr_buffer_pixx[pi_off + 76 * ccomps * dcomps];

            auto g_z_yyyyyy = contr_buffer_pixx[pi_off + 77 * ccomps * dcomps];

            auto g_z_yyyyyz = contr_buffer_pixx[pi_off + 78 * ccomps * dcomps];

            auto g_z_yyyyzz = contr_buffer_pixx[pi_off + 79 * ccomps * dcomps];

            auto g_z_yyyzzz = contr_buffer_pixx[pi_off + 80 * ccomps * dcomps];

            auto g_z_yyzzzz = contr_buffer_pixx[pi_off + 81 * ccomps * dcomps];

            auto g_z_yzzzzz = contr_buffer_pixx[pi_off + 82 * ccomps * dcomps];

            auto g_z_zzzzzz = contr_buffer_pixx[pi_off + 83 * ccomps * dcomps];

            /// Set up components of auxilary buffer : contr_buffer_pkxx

            const auto pk_off = i * dcomps + j;

            auto g_x_xxxxxxx = contr_buffer_pkxx[pk_off + 0 * ccomps * dcomps];

            auto g_x_xxxxxxy = contr_buffer_pkxx[pk_off + 1 * ccomps * dcomps];

            auto g_x_xxxxxxz = contr_buffer_pkxx[pk_off + 2 * ccomps * dcomps];

            auto g_x_xxxxxyy = contr_buffer_pkxx[pk_off + 3 * ccomps * dcomps];

            auto g_x_xxxxxyz = contr_buffer_pkxx[pk_off + 4 * ccomps * dcomps];

            auto g_x_xxxxxzz = contr_buffer_pkxx[pk_off + 5 * ccomps * dcomps];

            auto g_x_xxxxyyy = contr_buffer_pkxx[pk_off + 6 * ccomps * dcomps];

            auto g_x_xxxxyyz = contr_buffer_pkxx[pk_off + 7 * ccomps * dcomps];

            auto g_x_xxxxyzz = contr_buffer_pkxx[pk_off + 8 * ccomps * dcomps];

            auto g_x_xxxxzzz = contr_buffer_pkxx[pk_off + 9 * ccomps * dcomps];

            auto g_x_xxxyyyy = contr_buffer_pkxx[pk_off + 10 * ccomps * dcomps];

            auto g_x_xxxyyyz = contr_buffer_pkxx[pk_off + 11 * ccomps * dcomps];

            auto g_x_xxxyyzz = contr_buffer_pkxx[pk_off + 12 * ccomps * dcomps];

            auto g_x_xxxyzzz = contr_buffer_pkxx[pk_off + 13 * ccomps * dcomps];

            auto g_x_xxxzzzz = contr_buffer_pkxx[pk_off + 14 * ccomps * dcomps];

            auto g_x_xxyyyyy = contr_buffer_pkxx[pk_off + 15 * ccomps * dcomps];

            auto g_x_xxyyyyz = contr_buffer_pkxx[pk_off + 16 * ccomps * dcomps];

            auto g_x_xxyyyzz = contr_buffer_pkxx[pk_off + 17 * ccomps * dcomps];

            auto g_x_xxyyzzz = contr_buffer_pkxx[pk_off + 18 * ccomps * dcomps];

            auto g_x_xxyzzzz = contr_buffer_pkxx[pk_off + 19 * ccomps * dcomps];

            auto g_x_xxzzzzz = contr_buffer_pkxx[pk_off + 20 * ccomps * dcomps];

            auto g_x_xyyyyyy = contr_buffer_pkxx[pk_off + 21 * ccomps * dcomps];

            auto g_x_xyyyyyz = contr_buffer_pkxx[pk_off + 22 * ccomps * dcomps];

            auto g_x_xyyyyzz = contr_buffer_pkxx[pk_off + 23 * ccomps * dcomps];

            auto g_x_xyyyzzz = contr_buffer_pkxx[pk_off + 24 * ccomps * dcomps];

            auto g_x_xyyzzzz = contr_buffer_pkxx[pk_off + 25 * ccomps * dcomps];

            auto g_x_xyzzzzz = contr_buffer_pkxx[pk_off + 26 * ccomps * dcomps];

            auto g_x_xzzzzzz = contr_buffer_pkxx[pk_off + 27 * ccomps * dcomps];

            auto g_y_xxxxxxx = contr_buffer_pkxx[pk_off + 36 * ccomps * dcomps];

            auto g_y_xxxxxxy = contr_buffer_pkxx[pk_off + 37 * ccomps * dcomps];

            auto g_y_xxxxxxz = contr_buffer_pkxx[pk_off + 38 * ccomps * dcomps];

            auto g_y_xxxxxyy = contr_buffer_pkxx[pk_off + 39 * ccomps * dcomps];

            auto g_y_xxxxxyz = contr_buffer_pkxx[pk_off + 40 * ccomps * dcomps];

            auto g_y_xxxxxzz = contr_buffer_pkxx[pk_off + 41 * ccomps * dcomps];

            auto g_y_xxxxyyy = contr_buffer_pkxx[pk_off + 42 * ccomps * dcomps];

            auto g_y_xxxxyyz = contr_buffer_pkxx[pk_off + 43 * ccomps * dcomps];

            auto g_y_xxxxyzz = contr_buffer_pkxx[pk_off + 44 * ccomps * dcomps];

            auto g_y_xxxxzzz = contr_buffer_pkxx[pk_off + 45 * ccomps * dcomps];

            auto g_y_xxxyyyy = contr_buffer_pkxx[pk_off + 46 * ccomps * dcomps];

            auto g_y_xxxyyyz = contr_buffer_pkxx[pk_off + 47 * ccomps * dcomps];

            auto g_y_xxxyyzz = contr_buffer_pkxx[pk_off + 48 * ccomps * dcomps];

            auto g_y_xxxyzzz = contr_buffer_pkxx[pk_off + 49 * ccomps * dcomps];

            auto g_y_xxxzzzz = contr_buffer_pkxx[pk_off + 50 * ccomps * dcomps];

            auto g_y_xxyyyyy = contr_buffer_pkxx[pk_off + 51 * ccomps * dcomps];

            auto g_y_xxyyyyz = contr_buffer_pkxx[pk_off + 52 * ccomps * dcomps];

            auto g_y_xxyyyzz = contr_buffer_pkxx[pk_off + 53 * ccomps * dcomps];

            auto g_y_xxyyzzz = contr_buffer_pkxx[pk_off + 54 * ccomps * dcomps];

            auto g_y_xxyzzzz = contr_buffer_pkxx[pk_off + 55 * ccomps * dcomps];

            auto g_y_xxzzzzz = contr_buffer_pkxx[pk_off + 56 * ccomps * dcomps];

            auto g_y_xyyyyyy = contr_buffer_pkxx[pk_off + 57 * ccomps * dcomps];

            auto g_y_xyyyyyz = contr_buffer_pkxx[pk_off + 58 * ccomps * dcomps];

            auto g_y_xyyyyzz = contr_buffer_pkxx[pk_off + 59 * ccomps * dcomps];

            auto g_y_xyyyzzz = contr_buffer_pkxx[pk_off + 60 * ccomps * dcomps];

            auto g_y_xyyzzzz = contr_buffer_pkxx[pk_off + 61 * ccomps * dcomps];

            auto g_y_xyzzzzz = contr_buffer_pkxx[pk_off + 62 * ccomps * dcomps];

            auto g_y_xzzzzzz = contr_buffer_pkxx[pk_off + 63 * ccomps * dcomps];

            auto g_y_yyyyyyy = contr_buffer_pkxx[pk_off + 64 * ccomps * dcomps];

            auto g_y_yyyyyyz = contr_buffer_pkxx[pk_off + 65 * ccomps * dcomps];

            auto g_y_yyyyyzz = contr_buffer_pkxx[pk_off + 66 * ccomps * dcomps];

            auto g_y_yyyyzzz = contr_buffer_pkxx[pk_off + 67 * ccomps * dcomps];

            auto g_y_yyyzzzz = contr_buffer_pkxx[pk_off + 68 * ccomps * dcomps];

            auto g_y_yyzzzzz = contr_buffer_pkxx[pk_off + 69 * ccomps * dcomps];

            auto g_y_yzzzzzz = contr_buffer_pkxx[pk_off + 70 * ccomps * dcomps];

            auto g_z_xxxxxxx = contr_buffer_pkxx[pk_off + 72 * ccomps * dcomps];

            auto g_z_xxxxxxy = contr_buffer_pkxx[pk_off + 73 * ccomps * dcomps];

            auto g_z_xxxxxxz = contr_buffer_pkxx[pk_off + 74 * ccomps * dcomps];

            auto g_z_xxxxxyy = contr_buffer_pkxx[pk_off + 75 * ccomps * dcomps];

            auto g_z_xxxxxyz = contr_buffer_pkxx[pk_off + 76 * ccomps * dcomps];

            auto g_z_xxxxxzz = contr_buffer_pkxx[pk_off + 77 * ccomps * dcomps];

            auto g_z_xxxxyyy = contr_buffer_pkxx[pk_off + 78 * ccomps * dcomps];

            auto g_z_xxxxyyz = contr_buffer_pkxx[pk_off + 79 * ccomps * dcomps];

            auto g_z_xxxxyzz = contr_buffer_pkxx[pk_off + 80 * ccomps * dcomps];

            auto g_z_xxxxzzz = contr_buffer_pkxx[pk_off + 81 * ccomps * dcomps];

            auto g_z_xxxyyyy = contr_buffer_pkxx[pk_off + 82 * ccomps * dcomps];

            auto g_z_xxxyyyz = contr_buffer_pkxx[pk_off + 83 * ccomps * dcomps];

            auto g_z_xxxyyzz = contr_buffer_pkxx[pk_off + 84 * ccomps * dcomps];

            auto g_z_xxxyzzz = contr_buffer_pkxx[pk_off + 85 * ccomps * dcomps];

            auto g_z_xxxzzzz = contr_buffer_pkxx[pk_off + 86 * ccomps * dcomps];

            auto g_z_xxyyyyy = contr_buffer_pkxx[pk_off + 87 * ccomps * dcomps];

            auto g_z_xxyyyyz = contr_buffer_pkxx[pk_off + 88 * ccomps * dcomps];

            auto g_z_xxyyyzz = contr_buffer_pkxx[pk_off + 89 * ccomps * dcomps];

            auto g_z_xxyyzzz = contr_buffer_pkxx[pk_off + 90 * ccomps * dcomps];

            auto g_z_xxyzzzz = contr_buffer_pkxx[pk_off + 91 * ccomps * dcomps];

            auto g_z_xxzzzzz = contr_buffer_pkxx[pk_off + 92 * ccomps * dcomps];

            auto g_z_xyyyyyy = contr_buffer_pkxx[pk_off + 93 * ccomps * dcomps];

            auto g_z_xyyyyyz = contr_buffer_pkxx[pk_off + 94 * ccomps * dcomps];

            auto g_z_xyyyyzz = contr_buffer_pkxx[pk_off + 95 * ccomps * dcomps];

            auto g_z_xyyyzzz = contr_buffer_pkxx[pk_off + 96 * ccomps * dcomps];

            auto g_z_xyyzzzz = contr_buffer_pkxx[pk_off + 97 * ccomps * dcomps];

            auto g_z_xyzzzzz = contr_buffer_pkxx[pk_off + 98 * ccomps * dcomps];

            auto g_z_xzzzzzz = contr_buffer_pkxx[pk_off + 99 * ccomps * dcomps];

            auto g_z_yyyyyyy = contr_buffer_pkxx[pk_off + 100 * ccomps * dcomps];

            auto g_z_yyyyyyz = contr_buffer_pkxx[pk_off + 101 * ccomps * dcomps];

            auto g_z_yyyyyzz = contr_buffer_pkxx[pk_off + 102 * ccomps * dcomps];

            auto g_z_yyyyzzz = contr_buffer_pkxx[pk_off + 103 * ccomps * dcomps];

            auto g_z_yyyzzzz = contr_buffer_pkxx[pk_off + 104 * ccomps * dcomps];

            auto g_z_yyzzzzz = contr_buffer_pkxx[pk_off + 105 * ccomps * dcomps];

            auto g_z_yzzzzzz = contr_buffer_pkxx[pk_off + 106 * ccomps * dcomps];

            auto g_z_zzzzzzz = contr_buffer_pkxx[pk_off + 107 * ccomps * dcomps];

            /// set up bra offset for contr_buffer_dixx

            const auto di_off = i * dcomps + j;

            /// Set up 0-28 components of targeted buffer : contr_buffer_dixx

            auto g_xx_xxxxxx = contr_buffer_dixx[di_off + 0 * ccomps * dcomps];

            auto g_xx_xxxxxy = contr_buffer_dixx[di_off + 1 * ccomps * dcomps];

            auto g_xx_xxxxxz = contr_buffer_dixx[di_off + 2 * ccomps * dcomps];

            auto g_xx_xxxxyy = contr_buffer_dixx[di_off + 3 * ccomps * dcomps];

            auto g_xx_xxxxyz = contr_buffer_dixx[di_off + 4 * ccomps * dcomps];

            auto g_xx_xxxxzz = contr_buffer_dixx[di_off + 5 * ccomps * dcomps];

            auto g_xx_xxxyyy = contr_buffer_dixx[di_off + 6 * ccomps * dcomps];

            auto g_xx_xxxyyz = contr_buffer_dixx[di_off + 7 * ccomps * dcomps];

            auto g_xx_xxxyzz = contr_buffer_dixx[di_off + 8 * ccomps * dcomps];

            auto g_xx_xxxzzz = contr_buffer_dixx[di_off + 9 * ccomps * dcomps];

            auto g_xx_xxyyyy = contr_buffer_dixx[di_off + 10 * ccomps * dcomps];

            auto g_xx_xxyyyz = contr_buffer_dixx[di_off + 11 * ccomps * dcomps];

            auto g_xx_xxyyzz = contr_buffer_dixx[di_off + 12 * ccomps * dcomps];

            auto g_xx_xxyzzz = contr_buffer_dixx[di_off + 13 * ccomps * dcomps];

            auto g_xx_xxzzzz = contr_buffer_dixx[di_off + 14 * ccomps * dcomps];

            auto g_xx_xyyyyy = contr_buffer_dixx[di_off + 15 * ccomps * dcomps];

            auto g_xx_xyyyyz = contr_buffer_dixx[di_off + 16 * ccomps * dcomps];

            auto g_xx_xyyyzz = contr_buffer_dixx[di_off + 17 * ccomps * dcomps];

            auto g_xx_xyyzzz = contr_buffer_dixx[di_off + 18 * ccomps * dcomps];

            auto g_xx_xyzzzz = contr_buffer_dixx[di_off + 19 * ccomps * dcomps];

            auto g_xx_xzzzzz = contr_buffer_dixx[di_off + 20 * ccomps * dcomps];

            auto g_xx_yyyyyy = contr_buffer_dixx[di_off + 21 * ccomps * dcomps];

            auto g_xx_yyyyyz = contr_buffer_dixx[di_off + 22 * ccomps * dcomps];

            auto g_xx_yyyyzz = contr_buffer_dixx[di_off + 23 * ccomps * dcomps];

            auto g_xx_yyyzzz = contr_buffer_dixx[di_off + 24 * ccomps * dcomps];

            auto g_xx_yyzzzz = contr_buffer_dixx[di_off + 25 * ccomps * dcomps];

            auto g_xx_yzzzzz = contr_buffer_dixx[di_off + 26 * ccomps * dcomps];

            auto g_xx_zzzzzz = contr_buffer_dixx[di_off + 27 * ccomps * dcomps];

            #pragma omp simd aligned(g_x_xxxxxx, g_x_xxxxxxx, g_x_xxxxxxy, g_x_xxxxxxz, g_x_xxxxxy, g_x_xxxxxyy, g_x_xxxxxyz, g_x_xxxxxz, g_x_xxxxxzz, g_x_xxxxyy, g_x_xxxxyyy, g_x_xxxxyyz, g_x_xxxxyz, g_x_xxxxyzz, g_x_xxxxzz, g_x_xxxxzzz, g_x_xxxyyy, g_x_xxxyyyy, g_x_xxxyyyz, g_x_xxxyyz, g_x_xxxyyzz, g_x_xxxyzz, g_x_xxxyzzz, g_x_xxxzzz, g_x_xxxzzzz, g_x_xxyyyy, g_x_xxyyyyy, g_x_xxyyyyz, g_x_xxyyyz, g_x_xxyyyzz, g_x_xxyyzz, g_x_xxyyzzz, g_x_xxyzzz, g_x_xxyzzzz, g_x_xxzzzz, g_x_xxzzzzz, g_x_xyyyyy, g_x_xyyyyyy, g_x_xyyyyyz, g_x_xyyyyz, g_x_xyyyyzz, g_x_xyyyzz, g_x_xyyyzzz, g_x_xyyzzz, g_x_xyyzzzz, g_x_xyzzzz, g_x_xyzzzzz, g_x_xzzzzz, g_x_xzzzzzz, g_x_yyyyyy, g_x_yyyyyz, g_x_yyyyzz, g_x_yyyzzz, g_x_yyzzzz, g_x_yzzzzz, g_x_zzzzzz, g_xx_xxxxxx, g_xx_xxxxxy, g_xx_xxxxxz, g_xx_xxxxyy, g_xx_xxxxyz, g_xx_xxxxzz, g_xx_xxxyyy, g_xx_xxxyyz, g_xx_xxxyzz, g_xx_xxxzzz, g_xx_xxyyyy, g_xx_xxyyyz, g_xx_xxyyzz, g_xx_xxyzzz, g_xx_xxzzzz, g_xx_xyyyyy, g_xx_xyyyyz, g_xx_xyyyzz, g_xx_xyyzzz, g_xx_xyzzzz, g_xx_xzzzzz, g_xx_yyyyyy, g_xx_yyyyyz, g_xx_yyyyzz, g_xx_yyyzzz, g_xx_yyzzzz, g_xx_yzzzzz, g_xx_zzzzzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_xx_xxxxxx[k] = -g_x_xxxxxx[k] * ab_x + g_x_xxxxxxx[k];

                g_xx_xxxxxy[k] = -g_x_xxxxxy[k] * ab_x + g_x_xxxxxxy[k];

                g_xx_xxxxxz[k] = -g_x_xxxxxz[k] * ab_x + g_x_xxxxxxz[k];

                g_xx_xxxxyy[k] = -g_x_xxxxyy[k] * ab_x + g_x_xxxxxyy[k];

                g_xx_xxxxyz[k] = -g_x_xxxxyz[k] * ab_x + g_x_xxxxxyz[k];

                g_xx_xxxxzz[k] = -g_x_xxxxzz[k] * ab_x + g_x_xxxxxzz[k];

                g_xx_xxxyyy[k] = -g_x_xxxyyy[k] * ab_x + g_x_xxxxyyy[k];

                g_xx_xxxyyz[k] = -g_x_xxxyyz[k] * ab_x + g_x_xxxxyyz[k];

                g_xx_xxxyzz[k] = -g_x_xxxyzz[k] * ab_x + g_x_xxxxyzz[k];

                g_xx_xxxzzz[k] = -g_x_xxxzzz[k] * ab_x + g_x_xxxxzzz[k];

                g_xx_xxyyyy[k] = -g_x_xxyyyy[k] * ab_x + g_x_xxxyyyy[k];

                g_xx_xxyyyz[k] = -g_x_xxyyyz[k] * ab_x + g_x_xxxyyyz[k];

                g_xx_xxyyzz[k] = -g_x_xxyyzz[k] * ab_x + g_x_xxxyyzz[k];

                g_xx_xxyzzz[k] = -g_x_xxyzzz[k] * ab_x + g_x_xxxyzzz[k];

                g_xx_xxzzzz[k] = -g_x_xxzzzz[k] * ab_x + g_x_xxxzzzz[k];

                g_xx_xyyyyy[k] = -g_x_xyyyyy[k] * ab_x + g_x_xxyyyyy[k];

                g_xx_xyyyyz[k] = -g_x_xyyyyz[k] * ab_x + g_x_xxyyyyz[k];

                g_xx_xyyyzz[k] = -g_x_xyyyzz[k] * ab_x + g_x_xxyyyzz[k];

                g_xx_xyyzzz[k] = -g_x_xyyzzz[k] * ab_x + g_x_xxyyzzz[k];

                g_xx_xyzzzz[k] = -g_x_xyzzzz[k] * ab_x + g_x_xxyzzzz[k];

                g_xx_xzzzzz[k] = -g_x_xzzzzz[k] * ab_x + g_x_xxzzzzz[k];

                g_xx_yyyyyy[k] = -g_x_yyyyyy[k] * ab_x + g_x_xyyyyyy[k];

                g_xx_yyyyyz[k] = -g_x_yyyyyz[k] * ab_x + g_x_xyyyyyz[k];

                g_xx_yyyyzz[k] = -g_x_yyyyzz[k] * ab_x + g_x_xyyyyzz[k];

                g_xx_yyyzzz[k] = -g_x_yyyzzz[k] * ab_x + g_x_xyyyzzz[k];

                g_xx_yyzzzz[k] = -g_x_yyzzzz[k] * ab_x + g_x_xyyzzzz[k];

                g_xx_yzzzzz[k] = -g_x_yzzzzz[k] * ab_x + g_x_xyzzzzz[k];

                g_xx_zzzzzz[k] = -g_x_zzzzzz[k] * ab_x + g_x_xzzzzzz[k];
            }

            /// Set up 28-56 components of targeted buffer : contr_buffer_dixx

            auto g_xy_xxxxxx = contr_buffer_dixx[di_off + 28 * ccomps * dcomps];

            auto g_xy_xxxxxy = contr_buffer_dixx[di_off + 29 * ccomps * dcomps];

            auto g_xy_xxxxxz = contr_buffer_dixx[di_off + 30 * ccomps * dcomps];

            auto g_xy_xxxxyy = contr_buffer_dixx[di_off + 31 * ccomps * dcomps];

            auto g_xy_xxxxyz = contr_buffer_dixx[di_off + 32 * ccomps * dcomps];

            auto g_xy_xxxxzz = contr_buffer_dixx[di_off + 33 * ccomps * dcomps];

            auto g_xy_xxxyyy = contr_buffer_dixx[di_off + 34 * ccomps * dcomps];

            auto g_xy_xxxyyz = contr_buffer_dixx[di_off + 35 * ccomps * dcomps];

            auto g_xy_xxxyzz = contr_buffer_dixx[di_off + 36 * ccomps * dcomps];

            auto g_xy_xxxzzz = contr_buffer_dixx[di_off + 37 * ccomps * dcomps];

            auto g_xy_xxyyyy = contr_buffer_dixx[di_off + 38 * ccomps * dcomps];

            auto g_xy_xxyyyz = contr_buffer_dixx[di_off + 39 * ccomps * dcomps];

            auto g_xy_xxyyzz = contr_buffer_dixx[di_off + 40 * ccomps * dcomps];

            auto g_xy_xxyzzz = contr_buffer_dixx[di_off + 41 * ccomps * dcomps];

            auto g_xy_xxzzzz = contr_buffer_dixx[di_off + 42 * ccomps * dcomps];

            auto g_xy_xyyyyy = contr_buffer_dixx[di_off + 43 * ccomps * dcomps];

            auto g_xy_xyyyyz = contr_buffer_dixx[di_off + 44 * ccomps * dcomps];

            auto g_xy_xyyyzz = contr_buffer_dixx[di_off + 45 * ccomps * dcomps];

            auto g_xy_xyyzzz = contr_buffer_dixx[di_off + 46 * ccomps * dcomps];

            auto g_xy_xyzzzz = contr_buffer_dixx[di_off + 47 * ccomps * dcomps];

            auto g_xy_xzzzzz = contr_buffer_dixx[di_off + 48 * ccomps * dcomps];

            auto g_xy_yyyyyy = contr_buffer_dixx[di_off + 49 * ccomps * dcomps];

            auto g_xy_yyyyyz = contr_buffer_dixx[di_off + 50 * ccomps * dcomps];

            auto g_xy_yyyyzz = contr_buffer_dixx[di_off + 51 * ccomps * dcomps];

            auto g_xy_yyyzzz = contr_buffer_dixx[di_off + 52 * ccomps * dcomps];

            auto g_xy_yyzzzz = contr_buffer_dixx[di_off + 53 * ccomps * dcomps];

            auto g_xy_yzzzzz = contr_buffer_dixx[di_off + 54 * ccomps * dcomps];

            auto g_xy_zzzzzz = contr_buffer_dixx[di_off + 55 * ccomps * dcomps];

            #pragma omp simd aligned(g_xy_xxxxxx, g_xy_xxxxxy, g_xy_xxxxxz, g_xy_xxxxyy, g_xy_xxxxyz, g_xy_xxxxzz, g_xy_xxxyyy, g_xy_xxxyyz, g_xy_xxxyzz, g_xy_xxxzzz, g_xy_xxyyyy, g_xy_xxyyyz, g_xy_xxyyzz, g_xy_xxyzzz, g_xy_xxzzzz, g_xy_xyyyyy, g_xy_xyyyyz, g_xy_xyyyzz, g_xy_xyyzzz, g_xy_xyzzzz, g_xy_xzzzzz, g_xy_yyyyyy, g_xy_yyyyyz, g_xy_yyyyzz, g_xy_yyyzzz, g_xy_yyzzzz, g_xy_yzzzzz, g_xy_zzzzzz, g_y_xxxxxx, g_y_xxxxxxx, g_y_xxxxxxy, g_y_xxxxxxz, g_y_xxxxxy, g_y_xxxxxyy, g_y_xxxxxyz, g_y_xxxxxz, g_y_xxxxxzz, g_y_xxxxyy, g_y_xxxxyyy, g_y_xxxxyyz, g_y_xxxxyz, g_y_xxxxyzz, g_y_xxxxzz, g_y_xxxxzzz, g_y_xxxyyy, g_y_xxxyyyy, g_y_xxxyyyz, g_y_xxxyyz, g_y_xxxyyzz, g_y_xxxyzz, g_y_xxxyzzz, g_y_xxxzzz, g_y_xxxzzzz, g_y_xxyyyy, g_y_xxyyyyy, g_y_xxyyyyz, g_y_xxyyyz, g_y_xxyyyzz, g_y_xxyyzz, g_y_xxyyzzz, g_y_xxyzzz, g_y_xxyzzzz, g_y_xxzzzz, g_y_xxzzzzz, g_y_xyyyyy, g_y_xyyyyyy, g_y_xyyyyyz, g_y_xyyyyz, g_y_xyyyyzz, g_y_xyyyzz, g_y_xyyyzzz, g_y_xyyzzz, g_y_xyyzzzz, g_y_xyzzzz, g_y_xyzzzzz, g_y_xzzzzz, g_y_xzzzzzz, g_y_yyyyyy, g_y_yyyyyz, g_y_yyyyzz, g_y_yyyzzz, g_y_yyzzzz, g_y_yzzzzz, g_y_zzzzzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_xy_xxxxxx[k] = -g_y_xxxxxx[k] * ab_x + g_y_xxxxxxx[k];

                g_xy_xxxxxy[k] = -g_y_xxxxxy[k] * ab_x + g_y_xxxxxxy[k];

                g_xy_xxxxxz[k] = -g_y_xxxxxz[k] * ab_x + g_y_xxxxxxz[k];

                g_xy_xxxxyy[k] = -g_y_xxxxyy[k] * ab_x + g_y_xxxxxyy[k];

                g_xy_xxxxyz[k] = -g_y_xxxxyz[k] * ab_x + g_y_xxxxxyz[k];

                g_xy_xxxxzz[k] = -g_y_xxxxzz[k] * ab_x + g_y_xxxxxzz[k];

                g_xy_xxxyyy[k] = -g_y_xxxyyy[k] * ab_x + g_y_xxxxyyy[k];

                g_xy_xxxyyz[k] = -g_y_xxxyyz[k] * ab_x + g_y_xxxxyyz[k];

                g_xy_xxxyzz[k] = -g_y_xxxyzz[k] * ab_x + g_y_xxxxyzz[k];

                g_xy_xxxzzz[k] = -g_y_xxxzzz[k] * ab_x + g_y_xxxxzzz[k];

                g_xy_xxyyyy[k] = -g_y_xxyyyy[k] * ab_x + g_y_xxxyyyy[k];

                g_xy_xxyyyz[k] = -g_y_xxyyyz[k] * ab_x + g_y_xxxyyyz[k];

                g_xy_xxyyzz[k] = -g_y_xxyyzz[k] * ab_x + g_y_xxxyyzz[k];

                g_xy_xxyzzz[k] = -g_y_xxyzzz[k] * ab_x + g_y_xxxyzzz[k];

                g_xy_xxzzzz[k] = -g_y_xxzzzz[k] * ab_x + g_y_xxxzzzz[k];

                g_xy_xyyyyy[k] = -g_y_xyyyyy[k] * ab_x + g_y_xxyyyyy[k];

                g_xy_xyyyyz[k] = -g_y_xyyyyz[k] * ab_x + g_y_xxyyyyz[k];

                g_xy_xyyyzz[k] = -g_y_xyyyzz[k] * ab_x + g_y_xxyyyzz[k];

                g_xy_xyyzzz[k] = -g_y_xyyzzz[k] * ab_x + g_y_xxyyzzz[k];

                g_xy_xyzzzz[k] = -g_y_xyzzzz[k] * ab_x + g_y_xxyzzzz[k];

                g_xy_xzzzzz[k] = -g_y_xzzzzz[k] * ab_x + g_y_xxzzzzz[k];

                g_xy_yyyyyy[k] = -g_y_yyyyyy[k] * ab_x + g_y_xyyyyyy[k];

                g_xy_yyyyyz[k] = -g_y_yyyyyz[k] * ab_x + g_y_xyyyyyz[k];

                g_xy_yyyyzz[k] = -g_y_yyyyzz[k] * ab_x + g_y_xyyyyzz[k];

                g_xy_yyyzzz[k] = -g_y_yyyzzz[k] * ab_x + g_y_xyyyzzz[k];

                g_xy_yyzzzz[k] = -g_y_yyzzzz[k] * ab_x + g_y_xyyzzzz[k];

                g_xy_yzzzzz[k] = -g_y_yzzzzz[k] * ab_x + g_y_xyzzzzz[k];

                g_xy_zzzzzz[k] = -g_y_zzzzzz[k] * ab_x + g_y_xzzzzzz[k];
            }

            /// Set up 56-84 components of targeted buffer : contr_buffer_dixx

            auto g_xz_xxxxxx = contr_buffer_dixx[di_off + 56 * ccomps * dcomps];

            auto g_xz_xxxxxy = contr_buffer_dixx[di_off + 57 * ccomps * dcomps];

            auto g_xz_xxxxxz = contr_buffer_dixx[di_off + 58 * ccomps * dcomps];

            auto g_xz_xxxxyy = contr_buffer_dixx[di_off + 59 * ccomps * dcomps];

            auto g_xz_xxxxyz = contr_buffer_dixx[di_off + 60 * ccomps * dcomps];

            auto g_xz_xxxxzz = contr_buffer_dixx[di_off + 61 * ccomps * dcomps];

            auto g_xz_xxxyyy = contr_buffer_dixx[di_off + 62 * ccomps * dcomps];

            auto g_xz_xxxyyz = contr_buffer_dixx[di_off + 63 * ccomps * dcomps];

            auto g_xz_xxxyzz = contr_buffer_dixx[di_off + 64 * ccomps * dcomps];

            auto g_xz_xxxzzz = contr_buffer_dixx[di_off + 65 * ccomps * dcomps];

            auto g_xz_xxyyyy = contr_buffer_dixx[di_off + 66 * ccomps * dcomps];

            auto g_xz_xxyyyz = contr_buffer_dixx[di_off + 67 * ccomps * dcomps];

            auto g_xz_xxyyzz = contr_buffer_dixx[di_off + 68 * ccomps * dcomps];

            auto g_xz_xxyzzz = contr_buffer_dixx[di_off + 69 * ccomps * dcomps];

            auto g_xz_xxzzzz = contr_buffer_dixx[di_off + 70 * ccomps * dcomps];

            auto g_xz_xyyyyy = contr_buffer_dixx[di_off + 71 * ccomps * dcomps];

            auto g_xz_xyyyyz = contr_buffer_dixx[di_off + 72 * ccomps * dcomps];

            auto g_xz_xyyyzz = contr_buffer_dixx[di_off + 73 * ccomps * dcomps];

            auto g_xz_xyyzzz = contr_buffer_dixx[di_off + 74 * ccomps * dcomps];

            auto g_xz_xyzzzz = contr_buffer_dixx[di_off + 75 * ccomps * dcomps];

            auto g_xz_xzzzzz = contr_buffer_dixx[di_off + 76 * ccomps * dcomps];

            auto g_xz_yyyyyy = contr_buffer_dixx[di_off + 77 * ccomps * dcomps];

            auto g_xz_yyyyyz = contr_buffer_dixx[di_off + 78 * ccomps * dcomps];

            auto g_xz_yyyyzz = contr_buffer_dixx[di_off + 79 * ccomps * dcomps];

            auto g_xz_yyyzzz = contr_buffer_dixx[di_off + 80 * ccomps * dcomps];

            auto g_xz_yyzzzz = contr_buffer_dixx[di_off + 81 * ccomps * dcomps];

            auto g_xz_yzzzzz = contr_buffer_dixx[di_off + 82 * ccomps * dcomps];

            auto g_xz_zzzzzz = contr_buffer_dixx[di_off + 83 * ccomps * dcomps];

            #pragma omp simd aligned(g_xz_xxxxxx, g_xz_xxxxxy, g_xz_xxxxxz, g_xz_xxxxyy, g_xz_xxxxyz, g_xz_xxxxzz, g_xz_xxxyyy, g_xz_xxxyyz, g_xz_xxxyzz, g_xz_xxxzzz, g_xz_xxyyyy, g_xz_xxyyyz, g_xz_xxyyzz, g_xz_xxyzzz, g_xz_xxzzzz, g_xz_xyyyyy, g_xz_xyyyyz, g_xz_xyyyzz, g_xz_xyyzzz, g_xz_xyzzzz, g_xz_xzzzzz, g_xz_yyyyyy, g_xz_yyyyyz, g_xz_yyyyzz, g_xz_yyyzzz, g_xz_yyzzzz, g_xz_yzzzzz, g_xz_zzzzzz, g_z_xxxxxx, g_z_xxxxxxx, g_z_xxxxxxy, g_z_xxxxxxz, g_z_xxxxxy, g_z_xxxxxyy, g_z_xxxxxyz, g_z_xxxxxz, g_z_xxxxxzz, g_z_xxxxyy, g_z_xxxxyyy, g_z_xxxxyyz, g_z_xxxxyz, g_z_xxxxyzz, g_z_xxxxzz, g_z_xxxxzzz, g_z_xxxyyy, g_z_xxxyyyy, g_z_xxxyyyz, g_z_xxxyyz, g_z_xxxyyzz, g_z_xxxyzz, g_z_xxxyzzz, g_z_xxxzzz, g_z_xxxzzzz, g_z_xxyyyy, g_z_xxyyyyy, g_z_xxyyyyz, g_z_xxyyyz, g_z_xxyyyzz, g_z_xxyyzz, g_z_xxyyzzz, g_z_xxyzzz, g_z_xxyzzzz, g_z_xxzzzz, g_z_xxzzzzz, g_z_xyyyyy, g_z_xyyyyyy, g_z_xyyyyyz, g_z_xyyyyz, g_z_xyyyyzz, g_z_xyyyzz, g_z_xyyyzzz, g_z_xyyzzz, g_z_xyyzzzz, g_z_xyzzzz, g_z_xyzzzzz, g_z_xzzzzz, g_z_xzzzzzz, g_z_yyyyyy, g_z_yyyyyz, g_z_yyyyzz, g_z_yyyzzz, g_z_yyzzzz, g_z_yzzzzz, g_z_zzzzzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_xz_xxxxxx[k] = -g_z_xxxxxx[k] * ab_x + g_z_xxxxxxx[k];

                g_xz_xxxxxy[k] = -g_z_xxxxxy[k] * ab_x + g_z_xxxxxxy[k];

                g_xz_xxxxxz[k] = -g_z_xxxxxz[k] * ab_x + g_z_xxxxxxz[k];

                g_xz_xxxxyy[k] = -g_z_xxxxyy[k] * ab_x + g_z_xxxxxyy[k];

                g_xz_xxxxyz[k] = -g_z_xxxxyz[k] * ab_x + g_z_xxxxxyz[k];

                g_xz_xxxxzz[k] = -g_z_xxxxzz[k] * ab_x + g_z_xxxxxzz[k];

                g_xz_xxxyyy[k] = -g_z_xxxyyy[k] * ab_x + g_z_xxxxyyy[k];

                g_xz_xxxyyz[k] = -g_z_xxxyyz[k] * ab_x + g_z_xxxxyyz[k];

                g_xz_xxxyzz[k] = -g_z_xxxyzz[k] * ab_x + g_z_xxxxyzz[k];

                g_xz_xxxzzz[k] = -g_z_xxxzzz[k] * ab_x + g_z_xxxxzzz[k];

                g_xz_xxyyyy[k] = -g_z_xxyyyy[k] * ab_x + g_z_xxxyyyy[k];

                g_xz_xxyyyz[k] = -g_z_xxyyyz[k] * ab_x + g_z_xxxyyyz[k];

                g_xz_xxyyzz[k] = -g_z_xxyyzz[k] * ab_x + g_z_xxxyyzz[k];

                g_xz_xxyzzz[k] = -g_z_xxyzzz[k] * ab_x + g_z_xxxyzzz[k];

                g_xz_xxzzzz[k] = -g_z_xxzzzz[k] * ab_x + g_z_xxxzzzz[k];

                g_xz_xyyyyy[k] = -g_z_xyyyyy[k] * ab_x + g_z_xxyyyyy[k];

                g_xz_xyyyyz[k] = -g_z_xyyyyz[k] * ab_x + g_z_xxyyyyz[k];

                g_xz_xyyyzz[k] = -g_z_xyyyzz[k] * ab_x + g_z_xxyyyzz[k];

                g_xz_xyyzzz[k] = -g_z_xyyzzz[k] * ab_x + g_z_xxyyzzz[k];

                g_xz_xyzzzz[k] = -g_z_xyzzzz[k] * ab_x + g_z_xxyzzzz[k];

                g_xz_xzzzzz[k] = -g_z_xzzzzz[k] * ab_x + g_z_xxzzzzz[k];

                g_xz_yyyyyy[k] = -g_z_yyyyyy[k] * ab_x + g_z_xyyyyyy[k];

                g_xz_yyyyyz[k] = -g_z_yyyyyz[k] * ab_x + g_z_xyyyyyz[k];

                g_xz_yyyyzz[k] = -g_z_yyyyzz[k] * ab_x + g_z_xyyyyzz[k];

                g_xz_yyyzzz[k] = -g_z_yyyzzz[k] * ab_x + g_z_xyyyzzz[k];

                g_xz_yyzzzz[k] = -g_z_yyzzzz[k] * ab_x + g_z_xyyzzzz[k];

                g_xz_yzzzzz[k] = -g_z_yzzzzz[k] * ab_x + g_z_xyzzzzz[k];

                g_xz_zzzzzz[k] = -g_z_zzzzzz[k] * ab_x + g_z_xzzzzzz[k];
            }

            /// Set up 84-112 components of targeted buffer : contr_buffer_dixx

            auto g_yy_xxxxxx = contr_buffer_dixx[di_off + 84 * ccomps * dcomps];

            auto g_yy_xxxxxy = contr_buffer_dixx[di_off + 85 * ccomps * dcomps];

            auto g_yy_xxxxxz = contr_buffer_dixx[di_off + 86 * ccomps * dcomps];

            auto g_yy_xxxxyy = contr_buffer_dixx[di_off + 87 * ccomps * dcomps];

            auto g_yy_xxxxyz = contr_buffer_dixx[di_off + 88 * ccomps * dcomps];

            auto g_yy_xxxxzz = contr_buffer_dixx[di_off + 89 * ccomps * dcomps];

            auto g_yy_xxxyyy = contr_buffer_dixx[di_off + 90 * ccomps * dcomps];

            auto g_yy_xxxyyz = contr_buffer_dixx[di_off + 91 * ccomps * dcomps];

            auto g_yy_xxxyzz = contr_buffer_dixx[di_off + 92 * ccomps * dcomps];

            auto g_yy_xxxzzz = contr_buffer_dixx[di_off + 93 * ccomps * dcomps];

            auto g_yy_xxyyyy = contr_buffer_dixx[di_off + 94 * ccomps * dcomps];

            auto g_yy_xxyyyz = contr_buffer_dixx[di_off + 95 * ccomps * dcomps];

            auto g_yy_xxyyzz = contr_buffer_dixx[di_off + 96 * ccomps * dcomps];

            auto g_yy_xxyzzz = contr_buffer_dixx[di_off + 97 * ccomps * dcomps];

            auto g_yy_xxzzzz = contr_buffer_dixx[di_off + 98 * ccomps * dcomps];

            auto g_yy_xyyyyy = contr_buffer_dixx[di_off + 99 * ccomps * dcomps];

            auto g_yy_xyyyyz = contr_buffer_dixx[di_off + 100 * ccomps * dcomps];

            auto g_yy_xyyyzz = contr_buffer_dixx[di_off + 101 * ccomps * dcomps];

            auto g_yy_xyyzzz = contr_buffer_dixx[di_off + 102 * ccomps * dcomps];

            auto g_yy_xyzzzz = contr_buffer_dixx[di_off + 103 * ccomps * dcomps];

            auto g_yy_xzzzzz = contr_buffer_dixx[di_off + 104 * ccomps * dcomps];

            auto g_yy_yyyyyy = contr_buffer_dixx[di_off + 105 * ccomps * dcomps];

            auto g_yy_yyyyyz = contr_buffer_dixx[di_off + 106 * ccomps * dcomps];

            auto g_yy_yyyyzz = contr_buffer_dixx[di_off + 107 * ccomps * dcomps];

            auto g_yy_yyyzzz = contr_buffer_dixx[di_off + 108 * ccomps * dcomps];

            auto g_yy_yyzzzz = contr_buffer_dixx[di_off + 109 * ccomps * dcomps];

            auto g_yy_yzzzzz = contr_buffer_dixx[di_off + 110 * ccomps * dcomps];

            auto g_yy_zzzzzz = contr_buffer_dixx[di_off + 111 * ccomps * dcomps];

            #pragma omp simd aligned(g_y_xxxxxx, g_y_xxxxxxy, g_y_xxxxxy, g_y_xxxxxyy, g_y_xxxxxyz, g_y_xxxxxz, g_y_xxxxyy, g_y_xxxxyyy, g_y_xxxxyyz, g_y_xxxxyz, g_y_xxxxyzz, g_y_xxxxzz, g_y_xxxyyy, g_y_xxxyyyy, g_y_xxxyyyz, g_y_xxxyyz, g_y_xxxyyzz, g_y_xxxyzz, g_y_xxxyzzz, g_y_xxxzzz, g_y_xxyyyy, g_y_xxyyyyy, g_y_xxyyyyz, g_y_xxyyyz, g_y_xxyyyzz, g_y_xxyyzz, g_y_xxyyzzz, g_y_xxyzzz, g_y_xxyzzzz, g_y_xxzzzz, g_y_xyyyyy, g_y_xyyyyyy, g_y_xyyyyyz, g_y_xyyyyz, g_y_xyyyyzz, g_y_xyyyzz, g_y_xyyyzzz, g_y_xyyzzz, g_y_xyyzzzz, g_y_xyzzzz, g_y_xyzzzzz, g_y_xzzzzz, g_y_yyyyyy, g_y_yyyyyyy, g_y_yyyyyyz, g_y_yyyyyz, g_y_yyyyyzz, g_y_yyyyzz, g_y_yyyyzzz, g_y_yyyzzz, g_y_yyyzzzz, g_y_yyzzzz, g_y_yyzzzzz, g_y_yzzzzz, g_y_yzzzzzz, g_y_zzzzzz, g_yy_xxxxxx, g_yy_xxxxxy, g_yy_xxxxxz, g_yy_xxxxyy, g_yy_xxxxyz, g_yy_xxxxzz, g_yy_xxxyyy, g_yy_xxxyyz, g_yy_xxxyzz, g_yy_xxxzzz, g_yy_xxyyyy, g_yy_xxyyyz, g_yy_xxyyzz, g_yy_xxyzzz, g_yy_xxzzzz, g_yy_xyyyyy, g_yy_xyyyyz, g_yy_xyyyzz, g_yy_xyyzzz, g_yy_xyzzzz, g_yy_xzzzzz, g_yy_yyyyyy, g_yy_yyyyyz, g_yy_yyyyzz, g_yy_yyyzzz, g_yy_yyzzzz, g_yy_yzzzzz, g_yy_zzzzzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_yy_xxxxxx[k] = -g_y_xxxxxx[k] * ab_y + g_y_xxxxxxy[k];

                g_yy_xxxxxy[k] = -g_y_xxxxxy[k] * ab_y + g_y_xxxxxyy[k];

                g_yy_xxxxxz[k] = -g_y_xxxxxz[k] * ab_y + g_y_xxxxxyz[k];

                g_yy_xxxxyy[k] = -g_y_xxxxyy[k] * ab_y + g_y_xxxxyyy[k];

                g_yy_xxxxyz[k] = -g_y_xxxxyz[k] * ab_y + g_y_xxxxyyz[k];

                g_yy_xxxxzz[k] = -g_y_xxxxzz[k] * ab_y + g_y_xxxxyzz[k];

                g_yy_xxxyyy[k] = -g_y_xxxyyy[k] * ab_y + g_y_xxxyyyy[k];

                g_yy_xxxyyz[k] = -g_y_xxxyyz[k] * ab_y + g_y_xxxyyyz[k];

                g_yy_xxxyzz[k] = -g_y_xxxyzz[k] * ab_y + g_y_xxxyyzz[k];

                g_yy_xxxzzz[k] = -g_y_xxxzzz[k] * ab_y + g_y_xxxyzzz[k];

                g_yy_xxyyyy[k] = -g_y_xxyyyy[k] * ab_y + g_y_xxyyyyy[k];

                g_yy_xxyyyz[k] = -g_y_xxyyyz[k] * ab_y + g_y_xxyyyyz[k];

                g_yy_xxyyzz[k] = -g_y_xxyyzz[k] * ab_y + g_y_xxyyyzz[k];

                g_yy_xxyzzz[k] = -g_y_xxyzzz[k] * ab_y + g_y_xxyyzzz[k];

                g_yy_xxzzzz[k] = -g_y_xxzzzz[k] * ab_y + g_y_xxyzzzz[k];

                g_yy_xyyyyy[k] = -g_y_xyyyyy[k] * ab_y + g_y_xyyyyyy[k];

                g_yy_xyyyyz[k] = -g_y_xyyyyz[k] * ab_y + g_y_xyyyyyz[k];

                g_yy_xyyyzz[k] = -g_y_xyyyzz[k] * ab_y + g_y_xyyyyzz[k];

                g_yy_xyyzzz[k] = -g_y_xyyzzz[k] * ab_y + g_y_xyyyzzz[k];

                g_yy_xyzzzz[k] = -g_y_xyzzzz[k] * ab_y + g_y_xyyzzzz[k];

                g_yy_xzzzzz[k] = -g_y_xzzzzz[k] * ab_y + g_y_xyzzzzz[k];

                g_yy_yyyyyy[k] = -g_y_yyyyyy[k] * ab_y + g_y_yyyyyyy[k];

                g_yy_yyyyyz[k] = -g_y_yyyyyz[k] * ab_y + g_y_yyyyyyz[k];

                g_yy_yyyyzz[k] = -g_y_yyyyzz[k] * ab_y + g_y_yyyyyzz[k];

                g_yy_yyyzzz[k] = -g_y_yyyzzz[k] * ab_y + g_y_yyyyzzz[k];

                g_yy_yyzzzz[k] = -g_y_yyzzzz[k] * ab_y + g_y_yyyzzzz[k];

                g_yy_yzzzzz[k] = -g_y_yzzzzz[k] * ab_y + g_y_yyzzzzz[k];

                g_yy_zzzzzz[k] = -g_y_zzzzzz[k] * ab_y + g_y_yzzzzzz[k];
            }

            /// Set up 112-140 components of targeted buffer : contr_buffer_dixx

            auto g_yz_xxxxxx = contr_buffer_dixx[di_off + 112 * ccomps * dcomps];

            auto g_yz_xxxxxy = contr_buffer_dixx[di_off + 113 * ccomps * dcomps];

            auto g_yz_xxxxxz = contr_buffer_dixx[di_off + 114 * ccomps * dcomps];

            auto g_yz_xxxxyy = contr_buffer_dixx[di_off + 115 * ccomps * dcomps];

            auto g_yz_xxxxyz = contr_buffer_dixx[di_off + 116 * ccomps * dcomps];

            auto g_yz_xxxxzz = contr_buffer_dixx[di_off + 117 * ccomps * dcomps];

            auto g_yz_xxxyyy = contr_buffer_dixx[di_off + 118 * ccomps * dcomps];

            auto g_yz_xxxyyz = contr_buffer_dixx[di_off + 119 * ccomps * dcomps];

            auto g_yz_xxxyzz = contr_buffer_dixx[di_off + 120 * ccomps * dcomps];

            auto g_yz_xxxzzz = contr_buffer_dixx[di_off + 121 * ccomps * dcomps];

            auto g_yz_xxyyyy = contr_buffer_dixx[di_off + 122 * ccomps * dcomps];

            auto g_yz_xxyyyz = contr_buffer_dixx[di_off + 123 * ccomps * dcomps];

            auto g_yz_xxyyzz = contr_buffer_dixx[di_off + 124 * ccomps * dcomps];

            auto g_yz_xxyzzz = contr_buffer_dixx[di_off + 125 * ccomps * dcomps];

            auto g_yz_xxzzzz = contr_buffer_dixx[di_off + 126 * ccomps * dcomps];

            auto g_yz_xyyyyy = contr_buffer_dixx[di_off + 127 * ccomps * dcomps];

            auto g_yz_xyyyyz = contr_buffer_dixx[di_off + 128 * ccomps * dcomps];

            auto g_yz_xyyyzz = contr_buffer_dixx[di_off + 129 * ccomps * dcomps];

            auto g_yz_xyyzzz = contr_buffer_dixx[di_off + 130 * ccomps * dcomps];

            auto g_yz_xyzzzz = contr_buffer_dixx[di_off + 131 * ccomps * dcomps];

            auto g_yz_xzzzzz = contr_buffer_dixx[di_off + 132 * ccomps * dcomps];

            auto g_yz_yyyyyy = contr_buffer_dixx[di_off + 133 * ccomps * dcomps];

            auto g_yz_yyyyyz = contr_buffer_dixx[di_off + 134 * ccomps * dcomps];

            auto g_yz_yyyyzz = contr_buffer_dixx[di_off + 135 * ccomps * dcomps];

            auto g_yz_yyyzzz = contr_buffer_dixx[di_off + 136 * ccomps * dcomps];

            auto g_yz_yyzzzz = contr_buffer_dixx[di_off + 137 * ccomps * dcomps];

            auto g_yz_yzzzzz = contr_buffer_dixx[di_off + 138 * ccomps * dcomps];

            auto g_yz_zzzzzz = contr_buffer_dixx[di_off + 139 * ccomps * dcomps];

            #pragma omp simd aligned(g_yz_xxxxxx, g_yz_xxxxxy, g_yz_xxxxxz, g_yz_xxxxyy, g_yz_xxxxyz, g_yz_xxxxzz, g_yz_xxxyyy, g_yz_xxxyyz, g_yz_xxxyzz, g_yz_xxxzzz, g_yz_xxyyyy, g_yz_xxyyyz, g_yz_xxyyzz, g_yz_xxyzzz, g_yz_xxzzzz, g_yz_xyyyyy, g_yz_xyyyyz, g_yz_xyyyzz, g_yz_xyyzzz, g_yz_xyzzzz, g_yz_xzzzzz, g_yz_yyyyyy, g_yz_yyyyyz, g_yz_yyyyzz, g_yz_yyyzzz, g_yz_yyzzzz, g_yz_yzzzzz, g_yz_zzzzzz, g_z_xxxxxx, g_z_xxxxxxy, g_z_xxxxxy, g_z_xxxxxyy, g_z_xxxxxyz, g_z_xxxxxz, g_z_xxxxyy, g_z_xxxxyyy, g_z_xxxxyyz, g_z_xxxxyz, g_z_xxxxyzz, g_z_xxxxzz, g_z_xxxyyy, g_z_xxxyyyy, g_z_xxxyyyz, g_z_xxxyyz, g_z_xxxyyzz, g_z_xxxyzz, g_z_xxxyzzz, g_z_xxxzzz, g_z_xxyyyy, g_z_xxyyyyy, g_z_xxyyyyz, g_z_xxyyyz, g_z_xxyyyzz, g_z_xxyyzz, g_z_xxyyzzz, g_z_xxyzzz, g_z_xxyzzzz, g_z_xxzzzz, g_z_xyyyyy, g_z_xyyyyyy, g_z_xyyyyyz, g_z_xyyyyz, g_z_xyyyyzz, g_z_xyyyzz, g_z_xyyyzzz, g_z_xyyzzz, g_z_xyyzzzz, g_z_xyzzzz, g_z_xyzzzzz, g_z_xzzzzz, g_z_yyyyyy, g_z_yyyyyyy, g_z_yyyyyyz, g_z_yyyyyz, g_z_yyyyyzz, g_z_yyyyzz, g_z_yyyyzzz, g_z_yyyzzz, g_z_yyyzzzz, g_z_yyzzzz, g_z_yyzzzzz, g_z_yzzzzz, g_z_yzzzzzz, g_z_zzzzzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_yz_xxxxxx[k] = -g_z_xxxxxx[k] * ab_y + g_z_xxxxxxy[k];

                g_yz_xxxxxy[k] = -g_z_xxxxxy[k] * ab_y + g_z_xxxxxyy[k];

                g_yz_xxxxxz[k] = -g_z_xxxxxz[k] * ab_y + g_z_xxxxxyz[k];

                g_yz_xxxxyy[k] = -g_z_xxxxyy[k] * ab_y + g_z_xxxxyyy[k];

                g_yz_xxxxyz[k] = -g_z_xxxxyz[k] * ab_y + g_z_xxxxyyz[k];

                g_yz_xxxxzz[k] = -g_z_xxxxzz[k] * ab_y + g_z_xxxxyzz[k];

                g_yz_xxxyyy[k] = -g_z_xxxyyy[k] * ab_y + g_z_xxxyyyy[k];

                g_yz_xxxyyz[k] = -g_z_xxxyyz[k] * ab_y + g_z_xxxyyyz[k];

                g_yz_xxxyzz[k] = -g_z_xxxyzz[k] * ab_y + g_z_xxxyyzz[k];

                g_yz_xxxzzz[k] = -g_z_xxxzzz[k] * ab_y + g_z_xxxyzzz[k];

                g_yz_xxyyyy[k] = -g_z_xxyyyy[k] * ab_y + g_z_xxyyyyy[k];

                g_yz_xxyyyz[k] = -g_z_xxyyyz[k] * ab_y + g_z_xxyyyyz[k];

                g_yz_xxyyzz[k] = -g_z_xxyyzz[k] * ab_y + g_z_xxyyyzz[k];

                g_yz_xxyzzz[k] = -g_z_xxyzzz[k] * ab_y + g_z_xxyyzzz[k];

                g_yz_xxzzzz[k] = -g_z_xxzzzz[k] * ab_y + g_z_xxyzzzz[k];

                g_yz_xyyyyy[k] = -g_z_xyyyyy[k] * ab_y + g_z_xyyyyyy[k];

                g_yz_xyyyyz[k] = -g_z_xyyyyz[k] * ab_y + g_z_xyyyyyz[k];

                g_yz_xyyyzz[k] = -g_z_xyyyzz[k] * ab_y + g_z_xyyyyzz[k];

                g_yz_xyyzzz[k] = -g_z_xyyzzz[k] * ab_y + g_z_xyyyzzz[k];

                g_yz_xyzzzz[k] = -g_z_xyzzzz[k] * ab_y + g_z_xyyzzzz[k];

                g_yz_xzzzzz[k] = -g_z_xzzzzz[k] * ab_y + g_z_xyzzzzz[k];

                g_yz_yyyyyy[k] = -g_z_yyyyyy[k] * ab_y + g_z_yyyyyyy[k];

                g_yz_yyyyyz[k] = -g_z_yyyyyz[k] * ab_y + g_z_yyyyyyz[k];

                g_yz_yyyyzz[k] = -g_z_yyyyzz[k] * ab_y + g_z_yyyyyzz[k];

                g_yz_yyyzzz[k] = -g_z_yyyzzz[k] * ab_y + g_z_yyyyzzz[k];

                g_yz_yyzzzz[k] = -g_z_yyzzzz[k] * ab_y + g_z_yyyzzzz[k];

                g_yz_yzzzzz[k] = -g_z_yzzzzz[k] * ab_y + g_z_yyzzzzz[k];

                g_yz_zzzzzz[k] = -g_z_zzzzzz[k] * ab_y + g_z_yzzzzzz[k];
            }

            /// Set up 140-168 components of targeted buffer : contr_buffer_dixx

            auto g_zz_xxxxxx = contr_buffer_dixx[di_off + 140 * ccomps * dcomps];

            auto g_zz_xxxxxy = contr_buffer_dixx[di_off + 141 * ccomps * dcomps];

            auto g_zz_xxxxxz = contr_buffer_dixx[di_off + 142 * ccomps * dcomps];

            auto g_zz_xxxxyy = contr_buffer_dixx[di_off + 143 * ccomps * dcomps];

            auto g_zz_xxxxyz = contr_buffer_dixx[di_off + 144 * ccomps * dcomps];

            auto g_zz_xxxxzz = contr_buffer_dixx[di_off + 145 * ccomps * dcomps];

            auto g_zz_xxxyyy = contr_buffer_dixx[di_off + 146 * ccomps * dcomps];

            auto g_zz_xxxyyz = contr_buffer_dixx[di_off + 147 * ccomps * dcomps];

            auto g_zz_xxxyzz = contr_buffer_dixx[di_off + 148 * ccomps * dcomps];

            auto g_zz_xxxzzz = contr_buffer_dixx[di_off + 149 * ccomps * dcomps];

            auto g_zz_xxyyyy = contr_buffer_dixx[di_off + 150 * ccomps * dcomps];

            auto g_zz_xxyyyz = contr_buffer_dixx[di_off + 151 * ccomps * dcomps];

            auto g_zz_xxyyzz = contr_buffer_dixx[di_off + 152 * ccomps * dcomps];

            auto g_zz_xxyzzz = contr_buffer_dixx[di_off + 153 * ccomps * dcomps];

            auto g_zz_xxzzzz = contr_buffer_dixx[di_off + 154 * ccomps * dcomps];

            auto g_zz_xyyyyy = contr_buffer_dixx[di_off + 155 * ccomps * dcomps];

            auto g_zz_xyyyyz = contr_buffer_dixx[di_off + 156 * ccomps * dcomps];

            auto g_zz_xyyyzz = contr_buffer_dixx[di_off + 157 * ccomps * dcomps];

            auto g_zz_xyyzzz = contr_buffer_dixx[di_off + 158 * ccomps * dcomps];

            auto g_zz_xyzzzz = contr_buffer_dixx[di_off + 159 * ccomps * dcomps];

            auto g_zz_xzzzzz = contr_buffer_dixx[di_off + 160 * ccomps * dcomps];

            auto g_zz_yyyyyy = contr_buffer_dixx[di_off + 161 * ccomps * dcomps];

            auto g_zz_yyyyyz = contr_buffer_dixx[di_off + 162 * ccomps * dcomps];

            auto g_zz_yyyyzz = contr_buffer_dixx[di_off + 163 * ccomps * dcomps];

            auto g_zz_yyyzzz = contr_buffer_dixx[di_off + 164 * ccomps * dcomps];

            auto g_zz_yyzzzz = contr_buffer_dixx[di_off + 165 * ccomps * dcomps];

            auto g_zz_yzzzzz = contr_buffer_dixx[di_off + 166 * ccomps * dcomps];

            auto g_zz_zzzzzz = contr_buffer_dixx[di_off + 167 * ccomps * dcomps];

            #pragma omp simd aligned(g_z_xxxxxx, g_z_xxxxxxz, g_z_xxxxxy, g_z_xxxxxyz, g_z_xxxxxz, g_z_xxxxxzz, g_z_xxxxyy, g_z_xxxxyyz, g_z_xxxxyz, g_z_xxxxyzz, g_z_xxxxzz, g_z_xxxxzzz, g_z_xxxyyy, g_z_xxxyyyz, g_z_xxxyyz, g_z_xxxyyzz, g_z_xxxyzz, g_z_xxxyzzz, g_z_xxxzzz, g_z_xxxzzzz, g_z_xxyyyy, g_z_xxyyyyz, g_z_xxyyyz, g_z_xxyyyzz, g_z_xxyyzz, g_z_xxyyzzz, g_z_xxyzzz, g_z_xxyzzzz, g_z_xxzzzz, g_z_xxzzzzz, g_z_xyyyyy, g_z_xyyyyyz, g_z_xyyyyz, g_z_xyyyyzz, g_z_xyyyzz, g_z_xyyyzzz, g_z_xyyzzz, g_z_xyyzzzz, g_z_xyzzzz, g_z_xyzzzzz, g_z_xzzzzz, g_z_xzzzzzz, g_z_yyyyyy, g_z_yyyyyyz, g_z_yyyyyz, g_z_yyyyyzz, g_z_yyyyzz, g_z_yyyyzzz, g_z_yyyzzz, g_z_yyyzzzz, g_z_yyzzzz, g_z_yyzzzzz, g_z_yzzzzz, g_z_yzzzzzz, g_z_zzzzzz, g_z_zzzzzzz, g_zz_xxxxxx, g_zz_xxxxxy, g_zz_xxxxxz, g_zz_xxxxyy, g_zz_xxxxyz, g_zz_xxxxzz, g_zz_xxxyyy, g_zz_xxxyyz, g_zz_xxxyzz, g_zz_xxxzzz, g_zz_xxyyyy, g_zz_xxyyyz, g_zz_xxyyzz, g_zz_xxyzzz, g_zz_xxzzzz, g_zz_xyyyyy, g_zz_xyyyyz, g_zz_xyyyzz, g_zz_xyyzzz, g_zz_xyzzzz, g_zz_xzzzzz, g_zz_yyyyyy, g_zz_yyyyyz, g_zz_yyyyzz, g_zz_yyyzzz, g_zz_yyzzzz, g_zz_yzzzzz, g_zz_zzzzzz  : 64)
            for (int k = 0; k < ndims; k++)
            {
                g_zz_xxxxxx[k] = -g_z_xxxxxx[k] * ab_z + g_z_xxxxxxz[k];

                g_zz_xxxxxy[k] = -g_z_xxxxxy[k] * ab_z + g_z_xxxxxyz[k];

                g_zz_xxxxxz[k] = -g_z_xxxxxz[k] * ab_z + g_z_xxxxxzz[k];

                g_zz_xxxxyy[k] = -g_z_xxxxyy[k] * ab_z + g_z_xxxxyyz[k];

                g_zz_xxxxyz[k] = -g_z_xxxxyz[k] * ab_z + g_z_xxxxyzz[k];

                g_zz_xxxxzz[k] = -g_z_xxxxzz[k] * ab_z + g_z_xxxxzzz[k];

                g_zz_xxxyyy[k] = -g_z_xxxyyy[k] * ab_z + g_z_xxxyyyz[k];

                g_zz_xxxyyz[k] = -g_z_xxxyyz[k] * ab_z + g_z_xxxyyzz[k];

                g_zz_xxxyzz[k] = -g_z_xxxyzz[k] * ab_z + g_z_xxxyzzz[k];

                g_zz_xxxzzz[k] = -g_z_xxxzzz[k] * ab_z + g_z_xxxzzzz[k];

                g_zz_xxyyyy[k] = -g_z_xxyyyy[k] * ab_z + g_z_xxyyyyz[k];

                g_zz_xxyyyz[k] = -g_z_xxyyyz[k] * ab_z + g_z_xxyyyzz[k];

                g_zz_xxyyzz[k] = -g_z_xxyyzz[k] * ab_z + g_z_xxyyzzz[k];

                g_zz_xxyzzz[k] = -g_z_xxyzzz[k] * ab_z + g_z_xxyzzzz[k];

                g_zz_xxzzzz[k] = -g_z_xxzzzz[k] * ab_z + g_z_xxzzzzz[k];

                g_zz_xyyyyy[k] = -g_z_xyyyyy[k] * ab_z + g_z_xyyyyyz[k];

                g_zz_xyyyyz[k] = -g_z_xyyyyz[k] * ab_z + g_z_xyyyyzz[k];

                g_zz_xyyyzz[k] = -g_z_xyyyzz[k] * ab_z + g_z_xyyyzzz[k];

                g_zz_xyyzzz[k] = -g_z_xyyzzz[k] * ab_z + g_z_xyyzzzz[k];

                g_zz_xyzzzz[k] = -g_z_xyzzzz[k] * ab_z + g_z_xyzzzzz[k];

                g_zz_xzzzzz[k] = -g_z_xzzzzz[k] * ab_z + g_z_xzzzzzz[k];

                g_zz_yyyyyy[k] = -g_z_yyyyyy[k] * ab_z + g_z_yyyyyyz[k];

                g_zz_yyyyyz[k] = -g_z_yyyyyz[k] * ab_z + g_z_yyyyyzz[k];

                g_zz_yyyyzz[k] = -g_z_yyyyzz[k] * ab_z + g_z_yyyyzzz[k];

                g_zz_yyyzzz[k] = -g_z_yyyzzz[k] * ab_z + g_z_yyyzzzz[k];

                g_zz_yyzzzz[k] = -g_z_yyzzzz[k] * ab_z + g_z_yyzzzzz[k];

                g_zz_yzzzzz[k] = -g_z_yzzzzz[k] * ab_z + g_z_yzzzzzz[k];

                g_zz_zzzzzz[k] = -g_z_zzzzzz[k] * ab_z + g_z_zzzzzzz[k];
            }
        }
    }
}

} // erirec namespace
