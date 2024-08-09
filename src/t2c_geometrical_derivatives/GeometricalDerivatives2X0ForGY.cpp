#include "GeometricalDerivatives2X0ForGY.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_20_gx(CSimdArray<double>& pbuffer,
                        const size_t idx_op_geom_200_gs,
                        const size_t idx_op_ds,
                        const size_t idx_op_gs,
                        const size_t idx_op_is,
                        const size_t op_comps,
                        const size_t ket_comps,
                        const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    for (size_t i = 0; i < op_comps; i++)
    {
        for (size_t j = 0; j < ket_comps; j++)
        {
            // Set up components of auxiliary buffer : DS

            auto to_xx_0 = pbuffer.data(idx_op_ds + i * 6 * ket_comps + 0 * ket_comps + j);

            auto to_xy_0 = pbuffer.data(idx_op_ds + i * 6 * ket_comps + 1 * ket_comps + j);

            auto to_xz_0 = pbuffer.data(idx_op_ds + i * 6 * ket_comps + 2 * ket_comps + j);

            auto to_yy_0 = pbuffer.data(idx_op_ds + i * 6 * ket_comps + 3 * ket_comps + j);

            auto to_yz_0 = pbuffer.data(idx_op_ds + i * 6 * ket_comps + 4 * ket_comps + j);

            auto to_zz_0 = pbuffer.data(idx_op_ds + i * 6 * ket_comps + 5 * ket_comps + j);

            // Set up components of auxiliary buffer : GS

            auto to_xxxx_0 = pbuffer.data(idx_op_gs + i * 15 * ket_comps + 0 * ket_comps + j);

            auto to_xxxy_0 = pbuffer.data(idx_op_gs + i * 15 * ket_comps + 1 * ket_comps + j);

            auto to_xxxz_0 = pbuffer.data(idx_op_gs + i * 15 * ket_comps + 2 * ket_comps + j);

            auto to_xxyy_0 = pbuffer.data(idx_op_gs + i * 15 * ket_comps + 3 * ket_comps + j);

            auto to_xxyz_0 = pbuffer.data(idx_op_gs + i * 15 * ket_comps + 4 * ket_comps + j);

            auto to_xxzz_0 = pbuffer.data(idx_op_gs + i * 15 * ket_comps + 5 * ket_comps + j);

            auto to_xyyy_0 = pbuffer.data(idx_op_gs + i * 15 * ket_comps + 6 * ket_comps + j);

            auto to_xyyz_0 = pbuffer.data(idx_op_gs + i * 15 * ket_comps + 7 * ket_comps + j);

            auto to_xyzz_0 = pbuffer.data(idx_op_gs + i * 15 * ket_comps + 8 * ket_comps + j);

            auto to_xzzz_0 = pbuffer.data(idx_op_gs + i * 15 * ket_comps + 9 * ket_comps + j);

            auto to_yyyy_0 = pbuffer.data(idx_op_gs + i * 15 * ket_comps + 10 * ket_comps + j);

            auto to_yyyz_0 = pbuffer.data(idx_op_gs + i * 15 * ket_comps + 11 * ket_comps + j);

            auto to_yyzz_0 = pbuffer.data(idx_op_gs + i * 15 * ket_comps + 12 * ket_comps + j);

            auto to_yzzz_0 = pbuffer.data(idx_op_gs + i * 15 * ket_comps + 13 * ket_comps + j);

            auto to_zzzz_0 = pbuffer.data(idx_op_gs + i * 15 * ket_comps + 14 * ket_comps + j);

            // Set up components of auxiliary buffer : IS

            auto to_xxxxxx_0 = pbuffer.data(idx_op_is + i * 28 * ket_comps + 0 * ket_comps + j);

            auto to_xxxxxy_0 = pbuffer.data(idx_op_is + i * 28 * ket_comps + 1 * ket_comps + j);

            auto to_xxxxxz_0 = pbuffer.data(idx_op_is + i * 28 * ket_comps + 2 * ket_comps + j);

            auto to_xxxxyy_0 = pbuffer.data(idx_op_is + i * 28 * ket_comps + 3 * ket_comps + j);

            auto to_xxxxyz_0 = pbuffer.data(idx_op_is + i * 28 * ket_comps + 4 * ket_comps + j);

            auto to_xxxxzz_0 = pbuffer.data(idx_op_is + i * 28 * ket_comps + 5 * ket_comps + j);

            auto to_xxxyyy_0 = pbuffer.data(idx_op_is + i * 28 * ket_comps + 6 * ket_comps + j);

            auto to_xxxyyz_0 = pbuffer.data(idx_op_is + i * 28 * ket_comps + 7 * ket_comps + j);

            auto to_xxxyzz_0 = pbuffer.data(idx_op_is + i * 28 * ket_comps + 8 * ket_comps + j);

            auto to_xxxzzz_0 = pbuffer.data(idx_op_is + i * 28 * ket_comps + 9 * ket_comps + j);

            auto to_xxyyyy_0 = pbuffer.data(idx_op_is + i * 28 * ket_comps + 10 * ket_comps + j);

            auto to_xxyyyz_0 = pbuffer.data(idx_op_is + i * 28 * ket_comps + 11 * ket_comps + j);

            auto to_xxyyzz_0 = pbuffer.data(idx_op_is + i * 28 * ket_comps + 12 * ket_comps + j);

            auto to_xxyzzz_0 = pbuffer.data(idx_op_is + i * 28 * ket_comps + 13 * ket_comps + j);

            auto to_xxzzzz_0 = pbuffer.data(idx_op_is + i * 28 * ket_comps + 14 * ket_comps + j);

            auto to_xyyyyy_0 = pbuffer.data(idx_op_is + i * 28 * ket_comps + 15 * ket_comps + j);

            auto to_xyyyyz_0 = pbuffer.data(idx_op_is + i * 28 * ket_comps + 16 * ket_comps + j);

            auto to_xyyyzz_0 = pbuffer.data(idx_op_is + i * 28 * ket_comps + 17 * ket_comps + j);

            auto to_xyyzzz_0 = pbuffer.data(idx_op_is + i * 28 * ket_comps + 18 * ket_comps + j);

            auto to_xyzzzz_0 = pbuffer.data(idx_op_is + i * 28 * ket_comps + 19 * ket_comps + j);

            auto to_xzzzzz_0 = pbuffer.data(idx_op_is + i * 28 * ket_comps + 20 * ket_comps + j);

            auto to_yyyyyy_0 = pbuffer.data(idx_op_is + i * 28 * ket_comps + 21 * ket_comps + j);

            auto to_yyyyyz_0 = pbuffer.data(idx_op_is + i * 28 * ket_comps + 22 * ket_comps + j);

            auto to_yyyyzz_0 = pbuffer.data(idx_op_is + i * 28 * ket_comps + 23 * ket_comps + j);

            auto to_yyyzzz_0 = pbuffer.data(idx_op_is + i * 28 * ket_comps + 24 * ket_comps + j);

            auto to_yyzzzz_0 = pbuffer.data(idx_op_is + i * 28 * ket_comps + 25 * ket_comps + j);

            auto to_yzzzzz_0 = pbuffer.data(idx_op_is + i * 28 * ket_comps + 26 * ket_comps + j);

            auto to_zzzzzz_0 = pbuffer.data(idx_op_is + i * 28 * ket_comps + 27 * ket_comps + j);

            // Set up 0-15 components of targeted buffer : GS

            auto to_xx_0_xxxx_0 = pbuffer.data(idx_op_geom_200_gs + 0 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 0 * ket_comps + j);

            auto to_xx_0_xxxy_0 = pbuffer.data(idx_op_geom_200_gs + 0 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 1 * ket_comps + j);

            auto to_xx_0_xxxz_0 = pbuffer.data(idx_op_geom_200_gs + 0 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 2 * ket_comps + j);

            auto to_xx_0_xxyy_0 = pbuffer.data(idx_op_geom_200_gs + 0 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 3 * ket_comps + j);

            auto to_xx_0_xxyz_0 = pbuffer.data(idx_op_geom_200_gs + 0 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 4 * ket_comps + j);

            auto to_xx_0_xxzz_0 = pbuffer.data(idx_op_geom_200_gs + 0 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 5 * ket_comps + j);

            auto to_xx_0_xyyy_0 = pbuffer.data(idx_op_geom_200_gs + 0 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 6 * ket_comps + j);

            auto to_xx_0_xyyz_0 = pbuffer.data(idx_op_geom_200_gs + 0 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 7 * ket_comps + j);

            auto to_xx_0_xyzz_0 = pbuffer.data(idx_op_geom_200_gs + 0 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 8 * ket_comps + j);

            auto to_xx_0_xzzz_0 = pbuffer.data(idx_op_geom_200_gs + 0 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 9 * ket_comps + j);

            auto to_xx_0_yyyy_0 = pbuffer.data(idx_op_geom_200_gs + 0 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 10 * ket_comps + j);

            auto to_xx_0_yyyz_0 = pbuffer.data(idx_op_geom_200_gs + 0 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 11 * ket_comps + j);

            auto to_xx_0_yyzz_0 = pbuffer.data(idx_op_geom_200_gs + 0 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 12 * ket_comps + j);

            auto to_xx_0_yzzz_0 = pbuffer.data(idx_op_geom_200_gs + 0 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 13 * ket_comps + j);

            auto to_xx_0_zzzz_0 = pbuffer.data(idx_op_geom_200_gs + 0 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 14 * ket_comps + j);

            #pragma omp simd aligned(to_xx_0, to_xx_0_xxxx_0, to_xx_0_xxxy_0, to_xx_0_xxxz_0, to_xx_0_xxyy_0, to_xx_0_xxyz_0, to_xx_0_xxzz_0, to_xx_0_xyyy_0, to_xx_0_xyyz_0, to_xx_0_xyzz_0, to_xx_0_xzzz_0, to_xx_0_yyyy_0, to_xx_0_yyyz_0, to_xx_0_yyzz_0, to_xx_0_yzzz_0, to_xx_0_zzzz_0, to_xxxx_0, to_xxxxxx_0, to_xxxxxy_0, to_xxxxxz_0, to_xxxxyy_0, to_xxxxyz_0, to_xxxxzz_0, to_xxxy_0, to_xxxyyy_0, to_xxxyyz_0, to_xxxyzz_0, to_xxxz_0, to_xxxzzz_0, to_xxyy_0, to_xxyyyy_0, to_xxyyyz_0, to_xxyyzz_0, to_xxyz_0, to_xxyzzz_0, to_xxzz_0, to_xxzzzz_0, to_xy_0, to_xyyy_0, to_xyyz_0, to_xyzz_0, to_xz_0, to_xzzz_0, to_yy_0, to_yyyy_0, to_yyyz_0, to_yyzz_0, to_yz_0, to_yzzz_0, to_zz_0, to_zzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                const double tbe_0 = a_exp;

                to_xx_0_xxxx_0[k] = 12.0 * to_xx_0[k] - 18.0 * to_xxxx_0[k] * tbe_0 + 4.0 * to_xxxxxx_0[k] * tbe_0 * tbe_0;

                to_xx_0_xxxy_0[k] = 6.0 * to_xy_0[k] - 14.0 * to_xxxy_0[k] * tbe_0 + 4.0 * to_xxxxxy_0[k] * tbe_0 * tbe_0;

                to_xx_0_xxxz_0[k] = 6.0 * to_xz_0[k] - 14.0 * to_xxxz_0[k] * tbe_0 + 4.0 * to_xxxxxz_0[k] * tbe_0 * tbe_0;

                to_xx_0_xxyy_0[k] = 2.0 * to_yy_0[k] - 10.0 * to_xxyy_0[k] * tbe_0 + 4.0 * to_xxxxyy_0[k] * tbe_0 * tbe_0;

                to_xx_0_xxyz_0[k] = 2.0 * to_yz_0[k] - 10.0 * to_xxyz_0[k] * tbe_0 + 4.0 * to_xxxxyz_0[k] * tbe_0 * tbe_0;

                to_xx_0_xxzz_0[k] = 2.0 * to_zz_0[k] - 10.0 * to_xxzz_0[k] * tbe_0 + 4.0 * to_xxxxzz_0[k] * tbe_0 * tbe_0;

                to_xx_0_xyyy_0[k] = -6.0 * to_xyyy_0[k] * tbe_0 + 4.0 * to_xxxyyy_0[k] * tbe_0 * tbe_0;

                to_xx_0_xyyz_0[k] = -6.0 * to_xyyz_0[k] * tbe_0 + 4.0 * to_xxxyyz_0[k] * tbe_0 * tbe_0;

                to_xx_0_xyzz_0[k] = -6.0 * to_xyzz_0[k] * tbe_0 + 4.0 * to_xxxyzz_0[k] * tbe_0 * tbe_0;

                to_xx_0_xzzz_0[k] = -6.0 * to_xzzz_0[k] * tbe_0 + 4.0 * to_xxxzzz_0[k] * tbe_0 * tbe_0;

                to_xx_0_yyyy_0[k] = -2.0 * to_yyyy_0[k] * tbe_0 + 4.0 * to_xxyyyy_0[k] * tbe_0 * tbe_0;

                to_xx_0_yyyz_0[k] = -2.0 * to_yyyz_0[k] * tbe_0 + 4.0 * to_xxyyyz_0[k] * tbe_0 * tbe_0;

                to_xx_0_yyzz_0[k] = -2.0 * to_yyzz_0[k] * tbe_0 + 4.0 * to_xxyyzz_0[k] * tbe_0 * tbe_0;

                to_xx_0_yzzz_0[k] = -2.0 * to_yzzz_0[k] * tbe_0 + 4.0 * to_xxyzzz_0[k] * tbe_0 * tbe_0;

                to_xx_0_zzzz_0[k] = -2.0 * to_zzzz_0[k] * tbe_0 + 4.0 * to_xxzzzz_0[k] * tbe_0 * tbe_0;
            }

            // Set up 15-30 components of targeted buffer : GS

            auto to_xy_0_xxxx_0 = pbuffer.data(idx_op_geom_200_gs + 1 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 0 * ket_comps + j);

            auto to_xy_0_xxxy_0 = pbuffer.data(idx_op_geom_200_gs + 1 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 1 * ket_comps + j);

            auto to_xy_0_xxxz_0 = pbuffer.data(idx_op_geom_200_gs + 1 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 2 * ket_comps + j);

            auto to_xy_0_xxyy_0 = pbuffer.data(idx_op_geom_200_gs + 1 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 3 * ket_comps + j);

            auto to_xy_0_xxyz_0 = pbuffer.data(idx_op_geom_200_gs + 1 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 4 * ket_comps + j);

            auto to_xy_0_xxzz_0 = pbuffer.data(idx_op_geom_200_gs + 1 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 5 * ket_comps + j);

            auto to_xy_0_xyyy_0 = pbuffer.data(idx_op_geom_200_gs + 1 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 6 * ket_comps + j);

            auto to_xy_0_xyyz_0 = pbuffer.data(idx_op_geom_200_gs + 1 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 7 * ket_comps + j);

            auto to_xy_0_xyzz_0 = pbuffer.data(idx_op_geom_200_gs + 1 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 8 * ket_comps + j);

            auto to_xy_0_xzzz_0 = pbuffer.data(idx_op_geom_200_gs + 1 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 9 * ket_comps + j);

            auto to_xy_0_yyyy_0 = pbuffer.data(idx_op_geom_200_gs + 1 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 10 * ket_comps + j);

            auto to_xy_0_yyyz_0 = pbuffer.data(idx_op_geom_200_gs + 1 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 11 * ket_comps + j);

            auto to_xy_0_yyzz_0 = pbuffer.data(idx_op_geom_200_gs + 1 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 12 * ket_comps + j);

            auto to_xy_0_yzzz_0 = pbuffer.data(idx_op_geom_200_gs + 1 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 13 * ket_comps + j);

            auto to_xy_0_zzzz_0 = pbuffer.data(idx_op_geom_200_gs + 1 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 14 * ket_comps + j);

            #pragma omp simd aligned(to_xx_0, to_xxxx_0, to_xxxxxy_0, to_xxxxyy_0, to_xxxxyz_0, to_xxxy_0, to_xxxyyy_0, to_xxxyyz_0, to_xxxyzz_0, to_xxxz_0, to_xxyy_0, to_xxyyyy_0, to_xxyyyz_0, to_xxyyzz_0, to_xxyz_0, to_xxyzzz_0, to_xxzz_0, to_xy_0, to_xy_0_xxxx_0, to_xy_0_xxxy_0, to_xy_0_xxxz_0, to_xy_0_xxyy_0, to_xy_0_xxyz_0, to_xy_0_xxzz_0, to_xy_0_xyyy_0, to_xy_0_xyyz_0, to_xy_0_xyzz_0, to_xy_0_xzzz_0, to_xy_0_yyyy_0, to_xy_0_yyyz_0, to_xy_0_yyzz_0, to_xy_0_yzzz_0, to_xy_0_zzzz_0, to_xyyy_0, to_xyyyyy_0, to_xyyyyz_0, to_xyyyzz_0, to_xyyz_0, to_xyyzzz_0, to_xyzz_0, to_xyzzzz_0, to_xz_0, to_xzzz_0, to_yy_0, to_yyyy_0, to_yyyz_0, to_yyzz_0, to_yz_0, to_yzzz_0, to_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                const double tbe_0 = a_exp;

                to_xy_0_xxxx_0[k] = -8.0 * to_xxxy_0[k] * tbe_0 + 4.0 * to_xxxxxy_0[k] * tbe_0 * tbe_0;

                to_xy_0_xxxy_0[k] = 3.0 * to_xx_0[k] - 6.0 * to_xxyy_0[k] * tbe_0 - 2.0 * to_xxxx_0[k] * tbe_0 + 4.0 * to_xxxxyy_0[k] * tbe_0 * tbe_0;

                to_xy_0_xxxz_0[k] = -6.0 * to_xxyz_0[k] * tbe_0 + 4.0 * to_xxxxyz_0[k] * tbe_0 * tbe_0;

                to_xy_0_xxyy_0[k] = 4.0 * to_xy_0[k] - 4.0 * to_xyyy_0[k] * tbe_0 - 4.0 * to_xxxy_0[k] * tbe_0 + 4.0 * to_xxxyyy_0[k] * tbe_0 * tbe_0;

                to_xy_0_xxyz_0[k] = 2.0 * to_xz_0[k] - 4.0 * to_xyyz_0[k] * tbe_0 - 2.0 * to_xxxz_0[k] * tbe_0 + 4.0 * to_xxxyyz_0[k] * tbe_0 * tbe_0;

                to_xy_0_xxzz_0[k] = -4.0 * to_xyzz_0[k] * tbe_0 + 4.0 * to_xxxyzz_0[k] * tbe_0 * tbe_0;

                to_xy_0_xyyy_0[k] = 3.0 * to_yy_0[k] - 2.0 * to_yyyy_0[k] * tbe_0 - 6.0 * to_xxyy_0[k] * tbe_0 + 4.0 * to_xxyyyy_0[k] * tbe_0 * tbe_0;

                to_xy_0_xyyz_0[k] = 2.0 * to_yz_0[k] - 2.0 * to_yyyz_0[k] * tbe_0 - 4.0 * to_xxyz_0[k] * tbe_0 + 4.0 * to_xxyyyz_0[k] * tbe_0 * tbe_0;

                to_xy_0_xyzz_0[k] = to_zz_0[k] - 2.0 * to_yyzz_0[k] * tbe_0 - 2.0 * to_xxzz_0[k] * tbe_0 + 4.0 * to_xxyyzz_0[k] * tbe_0 * tbe_0;

                to_xy_0_xzzz_0[k] = -2.0 * to_yzzz_0[k] * tbe_0 + 4.0 * to_xxyzzz_0[k] * tbe_0 * tbe_0;

                to_xy_0_yyyy_0[k] = -8.0 * to_xyyy_0[k] * tbe_0 + 4.0 * to_xyyyyy_0[k] * tbe_0 * tbe_0;

                to_xy_0_yyyz_0[k] = -6.0 * to_xyyz_0[k] * tbe_0 + 4.0 * to_xyyyyz_0[k] * tbe_0 * tbe_0;

                to_xy_0_yyzz_0[k] = -4.0 * to_xyzz_0[k] * tbe_0 + 4.0 * to_xyyyzz_0[k] * tbe_0 * tbe_0;

                to_xy_0_yzzz_0[k] = -2.0 * to_xzzz_0[k] * tbe_0 + 4.0 * to_xyyzzz_0[k] * tbe_0 * tbe_0;

                to_xy_0_zzzz_0[k] = 4.0 * to_xyzzzz_0[k] * tbe_0 * tbe_0;
            }

            // Set up 30-45 components of targeted buffer : GS

            auto to_xz_0_xxxx_0 = pbuffer.data(idx_op_geom_200_gs + 2 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 0 * ket_comps + j);

            auto to_xz_0_xxxy_0 = pbuffer.data(idx_op_geom_200_gs + 2 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 1 * ket_comps + j);

            auto to_xz_0_xxxz_0 = pbuffer.data(idx_op_geom_200_gs + 2 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 2 * ket_comps + j);

            auto to_xz_0_xxyy_0 = pbuffer.data(idx_op_geom_200_gs + 2 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 3 * ket_comps + j);

            auto to_xz_0_xxyz_0 = pbuffer.data(idx_op_geom_200_gs + 2 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 4 * ket_comps + j);

            auto to_xz_0_xxzz_0 = pbuffer.data(idx_op_geom_200_gs + 2 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 5 * ket_comps + j);

            auto to_xz_0_xyyy_0 = pbuffer.data(idx_op_geom_200_gs + 2 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 6 * ket_comps + j);

            auto to_xz_0_xyyz_0 = pbuffer.data(idx_op_geom_200_gs + 2 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 7 * ket_comps + j);

            auto to_xz_0_xyzz_0 = pbuffer.data(idx_op_geom_200_gs + 2 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 8 * ket_comps + j);

            auto to_xz_0_xzzz_0 = pbuffer.data(idx_op_geom_200_gs + 2 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 9 * ket_comps + j);

            auto to_xz_0_yyyy_0 = pbuffer.data(idx_op_geom_200_gs + 2 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 10 * ket_comps + j);

            auto to_xz_0_yyyz_0 = pbuffer.data(idx_op_geom_200_gs + 2 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 11 * ket_comps + j);

            auto to_xz_0_yyzz_0 = pbuffer.data(idx_op_geom_200_gs + 2 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 12 * ket_comps + j);

            auto to_xz_0_yzzz_0 = pbuffer.data(idx_op_geom_200_gs + 2 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 13 * ket_comps + j);

            auto to_xz_0_zzzz_0 = pbuffer.data(idx_op_geom_200_gs + 2 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 14 * ket_comps + j);

            #pragma omp simd aligned(to_xx_0, to_xxxx_0, to_xxxxxz_0, to_xxxxyz_0, to_xxxxzz_0, to_xxxy_0, to_xxxyyz_0, to_xxxyzz_0, to_xxxz_0, to_xxxzzz_0, to_xxyy_0, to_xxyyyz_0, to_xxyyzz_0, to_xxyz_0, to_xxyzzz_0, to_xxzz_0, to_xxzzzz_0, to_xy_0, to_xyyy_0, to_xyyyyz_0, to_xyyyzz_0, to_xyyz_0, to_xyyzzz_0, to_xyzz_0, to_xyzzzz_0, to_xz_0, to_xz_0_xxxx_0, to_xz_0_xxxy_0, to_xz_0_xxxz_0, to_xz_0_xxyy_0, to_xz_0_xxyz_0, to_xz_0_xxzz_0, to_xz_0_xyyy_0, to_xz_0_xyyz_0, to_xz_0_xyzz_0, to_xz_0_xzzz_0, to_xz_0_yyyy_0, to_xz_0_yyyz_0, to_xz_0_yyzz_0, to_xz_0_yzzz_0, to_xz_0_zzzz_0, to_xzzz_0, to_xzzzzz_0, to_yy_0, to_yyyz_0, to_yyzz_0, to_yz_0, to_yzzz_0, to_zz_0, to_zzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                const double tbe_0 = a_exp;

                to_xz_0_xxxx_0[k] = -8.0 * to_xxxz_0[k] * tbe_0 + 4.0 * to_xxxxxz_0[k] * tbe_0 * tbe_0;

                to_xz_0_xxxy_0[k] = -6.0 * to_xxyz_0[k] * tbe_0 + 4.0 * to_xxxxyz_0[k] * tbe_0 * tbe_0;

                to_xz_0_xxxz_0[k] = 3.0 * to_xx_0[k] - 6.0 * to_xxzz_0[k] * tbe_0 - 2.0 * to_xxxx_0[k] * tbe_0 + 4.0 * to_xxxxzz_0[k] * tbe_0 * tbe_0;

                to_xz_0_xxyy_0[k] = -4.0 * to_xyyz_0[k] * tbe_0 + 4.0 * to_xxxyyz_0[k] * tbe_0 * tbe_0;

                to_xz_0_xxyz_0[k] = 2.0 * to_xy_0[k] - 4.0 * to_xyzz_0[k] * tbe_0 - 2.0 * to_xxxy_0[k] * tbe_0 + 4.0 * to_xxxyzz_0[k] * tbe_0 * tbe_0;

                to_xz_0_xxzz_0[k] = 4.0 * to_xz_0[k] - 4.0 * to_xzzz_0[k] * tbe_0 - 4.0 * to_xxxz_0[k] * tbe_0 + 4.0 * to_xxxzzz_0[k] * tbe_0 * tbe_0;

                to_xz_0_xyyy_0[k] = -2.0 * to_yyyz_0[k] * tbe_0 + 4.0 * to_xxyyyz_0[k] * tbe_0 * tbe_0;

                to_xz_0_xyyz_0[k] = to_yy_0[k] - 2.0 * to_yyzz_0[k] * tbe_0 - 2.0 * to_xxyy_0[k] * tbe_0 + 4.0 * to_xxyyzz_0[k] * tbe_0 * tbe_0;

                to_xz_0_xyzz_0[k] = 2.0 * to_yz_0[k] - 2.0 * to_yzzz_0[k] * tbe_0 - 4.0 * to_xxyz_0[k] * tbe_0 + 4.0 * to_xxyzzz_0[k] * tbe_0 * tbe_0;

                to_xz_0_xzzz_0[k] = 3.0 * to_zz_0[k] - 2.0 * to_zzzz_0[k] * tbe_0 - 6.0 * to_xxzz_0[k] * tbe_0 + 4.0 * to_xxzzzz_0[k] * tbe_0 * tbe_0;

                to_xz_0_yyyy_0[k] = 4.0 * to_xyyyyz_0[k] * tbe_0 * tbe_0;

                to_xz_0_yyyz_0[k] = -2.0 * to_xyyy_0[k] * tbe_0 + 4.0 * to_xyyyzz_0[k] * tbe_0 * tbe_0;

                to_xz_0_yyzz_0[k] = -4.0 * to_xyyz_0[k] * tbe_0 + 4.0 * to_xyyzzz_0[k] * tbe_0 * tbe_0;

                to_xz_0_yzzz_0[k] = -6.0 * to_xyzz_0[k] * tbe_0 + 4.0 * to_xyzzzz_0[k] * tbe_0 * tbe_0;

                to_xz_0_zzzz_0[k] = -8.0 * to_xzzz_0[k] * tbe_0 + 4.0 * to_xzzzzz_0[k] * tbe_0 * tbe_0;
            }

            // Set up 45-60 components of targeted buffer : GS

            auto to_yy_0_xxxx_0 = pbuffer.data(idx_op_geom_200_gs + 3 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 0 * ket_comps + j);

            auto to_yy_0_xxxy_0 = pbuffer.data(idx_op_geom_200_gs + 3 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 1 * ket_comps + j);

            auto to_yy_0_xxxz_0 = pbuffer.data(idx_op_geom_200_gs + 3 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 2 * ket_comps + j);

            auto to_yy_0_xxyy_0 = pbuffer.data(idx_op_geom_200_gs + 3 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 3 * ket_comps + j);

            auto to_yy_0_xxyz_0 = pbuffer.data(idx_op_geom_200_gs + 3 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 4 * ket_comps + j);

            auto to_yy_0_xxzz_0 = pbuffer.data(idx_op_geom_200_gs + 3 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 5 * ket_comps + j);

            auto to_yy_0_xyyy_0 = pbuffer.data(idx_op_geom_200_gs + 3 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 6 * ket_comps + j);

            auto to_yy_0_xyyz_0 = pbuffer.data(idx_op_geom_200_gs + 3 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 7 * ket_comps + j);

            auto to_yy_0_xyzz_0 = pbuffer.data(idx_op_geom_200_gs + 3 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 8 * ket_comps + j);

            auto to_yy_0_xzzz_0 = pbuffer.data(idx_op_geom_200_gs + 3 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 9 * ket_comps + j);

            auto to_yy_0_yyyy_0 = pbuffer.data(idx_op_geom_200_gs + 3 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 10 * ket_comps + j);

            auto to_yy_0_yyyz_0 = pbuffer.data(idx_op_geom_200_gs + 3 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 11 * ket_comps + j);

            auto to_yy_0_yyzz_0 = pbuffer.data(idx_op_geom_200_gs + 3 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 12 * ket_comps + j);

            auto to_yy_0_yzzz_0 = pbuffer.data(idx_op_geom_200_gs + 3 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 13 * ket_comps + j);

            auto to_yy_0_zzzz_0 = pbuffer.data(idx_op_geom_200_gs + 3 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 14 * ket_comps + j);

            #pragma omp simd aligned(to_xx_0, to_xxxx_0, to_xxxxyy_0, to_xxxy_0, to_xxxyyy_0, to_xxxyyz_0, to_xxxz_0, to_xxyy_0, to_xxyyyy_0, to_xxyyyz_0, to_xxyyzz_0, to_xxyz_0, to_xxzz_0, to_xy_0, to_xyyy_0, to_xyyyyy_0, to_xyyyyz_0, to_xyyyzz_0, to_xyyz_0, to_xyyzzz_0, to_xyzz_0, to_xz_0, to_xzzz_0, to_yy_0, to_yy_0_xxxx_0, to_yy_0_xxxy_0, to_yy_0_xxxz_0, to_yy_0_xxyy_0, to_yy_0_xxyz_0, to_yy_0_xxzz_0, to_yy_0_xyyy_0, to_yy_0_xyyz_0, to_yy_0_xyzz_0, to_yy_0_xzzz_0, to_yy_0_yyyy_0, to_yy_0_yyyz_0, to_yy_0_yyzz_0, to_yy_0_yzzz_0, to_yy_0_zzzz_0, to_yyyy_0, to_yyyyyy_0, to_yyyyyz_0, to_yyyyzz_0, to_yyyz_0, to_yyyzzz_0, to_yyzz_0, to_yyzzzz_0, to_yz_0, to_yzzz_0, to_zz_0, to_zzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                const double tbe_0 = a_exp;

                to_yy_0_xxxx_0[k] = -2.0 * to_xxxx_0[k] * tbe_0 + 4.0 * to_xxxxyy_0[k] * tbe_0 * tbe_0;

                to_yy_0_xxxy_0[k] = -6.0 * to_xxxy_0[k] * tbe_0 + 4.0 * to_xxxyyy_0[k] * tbe_0 * tbe_0;

                to_yy_0_xxxz_0[k] = -2.0 * to_xxxz_0[k] * tbe_0 + 4.0 * to_xxxyyz_0[k] * tbe_0 * tbe_0;

                to_yy_0_xxyy_0[k] = 2.0 * to_xx_0[k] - 10.0 * to_xxyy_0[k] * tbe_0 + 4.0 * to_xxyyyy_0[k] * tbe_0 * tbe_0;

                to_yy_0_xxyz_0[k] = -6.0 * to_xxyz_0[k] * tbe_0 + 4.0 * to_xxyyyz_0[k] * tbe_0 * tbe_0;

                to_yy_0_xxzz_0[k] = -2.0 * to_xxzz_0[k] * tbe_0 + 4.0 * to_xxyyzz_0[k] * tbe_0 * tbe_0;

                to_yy_0_xyyy_0[k] = 6.0 * to_xy_0[k] - 14.0 * to_xyyy_0[k] * tbe_0 + 4.0 * to_xyyyyy_0[k] * tbe_0 * tbe_0;

                to_yy_0_xyyz_0[k] = 2.0 * to_xz_0[k] - 10.0 * to_xyyz_0[k] * tbe_0 + 4.0 * to_xyyyyz_0[k] * tbe_0 * tbe_0;

                to_yy_0_xyzz_0[k] = -6.0 * to_xyzz_0[k] * tbe_0 + 4.0 * to_xyyyzz_0[k] * tbe_0 * tbe_0;

                to_yy_0_xzzz_0[k] = -2.0 * to_xzzz_0[k] * tbe_0 + 4.0 * to_xyyzzz_0[k] * tbe_0 * tbe_0;

                to_yy_0_yyyy_0[k] = 12.0 * to_yy_0[k] - 18.0 * to_yyyy_0[k] * tbe_0 + 4.0 * to_yyyyyy_0[k] * tbe_0 * tbe_0;

                to_yy_0_yyyz_0[k] = 6.0 * to_yz_0[k] - 14.0 * to_yyyz_0[k] * tbe_0 + 4.0 * to_yyyyyz_0[k] * tbe_0 * tbe_0;

                to_yy_0_yyzz_0[k] = 2.0 * to_zz_0[k] - 10.0 * to_yyzz_0[k] * tbe_0 + 4.0 * to_yyyyzz_0[k] * tbe_0 * tbe_0;

                to_yy_0_yzzz_0[k] = -6.0 * to_yzzz_0[k] * tbe_0 + 4.0 * to_yyyzzz_0[k] * tbe_0 * tbe_0;

                to_yy_0_zzzz_0[k] = -2.0 * to_zzzz_0[k] * tbe_0 + 4.0 * to_yyzzzz_0[k] * tbe_0 * tbe_0;
            }

            // Set up 60-75 components of targeted buffer : GS

            auto to_yz_0_xxxx_0 = pbuffer.data(idx_op_geom_200_gs + 4 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 0 * ket_comps + j);

            auto to_yz_0_xxxy_0 = pbuffer.data(idx_op_geom_200_gs + 4 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 1 * ket_comps + j);

            auto to_yz_0_xxxz_0 = pbuffer.data(idx_op_geom_200_gs + 4 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 2 * ket_comps + j);

            auto to_yz_0_xxyy_0 = pbuffer.data(idx_op_geom_200_gs + 4 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 3 * ket_comps + j);

            auto to_yz_0_xxyz_0 = pbuffer.data(idx_op_geom_200_gs + 4 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 4 * ket_comps + j);

            auto to_yz_0_xxzz_0 = pbuffer.data(idx_op_geom_200_gs + 4 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 5 * ket_comps + j);

            auto to_yz_0_xyyy_0 = pbuffer.data(idx_op_geom_200_gs + 4 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 6 * ket_comps + j);

            auto to_yz_0_xyyz_0 = pbuffer.data(idx_op_geom_200_gs + 4 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 7 * ket_comps + j);

            auto to_yz_0_xyzz_0 = pbuffer.data(idx_op_geom_200_gs + 4 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 8 * ket_comps + j);

            auto to_yz_0_xzzz_0 = pbuffer.data(idx_op_geom_200_gs + 4 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 9 * ket_comps + j);

            auto to_yz_0_yyyy_0 = pbuffer.data(idx_op_geom_200_gs + 4 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 10 * ket_comps + j);

            auto to_yz_0_yyyz_0 = pbuffer.data(idx_op_geom_200_gs + 4 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 11 * ket_comps + j);

            auto to_yz_0_yyzz_0 = pbuffer.data(idx_op_geom_200_gs + 4 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 12 * ket_comps + j);

            auto to_yz_0_yzzz_0 = pbuffer.data(idx_op_geom_200_gs + 4 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 13 * ket_comps + j);

            auto to_yz_0_zzzz_0 = pbuffer.data(idx_op_geom_200_gs + 4 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 14 * ket_comps + j);

            #pragma omp simd aligned(to_xx_0, to_xxxxyz_0, to_xxxy_0, to_xxxyyz_0, to_xxxyzz_0, to_xxxz_0, to_xxyy_0, to_xxyyyz_0, to_xxyyzz_0, to_xxyz_0, to_xxyzzz_0, to_xxzz_0, to_xy_0, to_xyyy_0, to_xyyyyz_0, to_xyyyzz_0, to_xyyz_0, to_xyyzzz_0, to_xyzz_0, to_xyzzzz_0, to_xz_0, to_xzzz_0, to_yy_0, to_yyyy_0, to_yyyyyz_0, to_yyyyzz_0, to_yyyz_0, to_yyyzzz_0, to_yyzz_0, to_yyzzzz_0, to_yz_0, to_yz_0_xxxx_0, to_yz_0_xxxy_0, to_yz_0_xxxz_0, to_yz_0_xxyy_0, to_yz_0_xxyz_0, to_yz_0_xxzz_0, to_yz_0_xyyy_0, to_yz_0_xyyz_0, to_yz_0_xyzz_0, to_yz_0_xzzz_0, to_yz_0_yyyy_0, to_yz_0_yyyz_0, to_yz_0_yyzz_0, to_yz_0_yzzz_0, to_yz_0_zzzz_0, to_yzzz_0, to_yzzzzz_0, to_zz_0, to_zzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                const double tbe_0 = a_exp;

                to_yz_0_xxxx_0[k] = 4.0 * to_xxxxyz_0[k] * tbe_0 * tbe_0;

                to_yz_0_xxxy_0[k] = -2.0 * to_xxxz_0[k] * tbe_0 + 4.0 * to_xxxyyz_0[k] * tbe_0 * tbe_0;

                to_yz_0_xxxz_0[k] = -2.0 * to_xxxy_0[k] * tbe_0 + 4.0 * to_xxxyzz_0[k] * tbe_0 * tbe_0;

                to_yz_0_xxyy_0[k] = -4.0 * to_xxyz_0[k] * tbe_0 + 4.0 * to_xxyyyz_0[k] * tbe_0 * tbe_0;

                to_yz_0_xxyz_0[k] = to_xx_0[k] - 2.0 * to_xxzz_0[k] * tbe_0 - 2.0 * to_xxyy_0[k] * tbe_0 + 4.0 * to_xxyyzz_0[k] * tbe_0 * tbe_0;

                to_yz_0_xxzz_0[k] = -4.0 * to_xxyz_0[k] * tbe_0 + 4.0 * to_xxyzzz_0[k] * tbe_0 * tbe_0;

                to_yz_0_xyyy_0[k] = -6.0 * to_xyyz_0[k] * tbe_0 + 4.0 * to_xyyyyz_0[k] * tbe_0 * tbe_0;

                to_yz_0_xyyz_0[k] = 2.0 * to_xy_0[k] - 4.0 * to_xyzz_0[k] * tbe_0 - 2.0 * to_xyyy_0[k] * tbe_0 + 4.0 * to_xyyyzz_0[k] * tbe_0 * tbe_0;

                to_yz_0_xyzz_0[k] = 2.0 * to_xz_0[k] - 2.0 * to_xzzz_0[k] * tbe_0 - 4.0 * to_xyyz_0[k] * tbe_0 + 4.0 * to_xyyzzz_0[k] * tbe_0 * tbe_0;

                to_yz_0_xzzz_0[k] = -6.0 * to_xyzz_0[k] * tbe_0 + 4.0 * to_xyzzzz_0[k] * tbe_0 * tbe_0;

                to_yz_0_yyyy_0[k] = -8.0 * to_yyyz_0[k] * tbe_0 + 4.0 * to_yyyyyz_0[k] * tbe_0 * tbe_0;

                to_yz_0_yyyz_0[k] = 3.0 * to_yy_0[k] - 6.0 * to_yyzz_0[k] * tbe_0 - 2.0 * to_yyyy_0[k] * tbe_0 + 4.0 * to_yyyyzz_0[k] * tbe_0 * tbe_0;

                to_yz_0_yyzz_0[k] = 4.0 * to_yz_0[k] - 4.0 * to_yzzz_0[k] * tbe_0 - 4.0 * to_yyyz_0[k] * tbe_0 + 4.0 * to_yyyzzz_0[k] * tbe_0 * tbe_0;

                to_yz_0_yzzz_0[k] = 3.0 * to_zz_0[k] - 2.0 * to_zzzz_0[k] * tbe_0 - 6.0 * to_yyzz_0[k] * tbe_0 + 4.0 * to_yyzzzz_0[k] * tbe_0 * tbe_0;

                to_yz_0_zzzz_0[k] = -8.0 * to_yzzz_0[k] * tbe_0 + 4.0 * to_yzzzzz_0[k] * tbe_0 * tbe_0;
            }

            // Set up 75-90 components of targeted buffer : GS

            auto to_zz_0_xxxx_0 = pbuffer.data(idx_op_geom_200_gs + 5 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 0 * ket_comps + j);

            auto to_zz_0_xxxy_0 = pbuffer.data(idx_op_geom_200_gs + 5 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 1 * ket_comps + j);

            auto to_zz_0_xxxz_0 = pbuffer.data(idx_op_geom_200_gs + 5 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 2 * ket_comps + j);

            auto to_zz_0_xxyy_0 = pbuffer.data(idx_op_geom_200_gs + 5 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 3 * ket_comps + j);

            auto to_zz_0_xxyz_0 = pbuffer.data(idx_op_geom_200_gs + 5 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 4 * ket_comps + j);

            auto to_zz_0_xxzz_0 = pbuffer.data(idx_op_geom_200_gs + 5 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 5 * ket_comps + j);

            auto to_zz_0_xyyy_0 = pbuffer.data(idx_op_geom_200_gs + 5 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 6 * ket_comps + j);

            auto to_zz_0_xyyz_0 = pbuffer.data(idx_op_geom_200_gs + 5 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 7 * ket_comps + j);

            auto to_zz_0_xyzz_0 = pbuffer.data(idx_op_geom_200_gs + 5 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 8 * ket_comps + j);

            auto to_zz_0_xzzz_0 = pbuffer.data(idx_op_geom_200_gs + 5 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 9 * ket_comps + j);

            auto to_zz_0_yyyy_0 = pbuffer.data(idx_op_geom_200_gs + 5 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 10 * ket_comps + j);

            auto to_zz_0_yyyz_0 = pbuffer.data(idx_op_geom_200_gs + 5 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 11 * ket_comps + j);

            auto to_zz_0_yyzz_0 = pbuffer.data(idx_op_geom_200_gs + 5 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 12 * ket_comps + j);

            auto to_zz_0_yzzz_0 = pbuffer.data(idx_op_geom_200_gs + 5 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 13 * ket_comps + j);

            auto to_zz_0_zzzz_0 = pbuffer.data(idx_op_geom_200_gs + 5 * op_comps * 15 * ket_comps + i * 15 * ket_comps + 14 * ket_comps + j);

            #pragma omp simd aligned(to_xx_0, to_xxxx_0, to_xxxxzz_0, to_xxxy_0, to_xxxyzz_0, to_xxxz_0, to_xxxzzz_0, to_xxyy_0, to_xxyyzz_0, to_xxyz_0, to_xxyzzz_0, to_xxzz_0, to_xxzzzz_0, to_xy_0, to_xyyy_0, to_xyyyzz_0, to_xyyz_0, to_xyyzzz_0, to_xyzz_0, to_xyzzzz_0, to_xz_0, to_xzzz_0, to_xzzzzz_0, to_yy_0, to_yyyy_0, to_yyyyzz_0, to_yyyz_0, to_yyyzzz_0, to_yyzz_0, to_yyzzzz_0, to_yz_0, to_yzzz_0, to_yzzzzz_0, to_zz_0, to_zz_0_xxxx_0, to_zz_0_xxxy_0, to_zz_0_xxxz_0, to_zz_0_xxyy_0, to_zz_0_xxyz_0, to_zz_0_xxzz_0, to_zz_0_xyyy_0, to_zz_0_xyyz_0, to_zz_0_xyzz_0, to_zz_0_xzzz_0, to_zz_0_yyyy_0, to_zz_0_yyyz_0, to_zz_0_yyzz_0, to_zz_0_yzzz_0, to_zz_0_zzzz_0, to_zzzz_0, to_zzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                const double tbe_0 = a_exp;

                to_zz_0_xxxx_0[k] = -2.0 * to_xxxx_0[k] * tbe_0 + 4.0 * to_xxxxzz_0[k] * tbe_0 * tbe_0;

                to_zz_0_xxxy_0[k] = -2.0 * to_xxxy_0[k] * tbe_0 + 4.0 * to_xxxyzz_0[k] * tbe_0 * tbe_0;

                to_zz_0_xxxz_0[k] = -6.0 * to_xxxz_0[k] * tbe_0 + 4.0 * to_xxxzzz_0[k] * tbe_0 * tbe_0;

                to_zz_0_xxyy_0[k] = -2.0 * to_xxyy_0[k] * tbe_0 + 4.0 * to_xxyyzz_0[k] * tbe_0 * tbe_0;

                to_zz_0_xxyz_0[k] = -6.0 * to_xxyz_0[k] * tbe_0 + 4.0 * to_xxyzzz_0[k] * tbe_0 * tbe_0;

                to_zz_0_xxzz_0[k] = 2.0 * to_xx_0[k] - 10.0 * to_xxzz_0[k] * tbe_0 + 4.0 * to_xxzzzz_0[k] * tbe_0 * tbe_0;

                to_zz_0_xyyy_0[k] = -2.0 * to_xyyy_0[k] * tbe_0 + 4.0 * to_xyyyzz_0[k] * tbe_0 * tbe_0;

                to_zz_0_xyyz_0[k] = -6.0 * to_xyyz_0[k] * tbe_0 + 4.0 * to_xyyzzz_0[k] * tbe_0 * tbe_0;

                to_zz_0_xyzz_0[k] = 2.0 * to_xy_0[k] - 10.0 * to_xyzz_0[k] * tbe_0 + 4.0 * to_xyzzzz_0[k] * tbe_0 * tbe_0;

                to_zz_0_xzzz_0[k] = 6.0 * to_xz_0[k] - 14.0 * to_xzzz_0[k] * tbe_0 + 4.0 * to_xzzzzz_0[k] * tbe_0 * tbe_0;

                to_zz_0_yyyy_0[k] = -2.0 * to_yyyy_0[k] * tbe_0 + 4.0 * to_yyyyzz_0[k] * tbe_0 * tbe_0;

                to_zz_0_yyyz_0[k] = -6.0 * to_yyyz_0[k] * tbe_0 + 4.0 * to_yyyzzz_0[k] * tbe_0 * tbe_0;

                to_zz_0_yyzz_0[k] = 2.0 * to_yy_0[k] - 10.0 * to_yyzz_0[k] * tbe_0 + 4.0 * to_yyzzzz_0[k] * tbe_0 * tbe_0;

                to_zz_0_yzzz_0[k] = 6.0 * to_yz_0[k] - 14.0 * to_yzzz_0[k] * tbe_0 + 4.0 * to_yzzzzz_0[k] * tbe_0 * tbe_0;

                to_zz_0_zzzz_0[k] = 12.0 * to_zz_0[k] - 18.0 * to_zzzz_0[k] * tbe_0 + 4.0 * to_zzzzzz_0[k] * tbe_0 * tbe_0;
            }

        }
    }

}

} // t2cgeom namespace

