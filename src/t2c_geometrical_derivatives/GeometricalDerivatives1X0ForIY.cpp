#include "GeometricalDerivatives1X0ForIY.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_10_ix(CSimdArray<double>& prim_buffer,
                        const int idx_op_geom_100_is,
                        const int idx_op_hs,
                        const int idx_op_ks,
                        const int op_comps,
                        const int ket_comps,
                        const double a_exp) -> void
{
    const auto nelems = prim_buffer.number_of_active_elements();

    for (size_t i = 0; i < op_comps; i++)
    {
        for (size_t j = 0; j < ket_comps; j++)
        {
            // Set up components of auxiliary buffer : HS

            auto to_xxxxx_0 = prim_buffer.data(idx_op_hs + i * 21 * ket_comps + 0 * ket_comps + j);

            auto to_xxxxy_0 = prim_buffer.data(idx_op_hs + i * 21 * ket_comps + 1 * ket_comps + j);

            auto to_xxxxz_0 = prim_buffer.data(idx_op_hs + i * 21 * ket_comps + 2 * ket_comps + j);

            auto to_xxxyy_0 = prim_buffer.data(idx_op_hs + i * 21 * ket_comps + 3 * ket_comps + j);

            auto to_xxxyz_0 = prim_buffer.data(idx_op_hs + i * 21 * ket_comps + 4 * ket_comps + j);

            auto to_xxxzz_0 = prim_buffer.data(idx_op_hs + i * 21 * ket_comps + 5 * ket_comps + j);

            auto to_xxyyy_0 = prim_buffer.data(idx_op_hs + i * 21 * ket_comps + 6 * ket_comps + j);

            auto to_xxyyz_0 = prim_buffer.data(idx_op_hs + i * 21 * ket_comps + 7 * ket_comps + j);

            auto to_xxyzz_0 = prim_buffer.data(idx_op_hs + i * 21 * ket_comps + 8 * ket_comps + j);

            auto to_xxzzz_0 = prim_buffer.data(idx_op_hs + i * 21 * ket_comps + 9 * ket_comps + j);

            auto to_xyyyy_0 = prim_buffer.data(idx_op_hs + i * 21 * ket_comps + 10 * ket_comps + j);

            auto to_xyyyz_0 = prim_buffer.data(idx_op_hs + i * 21 * ket_comps + 11 * ket_comps + j);

            auto to_xyyzz_0 = prim_buffer.data(idx_op_hs + i * 21 * ket_comps + 12 * ket_comps + j);

            auto to_xyzzz_0 = prim_buffer.data(idx_op_hs + i * 21 * ket_comps + 13 * ket_comps + j);

            auto to_xzzzz_0 = prim_buffer.data(idx_op_hs + i * 21 * ket_comps + 14 * ket_comps + j);

            auto to_yyyyy_0 = prim_buffer.data(idx_op_hs + i * 21 * ket_comps + 15 * ket_comps + j);

            auto to_yyyyz_0 = prim_buffer.data(idx_op_hs + i * 21 * ket_comps + 16 * ket_comps + j);

            auto to_yyyzz_0 = prim_buffer.data(idx_op_hs + i * 21 * ket_comps + 17 * ket_comps + j);

            auto to_yyzzz_0 = prim_buffer.data(idx_op_hs + i * 21 * ket_comps + 18 * ket_comps + j);

            auto to_yzzzz_0 = prim_buffer.data(idx_op_hs + i * 21 * ket_comps + 19 * ket_comps + j);

            auto to_zzzzz_0 = prim_buffer.data(idx_op_hs + i * 21 * ket_comps + 20 * ket_comps + j);

            // Set up components of auxiliary buffer : KS

            auto to_xxxxxxx_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 0 * ket_comps + j);

            auto to_xxxxxxy_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 1 * ket_comps + j);

            auto to_xxxxxxz_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 2 * ket_comps + j);

            auto to_xxxxxyy_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 3 * ket_comps + j);

            auto to_xxxxxyz_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 4 * ket_comps + j);

            auto to_xxxxxzz_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 5 * ket_comps + j);

            auto to_xxxxyyy_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 6 * ket_comps + j);

            auto to_xxxxyyz_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 7 * ket_comps + j);

            auto to_xxxxyzz_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 8 * ket_comps + j);

            auto to_xxxxzzz_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 9 * ket_comps + j);

            auto to_xxxyyyy_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 10 * ket_comps + j);

            auto to_xxxyyyz_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 11 * ket_comps + j);

            auto to_xxxyyzz_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 12 * ket_comps + j);

            auto to_xxxyzzz_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 13 * ket_comps + j);

            auto to_xxxzzzz_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 14 * ket_comps + j);

            auto to_xxyyyyy_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 15 * ket_comps + j);

            auto to_xxyyyyz_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 16 * ket_comps + j);

            auto to_xxyyyzz_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 17 * ket_comps + j);

            auto to_xxyyzzz_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 18 * ket_comps + j);

            auto to_xxyzzzz_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 19 * ket_comps + j);

            auto to_xxzzzzz_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 20 * ket_comps + j);

            auto to_xyyyyyy_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 21 * ket_comps + j);

            auto to_xyyyyyz_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 22 * ket_comps + j);

            auto to_xyyyyzz_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 23 * ket_comps + j);

            auto to_xyyyzzz_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 24 * ket_comps + j);

            auto to_xyyzzzz_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 25 * ket_comps + j);

            auto to_xyzzzzz_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 26 * ket_comps + j);

            auto to_xzzzzzz_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 27 * ket_comps + j);

            auto to_yyyyyyy_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 28 * ket_comps + j);

            auto to_yyyyyyz_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 29 * ket_comps + j);

            auto to_yyyyyzz_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 30 * ket_comps + j);

            auto to_yyyyzzz_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 31 * ket_comps + j);

            auto to_yyyzzzz_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 32 * ket_comps + j);

            auto to_yyzzzzz_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 33 * ket_comps + j);

            auto to_yzzzzzz_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 34 * ket_comps + j);

            auto to_zzzzzzz_0 = prim_buffer.data(idx_op_ks + i * 36 * ket_comps + 35 * ket_comps + j);

            // Set up 0-28 components of targeted buffer : IS

            auto to_x_0_xxxxxx_0 = prim_buffer.data(idx_op_geom_100_is + 0 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 0 * ket_comps + j);

            auto to_x_0_xxxxxy_0 = prim_buffer.data(idx_op_geom_100_is + 0 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 1 * ket_comps + j);

            auto to_x_0_xxxxxz_0 = prim_buffer.data(idx_op_geom_100_is + 0 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 2 * ket_comps + j);

            auto to_x_0_xxxxyy_0 = prim_buffer.data(idx_op_geom_100_is + 0 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 3 * ket_comps + j);

            auto to_x_0_xxxxyz_0 = prim_buffer.data(idx_op_geom_100_is + 0 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 4 * ket_comps + j);

            auto to_x_0_xxxxzz_0 = prim_buffer.data(idx_op_geom_100_is + 0 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 5 * ket_comps + j);

            auto to_x_0_xxxyyy_0 = prim_buffer.data(idx_op_geom_100_is + 0 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 6 * ket_comps + j);

            auto to_x_0_xxxyyz_0 = prim_buffer.data(idx_op_geom_100_is + 0 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 7 * ket_comps + j);

            auto to_x_0_xxxyzz_0 = prim_buffer.data(idx_op_geom_100_is + 0 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 8 * ket_comps + j);

            auto to_x_0_xxxzzz_0 = prim_buffer.data(idx_op_geom_100_is + 0 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 9 * ket_comps + j);

            auto to_x_0_xxyyyy_0 = prim_buffer.data(idx_op_geom_100_is + 0 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 10 * ket_comps + j);

            auto to_x_0_xxyyyz_0 = prim_buffer.data(idx_op_geom_100_is + 0 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 11 * ket_comps + j);

            auto to_x_0_xxyyzz_0 = prim_buffer.data(idx_op_geom_100_is + 0 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 12 * ket_comps + j);

            auto to_x_0_xxyzzz_0 = prim_buffer.data(idx_op_geom_100_is + 0 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 13 * ket_comps + j);

            auto to_x_0_xxzzzz_0 = prim_buffer.data(idx_op_geom_100_is + 0 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 14 * ket_comps + j);

            auto to_x_0_xyyyyy_0 = prim_buffer.data(idx_op_geom_100_is + 0 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 15 * ket_comps + j);

            auto to_x_0_xyyyyz_0 = prim_buffer.data(idx_op_geom_100_is + 0 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 16 * ket_comps + j);

            auto to_x_0_xyyyzz_0 = prim_buffer.data(idx_op_geom_100_is + 0 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 17 * ket_comps + j);

            auto to_x_0_xyyzzz_0 = prim_buffer.data(idx_op_geom_100_is + 0 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 18 * ket_comps + j);

            auto to_x_0_xyzzzz_0 = prim_buffer.data(idx_op_geom_100_is + 0 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 19 * ket_comps + j);

            auto to_x_0_xzzzzz_0 = prim_buffer.data(idx_op_geom_100_is + 0 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 20 * ket_comps + j);

            auto to_x_0_yyyyyy_0 = prim_buffer.data(idx_op_geom_100_is + 0 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 21 * ket_comps + j);

            auto to_x_0_yyyyyz_0 = prim_buffer.data(idx_op_geom_100_is + 0 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 22 * ket_comps + j);

            auto to_x_0_yyyyzz_0 = prim_buffer.data(idx_op_geom_100_is + 0 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 23 * ket_comps + j);

            auto to_x_0_yyyzzz_0 = prim_buffer.data(idx_op_geom_100_is + 0 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 24 * ket_comps + j);

            auto to_x_0_yyzzzz_0 = prim_buffer.data(idx_op_geom_100_is + 0 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 25 * ket_comps + j);

            auto to_x_0_yzzzzz_0 = prim_buffer.data(idx_op_geom_100_is + 0 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 26 * ket_comps + j);

            auto to_x_0_zzzzzz_0 = prim_buffer.data(idx_op_geom_100_is + 0 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 27 * ket_comps + j);

            #pragma omp simd aligned(to_x_0_xxxxxx_0, to_x_0_xxxxxy_0, to_x_0_xxxxxz_0, to_x_0_xxxxyy_0, to_x_0_xxxxyz_0, to_x_0_xxxxzz_0, to_x_0_xxxyyy_0, to_x_0_xxxyyz_0, to_x_0_xxxyzz_0, to_x_0_xxxzzz_0, to_x_0_xxyyyy_0, to_x_0_xxyyyz_0, to_x_0_xxyyzz_0, to_x_0_xxyzzz_0, to_x_0_xxzzzz_0, to_x_0_xyyyyy_0, to_x_0_xyyyyz_0, to_x_0_xyyyzz_0, to_x_0_xyyzzz_0, to_x_0_xyzzzz_0, to_x_0_xzzzzz_0, to_x_0_yyyyyy_0, to_x_0_yyyyyz_0, to_x_0_yyyyzz_0, to_x_0_yyyzzz_0, to_x_0_yyzzzz_0, to_x_0_yzzzzz_0, to_x_0_zzzzzz_0, to_xxxxx_0, to_xxxxxxx_0, to_xxxxxxy_0, to_xxxxxxz_0, to_xxxxxyy_0, to_xxxxxyz_0, to_xxxxxzz_0, to_xxxxy_0, to_xxxxyyy_0, to_xxxxyyz_0, to_xxxxyzz_0, to_xxxxz_0, to_xxxxzzz_0, to_xxxyy_0, to_xxxyyyy_0, to_xxxyyyz_0, to_xxxyyzz_0, to_xxxyz_0, to_xxxyzzz_0, to_xxxzz_0, to_xxxzzzz_0, to_xxyyy_0, to_xxyyyyy_0, to_xxyyyyz_0, to_xxyyyzz_0, to_xxyyz_0, to_xxyyzzz_0, to_xxyzz_0, to_xxyzzzz_0, to_xxzzz_0, to_xxzzzzz_0, to_xyyyy_0, to_xyyyyyy_0, to_xyyyyyz_0, to_xyyyyzz_0, to_xyyyz_0, to_xyyyzzz_0, to_xyyzz_0, to_xyyzzzz_0, to_xyzzz_0, to_xyzzzzz_0, to_xzzzz_0, to_xzzzzzz_0, to_yyyyy_0, to_yyyyz_0, to_yyyzz_0, to_yyzzz_0, to_yzzzz_0, to_zzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                const double tbe_0 = a_exp;

                to_x_0_xxxxxx_0[k] = -6.0 * to_xxxxx_0[k] + 2.0 * to_xxxxxxx_0[k] * tbe_0;

                to_x_0_xxxxxy_0[k] = -5.0 * to_xxxxy_0[k] + 2.0 * to_xxxxxxy_0[k] * tbe_0;

                to_x_0_xxxxxz_0[k] = -5.0 * to_xxxxz_0[k] + 2.0 * to_xxxxxxz_0[k] * tbe_0;

                to_x_0_xxxxyy_0[k] = -4.0 * to_xxxyy_0[k] + 2.0 * to_xxxxxyy_0[k] * tbe_0;

                to_x_0_xxxxyz_0[k] = -4.0 * to_xxxyz_0[k] + 2.0 * to_xxxxxyz_0[k] * tbe_0;

                to_x_0_xxxxzz_0[k] = -4.0 * to_xxxzz_0[k] + 2.0 * to_xxxxxzz_0[k] * tbe_0;

                to_x_0_xxxyyy_0[k] = -3.0 * to_xxyyy_0[k] + 2.0 * to_xxxxyyy_0[k] * tbe_0;

                to_x_0_xxxyyz_0[k] = -3.0 * to_xxyyz_0[k] + 2.0 * to_xxxxyyz_0[k] * tbe_0;

                to_x_0_xxxyzz_0[k] = -3.0 * to_xxyzz_0[k] + 2.0 * to_xxxxyzz_0[k] * tbe_0;

                to_x_0_xxxzzz_0[k] = -3.0 * to_xxzzz_0[k] + 2.0 * to_xxxxzzz_0[k] * tbe_0;

                to_x_0_xxyyyy_0[k] = -2.0 * to_xyyyy_0[k] + 2.0 * to_xxxyyyy_0[k] * tbe_0;

                to_x_0_xxyyyz_0[k] = -2.0 * to_xyyyz_0[k] + 2.0 * to_xxxyyyz_0[k] * tbe_0;

                to_x_0_xxyyzz_0[k] = -2.0 * to_xyyzz_0[k] + 2.0 * to_xxxyyzz_0[k] * tbe_0;

                to_x_0_xxyzzz_0[k] = -2.0 * to_xyzzz_0[k] + 2.0 * to_xxxyzzz_0[k] * tbe_0;

                to_x_0_xxzzzz_0[k] = -2.0 * to_xzzzz_0[k] + 2.0 * to_xxxzzzz_0[k] * tbe_0;

                to_x_0_xyyyyy_0[k] = -to_yyyyy_0[k] + 2.0 * to_xxyyyyy_0[k] * tbe_0;

                to_x_0_xyyyyz_0[k] = -to_yyyyz_0[k] + 2.0 * to_xxyyyyz_0[k] * tbe_0;

                to_x_0_xyyyzz_0[k] = -to_yyyzz_0[k] + 2.0 * to_xxyyyzz_0[k] * tbe_0;

                to_x_0_xyyzzz_0[k] = -to_yyzzz_0[k] + 2.0 * to_xxyyzzz_0[k] * tbe_0;

                to_x_0_xyzzzz_0[k] = -to_yzzzz_0[k] + 2.0 * to_xxyzzzz_0[k] * tbe_0;

                to_x_0_xzzzzz_0[k] = -to_zzzzz_0[k] + 2.0 * to_xxzzzzz_0[k] * tbe_0;

                to_x_0_yyyyyy_0[k] = 2.0 * to_xyyyyyy_0[k] * tbe_0;

                to_x_0_yyyyyz_0[k] = 2.0 * to_xyyyyyz_0[k] * tbe_0;

                to_x_0_yyyyzz_0[k] = 2.0 * to_xyyyyzz_0[k] * tbe_0;

                to_x_0_yyyzzz_0[k] = 2.0 * to_xyyyzzz_0[k] * tbe_0;

                to_x_0_yyzzzz_0[k] = 2.0 * to_xyyzzzz_0[k] * tbe_0;

                to_x_0_yzzzzz_0[k] = 2.0 * to_xyzzzzz_0[k] * tbe_0;

                to_x_0_zzzzzz_0[k] = 2.0 * to_xzzzzzz_0[k] * tbe_0;
            }

            // Set up 28-56 components of targeted buffer : IS

            auto to_y_0_xxxxxx_0 = prim_buffer.data(idx_op_geom_100_is + 1 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 0 * ket_comps + j);

            auto to_y_0_xxxxxy_0 = prim_buffer.data(idx_op_geom_100_is + 1 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 1 * ket_comps + j);

            auto to_y_0_xxxxxz_0 = prim_buffer.data(idx_op_geom_100_is + 1 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 2 * ket_comps + j);

            auto to_y_0_xxxxyy_0 = prim_buffer.data(idx_op_geom_100_is + 1 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 3 * ket_comps + j);

            auto to_y_0_xxxxyz_0 = prim_buffer.data(idx_op_geom_100_is + 1 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 4 * ket_comps + j);

            auto to_y_0_xxxxzz_0 = prim_buffer.data(idx_op_geom_100_is + 1 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 5 * ket_comps + j);

            auto to_y_0_xxxyyy_0 = prim_buffer.data(idx_op_geom_100_is + 1 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 6 * ket_comps + j);

            auto to_y_0_xxxyyz_0 = prim_buffer.data(idx_op_geom_100_is + 1 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 7 * ket_comps + j);

            auto to_y_0_xxxyzz_0 = prim_buffer.data(idx_op_geom_100_is + 1 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 8 * ket_comps + j);

            auto to_y_0_xxxzzz_0 = prim_buffer.data(idx_op_geom_100_is + 1 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 9 * ket_comps + j);

            auto to_y_0_xxyyyy_0 = prim_buffer.data(idx_op_geom_100_is + 1 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 10 * ket_comps + j);

            auto to_y_0_xxyyyz_0 = prim_buffer.data(idx_op_geom_100_is + 1 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 11 * ket_comps + j);

            auto to_y_0_xxyyzz_0 = prim_buffer.data(idx_op_geom_100_is + 1 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 12 * ket_comps + j);

            auto to_y_0_xxyzzz_0 = prim_buffer.data(idx_op_geom_100_is + 1 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 13 * ket_comps + j);

            auto to_y_0_xxzzzz_0 = prim_buffer.data(idx_op_geom_100_is + 1 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 14 * ket_comps + j);

            auto to_y_0_xyyyyy_0 = prim_buffer.data(idx_op_geom_100_is + 1 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 15 * ket_comps + j);

            auto to_y_0_xyyyyz_0 = prim_buffer.data(idx_op_geom_100_is + 1 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 16 * ket_comps + j);

            auto to_y_0_xyyyzz_0 = prim_buffer.data(idx_op_geom_100_is + 1 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 17 * ket_comps + j);

            auto to_y_0_xyyzzz_0 = prim_buffer.data(idx_op_geom_100_is + 1 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 18 * ket_comps + j);

            auto to_y_0_xyzzzz_0 = prim_buffer.data(idx_op_geom_100_is + 1 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 19 * ket_comps + j);

            auto to_y_0_xzzzzz_0 = prim_buffer.data(idx_op_geom_100_is + 1 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 20 * ket_comps + j);

            auto to_y_0_yyyyyy_0 = prim_buffer.data(idx_op_geom_100_is + 1 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 21 * ket_comps + j);

            auto to_y_0_yyyyyz_0 = prim_buffer.data(idx_op_geom_100_is + 1 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 22 * ket_comps + j);

            auto to_y_0_yyyyzz_0 = prim_buffer.data(idx_op_geom_100_is + 1 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 23 * ket_comps + j);

            auto to_y_0_yyyzzz_0 = prim_buffer.data(idx_op_geom_100_is + 1 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 24 * ket_comps + j);

            auto to_y_0_yyzzzz_0 = prim_buffer.data(idx_op_geom_100_is + 1 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 25 * ket_comps + j);

            auto to_y_0_yzzzzz_0 = prim_buffer.data(idx_op_geom_100_is + 1 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 26 * ket_comps + j);

            auto to_y_0_zzzzzz_0 = prim_buffer.data(idx_op_geom_100_is + 1 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 27 * ket_comps + j);

            #pragma omp simd aligned(to_xxxxx_0, to_xxxxxxy_0, to_xxxxxyy_0, to_xxxxxyz_0, to_xxxxy_0, to_xxxxyyy_0, to_xxxxyyz_0, to_xxxxyzz_0, to_xxxxz_0, to_xxxyy_0, to_xxxyyyy_0, to_xxxyyyz_0, to_xxxyyzz_0, to_xxxyz_0, to_xxxyzzz_0, to_xxxzz_0, to_xxyyy_0, to_xxyyyyy_0, to_xxyyyyz_0, to_xxyyyzz_0, to_xxyyz_0, to_xxyyzzz_0, to_xxyzz_0, to_xxyzzzz_0, to_xxzzz_0, to_xyyyy_0, to_xyyyyyy_0, to_xyyyyyz_0, to_xyyyyzz_0, to_xyyyz_0, to_xyyyzzz_0, to_xyyzz_0, to_xyyzzzz_0, to_xyzzz_0, to_xyzzzzz_0, to_xzzzz_0, to_y_0_xxxxxx_0, to_y_0_xxxxxy_0, to_y_0_xxxxxz_0, to_y_0_xxxxyy_0, to_y_0_xxxxyz_0, to_y_0_xxxxzz_0, to_y_0_xxxyyy_0, to_y_0_xxxyyz_0, to_y_0_xxxyzz_0, to_y_0_xxxzzz_0, to_y_0_xxyyyy_0, to_y_0_xxyyyz_0, to_y_0_xxyyzz_0, to_y_0_xxyzzz_0, to_y_0_xxzzzz_0, to_y_0_xyyyyy_0, to_y_0_xyyyyz_0, to_y_0_xyyyzz_0, to_y_0_xyyzzz_0, to_y_0_xyzzzz_0, to_y_0_xzzzzz_0, to_y_0_yyyyyy_0, to_y_0_yyyyyz_0, to_y_0_yyyyzz_0, to_y_0_yyyzzz_0, to_y_0_yyzzzz_0, to_y_0_yzzzzz_0, to_y_0_zzzzzz_0, to_yyyyy_0, to_yyyyyyy_0, to_yyyyyyz_0, to_yyyyyzz_0, to_yyyyz_0, to_yyyyzzz_0, to_yyyzz_0, to_yyyzzzz_0, to_yyzzz_0, to_yyzzzzz_0, to_yzzzz_0, to_yzzzzzz_0, to_zzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                const double tbe_0 = a_exp;

                to_y_0_xxxxxx_0[k] = 2.0 * to_xxxxxxy_0[k] * tbe_0;

                to_y_0_xxxxxy_0[k] = -to_xxxxx_0[k] + 2.0 * to_xxxxxyy_0[k] * tbe_0;

                to_y_0_xxxxxz_0[k] = 2.0 * to_xxxxxyz_0[k] * tbe_0;

                to_y_0_xxxxyy_0[k] = -2.0 * to_xxxxy_0[k] + 2.0 * to_xxxxyyy_0[k] * tbe_0;

                to_y_0_xxxxyz_0[k] = -to_xxxxz_0[k] + 2.0 * to_xxxxyyz_0[k] * tbe_0;

                to_y_0_xxxxzz_0[k] = 2.0 * to_xxxxyzz_0[k] * tbe_0;

                to_y_0_xxxyyy_0[k] = -3.0 * to_xxxyy_0[k] + 2.0 * to_xxxyyyy_0[k] * tbe_0;

                to_y_0_xxxyyz_0[k] = -2.0 * to_xxxyz_0[k] + 2.0 * to_xxxyyyz_0[k] * tbe_0;

                to_y_0_xxxyzz_0[k] = -to_xxxzz_0[k] + 2.0 * to_xxxyyzz_0[k] * tbe_0;

                to_y_0_xxxzzz_0[k] = 2.0 * to_xxxyzzz_0[k] * tbe_0;

                to_y_0_xxyyyy_0[k] = -4.0 * to_xxyyy_0[k] + 2.0 * to_xxyyyyy_0[k] * tbe_0;

                to_y_0_xxyyyz_0[k] = -3.0 * to_xxyyz_0[k] + 2.0 * to_xxyyyyz_0[k] * tbe_0;

                to_y_0_xxyyzz_0[k] = -2.0 * to_xxyzz_0[k] + 2.0 * to_xxyyyzz_0[k] * tbe_0;

                to_y_0_xxyzzz_0[k] = -to_xxzzz_0[k] + 2.0 * to_xxyyzzz_0[k] * tbe_0;

                to_y_0_xxzzzz_0[k] = 2.0 * to_xxyzzzz_0[k] * tbe_0;

                to_y_0_xyyyyy_0[k] = -5.0 * to_xyyyy_0[k] + 2.0 * to_xyyyyyy_0[k] * tbe_0;

                to_y_0_xyyyyz_0[k] = -4.0 * to_xyyyz_0[k] + 2.0 * to_xyyyyyz_0[k] * tbe_0;

                to_y_0_xyyyzz_0[k] = -3.0 * to_xyyzz_0[k] + 2.0 * to_xyyyyzz_0[k] * tbe_0;

                to_y_0_xyyzzz_0[k] = -2.0 * to_xyzzz_0[k] + 2.0 * to_xyyyzzz_0[k] * tbe_0;

                to_y_0_xyzzzz_0[k] = -to_xzzzz_0[k] + 2.0 * to_xyyzzzz_0[k] * tbe_0;

                to_y_0_xzzzzz_0[k] = 2.0 * to_xyzzzzz_0[k] * tbe_0;

                to_y_0_yyyyyy_0[k] = -6.0 * to_yyyyy_0[k] + 2.0 * to_yyyyyyy_0[k] * tbe_0;

                to_y_0_yyyyyz_0[k] = -5.0 * to_yyyyz_0[k] + 2.0 * to_yyyyyyz_0[k] * tbe_0;

                to_y_0_yyyyzz_0[k] = -4.0 * to_yyyzz_0[k] + 2.0 * to_yyyyyzz_0[k] * tbe_0;

                to_y_0_yyyzzz_0[k] = -3.0 * to_yyzzz_0[k] + 2.0 * to_yyyyzzz_0[k] * tbe_0;

                to_y_0_yyzzzz_0[k] = -2.0 * to_yzzzz_0[k] + 2.0 * to_yyyzzzz_0[k] * tbe_0;

                to_y_0_yzzzzz_0[k] = -to_zzzzz_0[k] + 2.0 * to_yyzzzzz_0[k] * tbe_0;

                to_y_0_zzzzzz_0[k] = 2.0 * to_yzzzzzz_0[k] * tbe_0;
            }

            // Set up 56-84 components of targeted buffer : IS

            auto to_z_0_xxxxxx_0 = prim_buffer.data(idx_op_geom_100_is + 2 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 0 * ket_comps + j);

            auto to_z_0_xxxxxy_0 = prim_buffer.data(idx_op_geom_100_is + 2 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 1 * ket_comps + j);

            auto to_z_0_xxxxxz_0 = prim_buffer.data(idx_op_geom_100_is + 2 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 2 * ket_comps + j);

            auto to_z_0_xxxxyy_0 = prim_buffer.data(idx_op_geom_100_is + 2 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 3 * ket_comps + j);

            auto to_z_0_xxxxyz_0 = prim_buffer.data(idx_op_geom_100_is + 2 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 4 * ket_comps + j);

            auto to_z_0_xxxxzz_0 = prim_buffer.data(idx_op_geom_100_is + 2 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 5 * ket_comps + j);

            auto to_z_0_xxxyyy_0 = prim_buffer.data(idx_op_geom_100_is + 2 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 6 * ket_comps + j);

            auto to_z_0_xxxyyz_0 = prim_buffer.data(idx_op_geom_100_is + 2 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 7 * ket_comps + j);

            auto to_z_0_xxxyzz_0 = prim_buffer.data(idx_op_geom_100_is + 2 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 8 * ket_comps + j);

            auto to_z_0_xxxzzz_0 = prim_buffer.data(idx_op_geom_100_is + 2 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 9 * ket_comps + j);

            auto to_z_0_xxyyyy_0 = prim_buffer.data(idx_op_geom_100_is + 2 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 10 * ket_comps + j);

            auto to_z_0_xxyyyz_0 = prim_buffer.data(idx_op_geom_100_is + 2 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 11 * ket_comps + j);

            auto to_z_0_xxyyzz_0 = prim_buffer.data(idx_op_geom_100_is + 2 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 12 * ket_comps + j);

            auto to_z_0_xxyzzz_0 = prim_buffer.data(idx_op_geom_100_is + 2 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 13 * ket_comps + j);

            auto to_z_0_xxzzzz_0 = prim_buffer.data(idx_op_geom_100_is + 2 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 14 * ket_comps + j);

            auto to_z_0_xyyyyy_0 = prim_buffer.data(idx_op_geom_100_is + 2 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 15 * ket_comps + j);

            auto to_z_0_xyyyyz_0 = prim_buffer.data(idx_op_geom_100_is + 2 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 16 * ket_comps + j);

            auto to_z_0_xyyyzz_0 = prim_buffer.data(idx_op_geom_100_is + 2 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 17 * ket_comps + j);

            auto to_z_0_xyyzzz_0 = prim_buffer.data(idx_op_geom_100_is + 2 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 18 * ket_comps + j);

            auto to_z_0_xyzzzz_0 = prim_buffer.data(idx_op_geom_100_is + 2 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 19 * ket_comps + j);

            auto to_z_0_xzzzzz_0 = prim_buffer.data(idx_op_geom_100_is + 2 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 20 * ket_comps + j);

            auto to_z_0_yyyyyy_0 = prim_buffer.data(idx_op_geom_100_is + 2 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 21 * ket_comps + j);

            auto to_z_0_yyyyyz_0 = prim_buffer.data(idx_op_geom_100_is + 2 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 22 * ket_comps + j);

            auto to_z_0_yyyyzz_0 = prim_buffer.data(idx_op_geom_100_is + 2 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 23 * ket_comps + j);

            auto to_z_0_yyyzzz_0 = prim_buffer.data(idx_op_geom_100_is + 2 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 24 * ket_comps + j);

            auto to_z_0_yyzzzz_0 = prim_buffer.data(idx_op_geom_100_is + 2 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 25 * ket_comps + j);

            auto to_z_0_yzzzzz_0 = prim_buffer.data(idx_op_geom_100_is + 2 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 26 * ket_comps + j);

            auto to_z_0_zzzzzz_0 = prim_buffer.data(idx_op_geom_100_is + 2 * op_comps * 28 * ket_comps + i * 28 * ket_comps + 27 * ket_comps + j);

            #pragma omp simd aligned(to_xxxxx_0, to_xxxxxxz_0, to_xxxxxyz_0, to_xxxxxzz_0, to_xxxxy_0, to_xxxxyyz_0, to_xxxxyzz_0, to_xxxxz_0, to_xxxxzzz_0, to_xxxyy_0, to_xxxyyyz_0, to_xxxyyzz_0, to_xxxyz_0, to_xxxyzzz_0, to_xxxzz_0, to_xxxzzzz_0, to_xxyyy_0, to_xxyyyyz_0, to_xxyyyzz_0, to_xxyyz_0, to_xxyyzzz_0, to_xxyzz_0, to_xxyzzzz_0, to_xxzzz_0, to_xxzzzzz_0, to_xyyyy_0, to_xyyyyyz_0, to_xyyyyzz_0, to_xyyyz_0, to_xyyyzzz_0, to_xyyzz_0, to_xyyzzzz_0, to_xyzzz_0, to_xyzzzzz_0, to_xzzzz_0, to_xzzzzzz_0, to_yyyyy_0, to_yyyyyyz_0, to_yyyyyzz_0, to_yyyyz_0, to_yyyyzzz_0, to_yyyzz_0, to_yyyzzzz_0, to_yyzzz_0, to_yyzzzzz_0, to_yzzzz_0, to_yzzzzzz_0, to_z_0_xxxxxx_0, to_z_0_xxxxxy_0, to_z_0_xxxxxz_0, to_z_0_xxxxyy_0, to_z_0_xxxxyz_0, to_z_0_xxxxzz_0, to_z_0_xxxyyy_0, to_z_0_xxxyyz_0, to_z_0_xxxyzz_0, to_z_0_xxxzzz_0, to_z_0_xxyyyy_0, to_z_0_xxyyyz_0, to_z_0_xxyyzz_0, to_z_0_xxyzzz_0, to_z_0_xxzzzz_0, to_z_0_xyyyyy_0, to_z_0_xyyyyz_0, to_z_0_xyyyzz_0, to_z_0_xyyzzz_0, to_z_0_xyzzzz_0, to_z_0_xzzzzz_0, to_z_0_yyyyyy_0, to_z_0_yyyyyz_0, to_z_0_yyyyzz_0, to_z_0_yyyzzz_0, to_z_0_yyzzzz_0, to_z_0_yzzzzz_0, to_z_0_zzzzzz_0, to_zzzzz_0, to_zzzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                const double tbe_0 = a_exp;

                to_z_0_xxxxxx_0[k] = 2.0 * to_xxxxxxz_0[k] * tbe_0;

                to_z_0_xxxxxy_0[k] = 2.0 * to_xxxxxyz_0[k] * tbe_0;

                to_z_0_xxxxxz_0[k] = -to_xxxxx_0[k] + 2.0 * to_xxxxxzz_0[k] * tbe_0;

                to_z_0_xxxxyy_0[k] = 2.0 * to_xxxxyyz_0[k] * tbe_0;

                to_z_0_xxxxyz_0[k] = -to_xxxxy_0[k] + 2.0 * to_xxxxyzz_0[k] * tbe_0;

                to_z_0_xxxxzz_0[k] = -2.0 * to_xxxxz_0[k] + 2.0 * to_xxxxzzz_0[k] * tbe_0;

                to_z_0_xxxyyy_0[k] = 2.0 * to_xxxyyyz_0[k] * tbe_0;

                to_z_0_xxxyyz_0[k] = -to_xxxyy_0[k] + 2.0 * to_xxxyyzz_0[k] * tbe_0;

                to_z_0_xxxyzz_0[k] = -2.0 * to_xxxyz_0[k] + 2.0 * to_xxxyzzz_0[k] * tbe_0;

                to_z_0_xxxzzz_0[k] = -3.0 * to_xxxzz_0[k] + 2.0 * to_xxxzzzz_0[k] * tbe_0;

                to_z_0_xxyyyy_0[k] = 2.0 * to_xxyyyyz_0[k] * tbe_0;

                to_z_0_xxyyyz_0[k] = -to_xxyyy_0[k] + 2.0 * to_xxyyyzz_0[k] * tbe_0;

                to_z_0_xxyyzz_0[k] = -2.0 * to_xxyyz_0[k] + 2.0 * to_xxyyzzz_0[k] * tbe_0;

                to_z_0_xxyzzz_0[k] = -3.0 * to_xxyzz_0[k] + 2.0 * to_xxyzzzz_0[k] * tbe_0;

                to_z_0_xxzzzz_0[k] = -4.0 * to_xxzzz_0[k] + 2.0 * to_xxzzzzz_0[k] * tbe_0;

                to_z_0_xyyyyy_0[k] = 2.0 * to_xyyyyyz_0[k] * tbe_0;

                to_z_0_xyyyyz_0[k] = -to_xyyyy_0[k] + 2.0 * to_xyyyyzz_0[k] * tbe_0;

                to_z_0_xyyyzz_0[k] = -2.0 * to_xyyyz_0[k] + 2.0 * to_xyyyzzz_0[k] * tbe_0;

                to_z_0_xyyzzz_0[k] = -3.0 * to_xyyzz_0[k] + 2.0 * to_xyyzzzz_0[k] * tbe_0;

                to_z_0_xyzzzz_0[k] = -4.0 * to_xyzzz_0[k] + 2.0 * to_xyzzzzz_0[k] * tbe_0;

                to_z_0_xzzzzz_0[k] = -5.0 * to_xzzzz_0[k] + 2.0 * to_xzzzzzz_0[k] * tbe_0;

                to_z_0_yyyyyy_0[k] = 2.0 * to_yyyyyyz_0[k] * tbe_0;

                to_z_0_yyyyyz_0[k] = -to_yyyyy_0[k] + 2.0 * to_yyyyyzz_0[k] * tbe_0;

                to_z_0_yyyyzz_0[k] = -2.0 * to_yyyyz_0[k] + 2.0 * to_yyyyzzz_0[k] * tbe_0;

                to_z_0_yyyzzz_0[k] = -3.0 * to_yyyzz_0[k] + 2.0 * to_yyyzzzz_0[k] * tbe_0;

                to_z_0_yyzzzz_0[k] = -4.0 * to_yyzzz_0[k] + 2.0 * to_yyzzzzz_0[k] * tbe_0;

                to_z_0_yzzzzz_0[k] = -5.0 * to_yzzzz_0[k] + 2.0 * to_yzzzzzz_0[k] * tbe_0;

                to_z_0_zzzzzz_0[k] = -6.0 * to_zzzzz_0[k] + 2.0 * to_zzzzzzz_0[k] * tbe_0;
            }

        }
    }

}

} // t2cgeom namespace

