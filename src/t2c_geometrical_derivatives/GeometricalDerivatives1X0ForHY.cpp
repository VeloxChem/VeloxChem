#include "GeometricalDerivatives1X0ForHY.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_10_hx(CSimdArray<double>& prim_buffer,
                        const int idx_op_geom_100_hs,
                        const int idx_op_gs,
                        const int idx_op_is,
                        const int op_comps,
                        const int ket_comps,
                        const double a_exp) -> void
{
    const auto nelems = prim_buffer.number_of_active_elements();

    for (size_t i = 0; i < op_comps; i++)
    {
        for (size_t j = 0; j < ket_comps; j++)
        {
            // Set up components of auxiliary buffer : GS

            auto to_xxxx_0 = prim_buffer.data(idx_op_gs + i * 15 * ket_comps + 0 * ket_comps + j);

            auto to_xxxy_0 = prim_buffer.data(idx_op_gs + i * 15 * ket_comps + 1 * ket_comps + j);

            auto to_xxxz_0 = prim_buffer.data(idx_op_gs + i * 15 * ket_comps + 2 * ket_comps + j);

            auto to_xxyy_0 = prim_buffer.data(idx_op_gs + i * 15 * ket_comps + 3 * ket_comps + j);

            auto to_xxyz_0 = prim_buffer.data(idx_op_gs + i * 15 * ket_comps + 4 * ket_comps + j);

            auto to_xxzz_0 = prim_buffer.data(idx_op_gs + i * 15 * ket_comps + 5 * ket_comps + j);

            auto to_xyyy_0 = prim_buffer.data(idx_op_gs + i * 15 * ket_comps + 6 * ket_comps + j);

            auto to_xyyz_0 = prim_buffer.data(idx_op_gs + i * 15 * ket_comps + 7 * ket_comps + j);

            auto to_xyzz_0 = prim_buffer.data(idx_op_gs + i * 15 * ket_comps + 8 * ket_comps + j);

            auto to_xzzz_0 = prim_buffer.data(idx_op_gs + i * 15 * ket_comps + 9 * ket_comps + j);

            auto to_yyyy_0 = prim_buffer.data(idx_op_gs + i * 15 * ket_comps + 10 * ket_comps + j);

            auto to_yyyz_0 = prim_buffer.data(idx_op_gs + i * 15 * ket_comps + 11 * ket_comps + j);

            auto to_yyzz_0 = prim_buffer.data(idx_op_gs + i * 15 * ket_comps + 12 * ket_comps + j);

            auto to_yzzz_0 = prim_buffer.data(idx_op_gs + i * 15 * ket_comps + 13 * ket_comps + j);

            auto to_zzzz_0 = prim_buffer.data(idx_op_gs + i * 15 * ket_comps + 14 * ket_comps + j);

            // Set up components of auxiliary buffer : IS

            auto to_xxxxxx_0 = prim_buffer.data(idx_op_is + i * 28 * ket_comps + 0 * ket_comps + j);

            auto to_xxxxxy_0 = prim_buffer.data(idx_op_is + i * 28 * ket_comps + 1 * ket_comps + j);

            auto to_xxxxxz_0 = prim_buffer.data(idx_op_is + i * 28 * ket_comps + 2 * ket_comps + j);

            auto to_xxxxyy_0 = prim_buffer.data(idx_op_is + i * 28 * ket_comps + 3 * ket_comps + j);

            auto to_xxxxyz_0 = prim_buffer.data(idx_op_is + i * 28 * ket_comps + 4 * ket_comps + j);

            auto to_xxxxzz_0 = prim_buffer.data(idx_op_is + i * 28 * ket_comps + 5 * ket_comps + j);

            auto to_xxxyyy_0 = prim_buffer.data(idx_op_is + i * 28 * ket_comps + 6 * ket_comps + j);

            auto to_xxxyyz_0 = prim_buffer.data(idx_op_is + i * 28 * ket_comps + 7 * ket_comps + j);

            auto to_xxxyzz_0 = prim_buffer.data(idx_op_is + i * 28 * ket_comps + 8 * ket_comps + j);

            auto to_xxxzzz_0 = prim_buffer.data(idx_op_is + i * 28 * ket_comps + 9 * ket_comps + j);

            auto to_xxyyyy_0 = prim_buffer.data(idx_op_is + i * 28 * ket_comps + 10 * ket_comps + j);

            auto to_xxyyyz_0 = prim_buffer.data(idx_op_is + i * 28 * ket_comps + 11 * ket_comps + j);

            auto to_xxyyzz_0 = prim_buffer.data(idx_op_is + i * 28 * ket_comps + 12 * ket_comps + j);

            auto to_xxyzzz_0 = prim_buffer.data(idx_op_is + i * 28 * ket_comps + 13 * ket_comps + j);

            auto to_xxzzzz_0 = prim_buffer.data(idx_op_is + i * 28 * ket_comps + 14 * ket_comps + j);

            auto to_xyyyyy_0 = prim_buffer.data(idx_op_is + i * 28 * ket_comps + 15 * ket_comps + j);

            auto to_xyyyyz_0 = prim_buffer.data(idx_op_is + i * 28 * ket_comps + 16 * ket_comps + j);

            auto to_xyyyzz_0 = prim_buffer.data(idx_op_is + i * 28 * ket_comps + 17 * ket_comps + j);

            auto to_xyyzzz_0 = prim_buffer.data(idx_op_is + i * 28 * ket_comps + 18 * ket_comps + j);

            auto to_xyzzzz_0 = prim_buffer.data(idx_op_is + i * 28 * ket_comps + 19 * ket_comps + j);

            auto to_xzzzzz_0 = prim_buffer.data(idx_op_is + i * 28 * ket_comps + 20 * ket_comps + j);

            auto to_yyyyyy_0 = prim_buffer.data(idx_op_is + i * 28 * ket_comps + 21 * ket_comps + j);

            auto to_yyyyyz_0 = prim_buffer.data(idx_op_is + i * 28 * ket_comps + 22 * ket_comps + j);

            auto to_yyyyzz_0 = prim_buffer.data(idx_op_is + i * 28 * ket_comps + 23 * ket_comps + j);

            auto to_yyyzzz_0 = prim_buffer.data(idx_op_is + i * 28 * ket_comps + 24 * ket_comps + j);

            auto to_yyzzzz_0 = prim_buffer.data(idx_op_is + i * 28 * ket_comps + 25 * ket_comps + j);

            auto to_yzzzzz_0 = prim_buffer.data(idx_op_is + i * 28 * ket_comps + 26 * ket_comps + j);

            auto to_zzzzzz_0 = prim_buffer.data(idx_op_is + i * 28 * ket_comps + 27 * ket_comps + j);

            // Set up 0-21 components of targeted buffer : HS

            auto to_x_0_xxxxx_0 = prim_buffer.data(idx_op_geom_100_hs + 0 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 0 * ket_comps + j);

            auto to_x_0_xxxxy_0 = prim_buffer.data(idx_op_geom_100_hs + 0 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 1 * ket_comps + j);

            auto to_x_0_xxxxz_0 = prim_buffer.data(idx_op_geom_100_hs + 0 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 2 * ket_comps + j);

            auto to_x_0_xxxyy_0 = prim_buffer.data(idx_op_geom_100_hs + 0 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 3 * ket_comps + j);

            auto to_x_0_xxxyz_0 = prim_buffer.data(idx_op_geom_100_hs + 0 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 4 * ket_comps + j);

            auto to_x_0_xxxzz_0 = prim_buffer.data(idx_op_geom_100_hs + 0 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 5 * ket_comps + j);

            auto to_x_0_xxyyy_0 = prim_buffer.data(idx_op_geom_100_hs + 0 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 6 * ket_comps + j);

            auto to_x_0_xxyyz_0 = prim_buffer.data(idx_op_geom_100_hs + 0 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 7 * ket_comps + j);

            auto to_x_0_xxyzz_0 = prim_buffer.data(idx_op_geom_100_hs + 0 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 8 * ket_comps + j);

            auto to_x_0_xxzzz_0 = prim_buffer.data(idx_op_geom_100_hs + 0 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 9 * ket_comps + j);

            auto to_x_0_xyyyy_0 = prim_buffer.data(idx_op_geom_100_hs + 0 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 10 * ket_comps + j);

            auto to_x_0_xyyyz_0 = prim_buffer.data(idx_op_geom_100_hs + 0 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 11 * ket_comps + j);

            auto to_x_0_xyyzz_0 = prim_buffer.data(idx_op_geom_100_hs + 0 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 12 * ket_comps + j);

            auto to_x_0_xyzzz_0 = prim_buffer.data(idx_op_geom_100_hs + 0 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 13 * ket_comps + j);

            auto to_x_0_xzzzz_0 = prim_buffer.data(idx_op_geom_100_hs + 0 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 14 * ket_comps + j);

            auto to_x_0_yyyyy_0 = prim_buffer.data(idx_op_geom_100_hs + 0 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 15 * ket_comps + j);

            auto to_x_0_yyyyz_0 = prim_buffer.data(idx_op_geom_100_hs + 0 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 16 * ket_comps + j);

            auto to_x_0_yyyzz_0 = prim_buffer.data(idx_op_geom_100_hs + 0 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 17 * ket_comps + j);

            auto to_x_0_yyzzz_0 = prim_buffer.data(idx_op_geom_100_hs + 0 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 18 * ket_comps + j);

            auto to_x_0_yzzzz_0 = prim_buffer.data(idx_op_geom_100_hs + 0 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 19 * ket_comps + j);

            auto to_x_0_zzzzz_0 = prim_buffer.data(idx_op_geom_100_hs + 0 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 20 * ket_comps + j);

            #pragma omp simd aligned(to_x_0_xxxxx_0, to_x_0_xxxxy_0, to_x_0_xxxxz_0, to_x_0_xxxyy_0, to_x_0_xxxyz_0, to_x_0_xxxzz_0, to_x_0_xxyyy_0, to_x_0_xxyyz_0, to_x_0_xxyzz_0, to_x_0_xxzzz_0, to_x_0_xyyyy_0, to_x_0_xyyyz_0, to_x_0_xyyzz_0, to_x_0_xyzzz_0, to_x_0_xzzzz_0, to_x_0_yyyyy_0, to_x_0_yyyyz_0, to_x_0_yyyzz_0, to_x_0_yyzzz_0, to_x_0_yzzzz_0, to_x_0_zzzzz_0, to_xxxx_0, to_xxxxxx_0, to_xxxxxy_0, to_xxxxxz_0, to_xxxxyy_0, to_xxxxyz_0, to_xxxxzz_0, to_xxxy_0, to_xxxyyy_0, to_xxxyyz_0, to_xxxyzz_0, to_xxxz_0, to_xxxzzz_0, to_xxyy_0, to_xxyyyy_0, to_xxyyyz_0, to_xxyyzz_0, to_xxyz_0, to_xxyzzz_0, to_xxzz_0, to_xxzzzz_0, to_xyyy_0, to_xyyyyy_0, to_xyyyyz_0, to_xyyyzz_0, to_xyyz_0, to_xyyzzz_0, to_xyzz_0, to_xyzzzz_0, to_xzzz_0, to_xzzzzz_0, to_yyyy_0, to_yyyz_0, to_yyzz_0, to_yzzz_0, to_zzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                const double tbe_0 = a_exp;

                to_x_0_xxxxx_0[k] = -5.0 * to_xxxx_0[k] + 2.0 * to_xxxxxx_0[k] * tbe_0;

                to_x_0_xxxxy_0[k] = -4.0 * to_xxxy_0[k] + 2.0 * to_xxxxxy_0[k] * tbe_0;

                to_x_0_xxxxz_0[k] = -4.0 * to_xxxz_0[k] + 2.0 * to_xxxxxz_0[k] * tbe_0;

                to_x_0_xxxyy_0[k] = -3.0 * to_xxyy_0[k] + 2.0 * to_xxxxyy_0[k] * tbe_0;

                to_x_0_xxxyz_0[k] = -3.0 * to_xxyz_0[k] + 2.0 * to_xxxxyz_0[k] * tbe_0;

                to_x_0_xxxzz_0[k] = -3.0 * to_xxzz_0[k] + 2.0 * to_xxxxzz_0[k] * tbe_0;

                to_x_0_xxyyy_0[k] = -2.0 * to_xyyy_0[k] + 2.0 * to_xxxyyy_0[k] * tbe_0;

                to_x_0_xxyyz_0[k] = -2.0 * to_xyyz_0[k] + 2.0 * to_xxxyyz_0[k] * tbe_0;

                to_x_0_xxyzz_0[k] = -2.0 * to_xyzz_0[k] + 2.0 * to_xxxyzz_0[k] * tbe_0;

                to_x_0_xxzzz_0[k] = -2.0 * to_xzzz_0[k] + 2.0 * to_xxxzzz_0[k] * tbe_0;

                to_x_0_xyyyy_0[k] = -to_yyyy_0[k] + 2.0 * to_xxyyyy_0[k] * tbe_0;

                to_x_0_xyyyz_0[k] = -to_yyyz_0[k] + 2.0 * to_xxyyyz_0[k] * tbe_0;

                to_x_0_xyyzz_0[k] = -to_yyzz_0[k] + 2.0 * to_xxyyzz_0[k] * tbe_0;

                to_x_0_xyzzz_0[k] = -to_yzzz_0[k] + 2.0 * to_xxyzzz_0[k] * tbe_0;

                to_x_0_xzzzz_0[k] = -to_zzzz_0[k] + 2.0 * to_xxzzzz_0[k] * tbe_0;

                to_x_0_yyyyy_0[k] = 2.0 * to_xyyyyy_0[k] * tbe_0;

                to_x_0_yyyyz_0[k] = 2.0 * to_xyyyyz_0[k] * tbe_0;

                to_x_0_yyyzz_0[k] = 2.0 * to_xyyyzz_0[k] * tbe_0;

                to_x_0_yyzzz_0[k] = 2.0 * to_xyyzzz_0[k] * tbe_0;

                to_x_0_yzzzz_0[k] = 2.0 * to_xyzzzz_0[k] * tbe_0;

                to_x_0_zzzzz_0[k] = 2.0 * to_xzzzzz_0[k] * tbe_0;
            }

            // Set up 21-42 components of targeted buffer : HS

            auto to_y_0_xxxxx_0 = prim_buffer.data(idx_op_geom_100_hs + 1 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 0 * ket_comps + j);

            auto to_y_0_xxxxy_0 = prim_buffer.data(idx_op_geom_100_hs + 1 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 1 * ket_comps + j);

            auto to_y_0_xxxxz_0 = prim_buffer.data(idx_op_geom_100_hs + 1 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 2 * ket_comps + j);

            auto to_y_0_xxxyy_0 = prim_buffer.data(idx_op_geom_100_hs + 1 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 3 * ket_comps + j);

            auto to_y_0_xxxyz_0 = prim_buffer.data(idx_op_geom_100_hs + 1 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 4 * ket_comps + j);

            auto to_y_0_xxxzz_0 = prim_buffer.data(idx_op_geom_100_hs + 1 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 5 * ket_comps + j);

            auto to_y_0_xxyyy_0 = prim_buffer.data(idx_op_geom_100_hs + 1 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 6 * ket_comps + j);

            auto to_y_0_xxyyz_0 = prim_buffer.data(idx_op_geom_100_hs + 1 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 7 * ket_comps + j);

            auto to_y_0_xxyzz_0 = prim_buffer.data(idx_op_geom_100_hs + 1 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 8 * ket_comps + j);

            auto to_y_0_xxzzz_0 = prim_buffer.data(idx_op_geom_100_hs + 1 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 9 * ket_comps + j);

            auto to_y_0_xyyyy_0 = prim_buffer.data(idx_op_geom_100_hs + 1 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 10 * ket_comps + j);

            auto to_y_0_xyyyz_0 = prim_buffer.data(idx_op_geom_100_hs + 1 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 11 * ket_comps + j);

            auto to_y_0_xyyzz_0 = prim_buffer.data(idx_op_geom_100_hs + 1 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 12 * ket_comps + j);

            auto to_y_0_xyzzz_0 = prim_buffer.data(idx_op_geom_100_hs + 1 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 13 * ket_comps + j);

            auto to_y_0_xzzzz_0 = prim_buffer.data(idx_op_geom_100_hs + 1 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 14 * ket_comps + j);

            auto to_y_0_yyyyy_0 = prim_buffer.data(idx_op_geom_100_hs + 1 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 15 * ket_comps + j);

            auto to_y_0_yyyyz_0 = prim_buffer.data(idx_op_geom_100_hs + 1 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 16 * ket_comps + j);

            auto to_y_0_yyyzz_0 = prim_buffer.data(idx_op_geom_100_hs + 1 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 17 * ket_comps + j);

            auto to_y_0_yyzzz_0 = prim_buffer.data(idx_op_geom_100_hs + 1 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 18 * ket_comps + j);

            auto to_y_0_yzzzz_0 = prim_buffer.data(idx_op_geom_100_hs + 1 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 19 * ket_comps + j);

            auto to_y_0_zzzzz_0 = prim_buffer.data(idx_op_geom_100_hs + 1 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 20 * ket_comps + j);

            #pragma omp simd aligned(to_xxxx_0, to_xxxxxy_0, to_xxxxyy_0, to_xxxxyz_0, to_xxxy_0, to_xxxyyy_0, to_xxxyyz_0, to_xxxyzz_0, to_xxxz_0, to_xxyy_0, to_xxyyyy_0, to_xxyyyz_0, to_xxyyzz_0, to_xxyz_0, to_xxyzzz_0, to_xxzz_0, to_xyyy_0, to_xyyyyy_0, to_xyyyyz_0, to_xyyyzz_0, to_xyyz_0, to_xyyzzz_0, to_xyzz_0, to_xyzzzz_0, to_xzzz_0, to_y_0_xxxxx_0, to_y_0_xxxxy_0, to_y_0_xxxxz_0, to_y_0_xxxyy_0, to_y_0_xxxyz_0, to_y_0_xxxzz_0, to_y_0_xxyyy_0, to_y_0_xxyyz_0, to_y_0_xxyzz_0, to_y_0_xxzzz_0, to_y_0_xyyyy_0, to_y_0_xyyyz_0, to_y_0_xyyzz_0, to_y_0_xyzzz_0, to_y_0_xzzzz_0, to_y_0_yyyyy_0, to_y_0_yyyyz_0, to_y_0_yyyzz_0, to_y_0_yyzzz_0, to_y_0_yzzzz_0, to_y_0_zzzzz_0, to_yyyy_0, to_yyyyyy_0, to_yyyyyz_0, to_yyyyzz_0, to_yyyz_0, to_yyyzzz_0, to_yyzz_0, to_yyzzzz_0, to_yzzz_0, to_yzzzzz_0, to_zzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                const double tbe_0 = a_exp;

                to_y_0_xxxxx_0[k] = 2.0 * to_xxxxxy_0[k] * tbe_0;

                to_y_0_xxxxy_0[k] = -to_xxxx_0[k] + 2.0 * to_xxxxyy_0[k] * tbe_0;

                to_y_0_xxxxz_0[k] = 2.0 * to_xxxxyz_0[k] * tbe_0;

                to_y_0_xxxyy_0[k] = -2.0 * to_xxxy_0[k] + 2.0 * to_xxxyyy_0[k] * tbe_0;

                to_y_0_xxxyz_0[k] = -to_xxxz_0[k] + 2.0 * to_xxxyyz_0[k] * tbe_0;

                to_y_0_xxxzz_0[k] = 2.0 * to_xxxyzz_0[k] * tbe_0;

                to_y_0_xxyyy_0[k] = -3.0 * to_xxyy_0[k] + 2.0 * to_xxyyyy_0[k] * tbe_0;

                to_y_0_xxyyz_0[k] = -2.0 * to_xxyz_0[k] + 2.0 * to_xxyyyz_0[k] * tbe_0;

                to_y_0_xxyzz_0[k] = -to_xxzz_0[k] + 2.0 * to_xxyyzz_0[k] * tbe_0;

                to_y_0_xxzzz_0[k] = 2.0 * to_xxyzzz_0[k] * tbe_0;

                to_y_0_xyyyy_0[k] = -4.0 * to_xyyy_0[k] + 2.0 * to_xyyyyy_0[k] * tbe_0;

                to_y_0_xyyyz_0[k] = -3.0 * to_xyyz_0[k] + 2.0 * to_xyyyyz_0[k] * tbe_0;

                to_y_0_xyyzz_0[k] = -2.0 * to_xyzz_0[k] + 2.0 * to_xyyyzz_0[k] * tbe_0;

                to_y_0_xyzzz_0[k] = -to_xzzz_0[k] + 2.0 * to_xyyzzz_0[k] * tbe_0;

                to_y_0_xzzzz_0[k] = 2.0 * to_xyzzzz_0[k] * tbe_0;

                to_y_0_yyyyy_0[k] = -5.0 * to_yyyy_0[k] + 2.0 * to_yyyyyy_0[k] * tbe_0;

                to_y_0_yyyyz_0[k] = -4.0 * to_yyyz_0[k] + 2.0 * to_yyyyyz_0[k] * tbe_0;

                to_y_0_yyyzz_0[k] = -3.0 * to_yyzz_0[k] + 2.0 * to_yyyyzz_0[k] * tbe_0;

                to_y_0_yyzzz_0[k] = -2.0 * to_yzzz_0[k] + 2.0 * to_yyyzzz_0[k] * tbe_0;

                to_y_0_yzzzz_0[k] = -to_zzzz_0[k] + 2.0 * to_yyzzzz_0[k] * tbe_0;

                to_y_0_zzzzz_0[k] = 2.0 * to_yzzzzz_0[k] * tbe_0;
            }

            // Set up 42-63 components of targeted buffer : HS

            auto to_z_0_xxxxx_0 = prim_buffer.data(idx_op_geom_100_hs + 2 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 0 * ket_comps + j);

            auto to_z_0_xxxxy_0 = prim_buffer.data(idx_op_geom_100_hs + 2 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 1 * ket_comps + j);

            auto to_z_0_xxxxz_0 = prim_buffer.data(idx_op_geom_100_hs + 2 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 2 * ket_comps + j);

            auto to_z_0_xxxyy_0 = prim_buffer.data(idx_op_geom_100_hs + 2 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 3 * ket_comps + j);

            auto to_z_0_xxxyz_0 = prim_buffer.data(idx_op_geom_100_hs + 2 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 4 * ket_comps + j);

            auto to_z_0_xxxzz_0 = prim_buffer.data(idx_op_geom_100_hs + 2 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 5 * ket_comps + j);

            auto to_z_0_xxyyy_0 = prim_buffer.data(idx_op_geom_100_hs + 2 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 6 * ket_comps + j);

            auto to_z_0_xxyyz_0 = prim_buffer.data(idx_op_geom_100_hs + 2 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 7 * ket_comps + j);

            auto to_z_0_xxyzz_0 = prim_buffer.data(idx_op_geom_100_hs + 2 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 8 * ket_comps + j);

            auto to_z_0_xxzzz_0 = prim_buffer.data(idx_op_geom_100_hs + 2 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 9 * ket_comps + j);

            auto to_z_0_xyyyy_0 = prim_buffer.data(idx_op_geom_100_hs + 2 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 10 * ket_comps + j);

            auto to_z_0_xyyyz_0 = prim_buffer.data(idx_op_geom_100_hs + 2 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 11 * ket_comps + j);

            auto to_z_0_xyyzz_0 = prim_buffer.data(idx_op_geom_100_hs + 2 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 12 * ket_comps + j);

            auto to_z_0_xyzzz_0 = prim_buffer.data(idx_op_geom_100_hs + 2 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 13 * ket_comps + j);

            auto to_z_0_xzzzz_0 = prim_buffer.data(idx_op_geom_100_hs + 2 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 14 * ket_comps + j);

            auto to_z_0_yyyyy_0 = prim_buffer.data(idx_op_geom_100_hs + 2 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 15 * ket_comps + j);

            auto to_z_0_yyyyz_0 = prim_buffer.data(idx_op_geom_100_hs + 2 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 16 * ket_comps + j);

            auto to_z_0_yyyzz_0 = prim_buffer.data(idx_op_geom_100_hs + 2 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 17 * ket_comps + j);

            auto to_z_0_yyzzz_0 = prim_buffer.data(idx_op_geom_100_hs + 2 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 18 * ket_comps + j);

            auto to_z_0_yzzzz_0 = prim_buffer.data(idx_op_geom_100_hs + 2 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 19 * ket_comps + j);

            auto to_z_0_zzzzz_0 = prim_buffer.data(idx_op_geom_100_hs + 2 * op_comps * 21 * ket_comps + i * 21 * ket_comps + 20 * ket_comps + j);

            #pragma omp simd aligned(to_xxxx_0, to_xxxxxz_0, to_xxxxyz_0, to_xxxxzz_0, to_xxxy_0, to_xxxyyz_0, to_xxxyzz_0, to_xxxz_0, to_xxxzzz_0, to_xxyy_0, to_xxyyyz_0, to_xxyyzz_0, to_xxyz_0, to_xxyzzz_0, to_xxzz_0, to_xxzzzz_0, to_xyyy_0, to_xyyyyz_0, to_xyyyzz_0, to_xyyz_0, to_xyyzzz_0, to_xyzz_0, to_xyzzzz_0, to_xzzz_0, to_xzzzzz_0, to_yyyy_0, to_yyyyyz_0, to_yyyyzz_0, to_yyyz_0, to_yyyzzz_0, to_yyzz_0, to_yyzzzz_0, to_yzzz_0, to_yzzzzz_0, to_z_0_xxxxx_0, to_z_0_xxxxy_0, to_z_0_xxxxz_0, to_z_0_xxxyy_0, to_z_0_xxxyz_0, to_z_0_xxxzz_0, to_z_0_xxyyy_0, to_z_0_xxyyz_0, to_z_0_xxyzz_0, to_z_0_xxzzz_0, to_z_0_xyyyy_0, to_z_0_xyyyz_0, to_z_0_xyyzz_0, to_z_0_xyzzz_0, to_z_0_xzzzz_0, to_z_0_yyyyy_0, to_z_0_yyyyz_0, to_z_0_yyyzz_0, to_z_0_yyzzz_0, to_z_0_yzzzz_0, to_z_0_zzzzz_0, to_zzzz_0, to_zzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                const double tbe_0 = a_exp;

                to_z_0_xxxxx_0[k] = 2.0 * to_xxxxxz_0[k] * tbe_0;

                to_z_0_xxxxy_0[k] = 2.0 * to_xxxxyz_0[k] * tbe_0;

                to_z_0_xxxxz_0[k] = -to_xxxx_0[k] + 2.0 * to_xxxxzz_0[k] * tbe_0;

                to_z_0_xxxyy_0[k] = 2.0 * to_xxxyyz_0[k] * tbe_0;

                to_z_0_xxxyz_0[k] = -to_xxxy_0[k] + 2.0 * to_xxxyzz_0[k] * tbe_0;

                to_z_0_xxxzz_0[k] = -2.0 * to_xxxz_0[k] + 2.0 * to_xxxzzz_0[k] * tbe_0;

                to_z_0_xxyyy_0[k] = 2.0 * to_xxyyyz_0[k] * tbe_0;

                to_z_0_xxyyz_0[k] = -to_xxyy_0[k] + 2.0 * to_xxyyzz_0[k] * tbe_0;

                to_z_0_xxyzz_0[k] = -2.0 * to_xxyz_0[k] + 2.0 * to_xxyzzz_0[k] * tbe_0;

                to_z_0_xxzzz_0[k] = -3.0 * to_xxzz_0[k] + 2.0 * to_xxzzzz_0[k] * tbe_0;

                to_z_0_xyyyy_0[k] = 2.0 * to_xyyyyz_0[k] * tbe_0;

                to_z_0_xyyyz_0[k] = -to_xyyy_0[k] + 2.0 * to_xyyyzz_0[k] * tbe_0;

                to_z_0_xyyzz_0[k] = -2.0 * to_xyyz_0[k] + 2.0 * to_xyyzzz_0[k] * tbe_0;

                to_z_0_xyzzz_0[k] = -3.0 * to_xyzz_0[k] + 2.0 * to_xyzzzz_0[k] * tbe_0;

                to_z_0_xzzzz_0[k] = -4.0 * to_xzzz_0[k] + 2.0 * to_xzzzzz_0[k] * tbe_0;

                to_z_0_yyyyy_0[k] = 2.0 * to_yyyyyz_0[k] * tbe_0;

                to_z_0_yyyyz_0[k] = -to_yyyy_0[k] + 2.0 * to_yyyyzz_0[k] * tbe_0;

                to_z_0_yyyzz_0[k] = -2.0 * to_yyyz_0[k] + 2.0 * to_yyyzzz_0[k] * tbe_0;

                to_z_0_yyzzz_0[k] = -3.0 * to_yyzz_0[k] + 2.0 * to_yyzzzz_0[k] * tbe_0;

                to_z_0_yzzzz_0[k] = -4.0 * to_yzzz_0[k] + 2.0 * to_yzzzzz_0[k] * tbe_0;

                to_z_0_zzzzz_0[k] = -5.0 * to_zzzz_0[k] + 2.0 * to_zzzzzz_0[k] * tbe_0;
            }

        }
    }

}

} // t2cgeom namespace

