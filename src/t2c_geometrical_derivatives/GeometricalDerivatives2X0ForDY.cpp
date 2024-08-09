#include "GeometricalDerivatives2X0ForDY.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_20_dx(CSimdArray<double>& pbuffer,
                        const size_t idx_op_geom_200_ds,
                        const size_t idx_op_ss,
                        const size_t idx_op_ds,
                        const size_t idx_op_gs,
                        const size_t op_comps,
                        const size_t ket_comps,
                        const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    for (size_t i = 0; i < op_comps; i++)
    {
        for (size_t j = 0; j < ket_comps; j++)
        {
            // Set up components of auxiliary buffer : SS

            auto to_0_0 = pbuffer.data(idx_op_ss + i * 1 * ket_comps + 0 * ket_comps + j);

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

            // Set up 0-6 components of targeted buffer : DS

            auto to_xx_0_xx_0 = pbuffer.data(idx_op_geom_200_ds + 0 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 0 * ket_comps + j);

            auto to_xx_0_xy_0 = pbuffer.data(idx_op_geom_200_ds + 0 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 1 * ket_comps + j);

            auto to_xx_0_xz_0 = pbuffer.data(idx_op_geom_200_ds + 0 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 2 * ket_comps + j);

            auto to_xx_0_yy_0 = pbuffer.data(idx_op_geom_200_ds + 0 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 3 * ket_comps + j);

            auto to_xx_0_yz_0 = pbuffer.data(idx_op_geom_200_ds + 0 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 4 * ket_comps + j);

            auto to_xx_0_zz_0 = pbuffer.data(idx_op_geom_200_ds + 0 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 5 * ket_comps + j);

            #pragma omp simd aligned(to_0_0, to_xx_0, to_xx_0_xx_0, to_xx_0_xy_0, to_xx_0_xz_0, to_xx_0_yy_0, to_xx_0_yz_0, to_xx_0_zz_0, to_xxxx_0, to_xxxy_0, to_xxxz_0, to_xxyy_0, to_xxyz_0, to_xxzz_0, to_xy_0, to_xz_0, to_yy_0, to_yz_0, to_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                const double tbe_0 = a_exp;

                to_xx_0_xx_0[k] = 2.0 * to_0_0[k] - 10.0 * to_xx_0[k] * tbe_0 + 4.0 * to_xxxx_0[k] * tbe_0 * tbe_0;

                to_xx_0_xy_0[k] = -6.0 * to_xy_0[k] * tbe_0 + 4.0 * to_xxxy_0[k] * tbe_0 * tbe_0;

                to_xx_0_xz_0[k] = -6.0 * to_xz_0[k] * tbe_0 + 4.0 * to_xxxz_0[k] * tbe_0 * tbe_0;

                to_xx_0_yy_0[k] = -2.0 * to_yy_0[k] * tbe_0 + 4.0 * to_xxyy_0[k] * tbe_0 * tbe_0;

                to_xx_0_yz_0[k] = -2.0 * to_yz_0[k] * tbe_0 + 4.0 * to_xxyz_0[k] * tbe_0 * tbe_0;

                to_xx_0_zz_0[k] = -2.0 * to_zz_0[k] * tbe_0 + 4.0 * to_xxzz_0[k] * tbe_0 * tbe_0;
            }

            // Set up 6-12 components of targeted buffer : DS

            auto to_xy_0_xx_0 = pbuffer.data(idx_op_geom_200_ds + 1 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 0 * ket_comps + j);

            auto to_xy_0_xy_0 = pbuffer.data(idx_op_geom_200_ds + 1 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 1 * ket_comps + j);

            auto to_xy_0_xz_0 = pbuffer.data(idx_op_geom_200_ds + 1 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 2 * ket_comps + j);

            auto to_xy_0_yy_0 = pbuffer.data(idx_op_geom_200_ds + 1 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 3 * ket_comps + j);

            auto to_xy_0_yz_0 = pbuffer.data(idx_op_geom_200_ds + 1 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 4 * ket_comps + j);

            auto to_xy_0_zz_0 = pbuffer.data(idx_op_geom_200_ds + 1 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 5 * ket_comps + j);

            #pragma omp simd aligned(to_0_0, to_xx_0, to_xxxy_0, to_xxyy_0, to_xxyz_0, to_xy_0, to_xy_0_xx_0, to_xy_0_xy_0, to_xy_0_xz_0, to_xy_0_yy_0, to_xy_0_yz_0, to_xy_0_zz_0, to_xyyy_0, to_xyyz_0, to_xyzz_0, to_xz_0, to_yy_0, to_yz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                const double tbe_0 = a_exp;

                to_xy_0_xx_0[k] = -4.0 * to_xy_0[k] * tbe_0 + 4.0 * to_xxxy_0[k] * tbe_0 * tbe_0;

                to_xy_0_xy_0[k] = to_0_0[k] - 2.0 * to_yy_0[k] * tbe_0 - 2.0 * to_xx_0[k] * tbe_0 + 4.0 * to_xxyy_0[k] * tbe_0 * tbe_0;

                to_xy_0_xz_0[k] = -2.0 * to_yz_0[k] * tbe_0 + 4.0 * to_xxyz_0[k] * tbe_0 * tbe_0;

                to_xy_0_yy_0[k] = -4.0 * to_xy_0[k] * tbe_0 + 4.0 * to_xyyy_0[k] * tbe_0 * tbe_0;

                to_xy_0_yz_0[k] = -2.0 * to_xz_0[k] * tbe_0 + 4.0 * to_xyyz_0[k] * tbe_0 * tbe_0;

                to_xy_0_zz_0[k] = 4.0 * to_xyzz_0[k] * tbe_0 * tbe_0;
            }

            // Set up 12-18 components of targeted buffer : DS

            auto to_xz_0_xx_0 = pbuffer.data(idx_op_geom_200_ds + 2 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 0 * ket_comps + j);

            auto to_xz_0_xy_0 = pbuffer.data(idx_op_geom_200_ds + 2 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 1 * ket_comps + j);

            auto to_xz_0_xz_0 = pbuffer.data(idx_op_geom_200_ds + 2 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 2 * ket_comps + j);

            auto to_xz_0_yy_0 = pbuffer.data(idx_op_geom_200_ds + 2 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 3 * ket_comps + j);

            auto to_xz_0_yz_0 = pbuffer.data(idx_op_geom_200_ds + 2 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 4 * ket_comps + j);

            auto to_xz_0_zz_0 = pbuffer.data(idx_op_geom_200_ds + 2 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 5 * ket_comps + j);

            #pragma omp simd aligned(to_0_0, to_xx_0, to_xxxz_0, to_xxyz_0, to_xxzz_0, to_xy_0, to_xyyz_0, to_xyzz_0, to_xz_0, to_xz_0_xx_0, to_xz_0_xy_0, to_xz_0_xz_0, to_xz_0_yy_0, to_xz_0_yz_0, to_xz_0_zz_0, to_xzzz_0, to_yz_0, to_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                const double tbe_0 = a_exp;

                to_xz_0_xx_0[k] = -4.0 * to_xz_0[k] * tbe_0 + 4.0 * to_xxxz_0[k] * tbe_0 * tbe_0;

                to_xz_0_xy_0[k] = -2.0 * to_yz_0[k] * tbe_0 + 4.0 * to_xxyz_0[k] * tbe_0 * tbe_0;

                to_xz_0_xz_0[k] = to_0_0[k] - 2.0 * to_zz_0[k] * tbe_0 - 2.0 * to_xx_0[k] * tbe_0 + 4.0 * to_xxzz_0[k] * tbe_0 * tbe_0;

                to_xz_0_yy_0[k] = 4.0 * to_xyyz_0[k] * tbe_0 * tbe_0;

                to_xz_0_yz_0[k] = -2.0 * to_xy_0[k] * tbe_0 + 4.0 * to_xyzz_0[k] * tbe_0 * tbe_0;

                to_xz_0_zz_0[k] = -4.0 * to_xz_0[k] * tbe_0 + 4.0 * to_xzzz_0[k] * tbe_0 * tbe_0;
            }

            // Set up 18-24 components of targeted buffer : DS

            auto to_yy_0_xx_0 = pbuffer.data(idx_op_geom_200_ds + 3 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 0 * ket_comps + j);

            auto to_yy_0_xy_0 = pbuffer.data(idx_op_geom_200_ds + 3 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 1 * ket_comps + j);

            auto to_yy_0_xz_0 = pbuffer.data(idx_op_geom_200_ds + 3 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 2 * ket_comps + j);

            auto to_yy_0_yy_0 = pbuffer.data(idx_op_geom_200_ds + 3 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 3 * ket_comps + j);

            auto to_yy_0_yz_0 = pbuffer.data(idx_op_geom_200_ds + 3 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 4 * ket_comps + j);

            auto to_yy_0_zz_0 = pbuffer.data(idx_op_geom_200_ds + 3 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 5 * ket_comps + j);

            #pragma omp simd aligned(to_0_0, to_xx_0, to_xxyy_0, to_xy_0, to_xyyy_0, to_xyyz_0, to_xz_0, to_yy_0, to_yy_0_xx_0, to_yy_0_xy_0, to_yy_0_xz_0, to_yy_0_yy_0, to_yy_0_yz_0, to_yy_0_zz_0, to_yyyy_0, to_yyyz_0, to_yyzz_0, to_yz_0, to_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                const double tbe_0 = a_exp;

                to_yy_0_xx_0[k] = -2.0 * to_xx_0[k] * tbe_0 + 4.0 * to_xxyy_0[k] * tbe_0 * tbe_0;

                to_yy_0_xy_0[k] = -6.0 * to_xy_0[k] * tbe_0 + 4.0 * to_xyyy_0[k] * tbe_0 * tbe_0;

                to_yy_0_xz_0[k] = -2.0 * to_xz_0[k] * tbe_0 + 4.0 * to_xyyz_0[k] * tbe_0 * tbe_0;

                to_yy_0_yy_0[k] = 2.0 * to_0_0[k] - 10.0 * to_yy_0[k] * tbe_0 + 4.0 * to_yyyy_0[k] * tbe_0 * tbe_0;

                to_yy_0_yz_0[k] = -6.0 * to_yz_0[k] * tbe_0 + 4.0 * to_yyyz_0[k] * tbe_0 * tbe_0;

                to_yy_0_zz_0[k] = -2.0 * to_zz_0[k] * tbe_0 + 4.0 * to_yyzz_0[k] * tbe_0 * tbe_0;
            }

            // Set up 24-30 components of targeted buffer : DS

            auto to_yz_0_xx_0 = pbuffer.data(idx_op_geom_200_ds + 4 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 0 * ket_comps + j);

            auto to_yz_0_xy_0 = pbuffer.data(idx_op_geom_200_ds + 4 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 1 * ket_comps + j);

            auto to_yz_0_xz_0 = pbuffer.data(idx_op_geom_200_ds + 4 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 2 * ket_comps + j);

            auto to_yz_0_yy_0 = pbuffer.data(idx_op_geom_200_ds + 4 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 3 * ket_comps + j);

            auto to_yz_0_yz_0 = pbuffer.data(idx_op_geom_200_ds + 4 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 4 * ket_comps + j);

            auto to_yz_0_zz_0 = pbuffer.data(idx_op_geom_200_ds + 4 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 5 * ket_comps + j);

            #pragma omp simd aligned(to_0_0, to_xxyz_0, to_xy_0, to_xyyz_0, to_xyzz_0, to_xz_0, to_yy_0, to_yyyz_0, to_yyzz_0, to_yz_0, to_yz_0_xx_0, to_yz_0_xy_0, to_yz_0_xz_0, to_yz_0_yy_0, to_yz_0_yz_0, to_yz_0_zz_0, to_yzzz_0, to_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                const double tbe_0 = a_exp;

                to_yz_0_xx_0[k] = 4.0 * to_xxyz_0[k] * tbe_0 * tbe_0;

                to_yz_0_xy_0[k] = -2.0 * to_xz_0[k] * tbe_0 + 4.0 * to_xyyz_0[k] * tbe_0 * tbe_0;

                to_yz_0_xz_0[k] = -2.0 * to_xy_0[k] * tbe_0 + 4.0 * to_xyzz_0[k] * tbe_0 * tbe_0;

                to_yz_0_yy_0[k] = -4.0 * to_yz_0[k] * tbe_0 + 4.0 * to_yyyz_0[k] * tbe_0 * tbe_0;

                to_yz_0_yz_0[k] = to_0_0[k] - 2.0 * to_zz_0[k] * tbe_0 - 2.0 * to_yy_0[k] * tbe_0 + 4.0 * to_yyzz_0[k] * tbe_0 * tbe_0;

                to_yz_0_zz_0[k] = -4.0 * to_yz_0[k] * tbe_0 + 4.0 * to_yzzz_0[k] * tbe_0 * tbe_0;
            }

            // Set up 30-36 components of targeted buffer : DS

            auto to_zz_0_xx_0 = pbuffer.data(idx_op_geom_200_ds + 5 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 0 * ket_comps + j);

            auto to_zz_0_xy_0 = pbuffer.data(idx_op_geom_200_ds + 5 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 1 * ket_comps + j);

            auto to_zz_0_xz_0 = pbuffer.data(idx_op_geom_200_ds + 5 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 2 * ket_comps + j);

            auto to_zz_0_yy_0 = pbuffer.data(idx_op_geom_200_ds + 5 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 3 * ket_comps + j);

            auto to_zz_0_yz_0 = pbuffer.data(idx_op_geom_200_ds + 5 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 4 * ket_comps + j);

            auto to_zz_0_zz_0 = pbuffer.data(idx_op_geom_200_ds + 5 * op_comps * 6 * ket_comps + i * 6 * ket_comps + 5 * ket_comps + j);

            #pragma omp simd aligned(to_0_0, to_xx_0, to_xxzz_0, to_xy_0, to_xyzz_0, to_xz_0, to_xzzz_0, to_yy_0, to_yyzz_0, to_yz_0, to_yzzz_0, to_zz_0, to_zz_0_xx_0, to_zz_0_xy_0, to_zz_0_xz_0, to_zz_0_yy_0, to_zz_0_yz_0, to_zz_0_zz_0, to_zzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                const double tbe_0 = a_exp;

                to_zz_0_xx_0[k] = -2.0 * to_xx_0[k] * tbe_0 + 4.0 * to_xxzz_0[k] * tbe_0 * tbe_0;

                to_zz_0_xy_0[k] = -2.0 * to_xy_0[k] * tbe_0 + 4.0 * to_xyzz_0[k] * tbe_0 * tbe_0;

                to_zz_0_xz_0[k] = -6.0 * to_xz_0[k] * tbe_0 + 4.0 * to_xzzz_0[k] * tbe_0 * tbe_0;

                to_zz_0_yy_0[k] = -2.0 * to_yy_0[k] * tbe_0 + 4.0 * to_yyzz_0[k] * tbe_0 * tbe_0;

                to_zz_0_yz_0[k] = -6.0 * to_yz_0[k] * tbe_0 + 4.0 * to_yzzz_0[k] * tbe_0 * tbe_0;

                to_zz_0_zz_0[k] = 2.0 * to_0_0[k] - 10.0 * to_zz_0[k] * tbe_0 + 4.0 * to_zzzz_0[k] * tbe_0 * tbe_0;
            }

        }
    }

}

} // t2cgeom namespace

