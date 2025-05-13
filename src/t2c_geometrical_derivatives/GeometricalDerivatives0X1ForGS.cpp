#include "GeometricalDerivatives0X1ForGS.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_geom_deriv_0x1_gs(CSimdArray<double>& pbuffer,
                       const int idx_op_geom_001_gs,
                       const int idx_op_gp,
                       const CSimdArray<double>& factors,
                       const int op_comps) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : GP

        auto to_xxxx_x = pbuffer.data(idx_op_gp + i * 45 + 0);

        auto to_xxxx_y = pbuffer.data(idx_op_gp + i * 45 + 1);

        auto to_xxxx_z = pbuffer.data(idx_op_gp + i * 45 + 2);

        auto to_xxxy_x = pbuffer.data(idx_op_gp + i * 45 + 3);

        auto to_xxxy_y = pbuffer.data(idx_op_gp + i * 45 + 4);

        auto to_xxxy_z = pbuffer.data(idx_op_gp + i * 45 + 5);

        auto to_xxxz_x = pbuffer.data(idx_op_gp + i * 45 + 6);

        auto to_xxxz_y = pbuffer.data(idx_op_gp + i * 45 + 7);

        auto to_xxxz_z = pbuffer.data(idx_op_gp + i * 45 + 8);

        auto to_xxyy_x = pbuffer.data(idx_op_gp + i * 45 + 9);

        auto to_xxyy_y = pbuffer.data(idx_op_gp + i * 45 + 10);

        auto to_xxyy_z = pbuffer.data(idx_op_gp + i * 45 + 11);

        auto to_xxyz_x = pbuffer.data(idx_op_gp + i * 45 + 12);

        auto to_xxyz_y = pbuffer.data(idx_op_gp + i * 45 + 13);

        auto to_xxyz_z = pbuffer.data(idx_op_gp + i * 45 + 14);

        auto to_xxzz_x = pbuffer.data(idx_op_gp + i * 45 + 15);

        auto to_xxzz_y = pbuffer.data(idx_op_gp + i * 45 + 16);

        auto to_xxzz_z = pbuffer.data(idx_op_gp + i * 45 + 17);

        auto to_xyyy_x = pbuffer.data(idx_op_gp + i * 45 + 18);

        auto to_xyyy_y = pbuffer.data(idx_op_gp + i * 45 + 19);

        auto to_xyyy_z = pbuffer.data(idx_op_gp + i * 45 + 20);

        auto to_xyyz_x = pbuffer.data(idx_op_gp + i * 45 + 21);

        auto to_xyyz_y = pbuffer.data(idx_op_gp + i * 45 + 22);

        auto to_xyyz_z = pbuffer.data(idx_op_gp + i * 45 + 23);

        auto to_xyzz_x = pbuffer.data(idx_op_gp + i * 45 + 24);

        auto to_xyzz_y = pbuffer.data(idx_op_gp + i * 45 + 25);

        auto to_xyzz_z = pbuffer.data(idx_op_gp + i * 45 + 26);

        auto to_xzzz_x = pbuffer.data(idx_op_gp + i * 45 + 27);

        auto to_xzzz_y = pbuffer.data(idx_op_gp + i * 45 + 28);

        auto to_xzzz_z = pbuffer.data(idx_op_gp + i * 45 + 29);

        auto to_yyyy_x = pbuffer.data(idx_op_gp + i * 45 + 30);

        auto to_yyyy_y = pbuffer.data(idx_op_gp + i * 45 + 31);

        auto to_yyyy_z = pbuffer.data(idx_op_gp + i * 45 + 32);

        auto to_yyyz_x = pbuffer.data(idx_op_gp + i * 45 + 33);

        auto to_yyyz_y = pbuffer.data(idx_op_gp + i * 45 + 34);

        auto to_yyyz_z = pbuffer.data(idx_op_gp + i * 45 + 35);

        auto to_yyzz_x = pbuffer.data(idx_op_gp + i * 45 + 36);

        auto to_yyzz_y = pbuffer.data(idx_op_gp + i * 45 + 37);

        auto to_yyzz_z = pbuffer.data(idx_op_gp + i * 45 + 38);

        auto to_yzzz_x = pbuffer.data(idx_op_gp + i * 45 + 39);

        auto to_yzzz_y = pbuffer.data(idx_op_gp + i * 45 + 40);

        auto to_yzzz_z = pbuffer.data(idx_op_gp + i * 45 + 41);

        auto to_zzzz_x = pbuffer.data(idx_op_gp + i * 45 + 42);

        auto to_zzzz_y = pbuffer.data(idx_op_gp + i * 45 + 43);

        auto to_zzzz_z = pbuffer.data(idx_op_gp + i * 45 + 44);

        // Set up 0-15 components of targeted buffer : GS

        auto to_0_x_xxxx_0 = pbuffer.data(idx_op_geom_001_gs + 0 * op_comps * 15 + i * 15 + 0);

        auto to_0_x_xxxy_0 = pbuffer.data(idx_op_geom_001_gs + 0 * op_comps * 15 + i * 15 + 1);

        auto to_0_x_xxxz_0 = pbuffer.data(idx_op_geom_001_gs + 0 * op_comps * 15 + i * 15 + 2);

        auto to_0_x_xxyy_0 = pbuffer.data(idx_op_geom_001_gs + 0 * op_comps * 15 + i * 15 + 3);

        auto to_0_x_xxyz_0 = pbuffer.data(idx_op_geom_001_gs + 0 * op_comps * 15 + i * 15 + 4);

        auto to_0_x_xxzz_0 = pbuffer.data(idx_op_geom_001_gs + 0 * op_comps * 15 + i * 15 + 5);

        auto to_0_x_xyyy_0 = pbuffer.data(idx_op_geom_001_gs + 0 * op_comps * 15 + i * 15 + 6);

        auto to_0_x_xyyz_0 = pbuffer.data(idx_op_geom_001_gs + 0 * op_comps * 15 + i * 15 + 7);

        auto to_0_x_xyzz_0 = pbuffer.data(idx_op_geom_001_gs + 0 * op_comps * 15 + i * 15 + 8);

        auto to_0_x_xzzz_0 = pbuffer.data(idx_op_geom_001_gs + 0 * op_comps * 15 + i * 15 + 9);

        auto to_0_x_yyyy_0 = pbuffer.data(idx_op_geom_001_gs + 0 * op_comps * 15 + i * 15 + 10);

        auto to_0_x_yyyz_0 = pbuffer.data(idx_op_geom_001_gs + 0 * op_comps * 15 + i * 15 + 11);

        auto to_0_x_yyzz_0 = pbuffer.data(idx_op_geom_001_gs + 0 * op_comps * 15 + i * 15 + 12);

        auto to_0_x_yzzz_0 = pbuffer.data(idx_op_geom_001_gs + 0 * op_comps * 15 + i * 15 + 13);

        auto to_0_x_zzzz_0 = pbuffer.data(idx_op_geom_001_gs + 0 * op_comps * 15 + i * 15 + 14);

        #pragma omp simd aligned(to_0_x_xxxx_0, to_0_x_xxxy_0, to_0_x_xxxz_0, to_0_x_xxyy_0, to_0_x_xxyz_0, to_0_x_xxzz_0, to_0_x_xyyy_0, to_0_x_xyyz_0, to_0_x_xyzz_0, to_0_x_xzzz_0, to_0_x_yyyy_0, to_0_x_yyyz_0, to_0_x_yyzz_0, to_0_x_yzzz_0, to_0_x_zzzz_0, to_xxxx_x, to_xxxy_x, to_xxxz_x, to_xxyy_x, to_xxyz_x, to_xxzz_x, to_xyyy_x, to_xyyz_x, to_xyzz_x, to_xzzz_x, to_yyyy_x, to_yyyz_x, to_yyzz_x, to_yzzz_x, to_zzzz_x, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxxx_0[k] = 2.0 * to_xxxx_x[k] * tke_0;

            to_0_x_xxxy_0[k] = 2.0 * to_xxxy_x[k] * tke_0;

            to_0_x_xxxz_0[k] = 2.0 * to_xxxz_x[k] * tke_0;

            to_0_x_xxyy_0[k] = 2.0 * to_xxyy_x[k] * tke_0;

            to_0_x_xxyz_0[k] = 2.0 * to_xxyz_x[k] * tke_0;

            to_0_x_xxzz_0[k] = 2.0 * to_xxzz_x[k] * tke_0;

            to_0_x_xyyy_0[k] = 2.0 * to_xyyy_x[k] * tke_0;

            to_0_x_xyyz_0[k] = 2.0 * to_xyyz_x[k] * tke_0;

            to_0_x_xyzz_0[k] = 2.0 * to_xyzz_x[k] * tke_0;

            to_0_x_xzzz_0[k] = 2.0 * to_xzzz_x[k] * tke_0;

            to_0_x_yyyy_0[k] = 2.0 * to_yyyy_x[k] * tke_0;

            to_0_x_yyyz_0[k] = 2.0 * to_yyyz_x[k] * tke_0;

            to_0_x_yyzz_0[k] = 2.0 * to_yyzz_x[k] * tke_0;

            to_0_x_yzzz_0[k] = 2.0 * to_yzzz_x[k] * tke_0;

            to_0_x_zzzz_0[k] = 2.0 * to_zzzz_x[k] * tke_0;
        }

        // Set up 15-30 components of targeted buffer : GS

        auto to_0_y_xxxx_0 = pbuffer.data(idx_op_geom_001_gs + 1 * op_comps * 15 + i * 15 + 0);

        auto to_0_y_xxxy_0 = pbuffer.data(idx_op_geom_001_gs + 1 * op_comps * 15 + i * 15 + 1);

        auto to_0_y_xxxz_0 = pbuffer.data(idx_op_geom_001_gs + 1 * op_comps * 15 + i * 15 + 2);

        auto to_0_y_xxyy_0 = pbuffer.data(idx_op_geom_001_gs + 1 * op_comps * 15 + i * 15 + 3);

        auto to_0_y_xxyz_0 = pbuffer.data(idx_op_geom_001_gs + 1 * op_comps * 15 + i * 15 + 4);

        auto to_0_y_xxzz_0 = pbuffer.data(idx_op_geom_001_gs + 1 * op_comps * 15 + i * 15 + 5);

        auto to_0_y_xyyy_0 = pbuffer.data(idx_op_geom_001_gs + 1 * op_comps * 15 + i * 15 + 6);

        auto to_0_y_xyyz_0 = pbuffer.data(idx_op_geom_001_gs + 1 * op_comps * 15 + i * 15 + 7);

        auto to_0_y_xyzz_0 = pbuffer.data(idx_op_geom_001_gs + 1 * op_comps * 15 + i * 15 + 8);

        auto to_0_y_xzzz_0 = pbuffer.data(idx_op_geom_001_gs + 1 * op_comps * 15 + i * 15 + 9);

        auto to_0_y_yyyy_0 = pbuffer.data(idx_op_geom_001_gs + 1 * op_comps * 15 + i * 15 + 10);

        auto to_0_y_yyyz_0 = pbuffer.data(idx_op_geom_001_gs + 1 * op_comps * 15 + i * 15 + 11);

        auto to_0_y_yyzz_0 = pbuffer.data(idx_op_geom_001_gs + 1 * op_comps * 15 + i * 15 + 12);

        auto to_0_y_yzzz_0 = pbuffer.data(idx_op_geom_001_gs + 1 * op_comps * 15 + i * 15 + 13);

        auto to_0_y_zzzz_0 = pbuffer.data(idx_op_geom_001_gs + 1 * op_comps * 15 + i * 15 + 14);

        #pragma omp simd aligned(to_0_y_xxxx_0, to_0_y_xxxy_0, to_0_y_xxxz_0, to_0_y_xxyy_0, to_0_y_xxyz_0, to_0_y_xxzz_0, to_0_y_xyyy_0, to_0_y_xyyz_0, to_0_y_xyzz_0, to_0_y_xzzz_0, to_0_y_yyyy_0, to_0_y_yyyz_0, to_0_y_yyzz_0, to_0_y_yzzz_0, to_0_y_zzzz_0, to_xxxx_y, to_xxxy_y, to_xxxz_y, to_xxyy_y, to_xxyz_y, to_xxzz_y, to_xyyy_y, to_xyyz_y, to_xyzz_y, to_xzzz_y, to_yyyy_y, to_yyyz_y, to_yyzz_y, to_yzzz_y, to_zzzz_y, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxxx_0[k] = 2.0 * to_xxxx_y[k] * tke_0;

            to_0_y_xxxy_0[k] = 2.0 * to_xxxy_y[k] * tke_0;

            to_0_y_xxxz_0[k] = 2.0 * to_xxxz_y[k] * tke_0;

            to_0_y_xxyy_0[k] = 2.0 * to_xxyy_y[k] * tke_0;

            to_0_y_xxyz_0[k] = 2.0 * to_xxyz_y[k] * tke_0;

            to_0_y_xxzz_0[k] = 2.0 * to_xxzz_y[k] * tke_0;

            to_0_y_xyyy_0[k] = 2.0 * to_xyyy_y[k] * tke_0;

            to_0_y_xyyz_0[k] = 2.0 * to_xyyz_y[k] * tke_0;

            to_0_y_xyzz_0[k] = 2.0 * to_xyzz_y[k] * tke_0;

            to_0_y_xzzz_0[k] = 2.0 * to_xzzz_y[k] * tke_0;

            to_0_y_yyyy_0[k] = 2.0 * to_yyyy_y[k] * tke_0;

            to_0_y_yyyz_0[k] = 2.0 * to_yyyz_y[k] * tke_0;

            to_0_y_yyzz_0[k] = 2.0 * to_yyzz_y[k] * tke_0;

            to_0_y_yzzz_0[k] = 2.0 * to_yzzz_y[k] * tke_0;

            to_0_y_zzzz_0[k] = 2.0 * to_zzzz_y[k] * tke_0;
        }

        // Set up 30-45 components of targeted buffer : GS

        auto to_0_z_xxxx_0 = pbuffer.data(idx_op_geom_001_gs + 2 * op_comps * 15 + i * 15 + 0);

        auto to_0_z_xxxy_0 = pbuffer.data(idx_op_geom_001_gs + 2 * op_comps * 15 + i * 15 + 1);

        auto to_0_z_xxxz_0 = pbuffer.data(idx_op_geom_001_gs + 2 * op_comps * 15 + i * 15 + 2);

        auto to_0_z_xxyy_0 = pbuffer.data(idx_op_geom_001_gs + 2 * op_comps * 15 + i * 15 + 3);

        auto to_0_z_xxyz_0 = pbuffer.data(idx_op_geom_001_gs + 2 * op_comps * 15 + i * 15 + 4);

        auto to_0_z_xxzz_0 = pbuffer.data(idx_op_geom_001_gs + 2 * op_comps * 15 + i * 15 + 5);

        auto to_0_z_xyyy_0 = pbuffer.data(idx_op_geom_001_gs + 2 * op_comps * 15 + i * 15 + 6);

        auto to_0_z_xyyz_0 = pbuffer.data(idx_op_geom_001_gs + 2 * op_comps * 15 + i * 15 + 7);

        auto to_0_z_xyzz_0 = pbuffer.data(idx_op_geom_001_gs + 2 * op_comps * 15 + i * 15 + 8);

        auto to_0_z_xzzz_0 = pbuffer.data(idx_op_geom_001_gs + 2 * op_comps * 15 + i * 15 + 9);

        auto to_0_z_yyyy_0 = pbuffer.data(idx_op_geom_001_gs + 2 * op_comps * 15 + i * 15 + 10);

        auto to_0_z_yyyz_0 = pbuffer.data(idx_op_geom_001_gs + 2 * op_comps * 15 + i * 15 + 11);

        auto to_0_z_yyzz_0 = pbuffer.data(idx_op_geom_001_gs + 2 * op_comps * 15 + i * 15 + 12);

        auto to_0_z_yzzz_0 = pbuffer.data(idx_op_geom_001_gs + 2 * op_comps * 15 + i * 15 + 13);

        auto to_0_z_zzzz_0 = pbuffer.data(idx_op_geom_001_gs + 2 * op_comps * 15 + i * 15 + 14);

        #pragma omp simd aligned(to_0_z_xxxx_0, to_0_z_xxxy_0, to_0_z_xxxz_0, to_0_z_xxyy_0, to_0_z_xxyz_0, to_0_z_xxzz_0, to_0_z_xyyy_0, to_0_z_xyyz_0, to_0_z_xyzz_0, to_0_z_xzzz_0, to_0_z_yyyy_0, to_0_z_yyyz_0, to_0_z_yyzz_0, to_0_z_yzzz_0, to_0_z_zzzz_0, to_xxxx_z, to_xxxy_z, to_xxxz_z, to_xxyy_z, to_xxyz_z, to_xxzz_z, to_xyyy_z, to_xyyz_z, to_xyzz_z, to_xzzz_z, to_yyyy_z, to_yyyz_z, to_yyzz_z, to_yzzz_z, to_zzzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxxx_0[k] = 2.0 * to_xxxx_z[k] * tke_0;

            to_0_z_xxxy_0[k] = 2.0 * to_xxxy_z[k] * tke_0;

            to_0_z_xxxz_0[k] = 2.0 * to_xxxz_z[k] * tke_0;

            to_0_z_xxyy_0[k] = 2.0 * to_xxyy_z[k] * tke_0;

            to_0_z_xxyz_0[k] = 2.0 * to_xxyz_z[k] * tke_0;

            to_0_z_xxzz_0[k] = 2.0 * to_xxzz_z[k] * tke_0;

            to_0_z_xyyy_0[k] = 2.0 * to_xyyy_z[k] * tke_0;

            to_0_z_xyyz_0[k] = 2.0 * to_xyyz_z[k] * tke_0;

            to_0_z_xyzz_0[k] = 2.0 * to_xyzz_z[k] * tke_0;

            to_0_z_xzzz_0[k] = 2.0 * to_xzzz_z[k] * tke_0;

            to_0_z_yyyy_0[k] = 2.0 * to_yyyy_z[k] * tke_0;

            to_0_z_yyyz_0[k] = 2.0 * to_yyyz_z[k] * tke_0;

            to_0_z_yyzz_0[k] = 2.0 * to_yyzz_z[k] * tke_0;

            to_0_z_yzzz_0[k] = 2.0 * to_yzzz_z[k] * tke_0;

            to_0_z_zzzz_0[k] = 2.0 * to_zzzz_z[k] * tke_0;
        }

    }

}

} // t2cgeom namespace

