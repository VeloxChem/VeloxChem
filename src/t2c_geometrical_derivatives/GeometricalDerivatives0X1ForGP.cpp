#include "GeometricalDerivatives0X1ForGP.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_geom_deriv_0x1_gp(CSimdArray<double>& pbuffer,
                       const int idx_op_geom_001_gp,
                       const int idx_op_gs,
                       const int idx_op_gd,
                       const CSimdArray<double>& factors,
                       const int op_comps) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : GS

        auto to_xxxx_0 = pbuffer.data(idx_op_gs + i * 15 + 0);

        auto to_xxxy_0 = pbuffer.data(idx_op_gs + i * 15 + 1);

        auto to_xxxz_0 = pbuffer.data(idx_op_gs + i * 15 + 2);

        auto to_xxyy_0 = pbuffer.data(idx_op_gs + i * 15 + 3);

        auto to_xxyz_0 = pbuffer.data(idx_op_gs + i * 15 + 4);

        auto to_xxzz_0 = pbuffer.data(idx_op_gs + i * 15 + 5);

        auto to_xyyy_0 = pbuffer.data(idx_op_gs + i * 15 + 6);

        auto to_xyyz_0 = pbuffer.data(idx_op_gs + i * 15 + 7);

        auto to_xyzz_0 = pbuffer.data(idx_op_gs + i * 15 + 8);

        auto to_xzzz_0 = pbuffer.data(idx_op_gs + i * 15 + 9);

        auto to_yyyy_0 = pbuffer.data(idx_op_gs + i * 15 + 10);

        auto to_yyyz_0 = pbuffer.data(idx_op_gs + i * 15 + 11);

        auto to_yyzz_0 = pbuffer.data(idx_op_gs + i * 15 + 12);

        auto to_yzzz_0 = pbuffer.data(idx_op_gs + i * 15 + 13);

        auto to_zzzz_0 = pbuffer.data(idx_op_gs + i * 15 + 14);

        // Set up components of auxiliary buffer : GD

        auto to_xxxx_xx = pbuffer.data(idx_op_gd + i * 90 + 0);

        auto to_xxxx_xy = pbuffer.data(idx_op_gd + i * 90 + 1);

        auto to_xxxx_xz = pbuffer.data(idx_op_gd + i * 90 + 2);

        auto to_xxxx_yy = pbuffer.data(idx_op_gd + i * 90 + 3);

        auto to_xxxx_yz = pbuffer.data(idx_op_gd + i * 90 + 4);

        auto to_xxxx_zz = pbuffer.data(idx_op_gd + i * 90 + 5);

        auto to_xxxy_xx = pbuffer.data(idx_op_gd + i * 90 + 6);

        auto to_xxxy_xy = pbuffer.data(idx_op_gd + i * 90 + 7);

        auto to_xxxy_xz = pbuffer.data(idx_op_gd + i * 90 + 8);

        auto to_xxxy_yy = pbuffer.data(idx_op_gd + i * 90 + 9);

        auto to_xxxy_yz = pbuffer.data(idx_op_gd + i * 90 + 10);

        auto to_xxxy_zz = pbuffer.data(idx_op_gd + i * 90 + 11);

        auto to_xxxz_xx = pbuffer.data(idx_op_gd + i * 90 + 12);

        auto to_xxxz_xy = pbuffer.data(idx_op_gd + i * 90 + 13);

        auto to_xxxz_xz = pbuffer.data(idx_op_gd + i * 90 + 14);

        auto to_xxxz_yy = pbuffer.data(idx_op_gd + i * 90 + 15);

        auto to_xxxz_yz = pbuffer.data(idx_op_gd + i * 90 + 16);

        auto to_xxxz_zz = pbuffer.data(idx_op_gd + i * 90 + 17);

        auto to_xxyy_xx = pbuffer.data(idx_op_gd + i * 90 + 18);

        auto to_xxyy_xy = pbuffer.data(idx_op_gd + i * 90 + 19);

        auto to_xxyy_xz = pbuffer.data(idx_op_gd + i * 90 + 20);

        auto to_xxyy_yy = pbuffer.data(idx_op_gd + i * 90 + 21);

        auto to_xxyy_yz = pbuffer.data(idx_op_gd + i * 90 + 22);

        auto to_xxyy_zz = pbuffer.data(idx_op_gd + i * 90 + 23);

        auto to_xxyz_xx = pbuffer.data(idx_op_gd + i * 90 + 24);

        auto to_xxyz_xy = pbuffer.data(idx_op_gd + i * 90 + 25);

        auto to_xxyz_xz = pbuffer.data(idx_op_gd + i * 90 + 26);

        auto to_xxyz_yy = pbuffer.data(idx_op_gd + i * 90 + 27);

        auto to_xxyz_yz = pbuffer.data(idx_op_gd + i * 90 + 28);

        auto to_xxyz_zz = pbuffer.data(idx_op_gd + i * 90 + 29);

        auto to_xxzz_xx = pbuffer.data(idx_op_gd + i * 90 + 30);

        auto to_xxzz_xy = pbuffer.data(idx_op_gd + i * 90 + 31);

        auto to_xxzz_xz = pbuffer.data(idx_op_gd + i * 90 + 32);

        auto to_xxzz_yy = pbuffer.data(idx_op_gd + i * 90 + 33);

        auto to_xxzz_yz = pbuffer.data(idx_op_gd + i * 90 + 34);

        auto to_xxzz_zz = pbuffer.data(idx_op_gd + i * 90 + 35);

        auto to_xyyy_xx = pbuffer.data(idx_op_gd + i * 90 + 36);

        auto to_xyyy_xy = pbuffer.data(idx_op_gd + i * 90 + 37);

        auto to_xyyy_xz = pbuffer.data(idx_op_gd + i * 90 + 38);

        auto to_xyyy_yy = pbuffer.data(idx_op_gd + i * 90 + 39);

        auto to_xyyy_yz = pbuffer.data(idx_op_gd + i * 90 + 40);

        auto to_xyyy_zz = pbuffer.data(idx_op_gd + i * 90 + 41);

        auto to_xyyz_xx = pbuffer.data(idx_op_gd + i * 90 + 42);

        auto to_xyyz_xy = pbuffer.data(idx_op_gd + i * 90 + 43);

        auto to_xyyz_xz = pbuffer.data(idx_op_gd + i * 90 + 44);

        auto to_xyyz_yy = pbuffer.data(idx_op_gd + i * 90 + 45);

        auto to_xyyz_yz = pbuffer.data(idx_op_gd + i * 90 + 46);

        auto to_xyyz_zz = pbuffer.data(idx_op_gd + i * 90 + 47);

        auto to_xyzz_xx = pbuffer.data(idx_op_gd + i * 90 + 48);

        auto to_xyzz_xy = pbuffer.data(idx_op_gd + i * 90 + 49);

        auto to_xyzz_xz = pbuffer.data(idx_op_gd + i * 90 + 50);

        auto to_xyzz_yy = pbuffer.data(idx_op_gd + i * 90 + 51);

        auto to_xyzz_yz = pbuffer.data(idx_op_gd + i * 90 + 52);

        auto to_xyzz_zz = pbuffer.data(idx_op_gd + i * 90 + 53);

        auto to_xzzz_xx = pbuffer.data(idx_op_gd + i * 90 + 54);

        auto to_xzzz_xy = pbuffer.data(idx_op_gd + i * 90 + 55);

        auto to_xzzz_xz = pbuffer.data(idx_op_gd + i * 90 + 56);

        auto to_xzzz_yy = pbuffer.data(idx_op_gd + i * 90 + 57);

        auto to_xzzz_yz = pbuffer.data(idx_op_gd + i * 90 + 58);

        auto to_xzzz_zz = pbuffer.data(idx_op_gd + i * 90 + 59);

        auto to_yyyy_xx = pbuffer.data(idx_op_gd + i * 90 + 60);

        auto to_yyyy_xy = pbuffer.data(idx_op_gd + i * 90 + 61);

        auto to_yyyy_xz = pbuffer.data(idx_op_gd + i * 90 + 62);

        auto to_yyyy_yy = pbuffer.data(idx_op_gd + i * 90 + 63);

        auto to_yyyy_yz = pbuffer.data(idx_op_gd + i * 90 + 64);

        auto to_yyyy_zz = pbuffer.data(idx_op_gd + i * 90 + 65);

        auto to_yyyz_xx = pbuffer.data(idx_op_gd + i * 90 + 66);

        auto to_yyyz_xy = pbuffer.data(idx_op_gd + i * 90 + 67);

        auto to_yyyz_xz = pbuffer.data(idx_op_gd + i * 90 + 68);

        auto to_yyyz_yy = pbuffer.data(idx_op_gd + i * 90 + 69);

        auto to_yyyz_yz = pbuffer.data(idx_op_gd + i * 90 + 70);

        auto to_yyyz_zz = pbuffer.data(idx_op_gd + i * 90 + 71);

        auto to_yyzz_xx = pbuffer.data(idx_op_gd + i * 90 + 72);

        auto to_yyzz_xy = pbuffer.data(idx_op_gd + i * 90 + 73);

        auto to_yyzz_xz = pbuffer.data(idx_op_gd + i * 90 + 74);

        auto to_yyzz_yy = pbuffer.data(idx_op_gd + i * 90 + 75);

        auto to_yyzz_yz = pbuffer.data(idx_op_gd + i * 90 + 76);

        auto to_yyzz_zz = pbuffer.data(idx_op_gd + i * 90 + 77);

        auto to_yzzz_xx = pbuffer.data(idx_op_gd + i * 90 + 78);

        auto to_yzzz_xy = pbuffer.data(idx_op_gd + i * 90 + 79);

        auto to_yzzz_xz = pbuffer.data(idx_op_gd + i * 90 + 80);

        auto to_yzzz_yy = pbuffer.data(idx_op_gd + i * 90 + 81);

        auto to_yzzz_yz = pbuffer.data(idx_op_gd + i * 90 + 82);

        auto to_yzzz_zz = pbuffer.data(idx_op_gd + i * 90 + 83);

        auto to_zzzz_xx = pbuffer.data(idx_op_gd + i * 90 + 84);

        auto to_zzzz_xy = pbuffer.data(idx_op_gd + i * 90 + 85);

        auto to_zzzz_xz = pbuffer.data(idx_op_gd + i * 90 + 86);

        auto to_zzzz_yy = pbuffer.data(idx_op_gd + i * 90 + 87);

        auto to_zzzz_yz = pbuffer.data(idx_op_gd + i * 90 + 88);

        auto to_zzzz_zz = pbuffer.data(idx_op_gd + i * 90 + 89);

        // Set up 0-3 components of targeted buffer : GP

        auto to_0_x_xxxx_x = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 0);

        auto to_0_x_xxxx_y = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 1);

        auto to_0_x_xxxx_z = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 2);

        #pragma omp simd aligned(to_0_x_xxxx_x, to_0_x_xxxx_y, to_0_x_xxxx_z, to_xxxx_0, to_xxxx_xx, to_xxxx_xy, to_xxxx_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxxx_x[k] = -to_xxxx_0[k] + 2.0 * to_xxxx_xx[k] * tke_0;

            to_0_x_xxxx_y[k] = 2.0 * to_xxxx_xy[k] * tke_0;

            to_0_x_xxxx_z[k] = 2.0 * to_xxxx_xz[k] * tke_0;
        }

        // Set up 3-6 components of targeted buffer : GP

        auto to_0_x_xxxy_x = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 3);

        auto to_0_x_xxxy_y = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 4);

        auto to_0_x_xxxy_z = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 5);

        #pragma omp simd aligned(to_0_x_xxxy_x, to_0_x_xxxy_y, to_0_x_xxxy_z, to_xxxy_0, to_xxxy_xx, to_xxxy_xy, to_xxxy_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxxy_x[k] = -to_xxxy_0[k] + 2.0 * to_xxxy_xx[k] * tke_0;

            to_0_x_xxxy_y[k] = 2.0 * to_xxxy_xy[k] * tke_0;

            to_0_x_xxxy_z[k] = 2.0 * to_xxxy_xz[k] * tke_0;
        }

        // Set up 6-9 components of targeted buffer : GP

        auto to_0_x_xxxz_x = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 6);

        auto to_0_x_xxxz_y = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 7);

        auto to_0_x_xxxz_z = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 8);

        #pragma omp simd aligned(to_0_x_xxxz_x, to_0_x_xxxz_y, to_0_x_xxxz_z, to_xxxz_0, to_xxxz_xx, to_xxxz_xy, to_xxxz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxxz_x[k] = -to_xxxz_0[k] + 2.0 * to_xxxz_xx[k] * tke_0;

            to_0_x_xxxz_y[k] = 2.0 * to_xxxz_xy[k] * tke_0;

            to_0_x_xxxz_z[k] = 2.0 * to_xxxz_xz[k] * tke_0;
        }

        // Set up 9-12 components of targeted buffer : GP

        auto to_0_x_xxyy_x = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 9);

        auto to_0_x_xxyy_y = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 10);

        auto to_0_x_xxyy_z = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 11);

        #pragma omp simd aligned(to_0_x_xxyy_x, to_0_x_xxyy_y, to_0_x_xxyy_z, to_xxyy_0, to_xxyy_xx, to_xxyy_xy, to_xxyy_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxyy_x[k] = -to_xxyy_0[k] + 2.0 * to_xxyy_xx[k] * tke_0;

            to_0_x_xxyy_y[k] = 2.0 * to_xxyy_xy[k] * tke_0;

            to_0_x_xxyy_z[k] = 2.0 * to_xxyy_xz[k] * tke_0;
        }

        // Set up 12-15 components of targeted buffer : GP

        auto to_0_x_xxyz_x = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 12);

        auto to_0_x_xxyz_y = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 13);

        auto to_0_x_xxyz_z = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 14);

        #pragma omp simd aligned(to_0_x_xxyz_x, to_0_x_xxyz_y, to_0_x_xxyz_z, to_xxyz_0, to_xxyz_xx, to_xxyz_xy, to_xxyz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxyz_x[k] = -to_xxyz_0[k] + 2.0 * to_xxyz_xx[k] * tke_0;

            to_0_x_xxyz_y[k] = 2.0 * to_xxyz_xy[k] * tke_0;

            to_0_x_xxyz_z[k] = 2.0 * to_xxyz_xz[k] * tke_0;
        }

        // Set up 15-18 components of targeted buffer : GP

        auto to_0_x_xxzz_x = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 15);

        auto to_0_x_xxzz_y = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 16);

        auto to_0_x_xxzz_z = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 17);

        #pragma omp simd aligned(to_0_x_xxzz_x, to_0_x_xxzz_y, to_0_x_xxzz_z, to_xxzz_0, to_xxzz_xx, to_xxzz_xy, to_xxzz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xxzz_x[k] = -to_xxzz_0[k] + 2.0 * to_xxzz_xx[k] * tke_0;

            to_0_x_xxzz_y[k] = 2.0 * to_xxzz_xy[k] * tke_0;

            to_0_x_xxzz_z[k] = 2.0 * to_xxzz_xz[k] * tke_0;
        }

        // Set up 18-21 components of targeted buffer : GP

        auto to_0_x_xyyy_x = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 18);

        auto to_0_x_xyyy_y = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 19);

        auto to_0_x_xyyy_z = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 20);

        #pragma omp simd aligned(to_0_x_xyyy_x, to_0_x_xyyy_y, to_0_x_xyyy_z, to_xyyy_0, to_xyyy_xx, to_xyyy_xy, to_xyyy_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xyyy_x[k] = -to_xyyy_0[k] + 2.0 * to_xyyy_xx[k] * tke_0;

            to_0_x_xyyy_y[k] = 2.0 * to_xyyy_xy[k] * tke_0;

            to_0_x_xyyy_z[k] = 2.0 * to_xyyy_xz[k] * tke_0;
        }

        // Set up 21-24 components of targeted buffer : GP

        auto to_0_x_xyyz_x = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 21);

        auto to_0_x_xyyz_y = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 22);

        auto to_0_x_xyyz_z = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 23);

        #pragma omp simd aligned(to_0_x_xyyz_x, to_0_x_xyyz_y, to_0_x_xyyz_z, to_xyyz_0, to_xyyz_xx, to_xyyz_xy, to_xyyz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xyyz_x[k] = -to_xyyz_0[k] + 2.0 * to_xyyz_xx[k] * tke_0;

            to_0_x_xyyz_y[k] = 2.0 * to_xyyz_xy[k] * tke_0;

            to_0_x_xyyz_z[k] = 2.0 * to_xyyz_xz[k] * tke_0;
        }

        // Set up 24-27 components of targeted buffer : GP

        auto to_0_x_xyzz_x = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 24);

        auto to_0_x_xyzz_y = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 25);

        auto to_0_x_xyzz_z = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 26);

        #pragma omp simd aligned(to_0_x_xyzz_x, to_0_x_xyzz_y, to_0_x_xyzz_z, to_xyzz_0, to_xyzz_xx, to_xyzz_xy, to_xyzz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xyzz_x[k] = -to_xyzz_0[k] + 2.0 * to_xyzz_xx[k] * tke_0;

            to_0_x_xyzz_y[k] = 2.0 * to_xyzz_xy[k] * tke_0;

            to_0_x_xyzz_z[k] = 2.0 * to_xyzz_xz[k] * tke_0;
        }

        // Set up 27-30 components of targeted buffer : GP

        auto to_0_x_xzzz_x = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 27);

        auto to_0_x_xzzz_y = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 28);

        auto to_0_x_xzzz_z = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 29);

        #pragma omp simd aligned(to_0_x_xzzz_x, to_0_x_xzzz_y, to_0_x_xzzz_z, to_xzzz_0, to_xzzz_xx, to_xzzz_xy, to_xzzz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xzzz_x[k] = -to_xzzz_0[k] + 2.0 * to_xzzz_xx[k] * tke_0;

            to_0_x_xzzz_y[k] = 2.0 * to_xzzz_xy[k] * tke_0;

            to_0_x_xzzz_z[k] = 2.0 * to_xzzz_xz[k] * tke_0;
        }

        // Set up 30-33 components of targeted buffer : GP

        auto to_0_x_yyyy_x = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 30);

        auto to_0_x_yyyy_y = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 31);

        auto to_0_x_yyyy_z = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 32);

        #pragma omp simd aligned(to_0_x_yyyy_x, to_0_x_yyyy_y, to_0_x_yyyy_z, to_yyyy_0, to_yyyy_xx, to_yyyy_xy, to_yyyy_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yyyy_x[k] = -to_yyyy_0[k] + 2.0 * to_yyyy_xx[k] * tke_0;

            to_0_x_yyyy_y[k] = 2.0 * to_yyyy_xy[k] * tke_0;

            to_0_x_yyyy_z[k] = 2.0 * to_yyyy_xz[k] * tke_0;
        }

        // Set up 33-36 components of targeted buffer : GP

        auto to_0_x_yyyz_x = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 33);

        auto to_0_x_yyyz_y = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 34);

        auto to_0_x_yyyz_z = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 35);

        #pragma omp simd aligned(to_0_x_yyyz_x, to_0_x_yyyz_y, to_0_x_yyyz_z, to_yyyz_0, to_yyyz_xx, to_yyyz_xy, to_yyyz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yyyz_x[k] = -to_yyyz_0[k] + 2.0 * to_yyyz_xx[k] * tke_0;

            to_0_x_yyyz_y[k] = 2.0 * to_yyyz_xy[k] * tke_0;

            to_0_x_yyyz_z[k] = 2.0 * to_yyyz_xz[k] * tke_0;
        }

        // Set up 36-39 components of targeted buffer : GP

        auto to_0_x_yyzz_x = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 36);

        auto to_0_x_yyzz_y = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 37);

        auto to_0_x_yyzz_z = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 38);

        #pragma omp simd aligned(to_0_x_yyzz_x, to_0_x_yyzz_y, to_0_x_yyzz_z, to_yyzz_0, to_yyzz_xx, to_yyzz_xy, to_yyzz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yyzz_x[k] = -to_yyzz_0[k] + 2.0 * to_yyzz_xx[k] * tke_0;

            to_0_x_yyzz_y[k] = 2.0 * to_yyzz_xy[k] * tke_0;

            to_0_x_yyzz_z[k] = 2.0 * to_yyzz_xz[k] * tke_0;
        }

        // Set up 39-42 components of targeted buffer : GP

        auto to_0_x_yzzz_x = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 39);

        auto to_0_x_yzzz_y = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 40);

        auto to_0_x_yzzz_z = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 41);

        #pragma omp simd aligned(to_0_x_yzzz_x, to_0_x_yzzz_y, to_0_x_yzzz_z, to_yzzz_0, to_yzzz_xx, to_yzzz_xy, to_yzzz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yzzz_x[k] = -to_yzzz_0[k] + 2.0 * to_yzzz_xx[k] * tke_0;

            to_0_x_yzzz_y[k] = 2.0 * to_yzzz_xy[k] * tke_0;

            to_0_x_yzzz_z[k] = 2.0 * to_yzzz_xz[k] * tke_0;
        }

        // Set up 42-45 components of targeted buffer : GP

        auto to_0_x_zzzz_x = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 42);

        auto to_0_x_zzzz_y = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 43);

        auto to_0_x_zzzz_z = pbuffer.data(idx_op_geom_001_gp + 0 * op_comps * 45 + i * 45 + 44);

        #pragma omp simd aligned(to_0_x_zzzz_x, to_0_x_zzzz_y, to_0_x_zzzz_z, to_zzzz_0, to_zzzz_xx, to_zzzz_xy, to_zzzz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_zzzz_x[k] = -to_zzzz_0[k] + 2.0 * to_zzzz_xx[k] * tke_0;

            to_0_x_zzzz_y[k] = 2.0 * to_zzzz_xy[k] * tke_0;

            to_0_x_zzzz_z[k] = 2.0 * to_zzzz_xz[k] * tke_0;
        }

        // Set up 45-48 components of targeted buffer : GP

        auto to_0_y_xxxx_x = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 0);

        auto to_0_y_xxxx_y = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 1);

        auto to_0_y_xxxx_z = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 2);

        #pragma omp simd aligned(to_0_y_xxxx_x, to_0_y_xxxx_y, to_0_y_xxxx_z, to_xxxx_0, to_xxxx_xy, to_xxxx_yy, to_xxxx_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxxx_x[k] = 2.0 * to_xxxx_xy[k] * tke_0;

            to_0_y_xxxx_y[k] = -to_xxxx_0[k] + 2.0 * to_xxxx_yy[k] * tke_0;

            to_0_y_xxxx_z[k] = 2.0 * to_xxxx_yz[k] * tke_0;
        }

        // Set up 48-51 components of targeted buffer : GP

        auto to_0_y_xxxy_x = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 3);

        auto to_0_y_xxxy_y = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 4);

        auto to_0_y_xxxy_z = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 5);

        #pragma omp simd aligned(to_0_y_xxxy_x, to_0_y_xxxy_y, to_0_y_xxxy_z, to_xxxy_0, to_xxxy_xy, to_xxxy_yy, to_xxxy_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxxy_x[k] = 2.0 * to_xxxy_xy[k] * tke_0;

            to_0_y_xxxy_y[k] = -to_xxxy_0[k] + 2.0 * to_xxxy_yy[k] * tke_0;

            to_0_y_xxxy_z[k] = 2.0 * to_xxxy_yz[k] * tke_0;
        }

        // Set up 51-54 components of targeted buffer : GP

        auto to_0_y_xxxz_x = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 6);

        auto to_0_y_xxxz_y = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 7);

        auto to_0_y_xxxz_z = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 8);

        #pragma omp simd aligned(to_0_y_xxxz_x, to_0_y_xxxz_y, to_0_y_xxxz_z, to_xxxz_0, to_xxxz_xy, to_xxxz_yy, to_xxxz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxxz_x[k] = 2.0 * to_xxxz_xy[k] * tke_0;

            to_0_y_xxxz_y[k] = -to_xxxz_0[k] + 2.0 * to_xxxz_yy[k] * tke_0;

            to_0_y_xxxz_z[k] = 2.0 * to_xxxz_yz[k] * tke_0;
        }

        // Set up 54-57 components of targeted buffer : GP

        auto to_0_y_xxyy_x = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 9);

        auto to_0_y_xxyy_y = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 10);

        auto to_0_y_xxyy_z = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 11);

        #pragma omp simd aligned(to_0_y_xxyy_x, to_0_y_xxyy_y, to_0_y_xxyy_z, to_xxyy_0, to_xxyy_xy, to_xxyy_yy, to_xxyy_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxyy_x[k] = 2.0 * to_xxyy_xy[k] * tke_0;

            to_0_y_xxyy_y[k] = -to_xxyy_0[k] + 2.0 * to_xxyy_yy[k] * tke_0;

            to_0_y_xxyy_z[k] = 2.0 * to_xxyy_yz[k] * tke_0;
        }

        // Set up 57-60 components of targeted buffer : GP

        auto to_0_y_xxyz_x = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 12);

        auto to_0_y_xxyz_y = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 13);

        auto to_0_y_xxyz_z = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 14);

        #pragma omp simd aligned(to_0_y_xxyz_x, to_0_y_xxyz_y, to_0_y_xxyz_z, to_xxyz_0, to_xxyz_xy, to_xxyz_yy, to_xxyz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxyz_x[k] = 2.0 * to_xxyz_xy[k] * tke_0;

            to_0_y_xxyz_y[k] = -to_xxyz_0[k] + 2.0 * to_xxyz_yy[k] * tke_0;

            to_0_y_xxyz_z[k] = 2.0 * to_xxyz_yz[k] * tke_0;
        }

        // Set up 60-63 components of targeted buffer : GP

        auto to_0_y_xxzz_x = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 15);

        auto to_0_y_xxzz_y = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 16);

        auto to_0_y_xxzz_z = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 17);

        #pragma omp simd aligned(to_0_y_xxzz_x, to_0_y_xxzz_y, to_0_y_xxzz_z, to_xxzz_0, to_xxzz_xy, to_xxzz_yy, to_xxzz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xxzz_x[k] = 2.0 * to_xxzz_xy[k] * tke_0;

            to_0_y_xxzz_y[k] = -to_xxzz_0[k] + 2.0 * to_xxzz_yy[k] * tke_0;

            to_0_y_xxzz_z[k] = 2.0 * to_xxzz_yz[k] * tke_0;
        }

        // Set up 63-66 components of targeted buffer : GP

        auto to_0_y_xyyy_x = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 18);

        auto to_0_y_xyyy_y = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 19);

        auto to_0_y_xyyy_z = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 20);

        #pragma omp simd aligned(to_0_y_xyyy_x, to_0_y_xyyy_y, to_0_y_xyyy_z, to_xyyy_0, to_xyyy_xy, to_xyyy_yy, to_xyyy_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xyyy_x[k] = 2.0 * to_xyyy_xy[k] * tke_0;

            to_0_y_xyyy_y[k] = -to_xyyy_0[k] + 2.0 * to_xyyy_yy[k] * tke_0;

            to_0_y_xyyy_z[k] = 2.0 * to_xyyy_yz[k] * tke_0;
        }

        // Set up 66-69 components of targeted buffer : GP

        auto to_0_y_xyyz_x = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 21);

        auto to_0_y_xyyz_y = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 22);

        auto to_0_y_xyyz_z = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 23);

        #pragma omp simd aligned(to_0_y_xyyz_x, to_0_y_xyyz_y, to_0_y_xyyz_z, to_xyyz_0, to_xyyz_xy, to_xyyz_yy, to_xyyz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xyyz_x[k] = 2.0 * to_xyyz_xy[k] * tke_0;

            to_0_y_xyyz_y[k] = -to_xyyz_0[k] + 2.0 * to_xyyz_yy[k] * tke_0;

            to_0_y_xyyz_z[k] = 2.0 * to_xyyz_yz[k] * tke_0;
        }

        // Set up 69-72 components of targeted buffer : GP

        auto to_0_y_xyzz_x = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 24);

        auto to_0_y_xyzz_y = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 25);

        auto to_0_y_xyzz_z = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 26);

        #pragma omp simd aligned(to_0_y_xyzz_x, to_0_y_xyzz_y, to_0_y_xyzz_z, to_xyzz_0, to_xyzz_xy, to_xyzz_yy, to_xyzz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xyzz_x[k] = 2.0 * to_xyzz_xy[k] * tke_0;

            to_0_y_xyzz_y[k] = -to_xyzz_0[k] + 2.0 * to_xyzz_yy[k] * tke_0;

            to_0_y_xyzz_z[k] = 2.0 * to_xyzz_yz[k] * tke_0;
        }

        // Set up 72-75 components of targeted buffer : GP

        auto to_0_y_xzzz_x = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 27);

        auto to_0_y_xzzz_y = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 28);

        auto to_0_y_xzzz_z = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 29);

        #pragma omp simd aligned(to_0_y_xzzz_x, to_0_y_xzzz_y, to_0_y_xzzz_z, to_xzzz_0, to_xzzz_xy, to_xzzz_yy, to_xzzz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xzzz_x[k] = 2.0 * to_xzzz_xy[k] * tke_0;

            to_0_y_xzzz_y[k] = -to_xzzz_0[k] + 2.0 * to_xzzz_yy[k] * tke_0;

            to_0_y_xzzz_z[k] = 2.0 * to_xzzz_yz[k] * tke_0;
        }

        // Set up 75-78 components of targeted buffer : GP

        auto to_0_y_yyyy_x = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 30);

        auto to_0_y_yyyy_y = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 31);

        auto to_0_y_yyyy_z = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 32);

        #pragma omp simd aligned(to_0_y_yyyy_x, to_0_y_yyyy_y, to_0_y_yyyy_z, to_yyyy_0, to_yyyy_xy, to_yyyy_yy, to_yyyy_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yyyy_x[k] = 2.0 * to_yyyy_xy[k] * tke_0;

            to_0_y_yyyy_y[k] = -to_yyyy_0[k] + 2.0 * to_yyyy_yy[k] * tke_0;

            to_0_y_yyyy_z[k] = 2.0 * to_yyyy_yz[k] * tke_0;
        }

        // Set up 78-81 components of targeted buffer : GP

        auto to_0_y_yyyz_x = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 33);

        auto to_0_y_yyyz_y = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 34);

        auto to_0_y_yyyz_z = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 35);

        #pragma omp simd aligned(to_0_y_yyyz_x, to_0_y_yyyz_y, to_0_y_yyyz_z, to_yyyz_0, to_yyyz_xy, to_yyyz_yy, to_yyyz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yyyz_x[k] = 2.0 * to_yyyz_xy[k] * tke_0;

            to_0_y_yyyz_y[k] = -to_yyyz_0[k] + 2.0 * to_yyyz_yy[k] * tke_0;

            to_0_y_yyyz_z[k] = 2.0 * to_yyyz_yz[k] * tke_0;
        }

        // Set up 81-84 components of targeted buffer : GP

        auto to_0_y_yyzz_x = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 36);

        auto to_0_y_yyzz_y = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 37);

        auto to_0_y_yyzz_z = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 38);

        #pragma omp simd aligned(to_0_y_yyzz_x, to_0_y_yyzz_y, to_0_y_yyzz_z, to_yyzz_0, to_yyzz_xy, to_yyzz_yy, to_yyzz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yyzz_x[k] = 2.0 * to_yyzz_xy[k] * tke_0;

            to_0_y_yyzz_y[k] = -to_yyzz_0[k] + 2.0 * to_yyzz_yy[k] * tke_0;

            to_0_y_yyzz_z[k] = 2.0 * to_yyzz_yz[k] * tke_0;
        }

        // Set up 84-87 components of targeted buffer : GP

        auto to_0_y_yzzz_x = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 39);

        auto to_0_y_yzzz_y = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 40);

        auto to_0_y_yzzz_z = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 41);

        #pragma omp simd aligned(to_0_y_yzzz_x, to_0_y_yzzz_y, to_0_y_yzzz_z, to_yzzz_0, to_yzzz_xy, to_yzzz_yy, to_yzzz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yzzz_x[k] = 2.0 * to_yzzz_xy[k] * tke_0;

            to_0_y_yzzz_y[k] = -to_yzzz_0[k] + 2.0 * to_yzzz_yy[k] * tke_0;

            to_0_y_yzzz_z[k] = 2.0 * to_yzzz_yz[k] * tke_0;
        }

        // Set up 87-90 components of targeted buffer : GP

        auto to_0_y_zzzz_x = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 42);

        auto to_0_y_zzzz_y = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 43);

        auto to_0_y_zzzz_z = pbuffer.data(idx_op_geom_001_gp + 1 * op_comps * 45 + i * 45 + 44);

        #pragma omp simd aligned(to_0_y_zzzz_x, to_0_y_zzzz_y, to_0_y_zzzz_z, to_zzzz_0, to_zzzz_xy, to_zzzz_yy, to_zzzz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_zzzz_x[k] = 2.0 * to_zzzz_xy[k] * tke_0;

            to_0_y_zzzz_y[k] = -to_zzzz_0[k] + 2.0 * to_zzzz_yy[k] * tke_0;

            to_0_y_zzzz_z[k] = 2.0 * to_zzzz_yz[k] * tke_0;
        }

        // Set up 90-93 components of targeted buffer : GP

        auto to_0_z_xxxx_x = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 0);

        auto to_0_z_xxxx_y = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 1);

        auto to_0_z_xxxx_z = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 2);

        #pragma omp simd aligned(to_0_z_xxxx_x, to_0_z_xxxx_y, to_0_z_xxxx_z, to_xxxx_0, to_xxxx_xz, to_xxxx_yz, to_xxxx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxxx_x[k] = 2.0 * to_xxxx_xz[k] * tke_0;

            to_0_z_xxxx_y[k] = 2.0 * to_xxxx_yz[k] * tke_0;

            to_0_z_xxxx_z[k] = -to_xxxx_0[k] + 2.0 * to_xxxx_zz[k] * tke_0;
        }

        // Set up 93-96 components of targeted buffer : GP

        auto to_0_z_xxxy_x = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 3);

        auto to_0_z_xxxy_y = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 4);

        auto to_0_z_xxxy_z = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 5);

        #pragma omp simd aligned(to_0_z_xxxy_x, to_0_z_xxxy_y, to_0_z_xxxy_z, to_xxxy_0, to_xxxy_xz, to_xxxy_yz, to_xxxy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxxy_x[k] = 2.0 * to_xxxy_xz[k] * tke_0;

            to_0_z_xxxy_y[k] = 2.0 * to_xxxy_yz[k] * tke_0;

            to_0_z_xxxy_z[k] = -to_xxxy_0[k] + 2.0 * to_xxxy_zz[k] * tke_0;
        }

        // Set up 96-99 components of targeted buffer : GP

        auto to_0_z_xxxz_x = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 6);

        auto to_0_z_xxxz_y = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 7);

        auto to_0_z_xxxz_z = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 8);

        #pragma omp simd aligned(to_0_z_xxxz_x, to_0_z_xxxz_y, to_0_z_xxxz_z, to_xxxz_0, to_xxxz_xz, to_xxxz_yz, to_xxxz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxxz_x[k] = 2.0 * to_xxxz_xz[k] * tke_0;

            to_0_z_xxxz_y[k] = 2.0 * to_xxxz_yz[k] * tke_0;

            to_0_z_xxxz_z[k] = -to_xxxz_0[k] + 2.0 * to_xxxz_zz[k] * tke_0;
        }

        // Set up 99-102 components of targeted buffer : GP

        auto to_0_z_xxyy_x = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 9);

        auto to_0_z_xxyy_y = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 10);

        auto to_0_z_xxyy_z = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 11);

        #pragma omp simd aligned(to_0_z_xxyy_x, to_0_z_xxyy_y, to_0_z_xxyy_z, to_xxyy_0, to_xxyy_xz, to_xxyy_yz, to_xxyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxyy_x[k] = 2.0 * to_xxyy_xz[k] * tke_0;

            to_0_z_xxyy_y[k] = 2.0 * to_xxyy_yz[k] * tke_0;

            to_0_z_xxyy_z[k] = -to_xxyy_0[k] + 2.0 * to_xxyy_zz[k] * tke_0;
        }

        // Set up 102-105 components of targeted buffer : GP

        auto to_0_z_xxyz_x = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 12);

        auto to_0_z_xxyz_y = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 13);

        auto to_0_z_xxyz_z = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 14);

        #pragma omp simd aligned(to_0_z_xxyz_x, to_0_z_xxyz_y, to_0_z_xxyz_z, to_xxyz_0, to_xxyz_xz, to_xxyz_yz, to_xxyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxyz_x[k] = 2.0 * to_xxyz_xz[k] * tke_0;

            to_0_z_xxyz_y[k] = 2.0 * to_xxyz_yz[k] * tke_0;

            to_0_z_xxyz_z[k] = -to_xxyz_0[k] + 2.0 * to_xxyz_zz[k] * tke_0;
        }

        // Set up 105-108 components of targeted buffer : GP

        auto to_0_z_xxzz_x = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 15);

        auto to_0_z_xxzz_y = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 16);

        auto to_0_z_xxzz_z = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 17);

        #pragma omp simd aligned(to_0_z_xxzz_x, to_0_z_xxzz_y, to_0_z_xxzz_z, to_xxzz_0, to_xxzz_xz, to_xxzz_yz, to_xxzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xxzz_x[k] = 2.0 * to_xxzz_xz[k] * tke_0;

            to_0_z_xxzz_y[k] = 2.0 * to_xxzz_yz[k] * tke_0;

            to_0_z_xxzz_z[k] = -to_xxzz_0[k] + 2.0 * to_xxzz_zz[k] * tke_0;
        }

        // Set up 108-111 components of targeted buffer : GP

        auto to_0_z_xyyy_x = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 18);

        auto to_0_z_xyyy_y = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 19);

        auto to_0_z_xyyy_z = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 20);

        #pragma omp simd aligned(to_0_z_xyyy_x, to_0_z_xyyy_y, to_0_z_xyyy_z, to_xyyy_0, to_xyyy_xz, to_xyyy_yz, to_xyyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xyyy_x[k] = 2.0 * to_xyyy_xz[k] * tke_0;

            to_0_z_xyyy_y[k] = 2.0 * to_xyyy_yz[k] * tke_0;

            to_0_z_xyyy_z[k] = -to_xyyy_0[k] + 2.0 * to_xyyy_zz[k] * tke_0;
        }

        // Set up 111-114 components of targeted buffer : GP

        auto to_0_z_xyyz_x = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 21);

        auto to_0_z_xyyz_y = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 22);

        auto to_0_z_xyyz_z = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 23);

        #pragma omp simd aligned(to_0_z_xyyz_x, to_0_z_xyyz_y, to_0_z_xyyz_z, to_xyyz_0, to_xyyz_xz, to_xyyz_yz, to_xyyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xyyz_x[k] = 2.0 * to_xyyz_xz[k] * tke_0;

            to_0_z_xyyz_y[k] = 2.0 * to_xyyz_yz[k] * tke_0;

            to_0_z_xyyz_z[k] = -to_xyyz_0[k] + 2.0 * to_xyyz_zz[k] * tke_0;
        }

        // Set up 114-117 components of targeted buffer : GP

        auto to_0_z_xyzz_x = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 24);

        auto to_0_z_xyzz_y = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 25);

        auto to_0_z_xyzz_z = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 26);

        #pragma omp simd aligned(to_0_z_xyzz_x, to_0_z_xyzz_y, to_0_z_xyzz_z, to_xyzz_0, to_xyzz_xz, to_xyzz_yz, to_xyzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xyzz_x[k] = 2.0 * to_xyzz_xz[k] * tke_0;

            to_0_z_xyzz_y[k] = 2.0 * to_xyzz_yz[k] * tke_0;

            to_0_z_xyzz_z[k] = -to_xyzz_0[k] + 2.0 * to_xyzz_zz[k] * tke_0;
        }

        // Set up 117-120 components of targeted buffer : GP

        auto to_0_z_xzzz_x = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 27);

        auto to_0_z_xzzz_y = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 28);

        auto to_0_z_xzzz_z = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 29);

        #pragma omp simd aligned(to_0_z_xzzz_x, to_0_z_xzzz_y, to_0_z_xzzz_z, to_xzzz_0, to_xzzz_xz, to_xzzz_yz, to_xzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xzzz_x[k] = 2.0 * to_xzzz_xz[k] * tke_0;

            to_0_z_xzzz_y[k] = 2.0 * to_xzzz_yz[k] * tke_0;

            to_0_z_xzzz_z[k] = -to_xzzz_0[k] + 2.0 * to_xzzz_zz[k] * tke_0;
        }

        // Set up 120-123 components of targeted buffer : GP

        auto to_0_z_yyyy_x = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 30);

        auto to_0_z_yyyy_y = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 31);

        auto to_0_z_yyyy_z = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 32);

        #pragma omp simd aligned(to_0_z_yyyy_x, to_0_z_yyyy_y, to_0_z_yyyy_z, to_yyyy_0, to_yyyy_xz, to_yyyy_yz, to_yyyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yyyy_x[k] = 2.0 * to_yyyy_xz[k] * tke_0;

            to_0_z_yyyy_y[k] = 2.0 * to_yyyy_yz[k] * tke_0;

            to_0_z_yyyy_z[k] = -to_yyyy_0[k] + 2.0 * to_yyyy_zz[k] * tke_0;
        }

        // Set up 123-126 components of targeted buffer : GP

        auto to_0_z_yyyz_x = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 33);

        auto to_0_z_yyyz_y = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 34);

        auto to_0_z_yyyz_z = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 35);

        #pragma omp simd aligned(to_0_z_yyyz_x, to_0_z_yyyz_y, to_0_z_yyyz_z, to_yyyz_0, to_yyyz_xz, to_yyyz_yz, to_yyyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yyyz_x[k] = 2.0 * to_yyyz_xz[k] * tke_0;

            to_0_z_yyyz_y[k] = 2.0 * to_yyyz_yz[k] * tke_0;

            to_0_z_yyyz_z[k] = -to_yyyz_0[k] + 2.0 * to_yyyz_zz[k] * tke_0;
        }

        // Set up 126-129 components of targeted buffer : GP

        auto to_0_z_yyzz_x = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 36);

        auto to_0_z_yyzz_y = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 37);

        auto to_0_z_yyzz_z = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 38);

        #pragma omp simd aligned(to_0_z_yyzz_x, to_0_z_yyzz_y, to_0_z_yyzz_z, to_yyzz_0, to_yyzz_xz, to_yyzz_yz, to_yyzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yyzz_x[k] = 2.0 * to_yyzz_xz[k] * tke_0;

            to_0_z_yyzz_y[k] = 2.0 * to_yyzz_yz[k] * tke_0;

            to_0_z_yyzz_z[k] = -to_yyzz_0[k] + 2.0 * to_yyzz_zz[k] * tke_0;
        }

        // Set up 129-132 components of targeted buffer : GP

        auto to_0_z_yzzz_x = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 39);

        auto to_0_z_yzzz_y = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 40);

        auto to_0_z_yzzz_z = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 41);

        #pragma omp simd aligned(to_0_z_yzzz_x, to_0_z_yzzz_y, to_0_z_yzzz_z, to_yzzz_0, to_yzzz_xz, to_yzzz_yz, to_yzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yzzz_x[k] = 2.0 * to_yzzz_xz[k] * tke_0;

            to_0_z_yzzz_y[k] = 2.0 * to_yzzz_yz[k] * tke_0;

            to_0_z_yzzz_z[k] = -to_yzzz_0[k] + 2.0 * to_yzzz_zz[k] * tke_0;
        }

        // Set up 132-135 components of targeted buffer : GP

        auto to_0_z_zzzz_x = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 42);

        auto to_0_z_zzzz_y = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 43);

        auto to_0_z_zzzz_z = pbuffer.data(idx_op_geom_001_gp + 2 * op_comps * 45 + i * 45 + 44);

        #pragma omp simd aligned(to_0_z_zzzz_x, to_0_z_zzzz_y, to_0_z_zzzz_z, to_zzzz_0, to_zzzz_xz, to_zzzz_yz, to_zzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_zzzz_x[k] = 2.0 * to_zzzz_xz[k] * tke_0;

            to_0_z_zzzz_y[k] = 2.0 * to_zzzz_yz[k] * tke_0;

            to_0_z_zzzz_z[k] = -to_zzzz_0[k] + 2.0 * to_zzzz_zz[k] * tke_0;
        }

    }

}

} // t2cgeom namespace

