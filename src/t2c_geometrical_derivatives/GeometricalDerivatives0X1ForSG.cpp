#include "GeometricalDerivatives0X1ForSG.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_geom_deriv_0x1_sg(CSimdArray<double>& pbuffer,
                       const int idx_op_geom_001_sg,
                       const int idx_op_sf,
                       const int idx_op_sh,
                       const CSimdArray<double>& factors,
                       const int op_comps) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : SF

        auto to_0_xxx = pbuffer.data(idx_op_sf + i * 10 + 0);

        auto to_0_xxy = pbuffer.data(idx_op_sf + i * 10 + 1);

        auto to_0_xxz = pbuffer.data(idx_op_sf + i * 10 + 2);

        auto to_0_xyy = pbuffer.data(idx_op_sf + i * 10 + 3);

        auto to_0_xyz = pbuffer.data(idx_op_sf + i * 10 + 4);

        auto to_0_xzz = pbuffer.data(idx_op_sf + i * 10 + 5);

        auto to_0_yyy = pbuffer.data(idx_op_sf + i * 10 + 6);

        auto to_0_yyz = pbuffer.data(idx_op_sf + i * 10 + 7);

        auto to_0_yzz = pbuffer.data(idx_op_sf + i * 10 + 8);

        auto to_0_zzz = pbuffer.data(idx_op_sf + i * 10 + 9);

        // Set up components of auxiliary buffer : SH

        auto to_0_xxxxx = pbuffer.data(idx_op_sh + i * 21 + 0);

        auto to_0_xxxxy = pbuffer.data(idx_op_sh + i * 21 + 1);

        auto to_0_xxxxz = pbuffer.data(idx_op_sh + i * 21 + 2);

        auto to_0_xxxyy = pbuffer.data(idx_op_sh + i * 21 + 3);

        auto to_0_xxxyz = pbuffer.data(idx_op_sh + i * 21 + 4);

        auto to_0_xxxzz = pbuffer.data(idx_op_sh + i * 21 + 5);

        auto to_0_xxyyy = pbuffer.data(idx_op_sh + i * 21 + 6);

        auto to_0_xxyyz = pbuffer.data(idx_op_sh + i * 21 + 7);

        auto to_0_xxyzz = pbuffer.data(idx_op_sh + i * 21 + 8);

        auto to_0_xxzzz = pbuffer.data(idx_op_sh + i * 21 + 9);

        auto to_0_xyyyy = pbuffer.data(idx_op_sh + i * 21 + 10);

        auto to_0_xyyyz = pbuffer.data(idx_op_sh + i * 21 + 11);

        auto to_0_xyyzz = pbuffer.data(idx_op_sh + i * 21 + 12);

        auto to_0_xyzzz = pbuffer.data(idx_op_sh + i * 21 + 13);

        auto to_0_xzzzz = pbuffer.data(idx_op_sh + i * 21 + 14);

        auto to_0_yyyyy = pbuffer.data(idx_op_sh + i * 21 + 15);

        auto to_0_yyyyz = pbuffer.data(idx_op_sh + i * 21 + 16);

        auto to_0_yyyzz = pbuffer.data(idx_op_sh + i * 21 + 17);

        auto to_0_yyzzz = pbuffer.data(idx_op_sh + i * 21 + 18);

        auto to_0_yzzzz = pbuffer.data(idx_op_sh + i * 21 + 19);

        auto to_0_zzzzz = pbuffer.data(idx_op_sh + i * 21 + 20);

        // Set up 0-15 components of targeted buffer : SG

        auto to_0_x_0_xxxx = pbuffer.data(idx_op_geom_001_sg + 0 * op_comps * 15 + i * 15 + 0);

        auto to_0_x_0_xxxy = pbuffer.data(idx_op_geom_001_sg + 0 * op_comps * 15 + i * 15 + 1);

        auto to_0_x_0_xxxz = pbuffer.data(idx_op_geom_001_sg + 0 * op_comps * 15 + i * 15 + 2);

        auto to_0_x_0_xxyy = pbuffer.data(idx_op_geom_001_sg + 0 * op_comps * 15 + i * 15 + 3);

        auto to_0_x_0_xxyz = pbuffer.data(idx_op_geom_001_sg + 0 * op_comps * 15 + i * 15 + 4);

        auto to_0_x_0_xxzz = pbuffer.data(idx_op_geom_001_sg + 0 * op_comps * 15 + i * 15 + 5);

        auto to_0_x_0_xyyy = pbuffer.data(idx_op_geom_001_sg + 0 * op_comps * 15 + i * 15 + 6);

        auto to_0_x_0_xyyz = pbuffer.data(idx_op_geom_001_sg + 0 * op_comps * 15 + i * 15 + 7);

        auto to_0_x_0_xyzz = pbuffer.data(idx_op_geom_001_sg + 0 * op_comps * 15 + i * 15 + 8);

        auto to_0_x_0_xzzz = pbuffer.data(idx_op_geom_001_sg + 0 * op_comps * 15 + i * 15 + 9);

        auto to_0_x_0_yyyy = pbuffer.data(idx_op_geom_001_sg + 0 * op_comps * 15 + i * 15 + 10);

        auto to_0_x_0_yyyz = pbuffer.data(idx_op_geom_001_sg + 0 * op_comps * 15 + i * 15 + 11);

        auto to_0_x_0_yyzz = pbuffer.data(idx_op_geom_001_sg + 0 * op_comps * 15 + i * 15 + 12);

        auto to_0_x_0_yzzz = pbuffer.data(idx_op_geom_001_sg + 0 * op_comps * 15 + i * 15 + 13);

        auto to_0_x_0_zzzz = pbuffer.data(idx_op_geom_001_sg + 0 * op_comps * 15 + i * 15 + 14);

        #pragma omp simd aligned(to_0_x_0_xxxx, to_0_x_0_xxxy, to_0_x_0_xxxz, to_0_x_0_xxyy, to_0_x_0_xxyz, to_0_x_0_xxzz, to_0_x_0_xyyy, to_0_x_0_xyyz, to_0_x_0_xyzz, to_0_x_0_xzzz, to_0_x_0_yyyy, to_0_x_0_yyyz, to_0_x_0_yyzz, to_0_x_0_yzzz, to_0_x_0_zzzz, to_0_xxx, to_0_xxxxx, to_0_xxxxy, to_0_xxxxz, to_0_xxxyy, to_0_xxxyz, to_0_xxxzz, to_0_xxy, to_0_xxyyy, to_0_xxyyz, to_0_xxyzz, to_0_xxz, to_0_xxzzz, to_0_xyy, to_0_xyyyy, to_0_xyyyz, to_0_xyyzz, to_0_xyz, to_0_xyzzz, to_0_xzz, to_0_xzzzz, to_0_yyy, to_0_yyz, to_0_yzz, to_0_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_0_xxxx[k] = -4.0 * to_0_xxx[k] + 2.0 * to_0_xxxxx[k] * tke_0;

            to_0_x_0_xxxy[k] = -3.0 * to_0_xxy[k] + 2.0 * to_0_xxxxy[k] * tke_0;

            to_0_x_0_xxxz[k] = -3.0 * to_0_xxz[k] + 2.0 * to_0_xxxxz[k] * tke_0;

            to_0_x_0_xxyy[k] = -2.0 * to_0_xyy[k] + 2.0 * to_0_xxxyy[k] * tke_0;

            to_0_x_0_xxyz[k] = -2.0 * to_0_xyz[k] + 2.0 * to_0_xxxyz[k] * tke_0;

            to_0_x_0_xxzz[k] = -2.0 * to_0_xzz[k] + 2.0 * to_0_xxxzz[k] * tke_0;

            to_0_x_0_xyyy[k] = -to_0_yyy[k] + 2.0 * to_0_xxyyy[k] * tke_0;

            to_0_x_0_xyyz[k] = -to_0_yyz[k] + 2.0 * to_0_xxyyz[k] * tke_0;

            to_0_x_0_xyzz[k] = -to_0_yzz[k] + 2.0 * to_0_xxyzz[k] * tke_0;

            to_0_x_0_xzzz[k] = -to_0_zzz[k] + 2.0 * to_0_xxzzz[k] * tke_0;

            to_0_x_0_yyyy[k] = 2.0 * to_0_xyyyy[k] * tke_0;

            to_0_x_0_yyyz[k] = 2.0 * to_0_xyyyz[k] * tke_0;

            to_0_x_0_yyzz[k] = 2.0 * to_0_xyyzz[k] * tke_0;

            to_0_x_0_yzzz[k] = 2.0 * to_0_xyzzz[k] * tke_0;

            to_0_x_0_zzzz[k] = 2.0 * to_0_xzzzz[k] * tke_0;
        }

        // Set up 15-30 components of targeted buffer : SG

        auto to_0_y_0_xxxx = pbuffer.data(idx_op_geom_001_sg + 1 * op_comps * 15 + i * 15 + 0);

        auto to_0_y_0_xxxy = pbuffer.data(idx_op_geom_001_sg + 1 * op_comps * 15 + i * 15 + 1);

        auto to_0_y_0_xxxz = pbuffer.data(idx_op_geom_001_sg + 1 * op_comps * 15 + i * 15 + 2);

        auto to_0_y_0_xxyy = pbuffer.data(idx_op_geom_001_sg + 1 * op_comps * 15 + i * 15 + 3);

        auto to_0_y_0_xxyz = pbuffer.data(idx_op_geom_001_sg + 1 * op_comps * 15 + i * 15 + 4);

        auto to_0_y_0_xxzz = pbuffer.data(idx_op_geom_001_sg + 1 * op_comps * 15 + i * 15 + 5);

        auto to_0_y_0_xyyy = pbuffer.data(idx_op_geom_001_sg + 1 * op_comps * 15 + i * 15 + 6);

        auto to_0_y_0_xyyz = pbuffer.data(idx_op_geom_001_sg + 1 * op_comps * 15 + i * 15 + 7);

        auto to_0_y_0_xyzz = pbuffer.data(idx_op_geom_001_sg + 1 * op_comps * 15 + i * 15 + 8);

        auto to_0_y_0_xzzz = pbuffer.data(idx_op_geom_001_sg + 1 * op_comps * 15 + i * 15 + 9);

        auto to_0_y_0_yyyy = pbuffer.data(idx_op_geom_001_sg + 1 * op_comps * 15 + i * 15 + 10);

        auto to_0_y_0_yyyz = pbuffer.data(idx_op_geom_001_sg + 1 * op_comps * 15 + i * 15 + 11);

        auto to_0_y_0_yyzz = pbuffer.data(idx_op_geom_001_sg + 1 * op_comps * 15 + i * 15 + 12);

        auto to_0_y_0_yzzz = pbuffer.data(idx_op_geom_001_sg + 1 * op_comps * 15 + i * 15 + 13);

        auto to_0_y_0_zzzz = pbuffer.data(idx_op_geom_001_sg + 1 * op_comps * 15 + i * 15 + 14);

        #pragma omp simd aligned(to_0_xxx, to_0_xxxxy, to_0_xxxyy, to_0_xxxyz, to_0_xxy, to_0_xxyyy, to_0_xxyyz, to_0_xxyzz, to_0_xxz, to_0_xyy, to_0_xyyyy, to_0_xyyyz, to_0_xyyzz, to_0_xyz, to_0_xyzzz, to_0_xzz, to_0_y_0_xxxx, to_0_y_0_xxxy, to_0_y_0_xxxz, to_0_y_0_xxyy, to_0_y_0_xxyz, to_0_y_0_xxzz, to_0_y_0_xyyy, to_0_y_0_xyyz, to_0_y_0_xyzz, to_0_y_0_xzzz, to_0_y_0_yyyy, to_0_y_0_yyyz, to_0_y_0_yyzz, to_0_y_0_yzzz, to_0_y_0_zzzz, to_0_yyy, to_0_yyyyy, to_0_yyyyz, to_0_yyyzz, to_0_yyz, to_0_yyzzz, to_0_yzz, to_0_yzzzz, to_0_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_0_xxxx[k] = 2.0 * to_0_xxxxy[k] * tke_0;

            to_0_y_0_xxxy[k] = -to_0_xxx[k] + 2.0 * to_0_xxxyy[k] * tke_0;

            to_0_y_0_xxxz[k] = 2.0 * to_0_xxxyz[k] * tke_0;

            to_0_y_0_xxyy[k] = -2.0 * to_0_xxy[k] + 2.0 * to_0_xxyyy[k] * tke_0;

            to_0_y_0_xxyz[k] = -to_0_xxz[k] + 2.0 * to_0_xxyyz[k] * tke_0;

            to_0_y_0_xxzz[k] = 2.0 * to_0_xxyzz[k] * tke_0;

            to_0_y_0_xyyy[k] = -3.0 * to_0_xyy[k] + 2.0 * to_0_xyyyy[k] * tke_0;

            to_0_y_0_xyyz[k] = -2.0 * to_0_xyz[k] + 2.0 * to_0_xyyyz[k] * tke_0;

            to_0_y_0_xyzz[k] = -to_0_xzz[k] + 2.0 * to_0_xyyzz[k] * tke_0;

            to_0_y_0_xzzz[k] = 2.0 * to_0_xyzzz[k] * tke_0;

            to_0_y_0_yyyy[k] = -4.0 * to_0_yyy[k] + 2.0 * to_0_yyyyy[k] * tke_0;

            to_0_y_0_yyyz[k] = -3.0 * to_0_yyz[k] + 2.0 * to_0_yyyyz[k] * tke_0;

            to_0_y_0_yyzz[k] = -2.0 * to_0_yzz[k] + 2.0 * to_0_yyyzz[k] * tke_0;

            to_0_y_0_yzzz[k] = -to_0_zzz[k] + 2.0 * to_0_yyzzz[k] * tke_0;

            to_0_y_0_zzzz[k] = 2.0 * to_0_yzzzz[k] * tke_0;
        }

        // Set up 30-45 components of targeted buffer : SG

        auto to_0_z_0_xxxx = pbuffer.data(idx_op_geom_001_sg + 2 * op_comps * 15 + i * 15 + 0);

        auto to_0_z_0_xxxy = pbuffer.data(idx_op_geom_001_sg + 2 * op_comps * 15 + i * 15 + 1);

        auto to_0_z_0_xxxz = pbuffer.data(idx_op_geom_001_sg + 2 * op_comps * 15 + i * 15 + 2);

        auto to_0_z_0_xxyy = pbuffer.data(idx_op_geom_001_sg + 2 * op_comps * 15 + i * 15 + 3);

        auto to_0_z_0_xxyz = pbuffer.data(idx_op_geom_001_sg + 2 * op_comps * 15 + i * 15 + 4);

        auto to_0_z_0_xxzz = pbuffer.data(idx_op_geom_001_sg + 2 * op_comps * 15 + i * 15 + 5);

        auto to_0_z_0_xyyy = pbuffer.data(idx_op_geom_001_sg + 2 * op_comps * 15 + i * 15 + 6);

        auto to_0_z_0_xyyz = pbuffer.data(idx_op_geom_001_sg + 2 * op_comps * 15 + i * 15 + 7);

        auto to_0_z_0_xyzz = pbuffer.data(idx_op_geom_001_sg + 2 * op_comps * 15 + i * 15 + 8);

        auto to_0_z_0_xzzz = pbuffer.data(idx_op_geom_001_sg + 2 * op_comps * 15 + i * 15 + 9);

        auto to_0_z_0_yyyy = pbuffer.data(idx_op_geom_001_sg + 2 * op_comps * 15 + i * 15 + 10);

        auto to_0_z_0_yyyz = pbuffer.data(idx_op_geom_001_sg + 2 * op_comps * 15 + i * 15 + 11);

        auto to_0_z_0_yyzz = pbuffer.data(idx_op_geom_001_sg + 2 * op_comps * 15 + i * 15 + 12);

        auto to_0_z_0_yzzz = pbuffer.data(idx_op_geom_001_sg + 2 * op_comps * 15 + i * 15 + 13);

        auto to_0_z_0_zzzz = pbuffer.data(idx_op_geom_001_sg + 2 * op_comps * 15 + i * 15 + 14);

        #pragma omp simd aligned(to_0_xxx, to_0_xxxxz, to_0_xxxyz, to_0_xxxzz, to_0_xxy, to_0_xxyyz, to_0_xxyzz, to_0_xxz, to_0_xxzzz, to_0_xyy, to_0_xyyyz, to_0_xyyzz, to_0_xyz, to_0_xyzzz, to_0_xzz, to_0_xzzzz, to_0_yyy, to_0_yyyyz, to_0_yyyzz, to_0_yyz, to_0_yyzzz, to_0_yzz, to_0_yzzzz, to_0_z_0_xxxx, to_0_z_0_xxxy, to_0_z_0_xxxz, to_0_z_0_xxyy, to_0_z_0_xxyz, to_0_z_0_xxzz, to_0_z_0_xyyy, to_0_z_0_xyyz, to_0_z_0_xyzz, to_0_z_0_xzzz, to_0_z_0_yyyy, to_0_z_0_yyyz, to_0_z_0_yyzz, to_0_z_0_yzzz, to_0_z_0_zzzz, to_0_zzz, to_0_zzzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_0_xxxx[k] = 2.0 * to_0_xxxxz[k] * tke_0;

            to_0_z_0_xxxy[k] = 2.0 * to_0_xxxyz[k] * tke_0;

            to_0_z_0_xxxz[k] = -to_0_xxx[k] + 2.0 * to_0_xxxzz[k] * tke_0;

            to_0_z_0_xxyy[k] = 2.0 * to_0_xxyyz[k] * tke_0;

            to_0_z_0_xxyz[k] = -to_0_xxy[k] + 2.0 * to_0_xxyzz[k] * tke_0;

            to_0_z_0_xxzz[k] = -2.0 * to_0_xxz[k] + 2.0 * to_0_xxzzz[k] * tke_0;

            to_0_z_0_xyyy[k] = 2.0 * to_0_xyyyz[k] * tke_0;

            to_0_z_0_xyyz[k] = -to_0_xyy[k] + 2.0 * to_0_xyyzz[k] * tke_0;

            to_0_z_0_xyzz[k] = -2.0 * to_0_xyz[k] + 2.0 * to_0_xyzzz[k] * tke_0;

            to_0_z_0_xzzz[k] = -3.0 * to_0_xzz[k] + 2.0 * to_0_xzzzz[k] * tke_0;

            to_0_z_0_yyyy[k] = 2.0 * to_0_yyyyz[k] * tke_0;

            to_0_z_0_yyyz[k] = -to_0_yyy[k] + 2.0 * to_0_yyyzz[k] * tke_0;

            to_0_z_0_yyzz[k] = -2.0 * to_0_yyz[k] + 2.0 * to_0_yyzzz[k] * tke_0;

            to_0_z_0_yzzz[k] = -3.0 * to_0_yzz[k] + 2.0 * to_0_yzzzz[k] * tke_0;

            to_0_z_0_zzzz[k] = -4.0 * to_0_zzz[k] + 2.0 * to_0_zzzzz[k] * tke_0;
        }

    }

}

} // t2cgeom namespace

