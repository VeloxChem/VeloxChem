#include "GeometricalDerivatives2X0ForFY.hpp"

namespace t2cgeom {  // t2cgeom namespace

auto
comp_prim_op_geom_20_fx(CSimdArray<double>& pbuffer,
                        const size_t        idx_op_geom_200_fs,
                        const size_t        idx_op_ps,
                        const size_t        idx_op_fs,
                        const size_t        idx_op_hs,
                        const size_t        op_comps,
                        const size_t        ket_comps,
                        const double        a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    for (size_t i = 0; i < op_comps; i++)
    {
        for (size_t j = 0; j < ket_comps; j++)
        {
            // Set up components of auxiliary buffer : PS

            auto to_x_0 = pbuffer.data(idx_op_ps + i * 3 * ket_comps + 0 * ket_comps + j);

            auto to_y_0 = pbuffer.data(idx_op_ps + i * 3 * ket_comps + 1 * ket_comps + j);

            auto to_z_0 = pbuffer.data(idx_op_ps + i * 3 * ket_comps + 2 * ket_comps + j);

            // Set up components of auxiliary buffer : FS

            auto to_xxx_0 = pbuffer.data(idx_op_fs + i * 10 * ket_comps + 0 * ket_comps + j);

            auto to_xxy_0 = pbuffer.data(idx_op_fs + i * 10 * ket_comps + 1 * ket_comps + j);

            auto to_xxz_0 = pbuffer.data(idx_op_fs + i * 10 * ket_comps + 2 * ket_comps + j);

            auto to_xyy_0 = pbuffer.data(idx_op_fs + i * 10 * ket_comps + 3 * ket_comps + j);

            auto to_xyz_0 = pbuffer.data(idx_op_fs + i * 10 * ket_comps + 4 * ket_comps + j);

            auto to_xzz_0 = pbuffer.data(idx_op_fs + i * 10 * ket_comps + 5 * ket_comps + j);

            auto to_yyy_0 = pbuffer.data(idx_op_fs + i * 10 * ket_comps + 6 * ket_comps + j);

            auto to_yyz_0 = pbuffer.data(idx_op_fs + i * 10 * ket_comps + 7 * ket_comps + j);

            auto to_yzz_0 = pbuffer.data(idx_op_fs + i * 10 * ket_comps + 8 * ket_comps + j);

            auto to_zzz_0 = pbuffer.data(idx_op_fs + i * 10 * ket_comps + 9 * ket_comps + j);

            // Set up components of auxiliary buffer : HS

            auto to_xxxxx_0 = pbuffer.data(idx_op_hs + i * 21 * ket_comps + 0 * ket_comps + j);

            auto to_xxxxy_0 = pbuffer.data(idx_op_hs + i * 21 * ket_comps + 1 * ket_comps + j);

            auto to_xxxxz_0 = pbuffer.data(idx_op_hs + i * 21 * ket_comps + 2 * ket_comps + j);

            auto to_xxxyy_0 = pbuffer.data(idx_op_hs + i * 21 * ket_comps + 3 * ket_comps + j);

            auto to_xxxyz_0 = pbuffer.data(idx_op_hs + i * 21 * ket_comps + 4 * ket_comps + j);

            auto to_xxxzz_0 = pbuffer.data(idx_op_hs + i * 21 * ket_comps + 5 * ket_comps + j);

            auto to_xxyyy_0 = pbuffer.data(idx_op_hs + i * 21 * ket_comps + 6 * ket_comps + j);

            auto to_xxyyz_0 = pbuffer.data(idx_op_hs + i * 21 * ket_comps + 7 * ket_comps + j);

            auto to_xxyzz_0 = pbuffer.data(idx_op_hs + i * 21 * ket_comps + 8 * ket_comps + j);

            auto to_xxzzz_0 = pbuffer.data(idx_op_hs + i * 21 * ket_comps + 9 * ket_comps + j);

            auto to_xyyyy_0 = pbuffer.data(idx_op_hs + i * 21 * ket_comps + 10 * ket_comps + j);

            auto to_xyyyz_0 = pbuffer.data(idx_op_hs + i * 21 * ket_comps + 11 * ket_comps + j);

            auto to_xyyzz_0 = pbuffer.data(idx_op_hs + i * 21 * ket_comps + 12 * ket_comps + j);

            auto to_xyzzz_0 = pbuffer.data(idx_op_hs + i * 21 * ket_comps + 13 * ket_comps + j);

            auto to_xzzzz_0 = pbuffer.data(idx_op_hs + i * 21 * ket_comps + 14 * ket_comps + j);

            auto to_yyyyy_0 = pbuffer.data(idx_op_hs + i * 21 * ket_comps + 15 * ket_comps + j);

            auto to_yyyyz_0 = pbuffer.data(idx_op_hs + i * 21 * ket_comps + 16 * ket_comps + j);

            auto to_yyyzz_0 = pbuffer.data(idx_op_hs + i * 21 * ket_comps + 17 * ket_comps + j);

            auto to_yyzzz_0 = pbuffer.data(idx_op_hs + i * 21 * ket_comps + 18 * ket_comps + j);

            auto to_yzzzz_0 = pbuffer.data(idx_op_hs + i * 21 * ket_comps + 19 * ket_comps + j);

            auto to_zzzzz_0 = pbuffer.data(idx_op_hs + i * 21 * ket_comps + 20 * ket_comps + j);

            // Set up 0-10 components of targeted buffer : FS

            auto to_xx_0_xxx_0 = pbuffer.data(idx_op_geom_200_fs + 0 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 0 * ket_comps + j);

            auto to_xx_0_xxy_0 = pbuffer.data(idx_op_geom_200_fs + 0 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 1 * ket_comps + j);

            auto to_xx_0_xxz_0 = pbuffer.data(idx_op_geom_200_fs + 0 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 2 * ket_comps + j);

            auto to_xx_0_xyy_0 = pbuffer.data(idx_op_geom_200_fs + 0 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 3 * ket_comps + j);

            auto to_xx_0_xyz_0 = pbuffer.data(idx_op_geom_200_fs + 0 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 4 * ket_comps + j);

            auto to_xx_0_xzz_0 = pbuffer.data(idx_op_geom_200_fs + 0 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 5 * ket_comps + j);

            auto to_xx_0_yyy_0 = pbuffer.data(idx_op_geom_200_fs + 0 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 6 * ket_comps + j);

            auto to_xx_0_yyz_0 = pbuffer.data(idx_op_geom_200_fs + 0 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 7 * ket_comps + j);

            auto to_xx_0_yzz_0 = pbuffer.data(idx_op_geom_200_fs + 0 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 8 * ket_comps + j);

            auto to_xx_0_zzz_0 = pbuffer.data(idx_op_geom_200_fs + 0 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 9 * ket_comps + j);

#pragma omp simd aligned(to_x_0,            \
                             to_xx_0_xxx_0, \
                             to_xx_0_xxy_0, \
                             to_xx_0_xxz_0, \
                             to_xx_0_xyy_0, \
                             to_xx_0_xyz_0, \
                             to_xx_0_xzz_0, \
                             to_xx_0_yyy_0, \
                             to_xx_0_yyz_0, \
                             to_xx_0_yzz_0, \
                             to_xx_0_zzz_0, \
                             to_xxx_0,      \
                             to_xxxxx_0,    \
                             to_xxxxy_0,    \
                             to_xxxxz_0,    \
                             to_xxxyy_0,    \
                             to_xxxyz_0,    \
                             to_xxxzz_0,    \
                             to_xxy_0,      \
                             to_xxyyy_0,    \
                             to_xxyyz_0,    \
                             to_xxyzz_0,    \
                             to_xxz_0,      \
                             to_xxzzz_0,    \
                             to_xyy_0,      \
                             to_xyz_0,      \
                             to_xzz_0,      \
                             to_y_0,        \
                             to_yyy_0,      \
                             to_yyz_0,      \
                             to_yzz_0,      \
                             to_z_0,        \
                             to_zzz_0 : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                const double tbe_0 = a_exp;

                to_xx_0_xxx_0[k] = 6.0 * to_x_0[k] - 14.0 * to_xxx_0[k] * tbe_0 + 4.0 * to_xxxxx_0[k] * tbe_0 * tbe_0;

                to_xx_0_xxy_0[k] = 2.0 * to_y_0[k] - 10.0 * to_xxy_0[k] * tbe_0 + 4.0 * to_xxxxy_0[k] * tbe_0 * tbe_0;

                to_xx_0_xxz_0[k] = 2.0 * to_z_0[k] - 10.0 * to_xxz_0[k] * tbe_0 + 4.0 * to_xxxxz_0[k] * tbe_0 * tbe_0;

                to_xx_0_xyy_0[k] = -6.0 * to_xyy_0[k] * tbe_0 + 4.0 * to_xxxyy_0[k] * tbe_0 * tbe_0;

                to_xx_0_xyz_0[k] = -6.0 * to_xyz_0[k] * tbe_0 + 4.0 * to_xxxyz_0[k] * tbe_0 * tbe_0;

                to_xx_0_xzz_0[k] = -6.0 * to_xzz_0[k] * tbe_0 + 4.0 * to_xxxzz_0[k] * tbe_0 * tbe_0;

                to_xx_0_yyy_0[k] = -2.0 * to_yyy_0[k] * tbe_0 + 4.0 * to_xxyyy_0[k] * tbe_0 * tbe_0;

                to_xx_0_yyz_0[k] = -2.0 * to_yyz_0[k] * tbe_0 + 4.0 * to_xxyyz_0[k] * tbe_0 * tbe_0;

                to_xx_0_yzz_0[k] = -2.0 * to_yzz_0[k] * tbe_0 + 4.0 * to_xxyzz_0[k] * tbe_0 * tbe_0;

                to_xx_0_zzz_0[k] = -2.0 * to_zzz_0[k] * tbe_0 + 4.0 * to_xxzzz_0[k] * tbe_0 * tbe_0;
            }

            // Set up 10-20 components of targeted buffer : FS

            auto to_xy_0_xxx_0 = pbuffer.data(idx_op_geom_200_fs + 1 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 0 * ket_comps + j);

            auto to_xy_0_xxy_0 = pbuffer.data(idx_op_geom_200_fs + 1 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 1 * ket_comps + j);

            auto to_xy_0_xxz_0 = pbuffer.data(idx_op_geom_200_fs + 1 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 2 * ket_comps + j);

            auto to_xy_0_xyy_0 = pbuffer.data(idx_op_geom_200_fs + 1 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 3 * ket_comps + j);

            auto to_xy_0_xyz_0 = pbuffer.data(idx_op_geom_200_fs + 1 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 4 * ket_comps + j);

            auto to_xy_0_xzz_0 = pbuffer.data(idx_op_geom_200_fs + 1 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 5 * ket_comps + j);

            auto to_xy_0_yyy_0 = pbuffer.data(idx_op_geom_200_fs + 1 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 6 * ket_comps + j);

            auto to_xy_0_yyz_0 = pbuffer.data(idx_op_geom_200_fs + 1 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 7 * ket_comps + j);

            auto to_xy_0_yzz_0 = pbuffer.data(idx_op_geom_200_fs + 1 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 8 * ket_comps + j);

            auto to_xy_0_zzz_0 = pbuffer.data(idx_op_geom_200_fs + 1 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 9 * ket_comps + j);

#pragma omp simd aligned(to_x_0,            \
                             to_xxx_0,      \
                             to_xxxxy_0,    \
                             to_xxxyy_0,    \
                             to_xxxyz_0,    \
                             to_xxy_0,      \
                             to_xxyyy_0,    \
                             to_xxyyz_0,    \
                             to_xxyzz_0,    \
                             to_xxz_0,      \
                             to_xy_0_xxx_0, \
                             to_xy_0_xxy_0, \
                             to_xy_0_xxz_0, \
                             to_xy_0_xyy_0, \
                             to_xy_0_xyz_0, \
                             to_xy_0_xzz_0, \
                             to_xy_0_yyy_0, \
                             to_xy_0_yyz_0, \
                             to_xy_0_yzz_0, \
                             to_xy_0_zzz_0, \
                             to_xyy_0,      \
                             to_xyyyy_0,    \
                             to_xyyyz_0,    \
                             to_xyyzz_0,    \
                             to_xyz_0,      \
                             to_xyzzz_0,    \
                             to_xzz_0,      \
                             to_y_0,        \
                             to_yyy_0,      \
                             to_yyz_0,      \
                             to_yzz_0,      \
                             to_z_0 : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                const double tbe_0 = a_exp;

                to_xy_0_xxx_0[k] = -6.0 * to_xxy_0[k] * tbe_0 + 4.0 * to_xxxxy_0[k] * tbe_0 * tbe_0;

                to_xy_0_xxy_0[k] = 2.0 * to_x_0[k] - 4.0 * to_xyy_0[k] * tbe_0 - 2.0 * to_xxx_0[k] * tbe_0 + 4.0 * to_xxxyy_0[k] * tbe_0 * tbe_0;

                to_xy_0_xxz_0[k] = -4.0 * to_xyz_0[k] * tbe_0 + 4.0 * to_xxxyz_0[k] * tbe_0 * tbe_0;

                to_xy_0_xyy_0[k] = 2.0 * to_y_0[k] - 2.0 * to_yyy_0[k] * tbe_0 - 4.0 * to_xxy_0[k] * tbe_0 + 4.0 * to_xxyyy_0[k] * tbe_0 * tbe_0;

                to_xy_0_xyz_0[k] = to_z_0[k] - 2.0 * to_yyz_0[k] * tbe_0 - 2.0 * to_xxz_0[k] * tbe_0 + 4.0 * to_xxyyz_0[k] * tbe_0 * tbe_0;

                to_xy_0_xzz_0[k] = -2.0 * to_yzz_0[k] * tbe_0 + 4.0 * to_xxyzz_0[k] * tbe_0 * tbe_0;

                to_xy_0_yyy_0[k] = -6.0 * to_xyy_0[k] * tbe_0 + 4.0 * to_xyyyy_0[k] * tbe_0 * tbe_0;

                to_xy_0_yyz_0[k] = -4.0 * to_xyz_0[k] * tbe_0 + 4.0 * to_xyyyz_0[k] * tbe_0 * tbe_0;

                to_xy_0_yzz_0[k] = -2.0 * to_xzz_0[k] * tbe_0 + 4.0 * to_xyyzz_0[k] * tbe_0 * tbe_0;

                to_xy_0_zzz_0[k] = 4.0 * to_xyzzz_0[k] * tbe_0 * tbe_0;
            }

            // Set up 20-30 components of targeted buffer : FS

            auto to_xz_0_xxx_0 = pbuffer.data(idx_op_geom_200_fs + 2 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 0 * ket_comps + j);

            auto to_xz_0_xxy_0 = pbuffer.data(idx_op_geom_200_fs + 2 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 1 * ket_comps + j);

            auto to_xz_0_xxz_0 = pbuffer.data(idx_op_geom_200_fs + 2 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 2 * ket_comps + j);

            auto to_xz_0_xyy_0 = pbuffer.data(idx_op_geom_200_fs + 2 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 3 * ket_comps + j);

            auto to_xz_0_xyz_0 = pbuffer.data(idx_op_geom_200_fs + 2 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 4 * ket_comps + j);

            auto to_xz_0_xzz_0 = pbuffer.data(idx_op_geom_200_fs + 2 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 5 * ket_comps + j);

            auto to_xz_0_yyy_0 = pbuffer.data(idx_op_geom_200_fs + 2 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 6 * ket_comps + j);

            auto to_xz_0_yyz_0 = pbuffer.data(idx_op_geom_200_fs + 2 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 7 * ket_comps + j);

            auto to_xz_0_yzz_0 = pbuffer.data(idx_op_geom_200_fs + 2 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 8 * ket_comps + j);

            auto to_xz_0_zzz_0 = pbuffer.data(idx_op_geom_200_fs + 2 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 9 * ket_comps + j);

#pragma omp simd aligned(to_x_0,            \
                             to_xxx_0,      \
                             to_xxxxz_0,    \
                             to_xxxyz_0,    \
                             to_xxxzz_0,    \
                             to_xxy_0,      \
                             to_xxyyz_0,    \
                             to_xxyzz_0,    \
                             to_xxz_0,      \
                             to_xxzzz_0,    \
                             to_xyy_0,      \
                             to_xyyyz_0,    \
                             to_xyyzz_0,    \
                             to_xyz_0,      \
                             to_xyzzz_0,    \
                             to_xz_0_xxx_0, \
                             to_xz_0_xxy_0, \
                             to_xz_0_xxz_0, \
                             to_xz_0_xyy_0, \
                             to_xz_0_xyz_0, \
                             to_xz_0_xzz_0, \
                             to_xz_0_yyy_0, \
                             to_xz_0_yyz_0, \
                             to_xz_0_yzz_0, \
                             to_xz_0_zzz_0, \
                             to_xzz_0,      \
                             to_xzzzz_0,    \
                             to_y_0,        \
                             to_yyz_0,      \
                             to_yzz_0,      \
                             to_z_0,        \
                             to_zzz_0 : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                const double tbe_0 = a_exp;

                to_xz_0_xxx_0[k] = -6.0 * to_xxz_0[k] * tbe_0 + 4.0 * to_xxxxz_0[k] * tbe_0 * tbe_0;

                to_xz_0_xxy_0[k] = -4.0 * to_xyz_0[k] * tbe_0 + 4.0 * to_xxxyz_0[k] * tbe_0 * tbe_0;

                to_xz_0_xxz_0[k] = 2.0 * to_x_0[k] - 4.0 * to_xzz_0[k] * tbe_0 - 2.0 * to_xxx_0[k] * tbe_0 + 4.0 * to_xxxzz_0[k] * tbe_0 * tbe_0;

                to_xz_0_xyy_0[k] = -2.0 * to_yyz_0[k] * tbe_0 + 4.0 * to_xxyyz_0[k] * tbe_0 * tbe_0;

                to_xz_0_xyz_0[k] = to_y_0[k] - 2.0 * to_yzz_0[k] * tbe_0 - 2.0 * to_xxy_0[k] * tbe_0 + 4.0 * to_xxyzz_0[k] * tbe_0 * tbe_0;

                to_xz_0_xzz_0[k] = 2.0 * to_z_0[k] - 2.0 * to_zzz_0[k] * tbe_0 - 4.0 * to_xxz_0[k] * tbe_0 + 4.0 * to_xxzzz_0[k] * tbe_0 * tbe_0;

                to_xz_0_yyy_0[k] = 4.0 * to_xyyyz_0[k] * tbe_0 * tbe_0;

                to_xz_0_yyz_0[k] = -2.0 * to_xyy_0[k] * tbe_0 + 4.0 * to_xyyzz_0[k] * tbe_0 * tbe_0;

                to_xz_0_yzz_0[k] = -4.0 * to_xyz_0[k] * tbe_0 + 4.0 * to_xyzzz_0[k] * tbe_0 * tbe_0;

                to_xz_0_zzz_0[k] = -6.0 * to_xzz_0[k] * tbe_0 + 4.0 * to_xzzzz_0[k] * tbe_0 * tbe_0;
            }

            // Set up 30-40 components of targeted buffer : FS

            auto to_yy_0_xxx_0 = pbuffer.data(idx_op_geom_200_fs + 3 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 0 * ket_comps + j);

            auto to_yy_0_xxy_0 = pbuffer.data(idx_op_geom_200_fs + 3 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 1 * ket_comps + j);

            auto to_yy_0_xxz_0 = pbuffer.data(idx_op_geom_200_fs + 3 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 2 * ket_comps + j);

            auto to_yy_0_xyy_0 = pbuffer.data(idx_op_geom_200_fs + 3 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 3 * ket_comps + j);

            auto to_yy_0_xyz_0 = pbuffer.data(idx_op_geom_200_fs + 3 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 4 * ket_comps + j);

            auto to_yy_0_xzz_0 = pbuffer.data(idx_op_geom_200_fs + 3 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 5 * ket_comps + j);

            auto to_yy_0_yyy_0 = pbuffer.data(idx_op_geom_200_fs + 3 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 6 * ket_comps + j);

            auto to_yy_0_yyz_0 = pbuffer.data(idx_op_geom_200_fs + 3 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 7 * ket_comps + j);

            auto to_yy_0_yzz_0 = pbuffer.data(idx_op_geom_200_fs + 3 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 8 * ket_comps + j);

            auto to_yy_0_zzz_0 = pbuffer.data(idx_op_geom_200_fs + 3 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 9 * ket_comps + j);

#pragma omp simd aligned(to_x_0,            \
                             to_xxx_0,      \
                             to_xxxyy_0,    \
                             to_xxy_0,      \
                             to_xxyyy_0,    \
                             to_xxyyz_0,    \
                             to_xxz_0,      \
                             to_xyy_0,      \
                             to_xyyyy_0,    \
                             to_xyyyz_0,    \
                             to_xyyzz_0,    \
                             to_xyz_0,      \
                             to_xzz_0,      \
                             to_y_0,        \
                             to_yy_0_xxx_0, \
                             to_yy_0_xxy_0, \
                             to_yy_0_xxz_0, \
                             to_yy_0_xyy_0, \
                             to_yy_0_xyz_0, \
                             to_yy_0_xzz_0, \
                             to_yy_0_yyy_0, \
                             to_yy_0_yyz_0, \
                             to_yy_0_yzz_0, \
                             to_yy_0_zzz_0, \
                             to_yyy_0,      \
                             to_yyyyy_0,    \
                             to_yyyyz_0,    \
                             to_yyyzz_0,    \
                             to_yyz_0,      \
                             to_yyzzz_0,    \
                             to_yzz_0,      \
                             to_z_0,        \
                             to_zzz_0 : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                const double tbe_0 = a_exp;

                to_yy_0_xxx_0[k] = -2.0 * to_xxx_0[k] * tbe_0 + 4.0 * to_xxxyy_0[k] * tbe_0 * tbe_0;

                to_yy_0_xxy_0[k] = -6.0 * to_xxy_0[k] * tbe_0 + 4.0 * to_xxyyy_0[k] * tbe_0 * tbe_0;

                to_yy_0_xxz_0[k] = -2.0 * to_xxz_0[k] * tbe_0 + 4.0 * to_xxyyz_0[k] * tbe_0 * tbe_0;

                to_yy_0_xyy_0[k] = 2.0 * to_x_0[k] - 10.0 * to_xyy_0[k] * tbe_0 + 4.0 * to_xyyyy_0[k] * tbe_0 * tbe_0;

                to_yy_0_xyz_0[k] = -6.0 * to_xyz_0[k] * tbe_0 + 4.0 * to_xyyyz_0[k] * tbe_0 * tbe_0;

                to_yy_0_xzz_0[k] = -2.0 * to_xzz_0[k] * tbe_0 + 4.0 * to_xyyzz_0[k] * tbe_0 * tbe_0;

                to_yy_0_yyy_0[k] = 6.0 * to_y_0[k] - 14.0 * to_yyy_0[k] * tbe_0 + 4.0 * to_yyyyy_0[k] * tbe_0 * tbe_0;

                to_yy_0_yyz_0[k] = 2.0 * to_z_0[k] - 10.0 * to_yyz_0[k] * tbe_0 + 4.0 * to_yyyyz_0[k] * tbe_0 * tbe_0;

                to_yy_0_yzz_0[k] = -6.0 * to_yzz_0[k] * tbe_0 + 4.0 * to_yyyzz_0[k] * tbe_0 * tbe_0;

                to_yy_0_zzz_0[k] = -2.0 * to_zzz_0[k] * tbe_0 + 4.0 * to_yyzzz_0[k] * tbe_0 * tbe_0;
            }

            // Set up 40-50 components of targeted buffer : FS

            auto to_yz_0_xxx_0 = pbuffer.data(idx_op_geom_200_fs + 4 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 0 * ket_comps + j);

            auto to_yz_0_xxy_0 = pbuffer.data(idx_op_geom_200_fs + 4 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 1 * ket_comps + j);

            auto to_yz_0_xxz_0 = pbuffer.data(idx_op_geom_200_fs + 4 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 2 * ket_comps + j);

            auto to_yz_0_xyy_0 = pbuffer.data(idx_op_geom_200_fs + 4 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 3 * ket_comps + j);

            auto to_yz_0_xyz_0 = pbuffer.data(idx_op_geom_200_fs + 4 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 4 * ket_comps + j);

            auto to_yz_0_xzz_0 = pbuffer.data(idx_op_geom_200_fs + 4 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 5 * ket_comps + j);

            auto to_yz_0_yyy_0 = pbuffer.data(idx_op_geom_200_fs + 4 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 6 * ket_comps + j);

            auto to_yz_0_yyz_0 = pbuffer.data(idx_op_geom_200_fs + 4 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 7 * ket_comps + j);

            auto to_yz_0_yzz_0 = pbuffer.data(idx_op_geom_200_fs + 4 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 8 * ket_comps + j);

            auto to_yz_0_zzz_0 = pbuffer.data(idx_op_geom_200_fs + 4 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 9 * ket_comps + j);

#pragma omp simd aligned(to_x_0,            \
                             to_xxxyz_0,    \
                             to_xxy_0,      \
                             to_xxyyz_0,    \
                             to_xxyzz_0,    \
                             to_xxz_0,      \
                             to_xyy_0,      \
                             to_xyyyz_0,    \
                             to_xyyzz_0,    \
                             to_xyz_0,      \
                             to_xyzzz_0,    \
                             to_xzz_0,      \
                             to_y_0,        \
                             to_yyy_0,      \
                             to_yyyyz_0,    \
                             to_yyyzz_0,    \
                             to_yyz_0,      \
                             to_yyzzz_0,    \
                             to_yz_0_xxx_0, \
                             to_yz_0_xxy_0, \
                             to_yz_0_xxz_0, \
                             to_yz_0_xyy_0, \
                             to_yz_0_xyz_0, \
                             to_yz_0_xzz_0, \
                             to_yz_0_yyy_0, \
                             to_yz_0_yyz_0, \
                             to_yz_0_yzz_0, \
                             to_yz_0_zzz_0, \
                             to_yzz_0,      \
                             to_yzzzz_0,    \
                             to_z_0,        \
                             to_zzz_0 : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                const double tbe_0 = a_exp;

                to_yz_0_xxx_0[k] = 4.0 * to_xxxyz_0[k] * tbe_0 * tbe_0;

                to_yz_0_xxy_0[k] = -2.0 * to_xxz_0[k] * tbe_0 + 4.0 * to_xxyyz_0[k] * tbe_0 * tbe_0;

                to_yz_0_xxz_0[k] = -2.0 * to_xxy_0[k] * tbe_0 + 4.0 * to_xxyzz_0[k] * tbe_0 * tbe_0;

                to_yz_0_xyy_0[k] = -4.0 * to_xyz_0[k] * tbe_0 + 4.0 * to_xyyyz_0[k] * tbe_0 * tbe_0;

                to_yz_0_xyz_0[k] = to_x_0[k] - 2.0 * to_xzz_0[k] * tbe_0 - 2.0 * to_xyy_0[k] * tbe_0 + 4.0 * to_xyyzz_0[k] * tbe_0 * tbe_0;

                to_yz_0_xzz_0[k] = -4.0 * to_xyz_0[k] * tbe_0 + 4.0 * to_xyzzz_0[k] * tbe_0 * tbe_0;

                to_yz_0_yyy_0[k] = -6.0 * to_yyz_0[k] * tbe_0 + 4.0 * to_yyyyz_0[k] * tbe_0 * tbe_0;

                to_yz_0_yyz_0[k] = 2.0 * to_y_0[k] - 4.0 * to_yzz_0[k] * tbe_0 - 2.0 * to_yyy_0[k] * tbe_0 + 4.0 * to_yyyzz_0[k] * tbe_0 * tbe_0;

                to_yz_0_yzz_0[k] = 2.0 * to_z_0[k] - 2.0 * to_zzz_0[k] * tbe_0 - 4.0 * to_yyz_0[k] * tbe_0 + 4.0 * to_yyzzz_0[k] * tbe_0 * tbe_0;

                to_yz_0_zzz_0[k] = -6.0 * to_yzz_0[k] * tbe_0 + 4.0 * to_yzzzz_0[k] * tbe_0 * tbe_0;
            }

            // Set up 50-60 components of targeted buffer : FS

            auto to_zz_0_xxx_0 = pbuffer.data(idx_op_geom_200_fs + 5 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 0 * ket_comps + j);

            auto to_zz_0_xxy_0 = pbuffer.data(idx_op_geom_200_fs + 5 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 1 * ket_comps + j);

            auto to_zz_0_xxz_0 = pbuffer.data(idx_op_geom_200_fs + 5 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 2 * ket_comps + j);

            auto to_zz_0_xyy_0 = pbuffer.data(idx_op_geom_200_fs + 5 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 3 * ket_comps + j);

            auto to_zz_0_xyz_0 = pbuffer.data(idx_op_geom_200_fs + 5 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 4 * ket_comps + j);

            auto to_zz_0_xzz_0 = pbuffer.data(idx_op_geom_200_fs + 5 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 5 * ket_comps + j);

            auto to_zz_0_yyy_0 = pbuffer.data(idx_op_geom_200_fs + 5 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 6 * ket_comps + j);

            auto to_zz_0_yyz_0 = pbuffer.data(idx_op_geom_200_fs + 5 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 7 * ket_comps + j);

            auto to_zz_0_yzz_0 = pbuffer.data(idx_op_geom_200_fs + 5 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 8 * ket_comps + j);

            auto to_zz_0_zzz_0 = pbuffer.data(idx_op_geom_200_fs + 5 * op_comps * 10 * ket_comps + i * 10 * ket_comps + 9 * ket_comps + j);

#pragma omp simd aligned(to_x_0,            \
                             to_xxx_0,      \
                             to_xxxzz_0,    \
                             to_xxy_0,      \
                             to_xxyzz_0,    \
                             to_xxz_0,      \
                             to_xxzzz_0,    \
                             to_xyy_0,      \
                             to_xyyzz_0,    \
                             to_xyz_0,      \
                             to_xyzzz_0,    \
                             to_xzz_0,      \
                             to_xzzzz_0,    \
                             to_y_0,        \
                             to_yyy_0,      \
                             to_yyyzz_0,    \
                             to_yyz_0,      \
                             to_yyzzz_0,    \
                             to_yzz_0,      \
                             to_yzzzz_0,    \
                             to_z_0,        \
                             to_zz_0_xxx_0, \
                             to_zz_0_xxy_0, \
                             to_zz_0_xxz_0, \
                             to_zz_0_xyy_0, \
                             to_zz_0_xyz_0, \
                             to_zz_0_xzz_0, \
                             to_zz_0_yyy_0, \
                             to_zz_0_yyz_0, \
                             to_zz_0_yzz_0, \
                             to_zz_0_zzz_0, \
                             to_zzz_0,      \
                             to_zzzzz_0 : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                const double tbe_0 = a_exp;

                to_zz_0_xxx_0[k] = -2.0 * to_xxx_0[k] * tbe_0 + 4.0 * to_xxxzz_0[k] * tbe_0 * tbe_0;

                to_zz_0_xxy_0[k] = -2.0 * to_xxy_0[k] * tbe_0 + 4.0 * to_xxyzz_0[k] * tbe_0 * tbe_0;

                to_zz_0_xxz_0[k] = -6.0 * to_xxz_0[k] * tbe_0 + 4.0 * to_xxzzz_0[k] * tbe_0 * tbe_0;

                to_zz_0_xyy_0[k] = -2.0 * to_xyy_0[k] * tbe_0 + 4.0 * to_xyyzz_0[k] * tbe_0 * tbe_0;

                to_zz_0_xyz_0[k] = -6.0 * to_xyz_0[k] * tbe_0 + 4.0 * to_xyzzz_0[k] * tbe_0 * tbe_0;

                to_zz_0_xzz_0[k] = 2.0 * to_x_0[k] - 10.0 * to_xzz_0[k] * tbe_0 + 4.0 * to_xzzzz_0[k] * tbe_0 * tbe_0;

                to_zz_0_yyy_0[k] = -2.0 * to_yyy_0[k] * tbe_0 + 4.0 * to_yyyzz_0[k] * tbe_0 * tbe_0;

                to_zz_0_yyz_0[k] = -6.0 * to_yyz_0[k] * tbe_0 + 4.0 * to_yyzzz_0[k] * tbe_0 * tbe_0;

                to_zz_0_yzz_0[k] = 2.0 * to_y_0[k] - 10.0 * to_yzz_0[k] * tbe_0 + 4.0 * to_yzzzz_0[k] * tbe_0 * tbe_0;

                to_zz_0_zzz_0[k] = 6.0 * to_z_0[k] - 14.0 * to_zzz_0[k] * tbe_0 + 4.0 * to_zzzzz_0[k] * tbe_0 * tbe_0;
            }
        }
    }
}

}  // namespace t2cgeom